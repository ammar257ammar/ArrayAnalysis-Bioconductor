# =============================================================================#
# ArrayAnalysis - affyAnalysisQC                                              #
# a tool for quality control and pre-processing of Affymetrix array data      #
#                                                                             #
# Copyright 2010-2011 BiGCaT Bioinformatics                                   #
#                                                                             #
# Licensed under the Apache License, Version 2.0 (the "License");             #
# you may not use this file except in compliance with the License.            #
# You may obtain a copy of the License at                                     #
#                                                                             #
# http://www.apache.org/licenses/LICENSE-2.0                                  #
#                                                                             #
# Unless required by applicable law or agreed to in writing, software         #
# distributed under the License is distributed on an "AS IS" BASIS,           #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #
# See the License for the specific language governing permissions and         #
# limitations under the License.                                              #
# =============================================================================#

#####################
## coverAndKeyPlot ##
#####################

#' Plot a cover sheet
#'
#' This function (from functions_imagesQC.R) plots a cover sheet,
#' and one or more key sheet indicating the links between CEL file names,
#' array names used in the plots, and to which experimental group they
#' belong based on the description file provided by the user, which has
#' been loaded into the description variable earlier in the main script,
#' passed as an argument (arrayGroup). For the key sheets, one sheet will
#' be created for every 35 arrays in the experiment.
#'
#' @param description (Status: required) The data.frame containing the
#' description file information (column 1: file names; column 2: names
#' to be used in the plots; column 3: experimental
#' groups the samples belong to)(datatype: data.frame)
#' @param refName (Status: optional, Default:"") dataset name. It is deduced
#' from the name of the zip file containing the CEL files when used from
#' arrayanalysis.org or as a GenePattern module.(datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width
#' (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height
#' (datatype: number)
#'
#' @return Two files, a file called ‘Cover_1’ that can be used as an opening
#' image. And, file(s) called ‘Description’, if needed followed by an index number,that
#' represent(s) the information as given in the input parameter.
#' @examples
#' # By default, the script will call:
#' # coverAndKeyPlot(description)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

coverAndKeyPlot <- function(description = NULL, refName = "", WIDTH = 1000,
                            HEIGHT = 1414) {
  if (is.null(description)) stop("The description parameters is required")

  # create a 'cover sheet'
  png(filename = "Cover_1.png", width = WIDTH, height = HEIGHT)
  plot(c(0, 2),
    type = "n", ann = FALSE, axes = FALSE,
    frame.plot = TRUE, xlim = c(0, 2), ylim = c(0, 2)
  )
  text(1, 1, "Quality Control & Pre-processing Evaluation\n\n\n\n", pos = 3, cex = 3)
  if (refName != "") {
    refName <- sub(
      "(_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}_\\d{2})", "",
      refName
    )
    text(1, 1, paste("of\n", refName, "\n\n"), pos = 3, cex = 2.5)
  }
  text(1, 1, "\n\n\n\n\n\n\n\n\nREPORT", pos = 3, cex = 5)
  dev.off()

  # create a page indicating the naming and grouping used (it corresponds to the
  # description file)
  linesPerPage <- 35
  noOfPages <- (dim(description)[1] %/% linesPerPage) +
    ((dim(description)[1] %% linesPerPage) > 0)
  for (j in 1:noOfPages) {
    png(
      filename = paste("Description", ifelse(refName != "", "_", ""), refName,
        ifelse(noOfPages > 1, paste("_", j, sep = ""), ""), ".png",
        sep = ""
      ),
      width = WIDTH, height = HEIGHT
    )
    par(oma = c(0, 0, 0, 0))
    plot(c(0, 2),
      type = "n", ann = FALSE, axes = FALSE,
      frame.plot = TRUE, xlim = c(0, 2), ylim = c(0, 2)
    )
    text(1, 1.90, "Array names and grouping", cex = 1.6, pos = 3)

    vertical_pos <- seq(1.8, 0.01, length.out = (linesPerPage + 1))
    line_coef <- (vertical_pos[1] - vertical_pos[2]) / 2
    k <- 1
    for (i in c(0, ((j - 1) * linesPerPage + 1):min(
      j * (linesPerPage),
      dim(description)[1]
    ))) {
      if (i == 0) {
        value <- colnames(description)
      } else {
        value <- description[i, ]
      }
      text(0.0, vertical_pos[k], value[1], pos = 4, cex = 1)
      text(0.9, vertical_pos[k], value[2], pos = 4, cex = 1)
      text(1.6, vertical_pos[k], value[3], pos = 4, cex = 1)
      abline(h = vertical_pos[k] - line_coef)
      k <- k + 1
    }
    dev.off()
  }
}


#################
## QCtablePlot ##
#################

#' Compute and plot a table of QC statistics
#'
#' This function (from functions_imagesQC.R) computes and plots a table
#' of QC statistics based on the qc (simpleaffy Bioconductor package)
#' and yaqc (yaqcaffy Bioconductor package) functions, which generally
#' only work for chiptypes with perfect match and mismatch probes, but
#' even not for all of those. As such, when these statistics are not
#' provided as parameters, trys are used in this function to compute
#' them internally. Values for which the try fails are not computed,
#' but the script will continue after giving a warning.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param quality (Status: optional, Default:NULL) object obtained by
#' calling the qc function (simpleaffy).When not provided, it
#' is computed within the function.(datatype: QCStats)
#' @param sprep (Status: optional, Default:NULL) A matrix of 3’probe
#' intensities for dap, thr, lys, and phe, taken from an object
#' of class YAQCStats (yaqc function, yaqcaffy). When not provided,
#' it is computed within the function. (datatype: matrix)
#' @param lys (Status: optional, Default:NULL) Matrix of A, M, P calls
#' for the 3’ probeset of Lys on each array, based on results from
#' the detection.p.val function (simpleaffy) . When not provided,
#' it is computed within the function. (datatype: matrix)
#' @param samplePrep (Status: optional, Default:TRUE) Does the table
#' have to contain sample preparation QC statistics? (datatype: logical)
#' @param ratio (Status: optional, Default:TRUE) Does the table have to contain 3’/5’ ratio statistics? (datatype: logical)
#' @param hybrid (Status: optional, Default:TRUE) Does the table have
#' to contain hybridisation QC statistics? (datatype: logical)
#' @param percPres (Status: optional, Default:TRUE) Does the table have
#' to contain percentage present QC statistics? (datatype: logical)
#' @param bgPlot (Status: optional, Default:TRUE) Does the table have
#' to contain background signal intensity QC statistics? (datatype: logical)
#' @param scaleFact (Status: optional, Default:TRUE) Does the table have
#' to contain scale factor QC statistics?(datatype: logical)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @return A PNG file called ‘QCtable’ with the requested statistics
#' @examples
#' # By default, the script will call:
#' # computeQCtable(rawData, quality, sprep, lys, samplePrep = samplePrep,
#' # ratio = ratio, hybrid = hybrid, percPres = percPres, bgPlot = bgPlot, scaleFact = scaleFact)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

QCtablePlot <- function(Data, quality = NULL, sprep = NULL, lys = NULL,
                        samplePrep = TRUE, ratio = TRUE, hybrid = TRUE, percPres = TRUE, bgPlot = TRUE,
                        scaleFact = TRUE, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24) {
  rows <- NULL
  indicator <- NULL
  indicCol <- NULL

  if (sum(ratio, hybrid, percPres, bgPlot, scaleFact) > 0 && is.null(quality)) {
    try(quality <- qc(Data), TRUE) # calculate Affymetrix quality data for PMMM (use MAS5)
    if (is.null(quality)) {
      warning("Plots based on the simpleaffy qc function cannot be created for this chip type")
    }
  }

  if (samplePrep) {
    if (is.null(sprep)) {
      # find the data
      try(yack <- yaqc(Data), TRUE)
      if (exists("yack")) {
        spnames <- rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", # only 3' !
          rownames(yack@morespikes),
          ignore.case = TRUE
        ), ])
        sprep <- t(yack@morespikes[spnames, ])
      } else {
        warning("Values based on the yaqc function cannot be computed for this chip type")
      }
    }


    if (is.null(lys)) {
      try(
        {
          calls <- detection.p.val(Data)$call
          lys <- calls[rownames(calls)[grep("lys.*3", rownames(calls), ignore.case = TRUE)], ]
          rm(calls)
        },
        TRUE
      )
      if (is.null(lys)) {
        warning("Values based on the detection.p.val function cannot be computed for this chip type")
      } else {
        if (length(lys) > length(sampleNames(Data))) {
          lys <- lys[1, ]
        }
      }
    }

    if (!is.null(sprep)) {
      print("indicators for sample prep controls")

      # test for Lys < Phe < Thr < Dap
      t1 <- (sprep[, 1] < sprep[, 2] & sprep[, 2] < sprep[, 3]
      & sprep[, 3] < sprep[, 4])
      Col1 <- t1
      Col1[t1 == TRUE] <- "blue"
      Col1[t1 == FALSE] <- "red"
      Text1 <- t1
      Text1[t1 == TRUE] <- "T"
      Text1[t1 == FALSE] <- "F"

      if (exists("rows")) {
        rows <- c(rows, "Sample Prep\nLys<Phe<Thr<Dap if\ncontrols were spiked-in")
      } else {
        rows <- c("Sample Prep\nLys<Phe<Thr<Dap if\ncontrols were spiked-in")
      }

      if (exists("indicator")) {
        indicator <- rbind(indicator, Text1)
      } else {
        indicator <- rbind(Text1)
      }
      if (exists("indicCol")) {
        indicCol <- rbind(indicCol, Col1)
      } else {
        indicCol <- rbind(Col1)
      }
    }

    if (!is.null(lys)) {
      # test for Lys presence
      t2 <- lys == "P"
      Col2 <- t2
      Col2[t2 == TRUE] <- "blue"
      Col2[t2 == FALSE] <- "red"

      if (exists("rows")) {
        rows <- c(rows, "Sample Prep\nLys=Present")
      } else {
        rows <- c("Sample Prep\nLys=Present")
      }
      if (exists("indicator")) {
        indicator <- rbind(indicator, lys)
      } else {
        indicator <- rbind(lys)
      }
      if (exists("indicCol")) {
        indicCol <- rbind(indicCol, Col2)
      } else {
        indicCol <- rbind(Col2)
      }
    }
  }

  if (!is.null(quality)) {
    if (ratio) {
      print("indicators for beta-actin & GAPDH 3'/5' ratio")

      ratio35_beta <- quality@qc.probes[, 1] / quality@qc.probes[, 3]
      ratio35_GAPDH <- quality@qc.probes[, 4] / quality@qc.probes[, 6]


      if (exists("rows")) {
        rows <- c(
          rows, "3'/5' beta-actin\n(cutoff=3)",
          "3'/5' GAPDH\n(cutoff=1.25)"
        )
      } else {
        rows <- c("3'/5' beta-actin\n(cutoff=3)", "3'/5' GAPDH\n(cutoff=1.25)")
      }

      # beta-actin
      t1 <- ratio35_beta <= 3
      Col1 <- t1
      Col1[t1 == TRUE] <- "blue"
      Col1[t1 == FALSE] <- "red"


      # GAPDH
      t2 <- ratio35_GAPDH <= 1.25
      Col2 <- t2
      Col2[t2 == TRUE] <- "blue"
      Col2[t2 == FALSE] <- "red"

      if (exists("indicator")) {
        indicator <- rbind(indicator, round(ratio35_beta, 2), round(ratio35_GAPDH, 2))
      } else {
        indicator <- rbind(round(ratio35_beta, 2), round(ratio35_GAPDH, 2))
      }
      if (exists("indicCol")) {
        indicCol <- rbind(indicCol, Col1, Col2)
      } else {
        indicCol <- rbind(Col1, Col2)
      }
    }

    if (hybrid) {
      print("indicators for spike-in hybridization controls")

      # test for BioB < BioC < BioD < CreX
      t1 <- (quality@spikes[, 1] < quality@spikes[, 2]
      & quality@spikes[, 2] < quality@spikes[, 3]
      & quality@spikes[, 3] < quality@spikes[, 4])
      Col1 <- t1
      Col1[t1 == TRUE] <- "blue"
      Col1[t1 == FALSE] <- "red"
      Text1 <- t1
      Text1[t1 == TRUE] <- "T"
      Text1[t1 == FALSE] <- "F"

      # test for BioB presence
      t2 <- quality@bioBCalls == "P"
      Col2 <- t2
      Col2[t2 == TRUE] <- "blue"
      Col2[t2 == FALSE] <- "red"

      if (exists("rows")) {
        rows <- c(
          rows, "Hybridization\nBioB<BioC<BioD<CreX",
          "Hybridization\nBioB=Present"
        )
      } else {
        rows <- c(
          "Hybridization\nBioB<BioC<BioD<CreX",
          "Hybridization\nBioB=Present"
        )
      }
      if (exists("indicator")) {
        indicator <- rbind(indicator, Text1, quality@bioBCalls)
      } else {
        indicator <- rbind(Text1, quality@bioBCalls)
      }
      if (exists("indicCol")) {
        indicCol <- rbind(indicCol, Col1, Col2)
      } else {
        indicCol <- rbind(Col1, Col2)
      }
    }

    if (percPres) {
      print("indicators for percent present")

      # Define the reference (it should minimize the number of outliers)
      ppMin <- min(quality@percent.present)
      ppMax <- max(quality@percent.present)
      ppMean <- mean(quality@percent.present)

      ref <- cbind(c(ppMean - 5, ppMin, ppMax - 10), c(ppMean + 5, ppMin + 10, ppMax))
      t1 <- (quality@percent.present >= (ppMean - 5) & quality@percent.present <= (ppMean + 5))
      t2 <- (quality@percent.present >= ppMin & quality@percent.present <= (ppMin + 10))
      t3 <- (quality@percent.present >= (ppMax - 10) & quality@percent.present <= ppMax)
      outliers <- c(length(t1[t1 == FALSE]), length(t2[t2 == FALSE]), length(t3[t3 == FALSE]))
      if (outliers[1] == 0) {
        ref <- ref[1, ]
      } else {
        ref <- ref[outliers == min(outliers), ] # less outliers
        if (length(ref) > 3) {
          ref <- ref[1, ]
        }
      }
      testMin <- ref[1]
      testMax <- ref[2]
      t1 <- (quality@percent.present >= testMin & quality@percent.present <= testMax)

      Col <- t1
      Col[t1 == TRUE] <- "blue"
      Col[t1 == FALSE] <- "red"

      if (exists("rows")) {
        rows <- c(rows, "Percent Present\nspread<=10%")
      } else {
        rows <- c("Percent Present\nspread<=10%")
      }
      if (exists("indicator")) {
        indicator <- rbind(indicator, paste(round(quality@percent.present, 0), "%"))
      } else {
        indicator <- paste(round(quality@percent.present, 0), "%")
      }
      if (exists("indicCol")) {
        indicCol <- rbind(indicCol, Col)
      } else {
        indicCol <- rbind(Col)
      }
    }

    if (bgPlot) {
      print("indicators for background intensities")

      # Define the reference for the grey rectangle (it should minimize the number of outliers)
      bgMean <- mean(quality@average.background)
      bgMin <- min(quality@average.background)
      bgMax <- max(quality@average.background)
      ref <- cbind(c(bgMean - 10, bgMin, bgMax - 20), c(bgMean + 10, bgMin + 20, bgMax))
      t1 <- (quality@average.background >= (bgMean - 10) & quality@average.background <= (bgMean + 10))
      t2 <- (quality@average.background >= bgMin & quality@average.background <= (bgMin + 20))
      t3 <- (quality@average.background >= (bgMax - 20) & quality@average.background <= bgMax)
      outliers <- c(length(t1[t1 == FALSE]), length(t2[t2 == FALSE]), length(t3[t3 == FALSE]))
      if (outliers[1] == 0) {
        ref <- ref[1, ]
      } else {
        ref <- ref[outliers == min(outliers), ] # less outliers
        if (length(ref) > 3) {
          ref <- ref[1, ]
        }
      }
      testMin <- ref[1]
      testMax <- ref[2]

      t1 <- (quality@average.background >= testMin & quality@average.background <= testMax)

      Col <- t1
      Col[t1 == TRUE] <- "blue"
      Col[t1 == FALSE] <- "red"

      if (exists("rows")) {
        rows <- c(rows, "Background\nspread<=20%")
      } else {
        rows <- c("Background\nspread<=20%")
      }
      if (exists("indicator")) {
        indicator <- rbind(indicator, round(quality@average.background, 0))
      } else {
        indicator <- round(quality@average.background, 0)
      }
      if (exists("indicCol")) {
        indicCol <- rbind(indicCol, Col)
      } else {
        indicCol <- rbind(Col)
      }
    }

    if (scaleFact) {
      print("indicators for scale factors")
      # Define the reference for the grey rectangle (it should minimize the number of outliers)
      sfMin <- min(log2(quality@scale.factors))
      sfMax <- max(log2(quality@scale.factors))
      sfMean <- mean(log2(quality@scale.factors))
      ref <- cbind(c(sfMean - 1.5, sfMin, sfMax - 3), c(sfMean + 1.5, sfMin + 3, sfMax))
      t1 <- (log2(quality@scale.factors) >= (sfMean - 1.5) & log2(quality@scale.factors) <= (sfMean + 1.5))
      t2 <- (log2(quality@scale.factors) >= sfMin & log2(quality@scale.factors) <= (sfMin + 3))
      t3 <- (log2(quality@scale.factors) >= (sfMax - 3) & log2(quality@scale.factors) <= sfMax)
      outliers <- c(length(t1[t1 == FALSE]), length(t2[t2 == FALSE]), length(t3[t3 == FALSE]))
      if (outliers[1] == 0) {
        ref <- ref[1, ]
      } else {
        ref <- ref[outliers == min(outliers), ] # less outliers
        if (length(ref) > 3) {
          ref <- ref[1, ]
        }
      }
      testMin <- ref[1]
      testMax <- ref[2]

      t1 <- (log2(quality@scale.factors) >= testMin & log2(quality@scale.factors) <= testMax)

      Col <- t1
      Col[t1 == TRUE] <- "blue"
      Col[t1 == FALSE] <- "red"

      #  if(sfMax - sfMin <= 3){
      #    Col <- rep("blue", length(sampleNames(Data)))
      #  }else{
      #    Col <- rep("red", length(sampleNames(Data)))
      #  }
      if (exists("rows")) {
        rows <- c(rows, "Log Scale Factor\nspread<=3")
      } else {
        rows <- c("Log Scale Factor\nspread<=3")
      }
      if (exists("indicator")) {
        indicator <- rbind(indicator, round(log2(quality@scale.factors), 2))
      } else {
        indicator <- round(log2(quality@scale.factors), 2)
      }
      if (exists("indicCol")) {
        indicCol <- rbind(indicCol, Col)
      } else {
        indicCol <- rbind(Col)
      }
    }
  }

  # par(parStart)
  if (!is.null(rows)) { # and as such indicator and indicCol
    png(filename = "QCtable.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    par(oma = c(2, 4.5, 0, 0))
    plot(c(0, length(rows)),
      type = "n", ann = FALSE, axes = FALSE,
      frame.plot = FALSE, ylim = c(0, length(sampleNames(Data))),
      xlim = c(0, length(rows))
    )
    grid(ny = length(sampleNames(Data)) + 1, nx = length(rows) + 1, equilogs = FALSE)
    for (i in 1:length(rows)) {
      for (j in 1:length(sampleNames(Data))) {
        if (length(rows) > 1 && length(sampleNames(Data)) > 0) {
          text(i, length(sampleNames(Data)) + 1 - j, as.character(indicator[i, j]), col = indicCol[i, j], cex = .6)
        } else {
          text(i, length(sampleNames(Data)) + 1 - j, as.character(indicator[ifelse(length(rows) > 1, j, i)]), col = indicCol[ifelse(length(rows) > 1, j, i)], cex = .6)
        }
      }
    }
    title(main = "Summary of raw data quality indicators", cex = 0.8)
    mtext("blue = \"within\" / red = \"out of\" recommended cut-off",
      side = 3,
      cex = 0.6
    )
    axis(1, at = 1:length(rows), las = 2, labels = rows, cex.axis = 0.5)
    axis(2,
      at = 1:length(sampleNames(Data)), las = 2, labels = rev(substr(sampleNames(Data), 1, 25)),
      cex.axis = .5
    )
    dev.off()
  }
}


####################
## samplePrepPlot ##
####################

#' Create an image of sample preparation controls
#'
#' This function (from functions_imagesQC.R) creates an image of
#' sample preparation controls based on the yaqc
#' (yaqcaffy Bioconductor package) function, which generally only
#' work for chiptypes with perfect match and mismatch probes, but
#' even not for all of those. As such, when these statistics are
#' not provided as parameters, trys are used in this function to
#' compute them internally. Values for which the try fails are not
#' computed, but the script will continue after giving a warning.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param sprep (Status: optional, Default:NULL) A matrix of 3’probe
#' intensities for dap, thr, lys, and phe, taken from an object
#' of class YAQCStats (yaqc function, yaqcaffy). When not provided,
#' it is computed within the function. (datatype: matrix)
#' @param lys (Status: optional, Default:NULL) Matrix of A, M, P calls
#' for the 3’ probeset of Lys on each array, based on results from
#' the detection.p.val function (simpleaffy) . When not provided,
#' it is computed within the function. (datatype: matrix)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image of sample preparation controls, called "RawDataSamplePrepControl"
#' @examples
#' # samplePrepPlot(rawData,sprep,lys,plotColors)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

samplePrepPlot <- function(Data, sprep = NULL, lys = NULL, plotColors = NULL,
                           WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")

  if (is.null(lys)) {
    try(
      {
        calls <- detection.p.val(Data)$call
        lys <- calls[rownames(calls)[grep("lys.*3", rownames(calls), ignore.case = TRUE)], ]
        rm(calls)
      },
      TRUE
    )
    if (is.null(lys)) {
      warning("Plots based on the detection.p.val function cannot be created for this chip type")
    } else {
      if (length(lys) > length(sampleNames(Data))) {
        lys <- lys[1, ]
      }
    }
  }

  if (!is.null(lys) && is.null(sprep)) {
    try(yack <- yaqc(Data), TRUE)
    if (exists("yack")) {
      spnames <- rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", # only 3' !
        rownames(yack@morespikes),
        ignore.case = TRUE
      ), ])
      sprep <- t(yack@morespikes[spnames, ])
    } else {
      warning("Values based on the yaqc function cannot be computed for this chip type")
    }
  }

  if (!is.null(sprep) && !is.null(lys)) {
    # par(parStart)
    lmin <- min(sprep)
    lmax <- max(sprep) + 50

    # create the plot

    png(filename = "RawDataSamplePrepControl.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    if (length(sampleNames(Data)) < MAXARRAY) {
      par(mfrow = c(1, 1), oma = c(13, 0, 0, 0))
    } else {
      par(mfrow = c(1, 1), oma = c(0.1, 0, 0, 0))
    }
    plot(c(lmin, lmax),
      type = "n", ann = FALSE, axes = FALSE,
      frame.plot = TRUE, ylim = c(lmin, lmax), xlim = c(1, 5)
    )

    dummy <- apply(
      cbind(sprep, plotColors, 1:length(sampleNames(Data))), 1,
      function(x) {
        par(new = TRUE)
        plot(as.numeric(x[1:(length(x) - 2)]),
          type = "b", ann = FALSE, axes = FALSE, pch = ".", cex = 4, lwd = 2,
          col = x[length(x) - 1], lty = as.numeric(x[length(x)]),
          ylim = c(lmin, lmax), xlim = c(1, 5)
        )
      }
    )

    title(
      main = "Spike-in Sample Prep controls intensities and calls",
      ylab = "Intensity", cex.lab = 0.8
    )

    label3 <- ""
    title1 <- "Intensity: OK (Lys < Phe < Thr < Dap for all arrays)"
    title2 <- "\nLys Present calls: OK (all Lys are called present)"

    # test for Lys < Phe < Thr < Dap
    bad1 <- colnames(t(sprep))[!(sprep[, 1] < sprep[, 2]
    & sprep[, 2] < sprep[, 3]
    & sprep[, 3] < sprep[, 4])]
    if (length(bad1) > 0) {
      title1 <- paste("Intensity: not OK (some array", ifelse(length(bad1) > 1, "s", ""),
        " do", ifelse(length(bad1) > 1, "", "es"), " not follow the rule Lys < Phe < Thr < Dap;",
        "\nmaybe no Sample Prep controls were spiked on these arrays.)",
        sep = ""
      )
    }

    # test for Lys presence
    bad2 <- names(lys[lys != "P"])
    bad2 <- gsub(".present", "", bad2)
    if (length(bad2) > 0) {
      title2 <- paste("Lys Present calls: ", length(bad2),
        " Lys not called present",
        sep = ""
      )
    }
    if (length(lys[lys == "P"]) == 0) {
      label3 <-
        "\nApparently no Sample Prep controls were spiked on these arrays."
    }

    title(xlab = c(title1, title2), cex.lab = 0.8)
    legend("topleft", paste(
      "Lys present calls = ",
      length(lys[lys == "P"]), "/", length(lys), label3
    ), bty = "n", cex = 0.7)
    if (length(sampleNames(Data)) < (MAXARRAY + 20)) {
      cexval <- 0.65
    } else {
      cexval <- 0.45
    }
    legend("topright", substr(sampleNames(Data), 1, 20),
      lwd = 2,
      col = plotColors, cex = cexval, bty = "n", lty = 1:length(sampleNames(Data))
    )
    par(cex.axis = 0.8)
    axis(2)
    axis(1, at = 1:4, labels = c("Lys", "Phe", "Thr", "Dap"))

    if (length(bad2) > 0) {
      text(c(rep(0.87, length(lys[lys != "P"]))), sprep[lys != "P", 1],
        lys[lys != "P"],
        pos = 4, offset = 0.2, cex = 0.7, col = "red"
      )
    }
    dev.off()
  } else {
    warning("Spike-in sample prep plot is cannot be computed for this chip type")
  }
}


###############
## ratioPlot ##
###############

#' Create an image of beta-actin and GAPDH 3'/5' ratios
#'
#' This function (from functions_imagesQC.R) creates an image of
#' beta-actin and GAPDH 3'/5' ratios based on the qc
#' (simpleaffy Bioconductor package) function, which generally
#' only work for chiptypes with perfect match and mismatch probes,
#' but even not for all of those. As such, when these statistics
#' are not provided as parameters, trys are used in this function
#' to compute them internally. Values for which the try fails are
#' not computed, but the script will continue after giving a warning.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param quality (Status: optional, Default:NULL) object obtained by
#' calling the qc function (simpleaffy).When not provided, it
#' is computed within the function.(datatype: QCStats)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return Two images. A PNG image of 3’/5’ratios for beta actin, called ‘RawData53ratioPlot_beta-actin’.
#' A PNG image of 3’/5’ ratios for GAPDH, called ‘RawData53ratioPlot_GAPDH’
#' @examples
#' # ratioPlot(rawData, quality=quality, experimentFactor, plotColors, legendColors)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

ratioPlot <- function(Data, quality = NULL, experimentFactor = NULL, plotColors = NULL, legendColors = NULL,
                      WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if (is.null(legendColors)) stop("the 'legendColors' parameter is required")
  if (is.null(quality)) try(quality <- qc(Data), TRUE)

  if (!is.null(quality)) {
    # par(parStart)

    plotFun <- function(i, j, cutoff, Cname) {
      ratio35 <- quality@qc.probes[, i] / quality@qc.probes[, j]
      ratioM <- quality@qc.probes[, i] / quality@qc.probes[, i + 1]

      cMin <- 0
      cMax <- max(max(cbind(ratio35, ratioM)) + 2, 5)
      if (Cname == "beta-actin") {
        symbol <- c(17, 2) # triangles for beta-actin
      } else {
        symbol <- c(19, 1) # circles for GAPDH (as in the simpleaffy QC report)
      }

      par(mfrow = c(1, 2), oma = c(17, 0.5, 0, 0.5), cex.axis = 0.8)
      plot(ratio35,
        type = "n", ann = FALSE, axes = FALSE,
        frame.plot = TRUE, pch = ".", cex = 10, ylim = c(cMin, cMax)
      )
      rect(0, 0, length(ratio35) + 1, cutoff, col = gray(0.9), border = FALSE)

      for (k in 0:cMax) {
        abline(h = k, lty = 2, col = gray(0.8))
      }

      par(new = TRUE)
      plot(ratio35,
        type = "h", ann = FALSE, axes = FALSE,
        frame.plot = TRUE, lwd = 3, col = plotColors, ylim = c(cMin, cMax)
      )
      par(new = TRUE)
      plot(ratio35,
        type = "p", ann = FALSE, axes = FALSE,
        frame.plot = TRUE, pch = symbol[1], col = plotColors, ylim = c(cMin, cMax)
      )
      par(new = TRUE)
      plot(ratioM,
        type = "p", ann = FALSE, axes = FALSE,
        frame.plot = TRUE, pch = symbol[2], col = "black", ylim = c(cMin, cMax)
      )

      title(main = paste("RNA degradation of", Cname), ylab = "3'/5' and 3'/M ratios")
      axis(2)
      par(cex.axis = 0.65)
      if (length(sampleNames(Data)) < (MAXARRAY / 2)) { # array names not reported if more than 20 arrays
        axis(1, at = 1:length(ratio35), las = 2, labels = sampleNames(Data))
      }
      if (length(levels(experimentFactor)) > 1) {
        legend("topright", levels(experimentFactor),
          col = legendColors, fill = legendColors, bty = "n", cex = 0.55
        )
      }
      legend("topleft", c(paste("3\'/5\' ratio (max=", round(max(ratio35), 2), ")",
        sep = ""
      ), paste("3'/M ratio (max=", round(max(ratioM), 2), ")", sep = "")),
      col = c(gray(0.3), "black"), pch = symbol, cex = 0.55, bty = "o"
      )
      t1 <- ratio35 <= cutoff
      if (length(t1[t1 == FALSE]) > 0 && length(sampleNames(Data)) >= (MAXARRAY / 2)) {
        if (length(t1[t1 == FALSE]) < 25) {
          textlab <- sampleNames(Data)[t1 == FALSE]
          legend("topleft", c(rep("", 4), "Outliers, from left to right", textlab),
            text.col = c(rep(1, 5), rep("red", length(textlab))), bty = "n", cex = 0.38
          )
        } else {
          mtext("Too many arrays and outliers; \nsee details on the QC table", side = 1, cex = 0.8, line = 3)
        }
      }
      mtext(paste("Ratios should stand within the grey rectangle [0,", cutoff, "]"),
        side = 4, font = 1, cex = 0.7
      )

      par(cex.axis = 0.8)
      boxplot(cbind(ratio35, ratioM), axes = FALSE, frame.plot = TRUE)
      title(main = paste("Boxplot of", Cname, "ratios"))
      par(cex.lab = 0.8)
      if (max(ratio35) < cutoff) {
        title(xlab = paste(Cname, " QC: OK (all 3'/5' ratios < ", cutoff, ")", sep = ""))
      } else {
        title(xlab = paste(Cname, " QC: not OK (some 3'/5' ratios >", cutoff,
          ")\nNote that the threshold of ", cutoff,
          " \nwas determined for Homo Sapiens.",
          sep = ""
        ), line = 4)
      }

      axis(1, at = 1:2, labels = c("3'/5' ratio", "3'/M ratio"))
      axis(2)
    }

    # Graph creation: two separated graphs
    png(filename = "RawData53ratioPlot_beta-actin.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    plotFun(1, 3, 3, "beta-actin")
    dev.off()
    png(filename = "RawData53ratioPlot_GAPDH.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    plotFun(4, 6, 1.25, "GAPDH")
    dev.off()
  } else {
    warning("3'/5' ratio plot is not computed for this chip type")
  }
}

################
## RNAdegPlot ##
################

#' Creates an image of overall RNA degradation for all arrays
#'
#' This function (from functions_imagesQC.R) creates an image of
#' overall RNA degradation for all arrays. It calls the function
#' AffyRNAdeg (\href{http://www.bioconductor.org/help/bioc-views/release/bioc/html/affy.html}{affy} Bioconductor package).
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param Data.rnadeg (Status: optional, Default:NULL) List as obtained
#' by calling the AffyRNAdeg function (affy). When not provided,
#' it is computed internally within this function. (datatype: list)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image of overall RNA degradation, called ‘RawDataRNAdegradationPlot’
#' @examples
#' # RNAdegPlot(rawData,plotColors=plotColors)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

RNAdegPlot <- function(Data, Data.rnadeg = NULL, plotColors = NULL,
                       WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if (is.null(Data.rnadeg)) Data.rnadeg <- AffyRNAdeg(Data)

  png(filename = "RawDataRNAdegradation.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)

  if (length(sampleNames(Data)) < MAXARRAY) {
    cexval <- 0.65
    par(lwd = 2, oma = c(13, 0, 0, 0))
  } else {
    cexval <- 0.45
    par(lwd = 2, oma = c(0.1, 0, 0, 0))
  }
  par(mar = c(4, 4, 4, 0), cex.axis = 0.6, cex.lab = 0.75)
  layout(matrix(c(1, 2), 1, 2, byrow = TRUE), c(2, 1), 1, FALSE)
  plotAffyRNAdeg(Data.rnadeg, cols = plotColors, lty = 1:length(sampleNames(Data)))
  par(mar = c(4, 0, 4, 0), cex.axis = 0.01, cex.lab = 0.01)
  plot(1, type = "n", ann = FALSE, axes = FALSE, frame.plot = FALSE)
  legend("topleft", sampleNames(Data), lwd = 2, col = plotColors, cex = cexval, bty = "n", lty = 1)
  dev.off()
}


################
## hybridPlot ##
################

#' Create an image of hybridisation controls
#'
#' This function (from functions_imagesQC.R) creates an image of
#' hybridisation controls based on the qc
#' (\href{http://bioinformatics.picr.man.ac.uk/research/software/simpleaffy/index.html}{simpleaffy} Bioconductor package)
#' function, which generally only work for chiptypes with perfect match
#' and mismatch probes, but even not for all of those. As such, when
#' these statistics are not provided as parameters, trys are used in
#' this function to compute them internally. Values for which the try
#' fails are not computed, but the script will continue after giving a warning.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param quality (Status: optional, Default:NULL) object obtained by
#' calling the qc function (simpleaffy).When not provided, it
#' is computed within the function.(datatype: QCStats)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return Two images. A PNG image of hybridisation controls, called "RawDataSpike-inPlot"
#' @examples
#' # hybridPlot(rawData,quality=quality,plotColors)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

hybridPlot <- function(Data, quality = NULL, plotColors = NULL,
                       WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if (is.null(quality)) try(quality <- qc(Data), TRUE)

  if (!is.null(quality)) {
    png(filename = "RawDataSpikeinHybridControl.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    if (length(sampleNames(Data)) < MAXARRAY) {
      par(mfrow = c(1, 1), oma = c(13, 0, 0, 0))
    } else {
      par(mfrow = c(1, 1), oma = c(0.1, 0, 0, 0))
    }
    plot(c(min(quality@spikes), max(quality@spikes)),
      type = "n", ann = FALSE, axes = FALSE, frame.plot = TRUE,
      ylim = c(min(quality@spikes), max(quality@spikes)),
      xlim = c(1, 5)
    )

    dummy <- apply(
      cbind(quality@spikes, plotColors, 1:length(sampleNames(Data))), 1,
      function(x) {
        par(new = TRUE)
        plot(as.numeric(x[1:(length(x) - 2)]),
          type = "b", ann = FALSE, axes = FALSE, pch = ".", cex = 4, lwd = 2,
          col = x[length(x) - 1], lty = as.numeric(x[length(x)]),
          ylim = c(min(quality@spikes), max(quality@spikes)), xlim = c(1, 5)
        )
      }
    )

    title(
      main = "Spike-in Hybridization controls intensities and calls",
      ylab = "Intensity", cex.lab = 0.8
    )
    label3 <- ""
    title1 <- "Intensities: OK (bioB < bioC < bioD < creX for all arrays)"
    title2 <- "BioB Present calls: OK (indeed all bioB are called present)"

    # test for bioB<bioC<bioD<creX
    bad1 <- colnames(t(quality@spikes))[!(quality@spikes[, 1] < quality@spikes[, 2]
    & quality@spikes[, 2] < quality@spikes[, 3]
    & quality@spikes[, 3] < quality@spikes[, 4])]
    if (length(bad1) > 0) {
      title1 <- paste("Intensity: not OK (some array", ifelse(length(bad1) > 1, "s", ""),
        " do", ifelse(length(bad1) > 1, "", "es"), " not follow the rule bioB < bioC < bioD < creX)",
        sep = ""
      )
    }

    # test for bioB presence
    bad2 <- gsub(".present", "", names(quality@bioBCalls[quality@bioBCalls != "P"]))
    if (length(bad2) > 0) {
      title2 <- paste(
        "BioB present calls: not OK (", length(bad2), "/",
        length(quality@bioBCalls), " bioB not called present)"
      )
    }
    if (length(quality@bioBCalls[quality@bioBCalls == "P"]) == 0) {
      label3 <-
        "\nApparently no Hybridization controls were spiked on these arrays."
    }

    title(xlab = c(title1, title2), cex.lab = 0.8)
    legend("topleft", paste(
      "bioB present calls = ",
      length(quality@bioBCalls[quality@bioBCalls == "P"]), "/",
      length(quality@bioBCalls), label3
    ), bty = "n", cex = 0.7)

    if (length(sampleNames(Data)) < (MAXARRAY + 20)) {
      cexval <- 0.65
    } else {
      cexval <- 0.45
    }

    legend("bottomright", substr(sampleNames(Data), 1, 20),
      lwd = 2,
      col = plotColors, cex = cexval, bty = "n", lty = 1:length(sampleNames(Data))
    )
    par(cex.axis = 0.8)
    axis(2)
    axis(1, at = 1:4, labels = names(colnames(quality@spikes)))

    if (length(bad2) > 0) {
      text(c(rep(0.89, length(quality@bioBCalls[quality@bioBCalls != "P"]))),
        quality@spikes[quality@bioBCalls != "P", 1],
        quality@bioBCalls[quality@bioBCalls != "P"],
        pos = 4, offset = 0.2,
        cex = 0.8, col = "red"
      )
    }
    dev.off()
  } else {
    warning("Spike-in hybridization plot is not computed for this chip type")
  }
}


##################
## percPresPlot ##
##################

#' Creates an image of the percentage present values of each array
#'
#' This function (from functions_imagesQC.R) creates an image of the
#' percentage present values of each array based on the qc
#' (\href{http://bioinformatics.picr.man.ac.uk/research/software/simpleaffy/index.html}{simpleaffy} Bioconductor package) function, which generally
#' only work for chiptypes with perfect match and mismatch probes,
#' but even not for all of those. As such, when these statistics
#' are not provided as parameters, trys are used in this function
#' to compute them internally. Values for which the try fails are
#' not computed, but the script will continue after giving a warning.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param quality (Status: optional, Default:NULL) object obtained by
#' calling the qc function (simpleaffy).When not provided, it
#' is computed within the function.(datatype: QCStats)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image of percent present values, called "RawDataPercentPresentPlot"
#' @examples
#' # percPresPlot(rawData, quality=quality, experimentFactor, plotColors, legendColors)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

percPresPlot <- function(Data, quality = NULL, experimentFactor = NULL, plotColors = NULL, legendColors = NULL,
                         WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if (is.null(legendColors)) stop("the 'legendColors' parameter is required")
  if (is.null(quality)) try(quality <- qc(Data), TRUE)

  if (!is.null(quality)) {
    # Define the reference for the grey rectangle (it should minimize the number of outliers)
    ppMin <- min(quality@percent.present)
    ppMax <- max(quality@percent.present)
    ppMean <- mean(quality@percent.present)

    ref <- cbind(c(ppMean - 5, ppMin, ppMax - 10), c(ppMean + 5, ppMin + 10, ppMax))
    reftext <- c("[mean-5% ; mean+5%]", "[min ; min+10%]", "[max-10% ; max]")
    t1 <- (quality@percent.present >= (ppMean - 5) & quality@percent.present <= (ppMean + 5))
    t2 <- (quality@percent.present >= ppMin & quality@percent.present <= (ppMin + 10))
    t3 <- (quality@percent.present >= (ppMax - 10) & quality@percent.present <= ppMax)
    outliers <- c(length(t1[t1 == FALSE]), length(t2[t2 == FALSE]), length(t3[t3 == FALSE]))
    if (outliers[1] == 0) {
      ref <- ref[1, ]
      reftext <- reftext[1]
    } else {
      ref <- ref[outliers == min(outliers), ] # less outliers
      if (length(ref) > 3) {
        ref <- ref[1, ]
      }
      reftext <- reftext[outliers == min(outliers)] # less outliers
      if (length(reftext) > 1) {
        reftext <- reftext[1]
      }
    }
    testMin <- ref[1]
    testMax <- ref[2]

    png(filename = "RawDataPercentPresent.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    par(mfrow = c(1, 2), oma = c(17, 0.5, 0, 0.5), cex.axis = 0.8)
    plot(quality@percent.present,
      type = "n", ann = FALSE, axes = FALSE,
      frame.plot = TRUE, pch = ".", cex = 10, ylim = c(0, 100)
    )
    rect(0, testMin, length(quality@scale.factors) + 1,
      testMax,
      col = gray(0.9), border = FALSE
    )

    for (i in seq(from = 0, to = 100, by = 20)) {
      abline(h = i, lty = 2, col = gray(0.8))
    }
    abline(h = 50, lty = 1, col = gray(0.5))

    par(new = TRUE)
    plot(quality@percent.present,
      type = "h", ann = FALSE, axes = FALSE,
      frame.plot = TRUE, pch = ".", lwd = 3, col = plotColors, ylim = c(0, 100)
    )
    par(new = TRUE)
    plot(quality@percent.present,
      type = "p", ann = FALSE, axes = FALSE,
      frame.plot = TRUE, pch = ".", cex = 10, col = plotColors, ylim = c(0, 100)
    )

    t1 <- (quality@percent.present >= testMin & quality@percent.present <= testMax)
    if (length(t1[t1 == FALSE]) > 0 && length(sampleNames(Data)) >= (MAXARRAY / 2)) {
      if (length(t1[t1 == FALSE]) < 25) {
        textlab <- sampleNames(Data)[t1 == FALSE]
        legend("topleft", c("Outliers, from left to right", textlab),
          text.col = c(1, rep("red", length(textlab))), bty = "n", cex = 0.38
        )
      } else {
        mtext("Too many arrays and outliers; \nsee details on the QC table", side = 1, cex = 0.8, line = 3)
      }
    }

    title(main = "Plot of percent present", xlab = "", ylab = "percentage")
    axis(2)
    par(cex.axis = 0.65)
    if (length(sampleNames(Data)) < (MAXARRAY / 2)) { # array names not reported if more than 20 arrays
      axis(1, at = 1:length(quality@percent.present), las = 2, labels = sampleNames(Data))
    }

    if (length(levels(experimentFactor)) > 1) {
      legend("topright", levels(experimentFactor),
        col = legendColors, fill = legendColors, cex = 0.55, bty = "n"
      )
    }

    legend("bottomright", c(
      paste("min = ", round(ppMin, 2), sep = ""),
      paste("max = ", round(ppMax, 2), sep = ""),
      paste("max-min = ", round(ppMax - ppMin, 2), sep = "")
    ), bty = "o", cex = 0.55)

    mtext(paste(
      "Data should be in the grey rectangle representing a spread of 10%"
      # "Data should stand within the grey rectangle",reftext
    ),
    side = 4, font = 1, cex = 0.7
    )
    par(cex.axis = 0.8, cex.lab = 0.8)
    boxplot(quality@percent.present)
    title(main = "Boxplot of percent present")
    if (ppMax - ppMin <= 10) {
      title(xlab = "Percent present QC: OK (spread <= 10%)")
    } else {
      title(xlab = "Percent present QC: not OK (spread > 10%)")
    }

    dev.off()
  } else {
    warning("Percent present plot is not computed for this chip type")
  }
}


#################
## PNdistrPlot ##
#################

#' Create an image with boxplots of positive and negative control intensities
#'
#' This function (from functions_imagesQC.R) creates an image with
#' boxplots of positive and negative control intensities for each array.
#' It calls the borderQC1 function (\href{http://www.bioconductor.org/help/bioc-views/release/bioc/html/affyQCReport.html}{affyQCReport} Bioconductor package).
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @return A PNG image of boxplots of positive and negative control intensities, called "RawDataPosNegDistribPlot"
#' @examples
#' # PNdistrPlot(rawData)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

PNdistrPlot <- function(Data, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24) {
  png(filename = "RawDataPosNegDistribution.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
  par(oma = c(17, 0, 0, 0), srt = 90)
  borderQC1(Data)
  mtext(paste("All distributions should be similar and", "extreme values should not be reached\n"), side = 3, cex = 0.7)
  dev.off()
}


############
## bgPlot ##
############

#' Create an image of the background intensities and deviations
#'
#' This function (from functions_imagesQC.R) creates an image of
#' the background intensities and deviations for each array based
#' on the qc (\href{http://bioinformatics.picr.man.ac.uk/research/software/simpleaffy/index.html}{simpleaffy} Bioconductor package) function, which
#' generally only work for chiptypes with perfect match and mismatch probes,
#' but even not for all of those. As such, when these statistics
#' are not provided as parameters, trys are used in this function
#' to compute them internally. Values for which the try fails are
#' not computed, but the script will continue after giving a warning.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param quality (Status: optional, Default:NULL) object obtained by
#' calling the qc function (simpleaffy).When not provided, it
#' is computed within the function.(datatype: QCStats)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image of background intensities, called "RawDataBackgroundPlot"
#' @examples
#' # backgroundPlot(rawData, quality=quality, experimentFactor, plotColors,legendColors)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

backgroundPlot <- function(Data, quality = NULL, experimentFactor = NULL, plotColors = NULL, legendColors = NULL,
                           WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if (is.null(legendColors)) stop("the 'legendColors' parameter is required")
  if (is.null(quality)) try(quality <- qc(Data), TRUE)

  if (!is.null(quality)) {

    # Define the reference for the grey rectangle (it should minimize the number of outliers)
    bgMean <- mean(quality@average.background)
    bgMin <- min(quality@average.background)
    bgMax <- max(quality@average.background)
    ref <- cbind(c(bgMean - 10, bgMin, bgMax - 20), c(bgMean + 10, bgMin + 20, bgMax))
    reftext <- c("[mean-10 ; mean+10]", "[min ; min+20]", "[max-20 ; max]")
    t1 <- (quality@average.background >= (bgMean - 10) & quality@average.background <= (bgMean + 10))
    t2 <- (quality@average.background >= bgMin & quality@average.background <= (bgMin + 20))
    t3 <- (quality@average.background >= (bgMax - 20) & quality@average.background <= bgMax)
    outliers <- c(length(t1[t1 == FALSE]), length(t2[t2 == FALSE]), length(t3[t3 == FALSE]))
    if (outliers[1] == 0) {
      ref <- ref[1, ]
      reftext <- reftext[1]
    } else {
      ref <- ref[outliers == min(outliers), ] # less outliers
      if (length(ref) > 3) {
        ref <- ref[1, ]
      }
      reftext <- reftext[outliers == min(outliers)] # less outliers
      if (length(reftext) > 1) {
        reftext <- reftext[1]
      }
    }
    testMin <- ref[1]
    testMax <- ref[2]

    ylimit <- c(0, max(
      max(quality@maximum.background),
      testMax
    ) + 5)

    png(filename = "RawDataBackground.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    par(mfrow = c(1, 2), oma = c(17, 0.5, 0, 0.5), cex.axis = 0.8)
    plot(c(quality@minimum.background, quality@maximum.background),
      type = "n", ann = FALSE, axes = FALSE, frame.plot = TRUE, pch = ".", cex = 10,
      xlim = c(1, length(quality@maximum.background + 1)),
      ylim = ylimit
    )
    rect(0, testMin,
      length(quality@average.background) + 1,
      testMax,
      col = gray(0.8), border = TRUE, density = 10, lty = 6
    )
    par(new = TRUE)
    plot(quality@maximum.background,
      type = "h", ann = FALSE, axes = FALSE, frame.plot = TRUE, pch = ".", lwd = 3,
      col = plotColors, ylim = ylimit
    )
    par(new = TRUE)
    plot(quality@minimum.background,
      type = "h", ann = FALSE, axes = FALSE, frame.plot = TRUE, pch = ".", lwd = 4,
      col = "white", ylim = ylimit
    )
    par(new = TRUE)
    plot(quality@maximum.background,
      type = "p", ann = FALSE, axes = FALSE, frame.plot = TRUE, pch = 6, cex = 0.5,
      ylim = ylimit
    )
    par(new = TRUE)
    plot(quality@minimum.background,
      type = "p", ann = FALSE, axes = FALSE, frame.plot = TRUE, pch = 2, cex = 0.5,
      ylim = ylimit
    )
    par(new = TRUE)
    plot(quality@average.background,
      type = "p", ann = FALSE, axes = FALSE, frame.plot = TRUE, pch = 3, cex = 1,
      ylim = ylimit
    )
    abline(h = testMin, col = gray(0.8), ylim = ylimit)
    abline(h = testMax, col = gray(0.8), ylim = ylimit)
    title(main = "Plot of background intensity", xlab = "", ylab = "background intensity")
    axis(2)
    par(cex.axis = 0.65)
    if (length(sampleNames(Data)) < (MAXARRAY / 2)) { # array names not reported if more than 20 arrays
      axis(1, at = 1:length(quality@average.background), las = 2, labels = sampleNames(Data))
    }
    if (length(levels(experimentFactor)) > 1) {
      legend("bottomright", c(
        levels(experimentFactor), "max bg", "average bg",
        "min bg"
      ),
      col = c(legendColors, 1, 1, 1),
      pch = c(rep(15, length(legendColors)), 6, 3, 2), bty = "n", cex = 0.55
      )
    } else {
      legend("bottomright", c(
        "max bg", "average bg",
        "min bg"
      ), col = c(1, 1, 1), pch = c(6, 3, 2), bty = "n", cex = 0.55)
    }
    t1 <- (quality@average.background >= testMin & quality@average.background <= testMax)
    if (length(t1[t1 == FALSE]) > 0 && length(sampleNames(Data)) >= (MAXARRAY / 2)) {
      if (length(t1[t1 == FALSE]) < 25) {
        textlab <- sampleNames(Data)[t1 == FALSE]
        legend("bottomleft", c("Outliers, from left to right", textlab),
          text.col = c(1, rep("red", length(textlab))), bty = "n", cex = 0.38
        )
      } else {
        mtext("Too many arrays and outliers; \nsee details on the QC table", side = 1, cex = 0.8, line = 3)
      }
    }

    mtext(paste(
      "Data should be in the grey rectangle representing a spread of 20"
      # "Data should stand within the dashed grey rectangle",reftext
    ), side = 4, font = 1, cex = 0.7)
    par(cex.axis = 0.8, cex.lab = 0.8)
    boxplot(quality@average.background)
    title(main = "Average background intensity")

    if (bgMax - bgMin <= 20) {
      title(xlab = "Background QC: OK (spread <= 20)")
    } else {
      title(xlab = "Background QC: not OK (spread > 20)")
    }

    legend("bottomleft", c(
      paste("min = ", round(bgMin, 2), sep = ""),
      paste("max = ", round(bgMax, 2), sep = ""),
      paste("max-min = ", round(bgMax - bgMin, 2), sep = "")
    ), cex = 0.7)

    dev.off()
  } else {
    warning("Background intensity plot is not computed for this chip type")
  }
}


###################
## scaleFactPlot ##
###################

#' Create an image of the scale factors of the array
#'
#' This function (from functions_imagesQC.R) creates an image of
#' the scale factors of the array based on the qc
#' (\href{http://bioinformatics.picr.man.ac.uk/research/software/simpleaffy/index.html}{simpleaffy} Bioconductor package) function, which generally
#' only work for chiptypes with perfect match and mismatch probes,
#' but even not for all of those. As such, when these statistics
#' are not provided as parameters, trys are used in this function
#' to compute them internally. Values for which the try fails are
#' not computed, but the script will continue after giving a warning.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param quality (Status: optional, Default:NULL) object obtained by
#' calling the qc function (simpleaffy).When not provided, it
#' is computed within the function.(datatype: QCStats)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image of scale factors, called "RawDataScaleFactorsPlot"
#' @examples
#' # scaleFactPlot(rawData, quality=quality, experimentFactor, plotColors,legendColors)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

scaleFactPlot <- function(Data, quality = NULL, experimentFactor = NULL, plotColors = NULL, legendColors = NULL,
                          WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if (is.null(legendColors)) stop("the 'legendColors' parameter is required")
  if (is.null(quality)) try(quality <- qc(Data), TRUE)

  if (!is.null(quality)) {

    # Define the reference for the grey rectangle (it should minimize the number of outliers)
    sfMin <- min(log2(quality@scale.factors))
    sfMax <- max(log2(quality@scale.factors))
    sfMean <- mean(log2(quality@scale.factors))
    ref <- cbind(c(sfMean - 1.5, sfMin, sfMax - 3), c(sfMean + 1.5, sfMin + 3, sfMax))
    t1 <- (log2(quality@scale.factors) >= (sfMean - 1.5) & log2(quality@scale.factors) <= (sfMean + 1.5))
    t2 <- (log2(quality@scale.factors) >= sfMin & log2(quality@scale.factors) <= (sfMin + 3))
    t3 <- (log2(quality@scale.factors) >= (sfMax - 3) & log2(quality@scale.factors) <= sfMax)
    outliers <- c(length(t1[t1 == FALSE]), length(t2[t2 == FALSE]), length(t3[t3 == FALSE]))
    if (outliers[1] == 0) {
      ref <- ref[1, ]
    } else {
      ref <- ref[outliers == min(outliers), ] # less outliers
      if (length(ref) > 3) {
        ref <- ref[1, ]
      }
    }
    testMin <- ref[1]
    testMax <- ref[2]

    ymin <- min(-3.5, sfMin - 1.5, testMin)
    ymax <- max(+3.5, sfMax + 1.5, testMax)
    ylimit <- c(ymin, ymax)

    png(filename = "RawDataScaleFactors.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)

    par(mfrow = c(1, 2), oma = c(17, 0.5, 0, 0.5), cex.axis = 0.8)
    plot(log2(quality@scale.factors),
      type = "n", axes = FALSE, frame.plot = TRUE,
      ann = FALSE, pch = ".", cex = 10, ylim = ylimit
    )
    rect(0, testMin,
      length(quality@scale.factors) + 1,
      testMax,
      col = gray(0.9), border = FALSE
    )

    for (i in seq(from = ceiling(ymin), to = floor(ymax), by = 1)) {
      abline(h = i, lty = 2, col = gray(0.8))
    }
    abline(h = 0, lty = 1, col = gray(0.5))

    par(new = TRUE)
    plot(log2(quality@scale.factors),
      type = "h", axes = FALSE, frame.plot = TRUE, ann = FALSE,
      pch = ".", lwd = 3, col = plotColors, ylim = ylimit
    )
    par(new = TRUE)
    plot(log2(quality@scale.factors),
      type = "p", axes = FALSE, frame.plot = TRUE, ann = FALSE, pch = ".",
      cex = 10, col = plotColors, ylim = ylimit
    )
    t1 <- (log2(quality@scale.factors) >= testMin & log2(quality@scale.factors) <= testMax)
    if (length(t1[t1 == FALSE]) > 0 && length(sampleNames(Data)) >= (MAXARRAY / 2)) {
      if (length(t1[t1 == FALSE]) < 16) {
        textlab <- sampleNames(Data)[t1 == FALSE]
        legend("bottomleft", c("Outliers, from left to right", textlab),
          text.col = c(1, rep("red", length(textlab))), bty = "n", cex = 0.38
        )
      } else {
        mtext("Too many arrays and outliers; \nsee details on the QC table", side = 1, cex = 0.8, line = 3)
      }
    }

    title(main = "Plot of Log scale factors", xlab = "", ylab = "Log2(scale factors)")
    axis(2)
    par(cex.axis = 0.65)
    if (length(sampleNames(Data)) < (MAXARRAY / 2)) { # array names not reported if more than 20 arrays
      axis(1, at = 1:length(quality@scale.factors), las = 2, labels = sampleNames(Data))
    }
    if (length(levels(experimentFactor)) > 1) {
      legend("topright", levels(experimentFactor),
        col = legendColors, fill = legendColors, cex = 0.55, bty = "n"
      )
    }
    legend("topleft", c(
      paste("min = ", round(sfMin, 2), sep = ""),
      paste("max = ", round(sfMax, 2), sep = ""),
      paste("max - min = ", round((sfMax - sfMin), 2), sep = "")
    ),
    bty = "n",
    cex = 0.55
    )

    mtext(paste(
      "Data should be in the grey rectangle representing 3-fold on a log scale"
    ), side = 4, font = 1, cex = 0.7)
    par(cex.axis = 0.8)
    boxplot(quality@scale.factors)
    title(main = "Boxplot of scale factors")
    mtext("(natural scale)\n", side = 3, font = 1, cex = 0.7)
    par(cex = 0.8)
    if ((sfMax - sfMin) < 3) {
      title(xlab = "Scale factors QC: OK (spread < 3-fold)")
    } else {
      title(xlab = "Scale factors QC: not OK (spread > 3-fold)")
    }
    dev.off()
  } else {
    warning("Scale factor plot is not computed for this chip type")
  }
}


##################
## controlPlots ##
##################

#' Create an image of the scale factors of the array
#'
#' This function (from functions_imagesQC.R) creates images of
#' the expression values of the affx controls, if present. One
#' image shows the expression profiles over all samples, for each
#' control separately. The other image shows the boxplots of all
#' controls (together), for each sample.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param affxplots (Status: optional, Default:TRUE) Does the AFFX expression profiles have to be plotted? (datatype: logical)
#' @param boxplots (Status: optional, Default:TRUE) Does the AFFX and other controls boxplots have to be plotted? (datatype: logical)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return Two images. A PNG image of the expression profiles of the affx controls, called "Profiles_affx_controls"
#' A PNG image of the boxplots of the affx controls, called "Boxplots_affx_controls"
#' @examples
#' # controlPlots(rawData,plotColors)
#' @import ArrayTools grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

controlPlots <-
  function(Data,
           plotColors = NULL,
           experimentFactor = NULL,
           legendColors = NULL,
           affxplots = TRUE,
           boxplots = TRUE,
           WIDTH = 1000,
           HEIGHT = 1414,
           POINTSIZE = 24,
           MAXARRAY = 41) {
    if (is.null(plotColors)) {
      stop("The 'plotColors' parameter must be provided")
    }
    if (is.null(experimentFactor)) {
      stop("the 'experimentFactor' parameter is required")
    }
    if (is.null(legendColors)) {
      stop("the 'legendColors' parameter is required")
    }

    affx1 <-
      unique(probeNames(Data)[grep("AFFX", probeNames(Data), fixed = TRUE)])

    if (requireNamespace("ArrayTools", quietly = TRUE)) {
      dataTable <-
        paste(substr(Data@annotation, 1, nchar(Data@annotation) - 2), "CONTROL",
          sep =
            ""
        )
      suppressWarnings(eval(parse(
        "", -1, paste("data(", dataTable, ")", sep = "")
      ))) # ArrayTools
      cntrl <- NULL
      try(cntrl <- get(dataTable), TRUE)

      cntrlVars <- NULL
      if (!is.null(cntrl)) {
        for (name in names(table(cntrl[, 2])[table(cntrl[, 2]) > 0])) {
          # often names have a format of xxx->yyy, we only want the last part
          name2 <-
            strsplit(name, "->")[[1]][length(strsplit(name, "->")[[1]])]
          # make sure that there are only letters in name2, to be used as variable name
          name2 <- gsub("[[:punct:]|[:digit:]]", "", name2)
          assign(name2, cntrl[grep(name, cntrl[, 2]), 1])
          assign(name2, get(name2)[get(name2) %in% geneNames(Data)])
          if (length(get(name2)) > 0) {
            cntrlVars <- c(cntrlVars, name2)
          }
        }
      }

      if (affxplots) {
        # plot affx controls

        if (length(affx1) > 0 | length(grep("affx", cntrlVars)) > 0) {
          if (length(grep("affx", cntrlVars)) > 0) {
            # affx2 cannot exists as automatically generated variable names are restricted to letters
            affx2 <- NULL
            for (a in grep("affx", cntrlVars)) {
              affx2 <- c(affx2, get(cntrlVars[a]))
            }
            if (length(affx1) > 0) {
              # both exist
              if (sort(affx2) != sort(affx1)) {
                # merge both
                affx2 <- unique(c(affx1, affx2))
              }
            }
          } else {
            # only length(affx1)>0
            affx2 <- affx1
          }
          # 		dimnCoef <- 1 # squared layout
          dimnCoef <- 2
          dimn <- ceiling(sqrt(length(affx2) / dimnCoef))
          number_empty <- dimnCoef * dimn * dimn - length(affx2)

          png(
            filename = "RawDataAFFXControlsProfiles.png",
            width = 250 * dimn + 250,
            height = dimnCoef * 250 * dimn,
            pointsize = POINTSIZE
          )
          if (number_empty >= dimn) {
            lineout <- dimn
            number_empty <- number_empty - dimn
          } else {
            lineout <- 0
          }
          legcol <- dimnCoef * dimn * dimn - lineout + 1
          layout(cbind(matrix(
            1:(dimnCoef * dimn * dimn - lineout),
            ncol = dimn,
            byrow = TRUE
          ), legcol))
          par(oma = c(5, 3, 4, 3))

          dummy <- sapply(
            affx2,
            function(y) {
              pm <- log(pm(Data, as.character(y)), 2)
              max_pm <- max(pm)
              min_pm <- min(pm)
              plot(
                0,
                type = "n",
                ylim = c(min_pm, max_pm),
                ylab = "expression",
                xlim = c(1, dim(pm)[1]),
                main = y,
                xlab = paste("\naverage: ", round(mean(pm))),
                cex.lab = 1.3,
                cex.main = 1.1
              )
              assign("col", 0)
              apply(
                pm,
                2,
                function(z) {
                  assign("col", get("col", 1) + 1)
                  points(z, type = "o", col = plotColors[get("col", 1)])
                }
              )
              points(
                apply(pm, 1, mean),
                type = "l",
                lwd = 2,
                col = "black"
              )
            }
          )


          for (p in 1:number_empty) {
            plot(
              0,
              type = "n",
              xaxt = "n",
              yaxt = "n",
              bty = "n",
              xlab = "",
              ylab = "",
              main = ""
            )
          }
          par(mar = c(0, 0, 0, 0))
          plot(
            0,
            type = "n",
            xaxt = "n",
            yaxt = "n",
            bty = "n",
            xlab = "",
            ylab = "",
            main = ""
          )
          legend(
            "topleft",
            c(sampleNames(Data), "", "average"),
            lwd = c(rep(1, length(
              sampleNames(Data)
            )), 0, 2),
            pch = c(rep(1, length(
              sampleNames(Data)
            )), -1, -1),
            col = c(plotColors, "white", "black"),
            cex = 1.5,
            bty = "n"
          )

          mtext(
            "affx control profiles",
            side = 3,
            outer = TRUE,
            font = 2,
            cex = 2
          )

          dev.off()
        } else {
          warning("plot of AFFX controls cannot be made, AFFX controls cannot be determined")
        }
      }
      ###################
      if (boxplots) {
        if (length(affx1) > 0) {
          # affx1 will not already be in controlVars as it contains a number
          cntrlVars <- c("affx1", cntrlVars)
        }
        n <- 0
        for (var in cntrlVars) {
          varName <- paste(var, "Data", sep = "")
          assign(varName, NULL)
          try(assign(varName, log(pm(Data, as.character(get(
            var
          ))), 2)), TRUE)
          if (!is.null(get(varName))) {
            n <- n + 1
          }
        }
        # make boxplots of controls
        if (n > 0) {
          for (i in 1:n) {
            if (!is.null(get(varName))) {
              png(
                paste(
                  "RawData",
                  toupper(cntrlVars[i]),
                  "ControlsBoxplot.png",
                  sep = ""
                ),
                width = WIDTH,
                height = HEIGHT,
                pointsize = POINTSIZE
              )
              par(oma = c(17, 0, 0, 0))
              boxplot(
                get(paste(cntrlVars, "Data", sep = "")[i]),
                main = paste(gsub("1", "",
                  cntrlVars[i],
                  fixed = TRUE
                ), "controls"),
                ylim = c(0, 16),
                col = plotColors,
                cex = 1,
                axes = FALSE
              )
              if (length(sampleNames(Data)) < MAXARRAY) {
                cexval <- 0.65
              } else {
                cexval <- 0.45
              }
              if (length(levels(experimentFactor)) > 1) {
                legend(
                  "topright",
                  levels(experimentFactor),
                  col = legendColors,
                  fill = legendColors,
                  cex = 0.7,
                  bg = "white",
                  bty = "o"
                )
              }
              axis(
                1,
                at = 1:length(sampleNames(Data)),
                las = 2,
                labels = sampleNames(Data),
                cex.axis = cexval
              )
              axis(2, cex.axis = 0.7)
              dev.off()
            }
          }
        } else {
          warning("no control boxplots can be made, control spots cannot be determined")
        }
      }
    }else{
      stop("ArrayTools package is not availabel", call. = FALSE)
    } # if requireNamespace()
  }


################
## boxplotFun ##
################

#' Create an image with intensity boxplots for all arrays
#'
#' This function (from functions_imagesQC.R) creates an image
#' with intensity boxplots for all arrays in the raw or normalized
#' dataset (depending on the object passed). It calls the boxplot
#' function applied to an AffyBatch object.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param normMeth (Status: required when Data is a normalized data object, Default:"")
#' String indicating the normalization method used (see
#' \link[ArrayAnalysis]{normalizeData} function for more information on the possible values). (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image of the raw or normalized boxplots of the arrays, called ‘DataBoxplot’
#' @examples
#' # By default, before the normalization the script will call:
#' # boxplotFun(Data=rawData, experimentFactor, plotColors, legendColors)
#' # and after normalization:
#' # boxplotFun(Data=normData, experimentFactor, plotColors, legendColors, normMeth=normMeth)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

boxplotFun <- function(Data, experimentFactor = NULL, plotColors = NULL, legendColors = NULL, normMeth = "",
                       WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if (is.null(legendColors)) stop("the 'legendColors' parameter is required")


  print(class(Data))
  #check affy or genefeatured (oligo)
  if(class(Data) == "AffyBatch" || class(Data) == "GeneFeatureSet"  || class(Data) == "ExonFeatureSet") {

    Type <- "Raw"
    tmain <- "Boxplot of raw intensities"
    tmtext2 <- "Raw log intensity\n\n\n"
  } else {
    if (normMeth == "") stop("When providing a normalised data object, the normMeth parameter is required")
    Type <- "Norm"
    tmain <- paste("Boxplot after ", normMeth, sep = "")
    tmtext2 <- "Normalized log intensity\n\n\n"
  }

  png(filename = paste(Type, "DataBoxplot.png", sep = ""), width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
  par(oma = c(17, 0, 0, 0), cex.axis = 1)

  if(class(Data) != "GeneFeatureSet"  && class(Data) != "ExonFeatureSet") {
    suppressWarnings(boxplot(Data, col=plotColors ,main=tmain, axes=FALSE))
  } else {
    suppressWarnings(boxplot(Data, which="all", col=plotColors ,main=tmain, axes=FALSE))
  }

  if (length(levels(experimentFactor)) > 1) {
    legend("topright", levels(experimentFactor),
      col = legendColors, fill = legendColors, cex = 0.7, bg = "white", bty = "o"
    )
  }
  if (length(sampleNames(Data)) < MAXARRAY) {
    cexval <- 0.65
  } else {
    cexval <- 0.45
  }
  axis(1, at = 1:length(sampleNames(Data)), las = 2, labels = sampleNames(Data), cex.axis = cexval)
  axis(2, cex.axis = 0.7)
  mtext(tmtext2, side = 2, cex = 0.8)
  mtext("Distributions should be comparable between arrays\n",
    side = 3, font = 1,
    cex = 0.7
  )
  dev.off()
}


################
## densityFun ##
################

#' Create an image with a density curve of the intensities
#'
#' This function (from functions_imagesQC.R) creates an image with
#' a density curve of the intensities for all arrays in the raw or
#' normalized dataset (depending on the object passed).
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param normMeth (Status: required when Data is a normalized data object, Default:"")
#' String indicating the normalization method used (see
#' \link[ArrayAnalysis]{normalizeData} function for more information on the possible values). (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image of the raw or normalized density plots of the arrays, called "DensityHistogram"
#' @examples
#' # By default, before the normalization the script will call:
#' # densityFun(Data=rawData, plotColors)
#' # and after normalization:
#' # densityFun(Data=normData, plotColors, normMeth=normMeth)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

densityFun <- function(Data, plotColors = NULL, normMeth = "",
                       WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")

  #check affy or genefeatured (oligo)
  Type <- ifelse(class(Data) == "AffyBatch" ||  class(Data) == "GeneFeatureSet" ||  class(Data) == "ExonFeatureSet" ,"Raw","Norm")

  png(filename = paste(Type,"DensityHistogram.png", sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)

  if (length(sampleNames(Data)) < MAXARRAY) {
    cexval <- 0.65
    par(oma = c(12, 0, 0, 0))
  } else {
    cexval <- 0.45
    par(oma = c(0.1, 0, 0, 0))
  }
  if (Type == "Raw") {
    hist(Data,
      lwd = 3, lt = 1:length(sampleNames(Data)), col = plotColors,
      which = "both", main = "Density histogram of raw intensities", cex.axis = 0.7, cex.lab = 0.8
    )
  } else {
    if (normMeth == "") stop("When providing a normalised data object, the normMeth parameter is required")
    hist(Data,
      lwd = 3, lt = 1:length(sampleNames(Data)), col = plotColors,
      main = paste("Density histogram after ", normMeth, "\n", sep = ""), cex.axis = 0.7, cex.lab = 0.8
    )
  }

  legend("topright", substr(sampleNames(Data), 1, 20),
    lwd = 3, lty = 1:length(sampleNames(Data)),
    col = plotColors, cex = cexval, bty = "n"
  )
  mtext("Curves should be comparable between arrays\n",
    side = 3, font = 1,
    cex = 0.7
  )
  dev.off()
}


##########################
## densityFunUnsmoothed ##
##########################

# alternative density function, not called by default by run_affyAnalysisQC.R
#----------------------------------------------------------------------------

#' Create an image with an unsmoothed density curve of the intensities
#'
#' This function (from functions_imagesQC.R) creates an image with
#' an unsmoothed density curve of the intensities for all arrays in
#' the raw or normalized dataset (depending on the object passed).
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param normMeth (Status: required when Data is a normalized data object, Default:"")
#' String indicating the normalization method used (see
#' \link[ArrayAnalysis]{normalizeData} function for more information on the possible values). (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image of the raw or normalized density plots of the arrays, called "DensityHistogram"
#' @examples
#' # By default, before the normalization the script will call:
#' # densityFunUnsmoothed(Data=rawData, plotColors)
#' # and after normalization:
#' # densityFunUnsmoothed(Data=normData, plotColors, normMeth=normMeth)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

densityFunUnsmoothed <- function(Data, plotColors = NULL, normMeth = "",
                                 WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")

  # plot(density(exprs(Data)))
  #check affy or genefeatured (oligo)
  Type <- ifelse(class(Data) == "AffyBatch" ||  class(Data) == "GeneFeatureSet" ,"Raw","Norm")

  if (Type == "Raw") {
    tmp_data <- log(exprs(Data), 2)
    main <- "Density histogram of raw intensities"
  } else { # normalized data
    if (normMeth == "") stop("When providing a normalised data object, the normMeth parameter is required")
    tmp_data <- exprs(Data)
    main <- paste("Density histogram after", normMeth)
  }

  png(filename = paste(Type, "DensityHistogramUnsmoothed.png", sep = ""), width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
  if (length(sampleNames(Data)) < MAXARRAY) {
    cexval <- 0.65
    par(oma = c(12, 0, 0, 0))
  } else {
    cexval <- 0.45
    par(oma = c(0.1, 0, 0, 0))
  }
  for (i in 1:length(sampleNames(Data))) {
    if (i > 1) par(new = TRUE)
    tmp_hist <- hist(tmp_data[, i], plot = FALSE, breaks = 64)
    plot(tmp_hist$mids, tmp_hist$density,
      type = "l", lty = i, col = plotColors[i], xlim = c(min(tmp_data), max(tmp_data)), ylim = c(0, 1),
      main = main, xlab = "log2 intensity", ylab = "density", cex.axis = 0.7, cex.lab = 0.8
    )
    text(tmp_hist$mids[order(tmp_hist$density)[length(tmp_hist$density)]], max(tmp_hist$density), sampleNames(Data)[i], col = i, pos = 3, offset = 0.2, cex = 0.6)
  }

  legend("topright", substr(sampleNames(Data), 1, 20),
    lwd = 3, lty = 1:length(sampleNames(Data)),
    col = plotColors, cex = cexval, bty = "n"
  )
  mtext("Curves should be comparable between arrays\n",
    side = 3, font = 1,
    cex = 0.7
  )

  dev.off()
}


###########
## maFun ##
###########

#' Create MA plots for each array versus the median array for the raw or normalized dataset
#'
#' This function (from functions_imagesQC.R) creates MA plots for each
#' array versus the median array for the raw or normalized dataset.
#' The median array is computed for the whole data set (if perGroup is FALSE)
#' of per experimental group (perGroup is TRUE). In the script this setting
#' will depend on the setting op the MAOption1 parameter, which can have
#' the values “dataset” or “group”.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch or ExpressionSet)
#' @param experimentFactor (Status: required when perGroup is TRUE, Default:NULL) The factor of groups. (datatype: factor)
#' @param perGroup (Status: optional, Default:FALSE) Are MA plots to be made for each experimental group separately or not? (datatype: logical)
#' @param normMeth (Status: required when Data is a normalized data object, Default:"")
#' String indicating the normalization method used (see
#' \link[ArrayAnalysis]{normalizeData} function for more information on the possible values). (datatype: character)
#' @param aType (Status: optional, Default:NULL) String indicating the type of
#' the current chip, either “PMMM” for chips with perfect match and
#' mismatch probes, or “PMonly” for chips with perfect match probes only.
#' Required when Data is a raw data object. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return Images of the MA plots of each array versus the median array, each file
#' contains MA plots for six arrays. The files contain the string ‘MAplot” and a
#' number if more than one are needed; in case of groupwise computation, the name
#' of the group is also included in the filename.
#' @examples
#' # By default, before the normalization the script will call:
#' # maFun(Data=rawData, experimentFactor, perGroup=(MAOption1=="group"), aType=aType)
#' # and after normalization:
#' # maFun(Data=normData, experimentFactor, perGroup=(MAOption1=="group"), normMeth=normMeth)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

maFun <- function(Data, experimentFactor = NULL, perGroup = FALSE, normMeth = "", aType = NULL,
                  WIDTH = 1000, HEIGHT = 1414, MAXARRAY = 41) {

  # verify whether raw or normalised data have been provided
  #check affy or genefeatured (oligo)
  if(class(Data) == "AffyBatch" ||  class(Data) == "GeneFeatureSet" ||  class(Data) == "ExonFeatureSet"){
    if (is.null(aType)) stop("When selecting MA plots of raw data, 'aType' must be provided")
    Type <- "Raw"
    tmain <- paste("MA plots of raw data", ifelse(perGroup, ", computed for group", ""), sep = "")
  } else {
    if (normMeth == "") stop("When selecting MA plots of normalized data, 'normMeth' must be provided")
    Type <- "Norm"
    tmain <- paste("MA plots after", normMeth, "normalization", ifelse(perGroup, ", computed for group ", ""), sep = "")
  }

  # check whether MA plots have to be created per experimental group or for the whole data set at once
  if (perGroup) {
    if (is.null(experimentFactor)) stop("When selecting MA plots per experimental group, 'experimentFactor' must be provided")
  } else {
    experimentFactor <- factor(rep("", length(sampleNames(Data))))
  }

  # plot the images
  for (k in (levels(experimentFactor))) {
    x2 <- Data[, experimentFactor == k]

    # amount of plots per page
    if (length(sampleNames(x2)) <= 12) {
      nPerPage <- 6
      nCol <- 2
      cexmain <- 1.9
      cexval <- 1.8
    } # max 3p of 6 plots
    if (length(sampleNames(x2)) > 12) {
      nPerPage <- 15
      nCol <- 3
      cexmain <- 1.3
      cexval <- 1.3
    }
    if (length(sampleNames(x2)) > 45) {
      nPerPage <- 24
      nCol <- 4
      cexmain <- 1.4
      cexval <- 1.2
    }

    # subselection of data
    if (length(sampleNames(x2)) < 2) {
      warning(
        "group ", ifelse(perGroup, k, ""), " consists of only one array: no MA plot can be made",
        ifelse(perGroup, ": consider selecting 'dataset' as option for the MA plots", "")
      )
    } else {
      nplots <- ceiling(length(sampleNames(x2)) / nPerPage)
      for (l in 1:nplots) {
        png(paste(Type, "DataMAplot", ifelse(k != "", "_", ""), gsub("[:/|\\?<>*\"]", "_", k),
          ifelse(nplots > 1, "_", ""), ifelse(nplots > 1, l, ""), ".png",
          sep = ""
        ), width = WIDTH, height = HEIGHT)
        par(mfrow = c((nPerPage / nCol), nCol), oma = c(0, 0, 6, 0), adj = 0)
        from <- nPerPage * l - (nPerPage - 1)
        to <- min(nPerPage * l, length(sampleNames(x2)))
        if (Type == "Raw" && aType != "PMMM") {
          if(class(Data) == "AffyBatch"){
            MAplot(x2, pairs = FALSE, which= from:to, plotFun = smoothScatter, lwd=3, type="pm", cex.main=cexmain,cex=cexval)
          } else {
            MAplot(x2, pairs = FALSE, which= from:to, plotFun = smoothScatter, lwd=3, what=exprs, cex.main=cexmain,cex=cexval)
          }
        } else {
          MAplot(x2, pairs = FALSE, which= from:to, plotFun = smoothScatter, lwd=3, cex.main=cexmain,cex=cexval)
        }
        mtext(paste(tmain, ifelse(perGroup, k, ""), ifelse(nplots > 1, paste(l, "/", nplots), "")), side = 3, outer = TRUE, font = 2, cex = 2)
        dev.off()
      }
    }
  }
}


#####################
## plotArrayLayout ##
#####################

#' Create an image of the layout of the current chiptype
#'
#' This function (from functions_imagesQC.R) creates an image of the layout
#' of the current chiptype. Thus, this plot does not plot any data, but shows
#' where control and regular perfect match (PM) and (if applicable)
#' mismatch (MM) probes are present on the array. Note: due to resolution issues,
#' banding may seem different from the real situation,
#' e.g. normally on classical chiptypes, PM and MM are present
#' in alternate lines, but patterns may appear due to image resolution.
#' This function tries to load annotation libraries, depending on the chiptype.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch or ExpressionSet)
#' @param aType (Status: optional, Default:NULL) String indicating the type of
#' the current chip, either “PMMM” for chips with perfect match and
#' mismatch probes, or “PMonly” for chips with perfect match probes only.
#' Required when Data is a raw data object. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @return A PNG image of the expression profiles of the affx controls, called "Array_layout_plot"
#' @examples
#' # plotArrayLayout(rawData,aType)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

plotArrayLayout <- function(Data,
                            aType = NULL,
                            WIDTH = 1000,
                            HEIGHT = 1414,
                            POINTSIZE = 24) {
  # note that this plot is based on the original Affymetrix cdf file
  # with updated annotations, some originale control probesets could be regular ones, many regular ones will be discarded

  # stop if array type is not provided
  if (is.null(aType)) {
    stop("The aType parameters must be provided")
  }

  # copy first array as a model
  data.tmp <- Data[, 1]

  # set every probe to intensity 1
  exprs(data.tmp)[] <- 1

  # determine which are the perfect match rows in exprs
  rn_pm <- sort(as.numeric(rownames(pm(data.tmp))))
  rn_mm <- NULL
  if (aType == "PMMM") {
    # determine which are the mismatch rows in exprs if present
    rn_mm <- sort(as.numeric(rownames(mm(data.tmp))))
    # set intensity of all mismatch probes if present to 2
    exprs(data.tmp)[rn_mm, ] <- 2
  }

  # determine control probes
  affx <- grep("AFFX", probeNames(data.tmp), fixed = TRUE)

  if (requireNamespace("ArrayTools", quietly = TRUE)) {
    library <-
      paste(substr(Data@annotation, 1, nchar(Data@annotation) - 2),
        "transcriptcluster.db",
        sep = ""
      )
    suppressWarnings(eval(parse(
      "",
      -1,
      paste("require(", library, ", quietly = TRUE)", sep = "")
    )))
    # get probe annotations
    all <- NULL
    eval(parse(
      "",
      -1,
      paste(
        "try(all2<-",
        substr(library, 1, nchar(library) - 3),
        "ACCNUM,TRUE)",
        sep = ""
      )
    ))
    if (exists("all2")) {
      eval(parse("", -1, paste(
        "try(all<-toTable(all2),TRUE)",
        sep = ""
      )))
    }

    dataTable <-
      paste(substr(Data@annotation, 1, nchar(Data@annotation) - 2), "CONTROL",
        sep =
          ""
      )
    suppressWarnings(eval(parse(
      "", -1, paste("data(", dataTable, ")", sep = "")
    ))) # ArrayTools
    cntrl <- NULL
    try(cntrl <- get(dataTable), TRUE)

    # controls<-all[all[,2]=="---",1]
    # sum(controls!=sort(cntrl[,1])) #should be 0, so all IDs are equal

    lab_add <- ""
    if (length(affx) > 0) {
      # set intensity of control probes to 5 for matches and to 8 for mismatches if present
      pm(data.tmp)[affx, ] <- 5
      if (aType == "PMMM") {
        mm(data.tmp)[affx, ] <- 8
      }
    }

    # if no cntrl object given: warn and skip plot
    if (!is.null(cntrl)) {
      # set control probes to intensity 2 (affx, intron, exon)
      # these normally are the affx, intron, and exon probes as ag controls are not present in pm
      pm(data.tmp)[probeNames(data.tmp) %in% cntrl[, 1], ] <- 5
      if (aType == "PMMM") {
        mm(data.tmp)[probeNames(data.tmp) %in% cntrl[, 1], ] <- 8
      }
    }

    # if given, set probes missing from list of annotated probes to NA
    if (!is.null(all)) {
      missingAnnotation <- !(probeNames(data.tmp) %in% all[, 1])
      if (sum(missingAnnotation) > 0) {
        pm(data.tmp)[missingAnnotation, ] <- NA
        if (aType == "PMMM") {
          mm(data.tmp)[missingAnnotation, ] <- NA
        }
        lab_add <- " - white: other unannotated probe"
      }
    }

    # set intensity of all probes in exprs, but not in pm (or if present mm) to 15 (all numbers chosen for the colors to work)
    exprs(data.tmp)[setdiff(1:dim(exprs(data.tmp))[1], c(rn_pm, rn_mm)), ] <-
      15

    # determine colors needed
    colors <-
      c("black", "#222222", "blue", "lightblue", "red")[c(TRUE, aType == "PMMM", TRUE, aType ==
        "PMMM", TRUE)]
    # set up label
    if (aType == "PMMM") {
      xlab <-
        paste(
          "black/gray: regular match/mismatch probe \n blue/light blue: ",
          "control match/mismatch probe \n red: unannotated probe (control region)\n",
          lab_add,
          sep = ""
        )
    } else {
      xlab <-
        paste(
          "black: regular probe \n blue: control probe \n red: unannotated probe (control region)\n",
          lab_add,
          sep = ""
        )
    }

    # plot the image
    png(
      "RawDataReferenceArrayLayout.png",
      width = WIDTH,
      height = HEIGHT,
      pointsize = POINTSIZE
    )
    par(oma = c(17, 2, 0, 3))
    image(
      data.tmp,
      col = colors,
      xlab = xlab,
      main = "Array reference layout",
      cex.lab = 0.8
    )
    dev.off()
  }else{
    stop("ArrayTools package is not availabel", call. = FALSE)
  } # if requireNamespace()
}


###################
## spatialImages ##
###################

#' Create an image per array, containing one to four spatial images
#'
#' This function (from functions_imagesQC.R) creates an image per array,
#' containing one to four spatial images. For any but the raw images,
#' an object obtained by calling fitPLM (affyPLM Bioconductor package)
#' is used. If this object is not provided, it will be computed internally
#' within this function.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch or ExpressionSet)
#' @param Data.pset (Status: optional, Default:NULL) An object obtained by
#' calling fitPLM (affyPLM), used for each but the raw plot. When not provided,
#' it is computed within the function.(datatype: PLMset)
#' @param Resid (Status: optional, Default:TRUE) Should a residual plot be made? (datatype: logical)
#' @param ResSign (Status: optional, Default:TRUE) Should a residual sign plot be made? (datatype: logical)
#' @param Raw (Status: optional, Default:TRUE) Should a raw plot be made? (datatype: logical)
#' @param Weight (Status: optional, Default:TRUE) Should a weight plot be made? (datatype: logical)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @return Images containing each of the requested plots, one file per array.
#' Naming is ‘virtual_image’ followed by the (sample)name of the array.
#' @examples
#' # There are two default calls by the script:
#' # 1/ Compute only the images showing the residuals of the PLM (if spatialImage parameter is TRUE):
#' # spatialImages(rawData, Data.pset=rawData.pset, Resid=TRUE, ResSign=FALSE, Raw=FALSE, Weight=FALSE)
#' # 2/ Compute the four images for all arrays (if PLMimage parameter is TRUE):
#' # spatialImages(rawData, Data.pset=rawData.pset)
#' # where rawData.pset has been constructed by calling:
#' # rawData.pset <- fitPLM(rawData)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

spatialImages <- function(Data, Data.pset = NULL, Resid = TRUE, ResSign = TRUE, Raw = TRUE, Weight = TRUE,
                          WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24) {
  if (is.null("Data.pset") & (Resid || ResSign || Weight)) {
    Data.pset <- fitPLM(Data)
  }

  splotFun <- function(typeVal) {
    if (length(sampleNames(Data)) <= 18) {
      nPerPage <- 6
      nCol <- 2
      cexval <- 1
    } # max 3p of 6 plots
    if (length(sampleNames(Data)) > 18) {
      nPerPage <- 12
      nCol <- 3
      cexval <- 0.8
    } # max 4p of 12 plots
    if (length(sampleNames(Data)) > 48) {
      nPerPage <- 20
      nCol <- 4
      cexval <- 0.8
    }

    nPages <- ceiling(length(sampleNames(Data)) / (nPerPage + 1))

    for (l in 1:nPages) {
      from <- nPerPage * l - (nPerPage - 1)
      to <- min(nPerPage * l, length(sampleNames(Data)))
      png(
        filename = paste("RawData2DVirtualImage_", typeVal, ifelse(nPages > 1, paste("_", l, sep = ""), ""), ".png", sep = ""),
        width = WIDTH, height = HEIGHT, pointsize = POINTSIZE
      )

      layout(matrix(1:nPerPage, ncol = nCol, byrow = TRUE))
      par(oma = c(0, 0, 3, 0), mar = c(1, 1, 3, 1))

      if (typeVal != "raw") {
        for (k in from:to) {
          image(Data.pset, which = k, type = typeVal, cex.main = cexval)
        }
        mtext(paste(
          "2D virtual PLM image for model characteristic:", typeVal,
          ifelse(nPages > 1, paste(l, "/", nPages), "")
        ),
        side = 3, outer = TRUE,
        font = 2, cex = 1.1
        )
      } else {
        for (k in from:to) {
          image(Data[, k], cex.main = cexval)
        }
        mtext(paste(
          "2D virtual image of raw data intensities",
          ifelse(nPages > 1, paste(l, "/", nPages), "")
        ),
        side = 3, outer = TRUE,
        font = 2, cex = 1.1
        )
      }

      dev.off()
    }
  }

  if (Resid) splotFun("resids")
  if (ResSign) splotFun("sign.resids")
  if (Raw) splotFun("raw")
  if (Weight) splotFun("weights")
}

#################
## array.image ##
#################

# alternative function to plot virtual arrays when PLM cannot be run (when > 6 arrays)
#-------------------------------------------------------------------------------------

#' Create several types of spatial images
#'
#' This function (from functions_imagesQC.R) creates several types of
#' spatial images. It does not make use of a PLMset object and is called
#' when the PLM images cannot be computed due to time or memory constaints.
#' Also, it offers more flexibility. By default, it creates
#' a relative intensity plot versus the median array on a blue to
#' red colour scale. Note that the median array will always be computed
#' based on the complete data set, also when plots for a subset of
#' arrays are requested. If the median array has to be computed on
#' subsets (e.g. on experimental group), a data object with only
#' those arrays can be provided to the function (without setting
#' the arrays parameter). The color ranges will saturate at pcut
#' percentage(s) of the data range. The color ranges can be modified
#' by tuning col.mod. By default, for the relative plots, the arrays
#' are first balanced for their overall intensity and a symmetrical
#' color range is used. When there is less than 6 arrays, the median
#' array is not used. Intensities are plotted using a virtual symetric
#' color scale, from blue to red.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch or ExpressionSet)
#' @param pcut (Status: optional, Default:NULL. This means that for relative plots a
#' value of 0.001 will be used, and for other plots values of c(0.01, 0.05).)
#' Either a numeric value or a vector of two numeric
#' values, both in the interval [0, 0.5], used a
#' color saturation limits as a percentage of the
#' data range. If one value is provided, this is
#' taken both a lower and upper saturation limit. (datatype: numeric)
#' @param relative (Status: optional, Default:TRUE) Is the plot to be created
#' a plot of each array relative to the median array, or not (i.e. an
#' absolute plot)? (datatype: logical)
#' @param symm (Status: optional, Default: TRUE if "relative" is TRUE,
#' FALSE otherwise.) Should a symmetric color scale be used? (datatype: logical)
#' @param balance (Status: optional, Default:Default: TRUE if "relative" is TRUE,
#' FALSE otherwise.) Should the arrays first be balanced for their
#' average intensities? (datatype: logical)
#' @param quantitative (Status: optional, Default:TRUE if "relative"is TRUE,
#' FALSE otherwise. Setting "quantitative" to TRUE has no effect if "relative" is FALSE,
#' a warning will be produced.) Should the plot be quantitative or qualitative
#' (i.e. only indicate the sign of the value)? (datatype: logical)
#' @param col.mod (Status: optional, Default:1) A numeric value used as a
#' modifier for the color range. A value of 1 means no modification
#' (linear), smaller leads to faster saturation, larger to slower saturation. (datatype: numeric)
#' @param postfix (Status: optional, Default:"") String to be attached to the file names produced. (datatype: character)
#' @param arrays (Status: optional, Default:NULL (in which case
#' all arrays are plotted)) Which arrays are to be plotted. NOTE: for
#' relative plot types, the median array is still computed using all arrays
#' in the dataset. (datatype: numeric)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @return Images containing each of the requested plots, one file per six arrays.
#' Naming is "virtual_array_plots" followed by the postfix given by the user.
#' @examples
#' # There are two default calls by the script:
#' # 1/ relative intensity plot versus the median array on a blue to red colour scale:
#' # array.image(rawData)
#' # 2/ absolute intensity plot when there are less than 6 arrays in the dataset:
#' # array.image(rawData,relative=FALSE,col.mod=4,symm=TRUE)
#' # Other calls could be made, such as:
#' # spatial intensity plot on a virtual colour scale (red to yellow)
#' # array.image(rawData,relative=FALSE)
#' # signs of the relative intensities versus the median array
#' # (e.g. lower or higher) with blue as lower and red as higher.
#' # array.image(rawData,quantitative=FALSE)
#' # similar to the default plot, but not balanced within arrays
#' # array.image(rawData,balance=FALSE)
#' # or: per experimental group
#' # for (i in levels(experimentFactor)) {
#' # array.image(rawData[, experimentFactor == i], postfix = paste("_",i), balance=FALSE)
#' # }
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

array.image <- function(Data, pcut = NULL, relative = TRUE, symm = relative,
                        balance = relative, quantitative = relative, col.mod = 1, postfix = "", arrays = NULL,
                        WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24) {

  # example calls
  #  array.image(rawData,postfix="_balanced")
  #  array.image(rawData,quantitative=FALSE,postfix="_updown")
  #  array.image(rawData,balance=FALSE)
  #  array.image(rawData,relative=FALSE,postfix="_Intensity")
  # or: per experiment group
  #  for (i in levels(experimentFactor)) {
  #    #remove bad characters because will become part of the file name
  #    array.image(rawData[, experimentFactor == i],
  #    postfix = gsub("[:/|\\?<>*\"]","_",i), balance=FALSE)
  #  }

  gc()
  if (!exists("Data")) {
    stop("ERROR on spatial graph: data parameter must be specified")
  }
  if (class(Data) != "AffyBatch") {
    stop("ERROR on spatial graph: data object must be of class 'AffyBatch'")
  }
  if (!relative & quantitative) {
    warning("When relative is set to FALSE, setting quantitative to TRUE has no effect")
  }
  if (is.null(pcut)) {
    if (relative) {
      pcut <- 0.001
    } else {
      pcut <- c(0, 0.1)
    }
  }
  if (length(pcut) == 1) pcut <- c(pcut, pcut)
  if (sum(pcut < 0 | pcut > .5) > 0) {
    stop("ERROR on spatial graph: pcut value(s) has/ve to be in the interval [0, 0.5]")
  }
  if (is.null(arrays)) arrays <- 1:dim(exprs(Data))[2]
  cat("---preprocessing data---\n")
  data2 <- Data
  if (relative) {
    cat("computing median array\n")
    av_array <- rowMedians(exprs(data2), na.rm = TRUE)
    cat("computing relative expression compared to median array\n")
    exprs(data2) <- exprs(data2) - av_array
  }
  # image function takes the log when operating on AffyBatch function!
  # so we visualise the log ratios as compared to the data set median
  # balance or not?
  if (balance) {
    cat("balancing for array-wide intensity differences\n")
    med <- median(exprs(data2), na.rm = TRUE)
    exprs(data2) <- apply(exprs(data2), 2, function(x) {
      x - median(x, na.rm = TRUE) + med
    })
  }

  if (relative & !quantitative) {
    exprs(data2) <- (exprs(data2) > 0) + 0
  } else {
    # symmetric or not?
    cat("computing color scale limits\n")
    if (!symm) {
      low.lim <- quantile(exprs(data2), pcut[1], na.rm = TRUE)
      up.lim <- quantile(exprs(data2), 1 - pcut[2], na.rm = TRUE)
    } else {
      low <- quantile(exprs(data2), pcut[1], na.rm = TRUE)
      up <- quantile(exprs(data2), 1 - pcut[2], na.rm = TRUE)
      low.lim <- min(low, -up)
      up.lim <- max(-low, up)
    }

    cat("replacing off-scale intensities by scale limits\n")
    exprs(data2)[exprs(data2) < low.lim] <- low.lim
    exprs(data2)[exprs(data2) > up.lim] <- up.lim
  }

  cat("---creating images---")
  col.range.part <- round(((seq(0, 127^col.mod, length.out = 128))^(1 / col.mod)))
  if (relative) {
    col.range <- c(
      rgb(col.range.part, col.range.part, 127, maxColorValue = 127),
      rgb(127, rev(col.range.part), rev(col.range.part), maxColorValue = 127)
    )
  } else {
    # col.range <- heat.colors(170)[1:128]
    col.range <- c(
      rgb(col.range.part, col.range.part, 127, maxColorValue = 127),
      rgb(127, rev(col.range.part), rev(col.range.part), maxColorValue = 127)
    )
    # col.range <- c(rgb(127, col.range.part, 0, max = 127))
  }

  exprs(data2) <- 2^exprs(data2) # image automatically takes log

  if (length(sampleNames(Data)) <= 18) {
    nPerPage <- 6
    nCol <- 2
    cexval <- 1
  } # max 3p of 6 plots
  if (length(sampleNames(Data)) > 18) {
    nPerPage <- 12
    nCol <- 3
    cexval <- 0.8
  } # max 4p of 12 plots
  if (length(sampleNames(Data)) > 48) {
    nPerPage <- 20
    nCol <- 4
    cexval <- 0.8
  }

  nPages <- ceiling(length(sampleNames(Data)) / (nPerPage + 1))

  for (l in 1:nPages) {
    from <- nPerPage * l - (nPerPage - 1)
    to <- min(nPerPage * l, length(sampleNames(Data)))
    png(filename = paste("RawDataArray.image", postfix,
      ifelse(nPages > 1, paste("_", l, sep = ""), ""), ".png",
      sep = ""
    ), width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    layout(matrix(1:nPerPage, ncol = nCol, byrow = TRUE))
    par(oma = c(0, 0, 3, 0), mar = c(1, 1, 3, 1), cex.main = cexval)
    cat("\nimage", l, sep = "")

    for (k in from:to) {
      cat(".")
      image(data2[, arrays][, k], col = col.range)
    }
    mtext(paste("2D virtual images", ifelse(nPages > 1, paste(l, "/", nPages), "")),
      side = 3, outer = TRUE, font = 2, cex = 1.1
    )
    dev.off()
  }
  gc()
  cat("\n")
}


###############
## PNposPlot ##
###############

#' Create an image showing the centres of intensity of the positive and negative border elements
#'
#' This function (from functions_imagesQC.R) creates an image showing the centres of intensity of the positive and negative border elements for each array. It calls the borderQC2 function ( affyQCReport Bioconductor package).
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @return An image of centres of intensity of positive and negative border
#' elements for each array, called "RawDataPosNegPositionPlot"
#' @examples
#' # PNposPlot(rawData)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

PNposPlot <- function(Data,
                      WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24) {
  png(filename = "RawDataPosNegPositions.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
  par(oma = c(17, 0, 0, 0), srt = 90)
  borderQC2(Data)
  mtext("COI absolute values should be < 0.5\n", side = 3, cex = 0.7)
  dev.off()
}

#############
## nuseFun ##
#############

#' Create an image with boxplots of normalized Unscaled Standard Errors (NUSE)
#'
#' This function (from functions_imagesQC.R) creates an image with boxplots of
#' normalized Unscaled Standard Errors (NUSE) for each array. It calls the NUSE
#' function (affyPLM Bioconductor package). An object obtained by calling
#' fitPLM (affyPLM) is needed. If this object is not provided, it will be
#' computed internally within this function.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param Data.pset (Status: optional, Default:NULL) An object obtained by
#' calling fitPLM (affyPLM), used for each but the raw plot. When not provided,
#' it is computed within the function.(datatype: PLMset)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image containing boxplots of the NUSE values per array, called "RawDataNUSEplot"
#' @examples
#' # By default, the script will call:
#' # nuseFun(rawData, Data.pset=rawData.pset, experimentFactor, plotColors, legendColors)
#' # where rawData.pset has been constructed by calling:
#' # rawData.pset <- fitPLM(rawData)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

nuseFun <- function(Data, Data.pset = NULL, experimentFactor = NULL, plotColors = NULL, legendColors = NULL,
                    WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
  if (is.null(plotColors)) stop("The 'plotColors' parameter must be specified")
  if (is.null(legendColors)) stop("The 'legendColors' parameter must be specified")

  if (is.null(Data.pset)) {
    Data.pset <- fitPLM(Data)
  }

  png(filename = "RawDataNUSE.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
  par(oma = c(17, 0, 0, 0), cex.axis = 1)
  NUSE(Data.pset, col = plotColors, main = "Normalized Unscaled Standard Errors (NUSE)", axes = FALSE)
  if (length(sampleNames(Data)) < MAXARRAY) {
    cexval <- 0.65
  } else {
    cexval <- 0.45
  }
  axis(1, at = 1:length(sampleNames(Data)), las = 2, labels = sampleNames(Data), cex.axis = cexval)
  axis(2, cex.axis = 0.7)
  if (length(levels(experimentFactor)) > 1) {
    legend("topright", levels(experimentFactor),
      col = legendColors,
      fill = legendColors, cex = 0.7, bg = "white", bty = "o"
    )
  }

  mtext("NUSE median value should be < 1.1\n", side = 3, font = 1, cex = 0.7)
  dev.off()
}


############
## rleFun ##
############

#' Create an image with boxplots of Relative Log Expression (RLE) values
#'
#' This function (from functions_imagesQC.R) creates an image with boxplots
#' of Relative Log Expression (RLE) values for each array. It calls the RLE
#' function (affyPLM Bioconductor package). An object obtained by calling
#' fitPLM (affyPLM) is needed. If this object is not provided, it will be
#' computed internally within this function.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param Data.pset (Status: optional, Default:NULL) An object obtained by
#' calling fitPLM (affyPLM), used for each but the raw plot. When not provided,
#' it is computed within the function.(datatype: PLMset)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image containing boxplots of the NUSE values per array, called "RawDataRLEplot"
#' @examples
#' # By default, the script will call:
#' # rleFun(rawData, Data.pset=rawData.pset, experimentFactor, plotColors, legendColors
#' # where rawData.pset has been constructed by calling:
#' # rawData.pset <- fitPLM(rawData)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

rleFun <- function(Data, Data.pset = NULL, experimentFactor = NULL, plotColors = NULL, legendColors = NULL,
                   WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
  if (is.null(plotColors)) stop("The 'plotColors' parameter must be specified")
  if (is.null(legendColors)) stop("The 'legendColors' parameter must be specified")

  if (is.null(Data.pset)) {
    Data.pset <- fitPLM(Data)
  }

  png(filename = "RawDataRLE.png", width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
  par(oma = c(17, 0, 0, 0), cex.axis = 1)
  Mbox(Data.pset, col = plotColors, main = "Relative Log Expression (RLE)", axes = FALSE)
  if (length(sampleNames(Data)) < MAXARRAY) {
    cexval <- 0.65
  } else {
    cexval <- 0.45
  }
  axis(1, at = 1:length(sampleNames(Data)), las = 2, labels = sampleNames(Data), cex.axis = cexval)
  axis(2, cex.axis = 0.7)
  if (length(levels(experimentFactor)) > 1) {
    legend("topright", levels(experimentFactor),
      col = legendColors,
      fill = legendColors, cex = 0.7, bg = "white", bty = "o"
    )
  }
  mtext("RLE distributions should be centered around 0\n", side = 3, font = 1, cex = 0.7)
  dev.off()
}


###############
## correlFun ##
###############

#' Create an image of the intensity correlation values
#'
#' This function (from functions_imagesQC.R) creates an image of the
#' intensity correlation values of the arrays in the raw or normalized dataset
#' (depending on the object passed). It calls correlationPlot
#' (affyQCReport Bioconductor package).
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param normMeth (Status: required when Data is a normalized data object, Default:"")
#' String indicating the normalization method used (see
#' \link[ArrayAnalysis]{normalizeData} function for more information on the possible values). (datatype: character)
#' @param clusterOption1 (Status: optional, Default:"Pearson") String indicating
#' the distance function to be used. Possible values are "pearson", “spearman”,
#' or “euclidean”. (datatype: character)
#' @param clusterOption2 (Status: optional, Default:"ward") String indicating
#' the hierarchical clustering function to be used. Possible values are "ward",
#' "single", "complete", "average", "mcquitty", "median" or "centroid". (datatype: character)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image of the intensity correlation values of the arrays, called "DataArrayCorrelationPlot"
#' @examples
#' # By default, before the normalization the script will call:
#' # correlFun(Data=rawData)
#' # and after normalization:
#' # correlFun(Data=normData, normMeth=normMeth)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

correlFun <- function(Data, clusterOption1 = "pearson", clusterOption2 = "ward.D2", normMeth = "",
                      experimentFactor = NULL, legendColors = NULL, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(experimentFactor)) stop("the 'exerimentFactor' parameter is required")
  if (is.null(legendColors)) stop("the 'legendColors' parameter is required")

  #check affy or genefeatured (oligo)
  if(class(Data) == "AffyBatch" || class(Data) == "GeneFeatureSet" || class(Data) == "ExonFeatureSet" ){
    Type <- "Raw"
    text1 <- "Raw data correlation plot"
  } else {
    if (normMeth == "") stop("When providing a normalized data object, the normMeth parameter is required")
    Type <- "Norm"
    text1 <- paste("Array correlation plot\nafter", normMeth, "normalization")
  }
  if (length(sampleNames(Data)) < 2) {
    warning("Only one array in dataset, no correlation plot made")
  } else {
    png(filename = paste(Type, "DataArrayCorrelation.png", sep = ""), width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    if (length(sampleNames(Data)) < MAXARRAY) {
      par(oma = c(17, 0, 0, 0), cex.axis = 0.7, cex.main = 0.8)
      # subval <- 10
    } else {
      par(oma = c(17, 0, 0, 0), srt = 90, las = 2, cex.axis = 0.5, cex.main = 0.8)
      # subval <- 16
    }

    # note: for computing array correlation, euclidean would not make sense
    # only use euclidean distance to compute the similarity of the correlation vectors for the arrays
    COpt1 <- "pearson"
    if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
    crp <- cor(exprs(Data), use = "complete.obs", method = COpt1)

    text1 <- paste(text1, "\ncorrelation method:", COpt1, "\ncluster method:", clusterOption2)

    switch(tolower(clusterOption1),
      "pearson" = {
        my.dist <- function(x) cor.dist(x, abs = FALSE)
      },
      "spearman" = {
        my.dist <- function(x) spearman.dist(x, abs = FALSE)
      },
      "euclidean" = {
        my.dist <- function(x) euc(x)
      }
    )

    my.hclust <- function(d) hclust(d, method = clusterOption2)

    # in order to create some space to put colored symbols as well
    # sampleNames(Data) <- paste(sampleNames(Data)," ")

    sideColors <- legendColors[as.numeric(experimentFactor)]

    heatmap.2(crp,
      distfun = my.dist, hclustfun = my.hclust, trace = "none", symm = TRUE, density.info = "density",
      main = text1, dendrogram = "row", ColSideColors = sideColors
    )

    # correlationPlot(Data)
    # axis(1,side=3,at=seq(from=0.5, to=(length(sampleNames(Data)))-0.5,by=1),
    #    labels=substr(as.character(sampleNames(Data)),1,subval),las=2)
    # par(srt=0)
    # plot(c(0,2), type = 'n', ann = FALSE, axes = FALSE,
    #    frame.plot = FALSE, xlim = c(0, 2), ylim = c(0,2))
    # text(1,1,text1,cex=1)
    dev.off()
  }
}


################
## clusterFun ##
################

#' Create a hierarchical clustering dendrogam
#'
#' This function (from functions_imagesQC.R) creates a hierarchical
#' clustering dendrogam of the arrays in the raw or normalized dataset
#' (depending on the object passed). When the dataset consists of less
#' than three arrays, no dendrogram is generated and a warning is given.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param clusterOption1 (Status: optional, Default:"Pearson") String indicating
#' the distance function to be used. Possible values are "pearson", “spearman”,
#' or “euclidean”. (datatype: character)
#' @param clusterOption2 (Status: optional, Default:"ward") String indicating
#' the hierarchical clustering function to be used. Possible values are "ward",
#' "single", "complete", "average", "mcquitty", "median" or "centroid". (datatype: character)
#' @param normMeth (Status: required when Data is a normalized data object, Default:"")
#' String indicating the normalization method used (see
#' \link[ArrayAnalysis]{normalizeData} function for more information on the possible values). (datatype: character)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param plotSymbols (Status: required, Default:NULL) symbol assigned to each array. (datatype: number)
#' @param legendSymbols (Status: required, Default:NULL) symbol assigned to each experimental group. (datatype: number)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @param MAXARRAY (Status: optional, Default:41) threshold to adapt the image to the number of arrays (datatype: number)
#' @return A PNG image with the clustering dendrogram of the arrays.
#' Naming is ‘DataCluster’ followed by the name of the distance function used
#' @examples
#' # By default, before the normalization the script will call:
#' # clusterFun(Data=rawData, clusterOption1=clusterOption1,
#' # clusterOption2=clusterOption2)
#' # and after normalization:
#' # clusterFun(Data=normData, clusterOption1=clusterOption1,
#' # clusterOption2=clusterOption2, normMeth=normMeth)
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

clusterFun <- function(Data, experimentFactor = NULL, clusterOption1 = "pearson", clusterOption2 = "ward.D2", normMeth = "",
                       plotColors = NULL, legendColors = NULL, plotSymbols = NULL, legendSymbols = NULL, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41) {
  if (is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if (is.null(legendColors)) stop("the 'legendColors' parameter is required")
  if (is.null(plotSymbols)) stop("the 'plotSymbols' parameter is required")
  if (is.null(legendSymbols)) stop("the 'legendSymbols' parameter is required")

  #check affy or genefeatured (oligo)
  if(class(Data)=="AffyBatch" || class(Data) == "GeneFeatureSet" || class(Data) == "ExonFeatureSet") {
    Type <- "Raw"
    main <- "Cluster dendrogram of raw data"
  } else {
    if (normMeth == "") stop("When providing a normalised data object, the normMeth parameter is required")
    Type <- "Norm"
    main <- paste("Cluster dendrogram of", normMeth, "normalized data")
  }
  if (length(sampleNames(Data)) < 3) {
    warning("Only ", length(sampleNames(Data)), " sample(s) in dataset, no clustering plot made")
  } else {
    switch(tolower(clusterOption1),
      "pearson" = {
        correl <- cor.dist(t(exprs(Data)), abs = FALSE)
      },
      "spearman" = {
        correl <- spearman.dist(t(exprs(Data)), abs = FALSE)
      },
      "euclidean" = {
        correl <- euc(t(exprs(Data)))
      }
    )
    if(tolower(clusterOption2)!="ward.d2" && tolower(clusterOption2)!="ward.d") {
      clust <- hclust(correl, method = tolower(clusterOption2))
    } else {
      if(tolower(clusterOption2)=="ward.d2") {
        clust <- hclust(correl, method = "ward.D2")
      } else {
        clust <- hclust(correl, method = "ward.D")
      }
    }
    png(filename = paste(Type, "DataCluster_", clusterOption1, ".png", sep = ""), width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
    if (length(sampleNames(Data)) < MAXARRAY) {
      cexval1 <- 0.75
      cexval2 <- 1.23
      cexval3 <- 0.55
    } else {
      cexval1 <- 0.55
      cexval2 <- 1.6
      cexval3 <- 0.41
    }
    par(cex = cexval1, oma = c(14, 1, 0, 0))
    par(cex.axis = cexval2, cex.lab = cexval2, cex.main = cexval2)
    plot(clust, hang = -1, main = main, xlab = paste("distance:", clusterOption1), sub = paste(" cluster method:", clusterOption2))
    points(1:length(clust$order), rep(0, length(clust$order)), pch = 15, col = "white", cex = 1.5)
    points(1:length(clust$order), rep(0, length(clust$order)), pch = plotSymbols[clust$order], col = plotColors[clust$order])
    if (length(levels(experimentFactor)) > 1) {
      legend("topright", levels(experimentFactor),
        pch = legendSymbols, col = legendColors
      )
    }
    par(cex = cexval3)
    dev.off()
  }
}


############
## pcaFun ##
############

#' Create a Principal Component Analysis (PCA) plot
#'
#' This function (from functions_imagesQC.R) creates a Principal Component
#' Analysis (PCA) plot of the arrays in the raw or normalized dataset
#' (depending on the object passed). When the dataset consists of less
#' than three arrays, no PCA plot is generated and a warning is given.
#' Before computing the PCA each probeset’s expression values are centred
#' on zero. If scaled_pca is TRUE, they will also be rescaled to unit variance.
#' When the maximum length of an array (sample)name is ten characters,
#' and there are no more than 16 samples, the array (sample)names are
#' put within the plot, otherwise they are put in the legend.
#' Since computing a PCA (using the prcomp function) can be memory intensive,
#' a try is used. Furthermore, in cases where scaling is not possible due
#' to loss of any variation, a second attempt is done using no scaling
#' (when scaled_pca had been set to TRUE), and a warning is given.
#' When no PCA can be computed the image is not created, and a warning is given.
#'
#' @param Data (Status: required) The raw data object (datatype: AffyBatch)
#' @param experimentFactor (Status: required, Default:NULL) The factor of groups. (datatype: factor)
#' @param normMeth (Status: required when Data is a normalized data object, Default:"")
#' String indicating the normalization method used (see
#' \link[ArrayAnalysis]{normalizeData} function for more information on the possible values). (datatype: character)
#' @param scaled_pca (Status: optional, Default:TRUE) Should each probeset’s
#' expression be scaled to unit variance before proceeding? Note that
#' the expression is centred on zero in any case. (datatype: logical)
#' @param plotColors (Status: required, Default:NULL) Vector of colors assigned to each array. (datatype: character)
#' @param legendColors (Status: required, Default:NULL) Vector of colors assigned to each experimental group. (datatype: character)
#' @param plotSymbols (Status: required, Default:NULL) symbol assigned to each array. (datatype: number)
#' @param legendSymbols (Status: required, Default:NULL) symbol assigned to each experimental group. (datatype: number)
#' @param namesInPlot (Status: optional, Default:FALSE) Should the array (sample)names
#' be put within the plot, or in the legend? (datatype: logical)
#' @param WIDTH (Status: optional, Default:1000) png image width (datatype: number)
#' @param HEIGHT (Status: optional, Default:1414) png image height (datatype: number)
#' @param POINTSIZE (Status: optional, Default:24) png image point size (datatype: number)
#' @return A PNG image with PCA plot of the arrays. Naming is either ‘Raw’ or
#' the normalization method, followed by "DataPCAanalysis"
#' @examples
#' # By default, before the normalization the script will call:
#' # pcaFun(Data=rawData, experimentFactor=experimentFactor,
#' # plotColors=plotColors, legendColors=legendColors,
#' # namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&&
#' #  (length(sampleNames(rawData))<=16))
#' # and after normalization:
#' # pcaFun(Data=normData, experimentFactor=experimentFactor,
#' # normMeth=normMeth, plotColors=plotColors,
#' # legendColors=legendColors,
#' # namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&&
#' #  (length(sampleNames(rawData))<=16))
#' @import grDevices graphics stats utils methods
#' @importFrom affyQCReport borderQC1 borderQC2
#' @importFrom simpleaffy qc detection.p.val
#' @importFrom yaqcaffy yaqc
#' @importFrom Biobase sampleNames rowMedians
#' @importFrom affy AffyRNAdeg plotAffyRNAdeg getCdfInfo exprs exprs<- probeNames geneNames pm pm<- mm mm<- MAplot rma mas5 Mbox
#' @importFrom bioDist cor.dist spearman.dist euc
#' @importFrom gplots heatmap.2
#' @importFrom biomaRt useMart getBM
#' @importFrom gcrma gcrma
#' @importFrom plier justPlier
#' @importFrom affyPLM fitPLM NUSE
#' @export

pcaFun <- function(Data, experimentFactor = NULL, normMeth = "", scaled_pca = TRUE, plotColors = NULL,
                   legendColors = NULL, plotSymbols = NULL, legendSymbols = NULL, namesInPlot = FALSE, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24) {
  # Scaled PCA by default
  if (is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
  if (is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if (is.null(legendColors)) stop("the 'legendColors' parameter is required")
  if (is.null(plotSymbols)) stop("the 'plotSymbols' parameter is required")
  if (is.null(legendSymbols)) stop("the 'legendSymbols' parameter is required")

  if (length(sampleNames(Data)) < 3) {
    warning("Only", length(sampleNames(Data)), "sample(s) in dataset, no PCA plot made")
  } else {
    #check affy or genefeatured (oligo)
    if(class(Data) == "AffyBatch" || class(Data) == "GeneFeatureSet" || class(Data) == "ExonFeatureSet" ){
      # raw data
      tmain <- "PCA analysis of Raw data"
      Type <- "Raw"
    } else {
      if (normMeth == "") stop("When providing a normalised data object, the normMeth parameter is required")
      tmain <- paste("PCA analysis after", normMeth, "normalization")
      Type <- "Norm"
    }

    pca1 <- NULL
    try(pca1 <- prcomp(t(exprs(Data)[apply(exprs(Data), 1, function(r) {
      sum(is.na(r)) == 0
    }), ]), retx = TRUE, center = TRUE, scale = scaled_pca), TRUE)
    if (is.null(pca1) & scaled_pca) {
      try(pca1 <- prcomp(t(exprs(Data)[apply(exprs(Data), 1, function(r) {
        sum(is.na(r)) == 0
      }), ]), retx = TRUE, center = TRUE, scale = FALSE), TRUE)
      if (!is.null(pca1)) warning("pca with scaling unsuccessful, successfully retried without scaling")
    }
    if (!is.null(pca1)) {
      perc_expl1 <- round(((pca1$sdev[1:3]^2) / sum(pca1$sdev^2)) * 100, 2)

      cex.circle <- 1.5
      cex.text <- 0.7
      cex.legend <- 0.75
      tcol <- "#444444"

      png(
        filename = paste(Type, "DataPCAanalysis.png", sep = ""), width = WIDTH + 200 * (!namesInPlot), height = HEIGHT + 283 * (!namesInPlot),
        pointsize = POINTSIZE
      )

      if (!namesInPlot) {
        layout(rbind(c(1, 1, 2, 2, 5), c(3, 3, 4, 4, 5)))
      } else {
        layout(rbind(c(1, 1, 2, 2), c(1, 1, 2, 2), c(3, 3, 4, 4), c(3, 3, 4, 4)))
      }
      par(oma = c(20, 0, 5, 0))
      plot(pca1$x[, 1], pca1$x[, 2],
        cex = cex.circle, pch = plotSymbols,
        col = plotColors, xlab = paste("PC1 (", perc_expl1[1], "%)", sep = ""),
        ylab = paste("PC2 (", perc_expl1[2], "%)", sep = "")
      )
      if (namesInPlot) text(pca1$x[, 1], pca1$x[, 2], sampleNames(Data), pos = 4, cex = cex.text, col = tcol)
      plot(pca1$x[, 1], pca1$x[, 3],
        cex = cex.circle, pch = plotSymbols,
        col = plotColors, xlab = paste("PC1 (", perc_expl1[1], "%)", sep = ""),
        ylab = paste("PC3 (", perc_expl1[3], "%)", sep = "")
      )
      if (namesInPlot) text(pca1$x[, 1], pca1$x[, 3], sampleNames(Data), pos = 4, cex = cex.text, col = tcol)
      plot(pca1$x[, 2], pca1$x[, 3],
        cex = cex.circle, pch = plotSymbols,
        col = plotColors, xlab = paste("PC2 (", perc_expl1[2], "%)", sep = ""),
        ylab = paste("PC3 (", perc_expl1[3], "%)", sep = "")
      )
      if (namesInPlot) text(pca1$x[, 2], pca1$x[, 3], sampleNames(Data), pos = 4, cex = cex.text, col = tcol)
      barplot((100 * pca1$sdev^2) / sum(pca1$sdev^2), xlab = "components", ylab = "% of total variance explained")

      if (namesInPlot) {
        if (length(levels(experimentFactor)) > 1) {
          legend("topright", levels(experimentFactor),
            pch = legendSymbols, col = legendColors, cex = cex.legend
          )
        }
      } else {
        par(mar = c(0, 0, 0, 0))
        plot(1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
        if (length(levels(experimentFactor)) > 1) {
          legend("topleft", c(levels(experimentFactor), "", sampleNames(Data)),
            #             pch=c(rep(20,length(unique(experimentFactor))+1),plotSymbols,
            pch = c(legendSymbols, 20, plotSymbols),
            col = c(legendColors, "white", plotColors), cex = (cex.legend + 0.1)
            #             ,fill=c(legendColors,rep("white",length(experimentFactor)+1)),
            #             border=c(legendColors,rep("white",length(experimentFactor)+1))
          )
        } else {
          legend("topleft", sampleNames(Data),
            pch = plotSymbols,
            col = plotColors, cex = 0.7, bty = "n"
          )
        }
      }

      mtext(tmain, side = 3, outer = TRUE, font = 2, cex = 1.2)
      dev.off()
    } else {
      warning("PCA on the", Type, "data set unsuccessful, image skipped")
    }
  }
}
