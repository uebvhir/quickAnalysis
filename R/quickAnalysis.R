#' quickAnalysis Function
#'
#' La función quickAnalysis realiza para una matriz de datos donde las filas corresponden a los genes y las columnas a las muestras un
#' análisis de expresión diferencial entre dos grupos concretos. Para el análisis usa dos métodos: un t.test para comparaciones multiples
#' (tabla resultante: reT) y un análisis limma básico (tabla resultante: topTable). Además, para aquellos genes diferencialmente expresados
#' realiza los boxplots para cada uno de los grupos.
#' @param expres a matrix-like data object containing log-expression or numeric values, with rows corresponding to genes and columns to samples.
#' @param groupingVar factor which codes the grouping to be tested. There must be 1 or 2 groups.The length of the factor needs to correspond to the sample size.
#' @param min minuend
#' @param sust subtrahend
#' @param aFName the name to the output file.
#' @param outputDir the path to the output file.
#' @param outputFType character string with the abbreviation for extension: "xls" or "html. Default value is NULL.
#' @param pvalThreshold setting the threshold p-value. Default value is 1.
#' @param useAdjP a logical value indicating use adjust p-value. Default value is TRUE.
#' @param plotSelected  A logical value indicating whether the output is a plot. Default value is TRUE.
#' @param plot2pdf A logical value indicating whether the plot is a pdf. Default value is FALSE.
#' @export quickAnalysis
#' @import beeswarm multtest limma xlsx hwriter
#' @author Alex Sánchez \email{alex.sanchez@@vhir.org} , Miriam Mota  \email{miriam.mota@@vhir.org}
#' @examples
#' par(mfrow=c(2,2))
#' quickAnalysis(expres = exprs(sampleSet),
#' groupingVar = pData(sampleSet)$group,
#' min = "Case", sust = "Control", pvalThreshold = 0.25)
#' @return genes: selected genes
#' @return resT: multtest
#' @return topTable: limma
#' @keywords quickAnalysis limma multtest plot differential expression microarrays


quickAnalysis <- function(expres, groupingVar,
                          min, sust,  #min = "Case", sust = "Control" ## en ordre Case - Control
                          aFName = "Results",
                          outputDir = ".", outputFType = NULL,
                          pvalThreshold = 1, useAdjP = TRUE,
                          plotSelected = TRUE, plot2pdf = FALSE)
{
  print("ANALYIS:...PROCESS STARTED")
  groupingVar <- factor(groupingVar, c(sust,min))


  ### ANALYSIS USING multtest library
  print(paste("ANALYSIS using multtest library"))
  resT <- mt.maxT(expres, classlabel = groupingVar)
  if (useAdjP) {
    resT.selected <- resT[resT$adjp <= pvalThreshold,]
  }else{
    resT.selected <- resT[resT$rawp <= pvalThreshold,]
  }
  gNames.multtest <- rownames(resT.selected)

  ### ANALYSIS USING limma
  print(paste("ANALYSIS using limma"))
  design <- model.matrix(~ -1 + groupingVar)
  if (all.equal(substr(colnames(design) ,
                       start = nchar(colnames(design) ) - nchar(levels(groupingVar)) + 1,
                       stop = nchar(colnames(design) )) ,
                levels(groupingVar))) {
    colnames(design) <- levels(groupingVar)
  }
  rownames(design) <- colnames(expres)

  print("Desing: ")
  print(design)

  contrastNames <- paste(levels(groupingVar)[2],levels(groupingVar)[1],sep = "-")
  contrastsMatrix <- matrix(c(-1,1),nrow = ncol(design))
  rownames(contrastsMatrix) <- colnames(design)
  colnames(contrastsMatrix) <- contrastNames
  print("Contrast matrix:")
  print(contrastsMatrix)


  fit <- lmFit(expres, design)
  fit.main <- contrasts.fit(fit, contrastsMatrix)
  fit.main <- eBayes(fit.main)

  top.Diff.all <- topTable(fit.main, n = nrow(expres), adjust = "fdr")

  top.Diff.all[,paste0("mean",min)] <-
    apply(expres,1, function(x) mean(x[groupingVar == min]))[rownames(top.Diff.all)]

  top.Diff.all[,paste0("mean",sust)] <-
    apply(expres,1, function(x) mean(x[groupingVar == sust]))[rownames(top.Diff.all)]

  if (useAdjP) {
    top.Diff <- top.Diff.all[top.Diff.all$adj.P.Val <= pvalThreshold,]
  }else{
    top.Diff <- top.Diff.all[top.Diff.all$P.Value <= pvalThreshold,]
  }

  gNames.limma <- rownames(top.Diff)
  selectedGenes <- union(gNames.multtest, gNames.limma)

  top.Diff.common <- top.Diff.all[selectedGenes, ]
  resT.common <- resT[selectedGenes, ]

  cat(paste("\nNumber of genes selected using permutations test (p <",pvalThreshold,"): ",
            nrow(resT.selected),sep = ""), "\n")
  cat(paste("Number of genes selected using linear model test (p <",pvalThreshold,"): ",
            nrow(top.Diff),sep = ""), "\n")
  cat(paste("Number of genes selected by both test (p <",pvalThreshold,")          : ",
            length(intersect(gNames.multtest,gNames.limma)),sep = ""), "\n")

  cat(paste("ANALYIS:", "... PROCESS COMPLETED", sep = " "), "\n")

  if (!is.null(outputFType)) {
    if (outputFType == "xls") {
      xlsFName <- file.path(outputDir,paste(aFName, "xls", sep = "."))
      write.xlsx(top.Diff, file = xlsFName, sheetName = "limma")
      write.xlsx(resT.selected, file = xlsFName, append = TRUE, sheetName = "multtest")
      cat(paste("RESULTS are in file",xlsFName, sep = " "), "\n")
    }
    if (outputFType == "html") {
      aFName1 <- file.path(outputDir,paste(aFName, "multtest","html", sep = "."))
      hwrite(resT.selected, page = aFName1)
      aFName2 <- file.path(outputDir,paste(aFName, "limma","html", sep = "."))
      hwrite(top.Diff, page = aFName2)
      cat(paste("RESULTS are in file",aFName1, aFName2, sep = " "), "\n")
    }
  }


  if (plotSelected && (length(selectedGenes) > 0)) {
    plotsFName <- file.path(outputDir,paste(aFName,"Plots", "pdf", sep = "."))
    if (plot2pdf) {pdf(file.path(outputDir, plotsFName))}
    myExpres <- as.matrix(expres[selectedGenes, ])
    varName <- "Gene Expression"
    for (i in 1:nrow(myExpres)) {
      desc <- paste("logFC=", round(top.Diff.common$logFC[i],3),
                    ", p-val=", round(top.Diff.common$P.Value[i],6),
                    ", Adj-p=", round(top.Diff.common$adj.P.Val[i],6), sep = "")
      beeswarm(myExpres[i,]~factor(groupingVar, c(min,sust)),
               ylab = "Expression", xlab = "Groups",
               main = paste(rownames(myExpres)[i],desc, sep = "\n"),
               labels = levels(factor(groupingVar, c(min,sust))))
      boxplot(myExpres[i,]~factor(groupingVar, c(min,sust)), add = T, names = c("",""), col = "#0000ff22")
      # Segons un post de: https://www.r-statistics.com/2011/03/beeswarm-boxplot-and-plotting-it-with-r/
    }
    if (plot2pdf) {
      dev.off()
      cat(paste("PLOTS are in file", plotsFName, sep = " "), "\n")
    }
  }
  return(list(genes = selectedGenes,
              resT = resT,
              topTable = top.Diff.all
  )
  )
}
