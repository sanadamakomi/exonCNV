#' @title Filter transcripts.
#' @param annoFile Path of annotation file.
#' @param outPath Path to write to.
#' @param transIdxFile Path of index file.
#' @param bedFile Path of bed file.
#' @export
batchFilterTrans <- function(annoFile, outPath, transIdxFile, bedFile) {
    if (is.null(annoFile)) stop('Input --annoFile')
    if (! file.exists(annoFile)) stop(paste0('File no exists', annoFile))
    if (is.null(outPath)) stop('Input --outPath')

    do <- filterTrans(file = annoFile, path = outPath,
                      index = transIdxFile, bed = bedFile)
}

#' @title Export exon bed file.
#' @param annoFile Path of annotation file.
#' @param outPath Path to write to.
#' @param gene A character string of gene symbols seperated by comma(,) to
#'   export.
#' @export
batchExonBed <- function(annoFile, outPath, gene) {
    if (is.null(annoFile)) stop('Input --annoFile')
    if (! file.exists(annoFile)) stop(paste0('File no exists:', annoFile))
    if (is.null(outPath)) stop('Input --outPath')
    if (! is.null(gene)) gene <- unlist(strsplit(gene, ','))
    do <- exonToBed(file = annoFile, path = outPath, gene = gene)
}

#' @title Calulate depth.
#' @param bamFiles A character string of BAM file paths seperated by comma(,).
#' @param bedFile A character string of BED file path.
#' @param outDir Path of directory to write to.
#' @param thread Integer, number of thread.
#' @param mapqFilter A numeric value, least map quality of reads to calculate.
#' @export
batchDepth <- function(bamFiles, bedFile, outDir, thread, mapqFilter) {
    if (is.null(bamFiles)) stop('Input --bamFiles')
    if (is.null(bedFile)) stop('Input --bedFile')
    if (is.null(outDir)) stop('Input --outDir')
    bamFiles <- unlist(strsplit(bamFiles, ','))
    checkBam(bamFiles)
    do <- performCreateCovFile(bamFiles = bamFiles,
                               bedFile = bedFile,
                               outDir = outDir,
                               thread = thread,
                               batch = 1000,
                               mapq.filter = mapqFilter)
}

#' @title Merge coverage files.
#' @param covFiles A character string of coverage file paths seperated by
#'   comma(,).
#' @param outPath Path to write to.
#' @export
batchMerge <- function (covFiles, outPath) {
    if (is.null(covFiles)) stop('Input --covFiles')
    if (is.null(outPath)) stop('Input --outPath')
    covFiles <- unlist(strsplit(covFiles, ','))
    for (f in covFiles) {
        if (! file.exists(f)) stop(paste0('File no exists: ', f))
    }
    do <- mergerCovFiles(covFiles, path = outPath)
}

#' @title Calculate metrics.
#' @param bamFiles A character string of BAM file paths seperated by comma(,).
#' @param covFiles A character string of coverage file paths seperated by
#'   comma(,).
#' @param outPath Path to write to.
#' @export
batchMetrics <- function (bamFiles, covFiles, outPath) {
    if (is.null(bamFiles)) stop('Input --bamFiles')
    if (is.null(covFiles)) stop('Input --covFiles')
    if (is.null(outPath)) stop('Input --outPath')
    bamFiles <- unlist(strsplit(bamFiles, ','))
    covFiles <- unlist(strsplit(covFiles, ','))
    checkBam(bamFiles)
    for (f in covFiles) {
        if (! file.exists(f)) stop(paste0('File no exists: ', f))
    }
    do <- getParameters(bamFiles = bamFiles, covFiles = covFiles, path = outPath)
}

#' @title Calculate metrics.
#' @param covFiles A character string of coverage file paths seperated by
#'   comma(,).
#' @param metricsFile Path of metrics file.
#' @param outPath Path to write to.
#' @param defaultDir Path of directory storing defaultCovFiles and
#'   defaultMetricsFile.
#' @param defaultCovFiles A character string of default coverage file paths
#'   seperated by comma(,).
#' @param defaultMetricsFile Path of default metrics file.
#' @param genderNum Integer, if counts of samples shares the same gender is less
#'   than genderNum, it will use default data as baseline.
#' @export
batchMergeDefault <- function(covFiles, metricsFile, outPath, defaultDir, defaultCovFiles, defaultMetricsFile, genderNum) {

    # option part 1
    if (is.null(covFiles)) stop('Input --covFiles')
    if (is.null(metricsFile)) stop('Input --metricsFile')
    if (is.null(outPath)) stop('Input --outPath')
    covFiles <- unlist(strsplit(covFiles, ','))
    for (f in c(covFiles, metricsFile)) {
        if (! file.exists(f)) stop(paste0("File no exists: ", f))
    }
    test.ids <- sapply(covFiles, getID, ptn = '.cov')
    ifMatchMetrics(x = test.ids, y = metricsFile)
     # gender counts
    test.dat <- read.table(metricsFile, header = TRUE, sep = "\t", quote = "",
                           comment.char = "#", na.strings = "NA",
                           fill = TRUE, stringsAsFactors = FALSE)
    test.gender <- as.data.frame(table(as.character(as.vector(test.dat[, 'gender']))))
    need.default.gender <- test.gender[which(test.gender[,'Freq'] < genderNum), 'Var1']
    need.default.gender <- as.character(as.vector(need.default.gender))
    if (length(need.default.gender) == 0) {
        write("No need to use default samples.", stdout())
        batchMerge(covFiles = paste(covFiles, collapse = ','), outPath = outPath)
        return()
    }
    write("Need to use default samples.", stdout())

    # option part 2
    if (is.null(defaultDir)) {
        if (is.null(defaultCovFiles) | is.null(defaultMetricsFile)) {
            stop("Input --defaultDir or --defaultCovFiles & --defaultMetricsFile")
        }
    } else {
        if (! file.exists(defaultDir)) stop(paste0("Directory no exists: ", defaultDir))
    }
    if (is.null(defaultCovFiles)) {
        defaultCovFiles <- list.files(path = defaultDir, pattern=".cov$",
            all.files = TRUE, full.names = TRUE, recursive = FALSE, include.dirs = FALSE)
    } else {
        defaultCovFiles <- unlist(strsplit(defaultCovFiles, ','))
        for (f in defaultCovFiles) {
            if (! file.exists(f)) stop(paste0("File no exists: ", f))
        }
    }
    if (is.null(defaultMetricsFile)) {
        defaultMetricsFile <- list.files(path = defaultDir,
            pattern = "metrics.txt", all.files = TRUE, full.names = TRUE, recursive = FALSE, include.dirs = FALSE)
    } else {
        if (! file.exists(defaultMetricsFile)) stop(paste0("File no exists: ", defaultMetricsFile))
    }
    default.ids <- sapply(defaultCovFiles, getID, ptn = '.cov')
    ifMatchMetrics(x = default.ids, y = defaultMetricsFile, n = length(default.ids))
    default.dat <- read.table(defaultMetricsFile, header = TRUE, sep = "\t", quote = "",
                              comment.char = "#", na.strings = "NA",
                              fill = TRUE, stringsAsFactors = FALSE)
    add.default.dat <- default.dat[which(as.character(as.vector(default.dat[, "gender"]))
                                         %in% need.default.gender),,drop = FALSE]
    add.covFiles <- defaultCovFiles[which(default.ids %in% add.default.dat$id)]
    update.test.dat <- rbind(test.dat, add.default.dat)
    write.table(update.test.dat, file = normalizePath(metricsFile), row.names = FALSE, col.names = TRUE,
                quote = FALSE, sep = "\t")
    do <- batchMerge(covFiles = paste(c(covFiles,add.covFiles), collapse = ','), outPath = outPath)

    write(paste0("Write to path: \n", normalizePath(outPath)), stdout())
}

#' @title Fit model to data.
#' @param mergedCovFile Path of merged coverage file.
#' @param metricsFile Path of metrics file.
#' @param outPath Path to write to.
#' @param lowDepth A numeric value, regions that avaerage depth less than this
#'   value will replaced by NA.
#' @export
batchFit <- function(mergedCovFile, metricsFile, outPath, lowDepth) {
    if (is.null(mergedCovFile)) stop('Input --mergedCovFile.')
    if (is.null(metricsFile)) stop('Input --metricsFile')
    if (is.null(outPath)) stop('Input --outPath')
    if (! file.exists(mergedCovFile)) stop(paste0('File no exists: ', mergedCovFile))
    if (! file.exists(metricsFile)) stop(paste0('File no exists: ', metricsFile))
    do <- performFitPoisson(mergedCovFile = mergedCovFile,
                            parameterFile = metricsFile,
                            path = outPath,
                            lowdepth = lowDepth)
}

#' @title Call exon CNV.
#' @param probFile Path of probability file.
#' @param mergedCovFile Path of merged coverage file.
#' @param outDir Path of directory to write to.
#' @param sampleId Sample ids to call CNV, seperated by comma(,).
#' @param cutoff A list of cutoff parameters to filter exon CNV, 'prob' is
#'   probability, 'pool.count' is the least counts of sample of which gene has
#'   CNV. 'baseline' is the least baseline depth to call CNV. 'lowdepth' is
#'   least depth in all samples to call CNV.
#' @export
batchCall <- function(probFile, mergedCovFile, outDir, sampleId, cutoff = list(prob = 1E-4, pool.count = 1, baseline = 50, lowdepth = 10)) {
    if (is.null(probFile)) stop('Input --probFile')
    if (is.null(mergedCovFile)) stop('Input --mergedCovFile')
    if (is.null(sampleId)) write('Default calling all samples.', stdout())
    if (! file.exists(probFile)) stop(paste0('File no exists: ', probFile))
    if (! file.exists(mergedCovFile)) stop(paste0('File no exists: ', mergedCovFile))
    if (! is.null(sampleId)) sampleId <- unlist(strsplit(sampleId, ','))
    do <- callExonCNV(probFile = probFile,
                      mergedCovFile = mergedCovFile,
                      path = outDir,
                      sample.id = sampleId,
                      cutoff = cutoff)
}

#' @title Export VCF file.
#' @param geneCNVFile Path of gene CNV result file.
#' @param outPath Path to write to.
#' @export
batchExportVCF <- function(geneCNVFile, outPath) {
    if (is.null(geneCNVFile)) stop('Input --geneCNVFile')
    if (is.null(outPath)) stop('Input --outPath')
    if (! file.exists(geneCNVFile)) stop(paste0('File no exists: ', geneCNVFile))
    do <- outputVcf(geneCNVFile, path = outPath)
}

#' @title Plot exon CNV result
#' @param geneCNVFile Path of gene CNV result file.
#' @param mergedCovFile Path of merged coverage file.
#' @param metricsFile Path of metrics file.
#' @param gene A character string of gene symbols to plot, seperated by comma(,).
#' @param outPath Path to write to.
#' @param draw Number of genes to plot.
#' @export
batchPlot <- function(geneCNVFile, mergedCovFile, metricsFile, gene, outPath, draw) {
    if (is.null(geneCNVFile)) stop('Input --geneCNVFile')
    if (is.null(mergedCovFile)) stop('Input --mergedCovFile')
    if (is.null(metricsFile)) stop('Input --metricsFile')
    if (is.null(outPath)) stop('Input --outPath')
    if (! is.null(gene)) gene <- unlist(strsplit(gene, ','))
    if (! file.exists(geneCNVFile)) stop(paste0('File no exists: ', geneCNVFile))
    if (! file.exists(mergedCovFile)) stop(paste0('File no exists: ', mergedCovFile))
    if (! file.exists(metricsFile)) stop(paste0('File no exists: ', metricsFile))
    do <- plotExonCNV(geneCNVFile = geneCNVFile,
                      mergedCovFile = mergedCovFile,
                      parameterFile = metricsFile,
                      gene = gene,
                      path = outPath,
                      draw = draw)
}
