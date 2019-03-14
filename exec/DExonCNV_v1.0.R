#!/usr/bin/Rscript
# Copyright(C) by Zhan-Ni Chen, version 1.0, 2019.3.14 compiled by R 3.5.1*

# dafault
VERSION="0.0.1"
file=NULL
bamFiles=NULL
covFiles=NULL
metricsFile=NULL
mergedCovFile=NULL
probFile=NULL
geneCNVFile=NULL
bedFile=NULL
annoFile=NULL
transIdxFile=NULL
sampleId=NULL
gene=NULL
outDir=NULL
outPath=NULL
defaultDir=NULL
defaultCovFiles=NULL
defaultMetricsFile=NULL
thread=1
lowDepth=10
probCutoff=1E-4
poolCount=3
baselineDepth=50
mapqFilter=30
draw=10
genderNum=3
mode="help"
modeHelp="help"

# options
library(GetoptLong)
GetoptLong(
    "modeHelp=s", "Please input --modeHelp name to check options in each mode",
    "mode=s", "'filterTrans', 'exonBed', 'depth', 'metrics', 'merge', 'mergeDefault', 'fit', 'call', 'export', 'format' or 'plot'",
    "file=s", "Path of file",
    "bamFiles=s", "BAM file or files' path, separated by comma(,)",
    "covFiles=s", "Coverage files' path, separated by comma(,)",
    "metricsFile=s", "Metrics file's path",
    "mergedCovFile=s", "Merged coverage file's path",
    "probFile=s", "Probability file's path",
    "geneCNVFile=s", "Gene CNV result file's path",
    "bedFile=s", "Bed file's path",
    "annoFile=s", "Annotation file's path",
    "transIdxFile=s", "A file that column1 is Gene and column 2 is transcript name",
    "sampleId=s", "Sample ids to call CNV, separated by comma(,)",
    "gene=s", "Gene symbols, separated by comma(,)",
    "outDir=s", "Output directory path",
    "outPath=s", "Output file path",
    "defaultDir=s", "Directory path of default data",
    "defaultCovFiles=s", "Default coverage files' path, separated by comma(,)",
    "defaultMetricsFile=s", "Default metrics file's path",
    "thread=i", "Number of thread [1]",
    "lowDepth=f", "Region of which depth less than lowdepth will not be called [10]",
    "probCutoff=f", "Cutoff of probability [0.0001]",
    "poolCount=i", "Cutoff of sample count in this pool [3]",
    "baselineDepth=f", "Cutoff of baseline depth [50]",
    "mapqFilter=i", "Reads of which mapq less than mapqFilter will not be calculated [30]",
    "draw=i", "Counts of gene to plot [10]",
    "genderNum=i", "Minimum sample counts to call CNV on sex chromosome [3]",
    head = 'An example to show how to use the script',
    foot = 'Please contact znchen@genokon.com for comments'
)

# mode help
if (modeHelp == "filterTrans") {
    write("Usage: Filter an annotation file by a index file or a BED file.", stdout())
    write("    --annoFile\tAnnotation file's path, details see http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz.", stdout())
    write("    --outPath\tOutput file path.", stdout())
    write("Optional:", stdout())
    write("    --transIdxFile\tA file that column1 is Gene and column 2 is transcript name.", stdout())
    write("    --bedFile\tBed file's path.", stdout())
    q()
} else if (modeHelp == "exonBed") {
    write("Usage: Extract exon infromation from annotation file.", stdout())
    write("    --annoFile\tAnnotation file's path, details see http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz.", stdout())
    write("    --outPath\tOutput file path.", stdout())
    write("Optional:", stdout())
    write("    --gene\tGene symbols in the annotation file, separated by comma(,).", stdout())
    q()
} else if (modeHelp == "depth") {
    write("Usage: Calculate coverage of depth.", stdout())
    write("    --bamFiles\tBAM file or files' path, separated by comma(,).", stdout())
    write("    --bedFile\tBed file's path.", stdout())
    write("    --outDir\tOutput directory path.", stdout())
    write("Optional:", stdout())
    write("    --thread\tNumber of thread [1]", stdout())
    write("    --mapqFilter\tReads of which mapq less than mapqFilter will not be calculated [30]", stdout())
    q()
} else if (modeHelp == "metrics") {
    write("Usage: Calculate metrics.", stdout())
    write("    --bamFiles\tBAM file or files' path, separated by comma(,).", stdout())
    write("    --covFiles\tCoverage files' path, separated by comma(,).", stdout())
    write("    --outPath\tOutput file path.", stdout())
    q()
} else if (modeHelp == "merge") {
    write("Usage: Merge coverage files.", stdout())
    write("    --covFiles\tCoverage files' path, separated by comma(,).", stdout())
    write("    --outPath\tOutput file path.", stdout())
    q()
} else if (modeHelp == "mergeDefault") {
    write("Usage: Fitting model to data.", stdout())
    write("    --covFiles\tCoverage files' path, separated by comma(,).", stdout())
    write("    --metricsFile\tMetrics file's path.", stdout())
    write("    --outPath\tOutput file path.", stdout())
    write("Optional:", stdout())
    write("    --defaultDir\tDirectory path of default data.", stdout())
    write("    --defaultCovFiles\tDefault coverage files' path, separated by comma(,).", stdout())
    write("    --defaultMetricsFile\tDefault metrics file's path.", stdout())
    write("    --genderNum\tMinimum sample counts to call CNV on sex chromosome [3].", stdout())
    q()
}else if (modeHelp == "fit") {
    write("Usage: Fitting model to data.", stdout())
    write("    --mergedCovFile\tMerged coverage file's path.", stdout())
    write("    --metricsFile\tMetrics file's path.", stdout())
    write("    --outPath\tOutput file path.", stdout())
    write("Optional:", stdout())
    write("    --lowDepth\tRegion of which depth less than lowdepth will not be called [10].", stdout())
    q()
} else if (modeHelp == "call") {
    write("Usage: Call exon CNV.", stdout())
    write("    --probFile\tProbability file's path.", stdout())
    write("    --mergedCovFile\tMerged coverage file's path.", stdout())
    write("    --outDir\tOutput directory's path.", stdout())
    write("Optional:", stdout())
    write("    --sampleId\tSample ids to call CNV, separated by comma(,).", stdout())
    write("    --lowDepth\tRegion of which depth less than lowdepth will not be called [10].", stdout())
    write("    --probCutoff\tCutoff of probability [0.0001].", stdout())
    write("    --poolCount\tCutoff of sample count in this pool [3].", stdout())
    write("    --baselineDepth\tCutoff of baseline depth [50].", stdout())
    q()
} else if (modeHelp == "export") {
    write("Usage: Input exon cnv tabular result file and export VCF file.", stdout())
    write("--geneCNVFile\tGene CNV result file's path.", stdout())
    write("--outPath\tOutput file path.", stdout())
    q()
} else if (modeHelp == "plot") {
    write("Usage: Plot exon cnv result.", stdout())
    write("    --geneCNVFile\tGene CNV result file's path.", stdout())
    write("    --mergedCovFile\tMerged coverage file's path.", stdout())
    write("    --metricsFile\tMetrics file's path.", stdout())
    write("    --outPath\tOutput file path.", stdout())
    write("Optional:", stdout())
    write("    --gene\tGene symbols in the annotation file, separated by comma(,).", stdout())
    write("    --draw\tCounts of gene in gene CNV result file to plot [10].", stdout())
    q()
} else if (modeHelp == "format") {
    write("Usage: Input annotation tabular result file and output reformat tabular file.", stdout())
    write("--file\tPath of *.hg19_multianno.txt file.", stdout())
    write("--outPath\tOutput file path.", stdout())
    q()
} else {
    if (! modeHelp == "help") stop(paste0('Error modeHelp "', modeHelp, '"'))
    if (mode == "help" & modeHelp == "help") {
        write("Please input --help to see details.", stdout())
        write("And input --modeHelp modename to see each mode's details.", stdout())
        q()
    }
}

# imporrt packages
suppressMessages(library(exonCNV))

# performing
cutoff <- list(prob = probCutoff, pool.count = poolCount, baseline = baselineDepth, lowdepth = lowDepth)
if (mode == "filterTrans") {

    batchFilterTrans(annoFile = annoFile, outPath = outPath,
                     transIdxFile = transIdxFile, bedFile = bedFile)

} else if (mode == "exonBed") {

    batchExonBed(annoFile = annoFile, outPath = outPath, gene = gene)

} else if (mode == "depth") {

    batchDepth(bamFiles = bamFiles, bedFile = bedFile,
               outDir = outDir, thread = thread, mapqFilter = mapqFilter)

} else if (mode == "metrics") {

    batchMetrics(bamFiles = bamFiles, covFiles = covFiles, outPath = outPath)

} else if (mode == "merge") {

    batchMerge(covFiles = covFiles, outPath = outPath)

} else if (mode == "mergeDefault") {

    batchMergeDefault(covFiles = covFiles, metricsFile = metricsFile, outPath = outPath,
                      defaultDir = defaultDir, defaultCovFiles = defaultCovFiles,
                      defaultMetricsFile = defaultMetricsFile, genderNum = genderNum)

} else if (mode == "fit") {

    batchFit(mergedCovFile = mergedCovFile, metricsFile = metricsFile,
             outPath= outPath, lowDepth = lowDepth)

} else if (mode == "call") {

    batchCall(probFile = probFile, mergedCovFile = mergedCovFile,
              outDir = outDir, sampleId = sampleId, cutoff = cutoff)

} else if (mode == "export") {

    batchExportVCF(geneCNVFile = geneCNVFile, outPath = outPath)

} else if (mode == "plot") {

    batchPlot(geneCNVFile = geneCNVFile, mergedCovFile = mergedCovFile,
              metricsFile = metricsFile, gene = gene, outPath = outPath, draw = draw)

} else if (mode == "format") {

    if (is.null(file)) stop("Input --file.")
    if (! file.exists(file)) stop(paste0("File no exists: ", file))
    batchChangeFormat(file = file, outPath = outPath)

} else {

    if (! mode == "help") stop(paste0('Error Mode "', mode, '"'))

}
