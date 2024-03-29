#' @title Filtering transcripts.
#' @description Filtering annotation file. If input an index file, which has two
#'   columns indicating \strong{gene symbol} and \emph{transcript name} without
#'   a header line, it will return transcripts based on names in the index file.
#'   If more than one transcript share the same name, it will return the
#'   transcript which has the most exons. If input a BED file, i.e. a
#'   target-capture region provided by the kit, it will filter out regions out
#'   of BED file.
#' @param file A character string of the annotation file's path.
#' @param path A character string of the file path to write to.
#' @param bed A character string of the BED file's path.
#' @param index A character string of the Index file's path.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomeInfoDb
#' @author Zhan-Ni Chen
filterTrans <- function(file, path = NULL, bed = NULL, index = NULL) {
    if (is.null(bed) & is.null(index)) {
        write('Nothing to do. Please input --bedFile or --transIdxFile to filter.', stdout())
        return()
    }
    if (! is.null(index)){
        if (! file.exists(index)) stop(paste0(index, 'no exists.'))
    }
    if (! is.null(bed)){
        if (! file.exists(bed)) stop(paste0(bed, 'no exists.'))
    }
    if (is.null(path)) path <- 'annotation.txt'

    write('Start to filter transcript...', stdout())

    annoDat <- readFile(file)
    if (!is.null(bed)) {

        write(paste0('Filter by BED file: ', normalizePath(bed)), stdout())

        bed.gr <- bedtoGRange(bed)
        chr <- as.character(as.vector(annoDat$V3))
        if (! grepl('chr', as.character(seqnames(bed.gr[1])))) {
            seqls <- paste0('chr', seqlevels(bed.gr))
            seqlevels(bed.gr) <- seqls
        }
        anno.gr <- GRanges(Rle(chr), ranges = IRanges(start = annoDat$V5, end = annoDat$V6))
        hit <- findOverlaps(anno.gr, bed.gr)
        annoDat <- annoDat[unique(queryHits(hit)),,drop = FALSE]
    }

    if (!is.null(index)) {

        write(paste0('Filter by index file: ', normalizePath(index)), stdout())

        indexDat <- readFile(index)
        annoDat <- annoDat[which(paste0(annoDat$V13, '_', annoDat$V2)
                                 %in% paste0(indexDat$V1,"_", indexDat$V2)),, drop=FALSE]
        annoDat.list <- split(annoDat, as.factor(annoDat$V13))
        annoDat <- lapply(annoDat.list, function(x){
            ec <- as.numeric(as.vector(x$V9))
            idx <- which(ec == max(ec)) # select the transcript which owns the most exon counts.
            idx <- idx[1]
            x[idx,, drop=FALSE]
        })
        annoDat <- do.call('rbind', annoDat)
    }

    cat(apply(annoDat, 1, function(x) {paste(x, collapse = '\t')}), sep = "\n", file = path)

    write(paste0("Write to path: \n", normalizePath(path)), stdout())
}

#' @title Bed of exon positions.
#' @description Extract exon informations from the annotation tabular file and
#'   output a BED file with 4 columns indicating chromosome, start and end
#'   position, exon id.
#' @param file A character string of the annotation file's path.
#' @param path A character string of the file path to write to.
#' @param gene A character string or vector of gene symbols to export.
#' @param expand An integer. If an exon's size is less than this value, its
#'   region will be expanded centrally to a width of this value.
#' @export
#' @importFrom methods as
#' @author Zhan-Ni Chen
exonToBed <- function(file, path = NULL, gene = NULL, expand=-1) {
    if (! file.exists(file)) stop(paste0(file, 'no exists.'))
    if (is.null(path)) path <- paste0(normalizePath('.'), '/exons.bed')
    indat <- readFile(file)
    allgene <- as.character(as.vector(indat[,13]))
    if (! is.null(gene)) {
        idx <- which(allgene %in% gene)
        if (length(gene) == 0) stop('Find no gene in file.')
        indat <- indat[idx,, drop = FALSE]
    }
    write("Start to creat exon BED file...", stdout())

    exon.grls <- sapply(1:nrow(indat), function(x){
        chr <- gsub('chr', '', indat[x, 3])
        exonS <- as.numeric(unlist(strsplit(indat[x, 10], ",")))
        exonE <- as.numeric(unlist(strsplit(indat[x, 11], ",")))
        strand <- as.character(indat[x, 4])
        if(strand == "+") {
            exon <- seq(1, length(exonS),1)
        }else{
            exon <- seq(length(exonS), 1, -1)
        }

        if (expand > 0) {
            r_idx <- which((exonE - exonS) < expand)
            if (length(r_idx) > 0) {
                e_size <- round((expand - (exonE[r_idx] - exonS[r_idx]))/2)
                exonS[r_idx] <- exonS[r_idx] - e_size
                exonE[r_idx] <- exonE[r_idx] + e_size
                exonS[which(exonS < 1)] <- 1
            }
        }
        GRanges(Rle(rep(chr, length(exonS))),
                ranges = IRanges(start = exonS, end = exonE),
                id = paste0(indat[x, 13], '(', indat[x, 2], ')' , '_', exon))
    })
    exon.grls <- as(exon.grls, "GRangesList")
    exon.grls <- unlist(exon.grls)
    exon.grls <- as.data.frame(exon.grls)
    exon.grls$seqnames <- factor(exon.grls$seqnames,
                                 levels = c(as.character(seq(1, 22, 1)), 'X', 'Y'))
    exon.grls <- exon.grls[with(exon.grls, order(exon.grls$seqnames, exon.grls$start, exon.grls$end)), ]
    exon.grls$id <- paste0(exon.grls$id, '_', seq(1, nrow(exon.grls), 1))
    cat(apply(exon.grls[,c("seqnames", "start", "end", "id"), drop = FALSE], 1, function(x)
        {paste(x, collapse = '\t')}), sep = "\n", file = path)

    write(paste0("Write to path: \n", normalizePath(path)), stdout())
}

#' @title Check BAM file.
#' @description It will stop if BAM file is illegal.
#' @param x A character string or vector of BAM File path.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom Rsamtools BamFile
#' @author Zhan-Ni Chen ####### Check BAM file #######
checkBam <- function(x) {
    a <- sapply(x, function(bam) {
        if (! file.exists(bam)) stop(paste0(bam, ' file is missing.'))
        bai <- BamFile(bam)$index
        if (is.na(bai)) stop(paste0(bam, '.bai index file is missing.'))
        bam_info <- file.info(bam)
        bai_info <- file.info(bai)
        dt <- difftime(bai_info$mtime, bam_info$mtime, units = 'secs')
        dt <- as.numeric(dt)
        if (dt < 0) stop(paste0(bam, ' index file is older than BAM file.'))
    })
}

#' @title Count total reads.
#' @param x A character string or vector of BAM File path.
#' @return A numeric value.
#' @importFrom Rsamtools idxstatsBam
countTR <- function(x) {
    checkBam(x)
    df <- idxstatsBam(x)
    sum(df[which(df$seqnames %in%
                     c(as.character(seq(1, 22, 1)), 'X', 'Y')),3])
}

#' @title Call gender from coverage file.
#' @description Call gender by read depth of chromosome X and Y. Compare the
#'   chi-squared values obtained to infer whether the male or female assumption
#'   fits read depth better.
#' @param x A character string of coverage file (<fileName>.cov) path.
#' @return A character string, \emph{Unknow}, \emph{Female} or \emph{Male}.
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom stats chisq.test
#' @author Zhan-Ni Chen
callGenderByCov <- function (x) {
    gr <- readCovFile(x)
    sry <- grep('^SRY\\(', gr$id, perl = TRUE)
    if (length(sry) > 0) {
        sry_depth <- gr[sry]$depth
        sry_depth <- mean(sry_depth, na.rm = TRUE)
        if (sry_depth >= 10) {
            return('Male')
        } else {
            return('Female')
        }
    }
    seqn <- as.character(as.vector(runValue(seqnames(gr))))
    sexn <- lapply(c('X', 'Y'), function(x) {seqn[grep(x, seqn)]})
    sexn <- unlist(sexn)
    if (length(sexn) == 0) return('Unknow')
    sexgr <- split(gr, seqnames(gr))
    sexgr <- sexgr[sexn]
    ratio <- sapply(sexgr, function(x) {mean(x$depth, na.rm = TRUE)})
    ratio_a <- mean(gr$depth, na.rm = TRUE )
    ratio_x <- 0
    ratio_y <- 0
    if (length(grep('X', names(ratio))) == 1) ratio_x <- ratio[grep('X', names(ratio))]
    if (length(grep('Y', names(ratio))) == 1) ratio_y <- ratio[grep('Y', names(ratio))]
    f_a_x <- chisq.test(c(ratio_a, ratio_x), p = c(2/4, 2/4))$p.value
    m_a_y <- chisq.test(c(ratio_a, ratio_y), p = c(2/3, 1/3))$p.value
    if (ratio_x == 0 & ratio_y == 0) return('Unknow')
    if (ratio_x == 0) {
        if (m_a_y > 5E-2) return('Male')
        return('Female')
    } else if (ratio_y == 0) {
        if (f_a_x > 5E-2) return('Female')
        return('Male')
    } else {
        if (f_a_x > m_a_y) return('Female')
        if (f_a_x < m_a_y) return('Male')
        return('Unknow')
    }
}

#' @title Get bam file parameters.
#' @param bamFiles Paths of BAM files.
#' @param covFiles Paths of coverage files.
#' @param path Path to write to.
#' @importFrom utils read.table
#' @importFrom utils write.table
getParameters <- function(bamFiles, covFiles, path = NULL) {
    write("Start to calculat metrics...", stdout())

    dat1 <- data.frame(
        id=sapply(bamFiles, getID),
        total.read=sapply(bamFiles, countTR)
    )
    dat2 <- data.frame(
        id=sapply(covFiles, getID),
        gender=sapply(covFiles, callGenderByCov)
    )
    dat <- merge(dat1, dat2, by = 'id')
    if (is.null(path)) path <- 'parameter.txt'
    write.table(dat, file = path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    write(paste0("Write to path: \n", normalizePath(path)), stdout())
}

ifMatchMetrics <- function(x, y, n = 0) {
    dat <- read.table(y, header = TRUE, sep = "\t", quote = "",
                      comment.char = "#", na.strings = "NA",
                      fill = TRUE, stringsAsFactors = FALSE)
    m <- match(x, as.character(as.vector(dat[, 'id'])))
    no.match <- x[which(is.na(m))]
    if (length(no.match) == length(x)) stop(paste0("Cannot find sample in file ", y))
    if (length(no.match) > n) write(paste0("Warning: cannot find ",
                                           paste(no.match, collapse = ', '), " in file ", y), stdout())
}

getID <- function(file, ptn = NULL) {
    sampleid  <- rev(unlist(strsplit(file, "/")))[1]
    if(is.null(ptn)) {
        sampleid  <- unlist(strsplit(sampleid, '\\.'))[1]
    }else{
        sampleid  <- unlist(strsplit(sampleid, ptn))[1]
    }
    sampleid
}
