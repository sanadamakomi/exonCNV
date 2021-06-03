readFile <- function(file) {
    header <- readLines(file, 1)
    if (grepl('^#', header, perl = TRUE)) {
        indat <- read.table(file, header = TRUE, sep = "\t", quote = "", comment.char = "", na.strings = "NA",fill = TRUE, stringsAsFactors = FALSE)
    } else {
        indat <- read.table(file, header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",fill = TRUE, stringsAsFactors = FALSE)
    }
    indat
}

bedtoGRange <- function(file, genomeSeqinfo = NULL) {
    indat <- read.table(file, header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",
                        fill = TRUE, stringsAsFactors = FALSE)
    if (nrow(indat) < 3) stop(paste0('Error format: ', file))
    chr <- as.character(as.vector(indat[,1]))
    chr <- gsub('chr', '', chr)
    idx <- which(chr %in% c(as.character(seq(1, 22, 1)), 'X', 'Y'))
    if (! is.null(genomeSeqinfo)) {
        if(grepl('chr', seqnames(genomeSeqinfo)[1])) chr <- paste0('chr', chr)
    }
    if (nrow(indat) == 3) {
        gr <- GRanges(Rle(chr[idx]), IRanges(start = as.numeric(as.vector(indat[idx, 2])),
                                   end = as.numeric(as.vector(indat[idx, 3]))))
    } else {
        gr <- GRanges(Rle(chr[idx]), IRanges(start = as.numeric(as.vector(indat[idx, 2])),
                                             end = as.numeric(as.vector(indat[idx, 3]))),
                      id = indat[idx, 4])
    }
}

readCovFile <- function(file) {
    indat <- read.table(file, header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",
                        fill = TRUE, stringsAsFactors = FALSE)
    sort( GRanges(Rle(as.character(as.vector(indat[,1]))),
                  IRanges(start = as.numeric(as.vector(indat[,2])), end = as.numeric(as.vector(indat[,3]))),
                  depth = as.numeric(as.vector(indat[,5])),
                  id = as.character(as.vector(indat[,4]))))
}

writeCovFile <- function(gr, path) {
    df <- as.data.frame(gr)
    cat(paste(as.character(as.vector(df[, "seqnames"])),
              as.numeric(as.vector(df[, "start"])), as.numeric(as.vector(df[, "end"])),
              as.character(as.vector(df[, "id"])), as.numeric(as.vector(df[, "depth"])),
              sep = '\t'), sep = '\n', file = path)
}

outputVcf <- function(x, path) {
    if (! file.exists(x)) stop(paste0(x, 'no exists.'))
    id <- getID(x)
    dat <- read.table(x, header = TRUE, sep = "\t", quote = "",
                      comment.char = "#", na.strings = "NA",
                      fill = TRUE, stringsAsFactors = FALSE)
    dat$svlen <- as.numeric(as.vector(dat[,'end'])) - as.numeric(as.vector(dat[,'start'])) + 1
    header <- makeVcfHeader()
    vcf.matrix <- makeVcfMatrix(dat)
    cat(header, file=path, sep="\n")
    cat(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", id), collapse = "\t"),
        file = path, sep = "\n", append = TRUE)
    cat(apply(vcf.matrix, 1, paste, collapse = "\t"), file = path, sep = "\n", append = TRUE)

    write(paste0("Write to path: \n", normalizePath(path)), stdout())
}


makeVcfMatrix <- function(x) {
    CHROM <- as.character(as.vector(x[, "chr"]))
    POS <- as.numeric(as.vector(x[, "start"]))
    if ("id" %in% names(x)) {
        ID <- as.character(as.vector(x[, "id"]))
    } else {
        ID <- rep(".", nrow(x))
    }
    REF <- rep("N", nrow(x))
    ALT <- paste0("<", as.character(as.vector(x[, "svtype"])), ">")
    if ("qual" %in% names(x)) {
        QUAL <- as.character(as.vector(x[, "qual"]))
    } else {
        QUAL <- rep(".", nrow(x))
    }
    if ("filter" %in% names(x)) {
        FILTER <- as.character(as.vector(x[, "filter"]))
    } else {
        FILTER <- rep(".", nrow(x))
    }
    info.vec <- c('end', 'svlen', 'svtype', 'gene', 'trans', 'exon.num', 'exon.sv.id',
                  'exon.sv.num', 'cn', 'prob', 'depth','baseline', 'pool.count', 'past.count', 'exon.low.num')
    new.name <- c("END", "SVLEN", "SVTYPE", 'GENE', "TRANS", "EXONC", "EXONSV",
                  "EXONSVC", "CN", "PROB", "DEPTH", "BASELINE", "POOLC", "PASTC", 'EXONLOWC')
    info.df <- x[, info.vec, drop = FALSE]
    colnames(info.df) <- new.name
    INFO <- makeVcfInfo(info.df)
    x$gt <- ifelse(as.numeric(as.vector(x[,'cn'])) <= 0.5, "1/1", "0/1")
    if(! "gq" %in% names(x)) {
        x$gq <- rep('.', nrow(x))
    }
    FORMAT <- makeVcfFormat(x[, c("gt", "gq", "cn"), drop=FALSE])
    cbind(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, rep("GT:GQ:CN", nrow(x)), FORMAT)
}

makeVcfHeader <- function() {
    fileformat <- "##fileformat=VCFv4.2"
    fileDate <- paste0("##fileDate=", format(Sys.time(), "%Y%m%d"))
    alt <- c('##ALT=<ID=DEL,Description="Deletion">',
             '##ALT=<ID=DUP,Description="Duplication">')
    info <- c('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record.">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles.">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">',
    '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol.">',
    '##INFO=<ID=TRANS,Number=1,Type=String,Description="Transcript name.">',
    '##INFO=<ID=EXONC,Number=1,Type=Integer,Description="Exon counts of provided transcript.">',
    '##INFO=<ID=EXONSV,Number=1,Type=String,Description="Variant exon index according to provided transcript.">',
    '##INFO=<ID=EXONSVC,Number=1,Type=Integer,Description="Variant exon counts.">',
    '##INFO=<ID=CN,Number=1,Type=Float,Description="Copy number.">',
    '##INFO=<ID=PROB,Number=1,Type=String,Description="Probability of cnv.">',
    '##INFO=<ID=DEPTH,Number=1,Type=Integer,Description="Weighted mean of normalized read depths across all this gene bins.">',
    '##INFO=<ID=BASELINE,Number=1,Type=Integer,Description="Weighted mean of normalized read depths across control samples.">',
    '##INFO=<ID=POOLC,Number=1,Type=Integer,Description="Counts of samples sharing this CNV in the same pool.">',
    '##INFO=<ID=PASTC,Number=1,Type=String,Description="Counts of samples sharing this CNV in database.">',
    '##INFO=<ID=EXONLOWC,Number=1,Type=Integer,Description="Counts of samples low-depth exon.">')
    format <- c('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype.">',
                '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality.">',
                '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events.">')
    c(fileformat, fileDate, alt, info, format)
}

makeVcfInfo <- function(x) {
    n_col <- ncol(x)
    n_row <- nrow(x)
    col_n <- colnames(x)
    x <- apply(x, 2, function(y) { gsub(" ", "", as.character(as.vector(y))) })
    x <- matrix(data=x, ncol=n_col)
    sapply(1:nrow(x), function(i) {
        paste(paste0(toupper(col_n), "=", as.character(as.vector(x[i,,drop=FALSE]))), collapse=";")
    })
}

makeVcfFormat <- function(x) {
    n_col <- ncol(x)
    n_row <- nrow(x)
    col_n <- colnames(x)
    x <- apply(x, 2, function(y) { gsub(" ", "", as.character(as.vector(y))) })
    x <- matrix(data=x, ncol=n_col)
    sapply(1:nrow(x), function(i) {
        paste(as.character(as.vector(x[i,,drop=FALSE])), collapse=":")
    })
}


