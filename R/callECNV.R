#' @title Fit data
#' @param mergedCovFile Path of merged coverage file.
#' @param parameterFile Path of metrics file.
#' @param path Path to write to.
#' @param lowdepth A numeric value, regions that avaerage depth less than this
#'   value will replaced by NA.
#' @export
#' @author Zhan-Ni Chen
performFitPoisson <- function(mergedCovFile, parameterFile, path = NULL, lowdepth = 10) {
    if (! file.exists(mergedCovFile)) stop(paste0(mergedCovFile, 'no exists.'))
    if (! file.exists(parameterFile)) stop(paste0(parameterFile, 'no exists.'))
    if (is.null(path)) path <- paste0(normalizePath('.'), 'probability.txt')

    write(paste0("Start to fit data...\nCutoff: lowdepth=", lowdepth), stdout())

    covDat <- read.table(mergedCovFile, header = TRUE, sep = "\t", quote = "",
                         comment.char = "#", na.strings = "NA",
                         fill = TRUE, stringsAsFactors = FALSE)
    paraDat <- read.table(parameterFile, header = TRUE, sep = "\t", quote = "",
                          comment.char = "#", na.strings = "NA",
                          fill = TRUE, stringsAsFactors = FALSE)

    # 1. get ids in parameter file
    x.ids <- setdiff(colnames(covDat), c('chr', 'start', 'end', 'id'))
    x.nopara <- x.ids[which(! x.ids %in% paraDat$id)]
    if (length(x.nopara) > 0) stop(paste0('Error: ', paste(x.nopara, collapse = ', '), ' no exists in ', parameterFile))
    paraDat <- paraDat[which(paraDat$id %in% x.ids), ]
    x.ids.m <- paraDat[which(paraDat[, 'gender'] == 'Male'), 'id']
    x.ids.f <- paraDat[which(paraDat[, 'gender'] == 'Female'), 'id']

    # 2. normalize coverage and
    ratio <- min(paraDat$total.read) / paraDat$total.read
    names(ratio) <- paraDat$id
    for(i in paraDat$id) {
        covDat[,i] <- as.integer(covDat[,i] * ratio[i])
    }

    # 3. fit poisson
    ppmatrix <- matrix(data = 1, nrow = nrow(covDat), ncol = length(x.ids))
    colnames(ppmatrix) <- x.ids
    for (i in 1:nrow(covDat)){
        if (covDat[i, 'chr'] %in% c('X', 'Y')) {
            if (length(which(covDat[i, x.ids.m] > lowdepth)) > length(x.ids.m) * 1/2) {
                ppmatrix[i, x.ids.m] <- doPois(covDat[i, x.ids.m])
            }
            if (length(which(covDat[i, x.ids.f] > lowdepth)) > length(x.ids.f) * 1/2) {
                ppmatrix[i, x.ids.f] <- doPois(covDat[i, x.ids.f])
            }
        } else {
            if (length(which(covDat[i, x.ids] > lowdepth)) > length(x.ids) * 1/2) {
                ppmatrix[i, x.ids] <- doPois(covDat[i, x.ids])
            }
        }
    }
    ppmatrix <- cbind(covDat[,c('chr', 'start', 'end', 'id')], ppmatrix)
    write.table(ppmatrix, file = path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    write(paste0("Write to path: \n", normalizePath(path)), stdout())
}

#' @title Call exon CNVs.
#' @param probFile Path of probability file.
#' @param mergedCovFile Path of merged coverage file.
#' @param path Directory path to write to.
#' @param sample.id Sample ids to call CNV. By default it will run in all
#'   samples.
#' @param cutoff A list of cutoff parameters to filter exon CNV, 'prob' is
#'   probability, 'pool.count' is the least counts of sample of which gene has
#'   CNV. 'baseline' is the least baseline depth to call CNV. 'lowdepth' is
#'   least depth in all samples to call CNV.
#' @export
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom stringr str_match
#' @importFrom stats dpois
#' @importFrom stats nlm
#' @author Zhan-Ni Chen
callExonCNV <- function(probFile, mergedCovFile, path = NULL, sample.id = NULL, cutoff = list(prob = 1E-4, pool.count = 1, baseline = 50, lowdepth = 10)) {
    if (! file.exists(probFile)) stop(paste0(probFile, 'no exists.'))
    if (! file.exists(mergedCovFile)) stop(paste0(mergedCovFile, 'no exists.'))
    if (is.null(path)) path <- '.'
    if (! dir.exists(path)) dir.create(path)

    write("Start to call exon CNV...", stdout())
    write(paste0("Cutoff: ", paste(paste0(names(cutoff), "=", cutoff[names(cutoff)]), collapse = ', ')), stdout())

    path <- normalizePath(path)
    dat <- read.table(probFile, header = TRUE, sep = "\t", quote = "",
                      comment.char = "#", na.strings = "NA",
                      fill = TRUE, stringsAsFactors = FALSE)
    depth.dat <- read.table(mergedCovFile, header = TRUE, sep = "\t", quote = "",
                            comment.char = "#", na.strings = "NA",
                        fill = TRUE, stringsAsFactors = FALSE)
    all.id <- setdiff(colnames(dat), c('chr', 'start', 'end', 'id'))
    if (is.null(sample.id)) {
        sample.id <- all.id
    } else {
        sample.id <- intersect(sample.id, all.id)
        if (length(sample.id) == 0) stop(paste0(paste(sample.id, collapse = ' '), 'no exists.'))
    }
    counts <- apply(dat[, all.id], 1,
                   function(x) {length(which(x < cutoff$prob))})
    infos <- str_match(dat[,'id'], "^(.+)\\((.+?)\\)_([0-9]+)_[0-9]+$")
    gene.table <- table(infos[,2])
    exon.num <- as.data.frame(gene.table[infos[,2]])
    dosample <- sapply(sample.id, function(i){

        write(paste0('Perform ', i), stdout())

        raw.path <- paste0(path, '/', i, '.raw.cnv.txt')
        positive.path <- paste0(path, '/', i, '.positive.cnv.txt')
        filter.path <- paste0(path, '/', i, '.gene.cnv.txt')

        baseline <- apply(depth.dat[, setdiff(colnames(depth.dat),
                          c('chr', 'start', 'end', 'id', i)), drop = FALSE],
                          1, median, na.rm = TRUE)
        log2 <- calculateLog2ratio(x = depth.dat[,i], baseline = baseline, badDepth = cutoff$lowdepth)
        cn <- round(2 * 2 ^ log2, 1)
        svType <- rep('.', nrow(dat))
        dup.idx <- which(dat[,i] < cutoff$prob & log2 > 0)
        del.idx <- which(dat[,i] < cutoff$prob & log2 < 0)
        svType[dup.idx] <- 'DUP'
        svType[del.idx] <- 'DEL'
        df <- data.frame(
            chr=dat[,'chr'],
            start=dat[,'start'],
            end=dat[,'end'],
            gene=infos[,2],
            transcript=infos[,3],
            exon.num=exon.num$Freq,
            exon=infos[,4],
            exon.len=dat[,'end']-dat[,'start'],
            svtype=svType,
            cn=cn,
            prob=format(dat[,i], digits = 3, scientific = TRUE),
            depth=round(depth.dat[,i], 0),
            baseline=round(baseline, 0),
            pool.count=counts,
            past.count=rep(NA, nrow(dat))
        )
        write.table(df, file = raw.path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        write.table(df[which(df$svtype %in% c('DEL', 'DUP')),, drop = FALSE],
                    file = positive.path, row.names = FALSE,
                    col.names = TRUE, quote = FALSE, sep = "\t")

        # filter
        gene.df <- filterExonCNV(df, cutoff = cutoff)
        count1 <- as.numeric(as.vector(gene.df[, ncol(gene.df)]))
        bad.index <- which(count1 / gene.df$exon.num >= 0.2)
        gene.df.good <- gene.df[setdiff(1:nrow(gene.df), bad.index),,drop = FALSE]
        gene.df.bad <- gene.df[bad.index,,drop = FALSE]
        rank.coefficient <- -log10(gene.df.good$prob) + gene.df.good$exon.sv.num * 2
        gene.df.good <- gene.df.good[order(rank.coefficient, decreasing = TRUE),,drop = FALSE]
        write.table(rbind(gene.df.good, gene.df.bad),
                    file = filter.path, row.names = FALSE,
                    col.names = TRUE, quote = FALSE, sep = "\t")

        write(paste0('Write to path:\n', normalizePath(filter.path)), stdout())
    })
}

filterExonCNV <- function(x, cutoff = list(prob = 1E-4, pool.count = 3, baseline = 50, lowdepth = 10)) {

    # filter cutoff
    filter.x <- x[which(as.numeric(as.vector(x[,'prob'])) < cutoff$prob &
      as.numeric(as.vector(x[,'baseline'])) >= cutoff$baseline ),, drop = FALSE]

    # merge exons sharing the same gene symbol
    filter.gene <- unique(as.character(as.vector(filter.x[,'gene'])))
    result <- lapply(filter.gene, function(gene) {
        idx <- which(as.character(as.vector(x[,'gene'])) %in% gene)
        gene.whole.dat <- x[idx,,drop = FALSE]
        count1 <- length(which(as.numeric(as.vector(gene.whole.dat[,'prob'])) == 1))
        gene.dat <- gene.whole.dat[which(gene.whole.dat[,'svtype'] %in% c('DEL', 'DUP')),,drop = FALSE]

        # call by sv type
        svtype.class <- unique(as.character(as.vector(gene.dat[,'svtype'])))
        out.dat <- list()

        for (sv in svtype.class) {
            sub.dat <- gene.dat[which(as.character(as.vector(gene.dat[,'svtype'])) %in% sv),,drop = FALSE]
            exon.sv.id <- continuousInteger(sort(as.numeric(as.vector(sub.dat[, 'exon']))))
            past.count <- na.omit(as.numeric(as.vector(sub.dat[, 'past.count'])))
            if (length(past.count) > 0) {
                past.count <- min(past.count)
            } else {
                past.count <- NA
            }
            out.dat[[sv]] <- data.frame(chr=sub.dat[1, 'chr'],
                start=min(as.numeric(as.vector(sub.dat[, 'start']))),
                end=max(as.numeric(as.vector(sub.dat[, 'end']))),
                gene=gene,
                trans=sub.dat[1, 'transcript'],
                exon.num=sub.dat[1, 'exon.num'],
                exon.sv.id=paste(exon.sv.id, collapse = ','),
                exon.sv.num=nrow(sub.dat),
                svtype=sv,
                cn=median(as.numeric(as.vector(sub.dat[, 'cn'])), na.rm = TRUE),
                prob=min(as.numeric(as.vector(sub.dat[, 'prob'])), na.rm = TRUE),
                depth=median(as.numeric(as.vector(sub.dat[, 'depth'])), na.rm = TRUE),
                baseline=median(as.numeric(as.vector(sub.dat[, 'baseline'])), na.rm = TRUE),
                pool.count=min(as.numeric(as.vector(sub.dat[, 'pool.count'])), na.rm = TRUE),
                past.count=past.count,
                exon.low.num=count1)
        }

        out.dat <- do.call('rbind', out.dat)
        rownames(out.dat) <- NULL
        out.dat
    })

    result <- do.call('rbind', result)
    result[which(as.numeric(as.vector(result[,'pool.count'])) <= cutoff$pool.count),,drop = FALSE]
}

doPois <- function(depth, ...){
    depth <- as.integer(depth)
    negpois.LL <- function(par, ...){
        mu <- par[1]
        LL <- sum(dpois(depth, lambda = mu, log=T))
        -LL
    }
    out.pois <- nlm(negpois.LL, 2)
    dpois(depth, lambda = out.pois$estimate[1])
}

continuousInteger <- function(x) {
    t <- as.data.frame(reduce(IRanges(x,x)))[,1:2]
    apply(t, 1, function(x) paste(unique(x), collapse = "~"))
}

revContinuousInteger <- function(x) {
    s <- unlist(strsplit(x, ','))
    out <- c()
    for (i in s){
        if (grepl('~', i)) {
            s1 <- as.numeric(unlist(strsplit(i, '~')))
            out <- c(out, seq(s1[1], s1[2], 1))
        } else{
            out <- c(out, as.numeric(i))

        }

    }
    sort(out)
}

#' @title Compute the log2 ratio between sample depth and baseline.
#' @description Input vectors of sample depth and baseline to compute the log2
#'   ratio.
#' @param x A numeric vector of sample depth.
#' @param baseline A numeric vector contains baseline depth of region with the
#'   same order in sample depth vector.
#' @param badDepth A numeric value; a baseline depth low than badDepth will be
#'   set to \code{NA}.
#' @return A numeric vector of log2 ratio.
#' @author Zhan-Ni Chen
calculateLog2ratio <- function(x, baseline, badDepth) {
    x <- as.numeric(as.vector(x))
    baseline <- as.numeric(as.vector(baseline))
    x[which(x == 0)] <- 0.01
    baseline[which(baseline < badDepth)] <- NA
    log2ratio <- log2(x / baseline)
    log2ratio[which(is.na(log2ratio))] <- NA
    log2ratio
}
