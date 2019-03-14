#' @title Plot the unifomity of coverage among samples.
#' @description Compute correlation between columns of a numeric matrix or data
#'   frame and plot the distribution of correlations.
#' @param geneCNVFile Path of gene cnv result file.
#' @param mergedCovFile Path of merged coverage file.
#' @param parameterFile Path of metrics file.
#' @param gene A character string or vector of gene symbol.
#' @param path Path to write to.
#' @param draw An integer of gene counts to plot.
#' @export
#' @import ggplot2
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom stringr str_match
#' @importFrom reshape2 melt
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @author Zhan-Ni Chen
plotExonCNV <- function(geneCNVFile, mergedCovFile, parameterFile, gene = NULL, path = NULL, draw = 10){
    if (! file.exists(geneCNVFile)) stop(paste0(geneCNVFile, 'no exists.'))
    if (! file.exists(mergedCovFile)) stop(paste0(mergedCovFile, 'no exists.'))
    if (! file.exists(parameterFile)) stop(paste0(parameterFile, 'no exists.'))
    if (is.null(path)) path <- paste0(normalizePath('.'), '/', sample.id, '.exon.cnv.pdf')

    write("Start to plot exon CNV...", stdout())

    path <- normalizePath(path)
    sample.id <- getID(geneCNVFile)
    gene.dat <- read.table(geneCNVFile, header = TRUE, sep = "\t", quote = "",
                           comment.char = "#", na.strings = "NA",
                           fill = TRUE, stringsAsFactors = FALSE)
    if (is.null(gene)) {
        n.draw <- nrow(gene.dat)
        if (n.draw >= draw) n.draw <- draw
        gene <- unique(as.character(as.vector(gene.dat[1:n.draw, 'gene'])))
    } else {
        gene <- intersect(gene, as.character(as.vector(gene.dat[, 'gene'])))
        if (length(gene) == 0) stop('Input gene no exists in gene cnv result file.')
    }

    gene.dat <- gene.dat[which(gene.dat[,'gene'] %in% gene),, drop = FALSE]
    cov.dat <- read.table(mergedCovFile, header = TRUE, sep = "\t", quote = "",
                          comment.char = "#", na.strings = "NA",
                          fill = TRUE, stringsAsFactors = FALSE)
    infos <- str_match(cov.dat[,'id'], "^(.+)\\((.+?)\\)_([0-9]+)_[0-9]+$")
    cov.dat <- cov.dat[which(as.character(as.vector(infos[,2])) %in% gene),,drop = FALSE]
    infos <- str_match(cov.dat[,'id'], "^(.+)\\((.+?)\\)_([0-9]+)_[0-9]+$")
    para.dat <- read.table(parameterFile, header = TRUE, sep = "\t", quote = "",
                          comment.char = "#", na.strings = "NA",
                          fill = TRUE, stringsAsFactors = FALSE)

    ## plot by gene
    pdf(file = path, width = 12, height = 4.5)

    dogene <- sapply(1:nrow(gene.dat), function(i) {
        cnv.info <- gene.dat[i,,drop = FALSE]

        write(paste0('Perform gene: ', cnv.info[, 'gene'], ", type: ", cnv.info[,'svtype']), stdout())

        sample.gender <- para.dat[which(para.dat[, 'id'] %in% sample.id), 'gender']
        sub.dat <- cov.dat[which(infos[, 2] %in% cnv.info[,'gene']),, drop = FALSE]
        baseline.id <-setdiff(colnames(sub.dat), c('chr', 'start', 'end', 'id', sample.id))
        if (sub.dat[1, 'chr'] %in% c('X, Y')) {
            baseline.id <- intersect(baseline.id,
                                     para.dat[which(para.dat[, 'gender'] %in% sample.gender), 'id'])
            sub.dat <- sub.dat[, c('chr', 'start', 'end', 'id', sample.id, baseline.id)]
        }
        exon.infos <- str_match(sub.dat[,'id'], "^(.+)\\((.+?)\\)_([0-9]+)_[0-9]+$")
        sub.dat$exon <- exon.infos[,4]
        gene.start <- sub.dat[1,'start']
        gene.end <- sub.dat[1,'end']
        if(nrow(sub.dat) > 1 ) {
            gene.start <- min(as.numeric(as.vector(sub.dat[,'start'])), na.rm = TRUE)
            gene.end <- max(as.numeric(as.vector(sub.dat[,'end'])), na.rm = TRUE)
        }

        ## plot
        ## labels
        positive.exon <- revContinuousInteger(cnv.info[,'exon.sv.id'])
        sv.type <- ifelse(cnv.info[,'svtype'] == 'DUP', 'duplicates', 'deletes')
        sv.type2 <- ifelse(cnv.info[,'svtype'] == 'DUP', 'Duplicated', 'Deleted')
        x.bar.label <- as.numeric(as.vector(sub.dat$exon))
        x.bar.label[which(! x.bar.label %in% positive.exon)] <- " "
        title <- paste0(sample.id, "(", sample.gender, ")'s ",
                                 cnv.info[,'gene'], " gene ", sv.type,
                                 " in exon ", cnv.info[,'exon.sv.id'], '.')
        sub.title <- paste0(sv.type2, ' exon region is chr', cnv.info[1, 'chr'], ":",
                            cnv.info[1, 'start'], '-', cnv.info[1, 'end'],
                            " according to tanscript ", exon.infos[1,3], ".")
        x.axis.label <- paste0(cnv.info[,'gene'], " Exon Rank(chr",
                               sub.dat[1,'chr'], ":", gene.start, "-", gene.end, ")")
        y.axis.label <- "Depth of Exon"

        ## boxplot
        plot.df <- melt(sub.dat[, c('exon', sample.id, baseline.id)], id.vars = "exon")
        names(plot.df) <- c('exon', 'sampleID', 'depth')
        plot.df$sampleID <- as.character(as.vector(plot.df$sampleID))
        plot.df[which(! plot.df[,'sampleID'] %in% sample.id), 'sampleID'] <- rep("control samples",
             length(which(!plot.df[,'sampleID'] %in% sample.id)))
        plot.df$sampleID <- as.factor(plot.df$sampleID)
        plot.df$exon <- as.integer(plot.df$exon)

        y.time <- nchar(as.character(ceiling(max(plot.df$depth))))
        y.lim <- ceiling(max(plot.df$depth)/((y.time-1)*10))*((y.time-1)*10)
        if (is.na(y.lim)) y.lim <- 100

        exon <- depth <- sampleID <- sampleLable <- NULL
        xmax <- xmin <- ymax <- ymin <- NULL
        p.box <- ggplot(plot.df, aes(x = as.factor(exon), y = depth)) +
            geom_boxplot(notch = FALSE, outlier.colour = NA) +
            geom_jitter(aes_string(colour = "sampleID"), shape = 16, position = position_jitter(0.2)) +
            scale_color_manual(breaks = c(sample.id, "control samples"), values=c("gray", "red")) +
            scale_y_continuous(breaks = seq(0, y.lim, 20), labels=seq(0, y.lim, 20)) +
            labs(title = title,
                  subtitle  = sub.title,
                  x = x.axis.label,
                  y = y.axis.label)+
            theme(axis.text.x = element_text(size = 6, angle = 0), legend.position = "bottom")

        rect.df <- data.frame(xmin = positive.exon - 0.5,
                              xmax = positive.exon + 0.5,
                              ymin = rep(-Inf,length(positive.exon)),
                              ymax = rep(Inf, length(positive.exon)))
        p.box <- p.box + geom_rect(data = rect.df,
                   aes(xmin = xmin, xmax = xmax, ymin = ymin,
                   ymax = ymax, fill = "red"), alpha = 0.2, inherit.aes = FALSE)+
                   guides(fill=FALSE)

        # bar box
        plot.df2 <- data.frame(exon = as.integer(sub.dat[,'exon']),
                               sample = sub.dat[, sample.id],
                               control = apply(sub.dat[, baseline.id], 1, mean, na.rm = TRUE))
        names(plot.df2) <- c('exon', sample.id, 'control')
        plot.df2 <- melt(plot.df2, id.vars = "exon")
        names(plot.df2) <- c('exon', 'sampleID', 'depth')
        plot.df2$sampleLable <- rep(' ', nrow(plot.df2))
        plot.df2[which(plot.df2[,'sampleID'] %in% sample.id), 'sampleLable'] <- x.bar.label

        p.bar <- ggplot(plot.df2, aes(x = exon, y = depth, fill = sampleID))+
            geom_bar(stat = "identity", alpha = 2/3, aes(fill = factor(sampleID)), position = "dodge") +
            geom_text(aes(label = sampleLable), vjust = - 0.5, hjust = 1.2, color = "black", size = 2)+
            theme_classic()+
            guides(fill = guide_legend(title = NULL))+
            scale_fill_discrete(guide = guide_legend(reverse = TRUE))+
            labs(title = title,
                 subtitle  = sub.title,
                 x = x.axis.label,
                 y = y.axis.label)+
            theme(legend.title = element_text(size = rel(2.5)), legend.position = 'bottom',
                  axis.text.x = element_text(size = 6, angle = 0))

        print(p.bar)
        print(p.box)
    })
    dev.off()

    write(paste0("Write to path: \n", normalizePath(path)), stdout())
}
