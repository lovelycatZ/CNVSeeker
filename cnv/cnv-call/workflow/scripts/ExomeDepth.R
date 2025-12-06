library(ExomeDepth)
library(GenomicRanges)


args <- commandArgs(TRUE)

print(paste0("参考基因组: ", args[1]))
print(paste0("bed文件: ", args[2]))
print(paste0("结果目录: ", args[3]))
print(paste0("样本: ", tail(args, -3)))


ref = args[1]
bed = args[2]
outdir = args[3]
bam_files = tail(args, -3)


bed <- read.table(bed, header = FALSE, stringsAsFactors = FALSE)
if(ncol(bed) == 3) {
  colnames(bed) <- c("chromosome", "start", "end")
} else {
  colnames(bed) <- c("chromosome", "start", "end", "names")
}


# Get counts
bam_counts_gr <- GRanges(getBamCounts(bed.frame = bed,
    bam.files = bam_files,
    include.chr = FALSE,
    referenceFasta = ref
))
my.count.dafr <- data.frame(bam_counts_gr)
my.count.dafr$names <- c(1:nrow(bed))
colnames(my.count.dafr) <- sub("[.]","-",sub("[.]","-",colnames(my.count.dafr)))

ExomeCount.mat <- as.matrix(my.count.dafr[,grep(names(my.count.dafr), pattern = '*.bam')]) 



# 循环运行
bam_names <- basename(bam_files)
nsamples <- ncol(ExomeCount.mat) 
for (i in 1:nsamples) 
{ 
    my.choice <- select.reference.set (
        test.counts = ExomeCount.mat[,i],
        reference.counts = ExomeCount.mat[,-i], 
        bin.length = (my.count.dafr$end - my.count.dafr$start)/1000, 
        n.bins.reduced = 10000)

    my.reference.selected <- apply(
        X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE], 
        MAR = 1, 
        FUN = sum) 

    all.exons <- new(
        'ExomeDepth', 
        test = ExomeCount.mat[,i], 
        reference = my.reference.selected, 
        formula = 'cbind(test, reference) ~ 1') 
    
    all.exons <- CallCNVs(
        x = all.exons, 
        transition.probability = 10^-4, 
        chromosome = my.count.dafr$seqnames, 
        start = my.count.dafr$start, 
        end = my.count.dafr$end,
        name = my.count.dafr$names
    )
    
    output.file <-  paste(outdir, paste(sub("[.].*", "", bam_names[i]), "tsv", sep="."), sep="")
    write.table(all.exons@CNV.calls,
            output.file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
}
