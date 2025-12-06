library(cn.mops)
library(magrittr)


args <- commandArgs(TRUE)
NGS_DATA_TYPE = args[1]
chrom_style = as.character(args[2])
window_size = as.integer(args[3])
exome_bed = args[4]
threads = as.integer(args[5])
cnvs = args[6]
cnvr = args[7]
bam_files = tail(args, -7)


if (NGS_DATA_TYPE == "WES"){

    segments <- read.table(exome_bed, sep="\t", as.is = TRUE)
    gr <- GRanges(segments[,1], IRanges(segments[,2], segments[,3]))
    bam_data_ranges <- getSegmentReadCountsFromBAM(bam_files, GR = gr, parallel = threads)
    results <- exomecn.mops(bam_data_ranges) %>% calcIntegerCopyNumbers()
    
} else {

    if (chrom_style == "l"){
        seq_names <- paste("chr", c(as.character(seq(22)), "X", "Y"), sep = "")
    } else {
        seq_names <- c(as.character(seq(22)), "X", "Y")
    }
    
    bam_data_ranges <- getReadCountsFromBAM(bam_files, refSeqNames = seq_names, WL = window_size, parallel = threads)
    results <- cn.mops(bam_data_ranges) %>% calcIntegerCopyNumbers()
}


CNVs <- cnvs(results)
CNVRegions <- cnvr(results)
write.table(CNVs, file = cnvs, sep = "\t", quote = F)
write.table(CNVRegions, file = cnvr, sep = "\t", quote = F)
