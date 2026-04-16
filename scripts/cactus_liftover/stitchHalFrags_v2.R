suppressMessages(suppressWarnings(library(GenomicFeatures)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(dplyr)))

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("At least two argument must be supplied (original bed file) (hal lift bed file) (output bed file name) (minMatch) (max fraction).", call.=FALSE)
}


orig.bed <- args[1]
lift.bed <- args[2]
if(!file.exists(orig.bed)){
    stop(paste(orig.bed,"does not exist."))
}

if(!file.exists(lift.bed)){
    stop(paste(lift.bed,"does not exist."))
}

output.bed <- args[3]
minMatch <- as.numeric(args[4])
if (args[4] == "") {
  minMatch <- 0.5
}

max.frac <- as.numeric(args[5])
if (args[5] == "") {
  max.frac <- 1.5
}


orig.bed <- read.delim(orig.bed, sep ='\t', header = F)
orig.bed <- orig.bed %>% mutate(width = V3 - V2) %>% mutate(width_cutoff = width * minMatch, width_frac = width * max.frac)
lift.bed <- read.delim(lift.bed, sep = '\t', header = F)
filter.bed <- lift.bed %>% group_by(V4) %>% summarize(width = sum(V3 - V2))
mm <- match(filter.bed$V4, orig.bed$V4)
filter.bed$cutoff <- orig.bed$width_cutoff[mm]
filter.bed <- filter.bed %>%
              filter(width >= cutoff)
lift.bed <- lift.bed %>% filter(V4%in%filter.bed$V4)
print(paste("Filtered", nrow(orig.bed) - nrow(filter.bed), "interval(s) from original bed file with minMatch =", minMatch))

### consolidate fragments by chromosome and name
hal.join <- lift.bed %>% group_by(V1,V4) %>% summarise(start = min(V2), end = max(V3)) %>% mutate(width = end - start) 
hal.filter <- hal.join %>% group_by(V4) %>% summarise(width = sum(width))
mm <- match(hal.filter$V4, orig.bed$V4)
hal.filter$cutoff <- orig.bed$width_frac[mm]
hal.filter <- hal.filter %>%  filter(width <= cutoff) 
hal.join <- hal.join %>% filter(V4%in%hal.filter$V4)

print(paste("Filtered", nrow(filter.bed) - nrow(hal.join), "interval(s) from original bed file with and maxFrac =",max.frac))
print(paste("Total interval(s) filtered =", nrow(orig.bed) - nrow(hal.join)))

if(nrow(hal.join) > 0){
    colnames(hal.join) <- c('seqnames','name','start','end','width')
    hal.join <- hal.join[,c('seqnames','start','end','name')]
    hal.join <- GRanges(hal.join)    
    hal.join <- data.frame(hal.join)[,-4] ## remove width column
    hal.join <- data.frame(hal.join)[,-4] ## remove strand column
    write.table(hal.join, output.bed, sep = '\t', col.names = F, row.names = F, quote = F)
} else{
    write.table(hal.join, output.bed, sep = '\t', col.names = F, row.names = F, quote = F)
    print(paste("no regions recovered in",lift.bed))
}
