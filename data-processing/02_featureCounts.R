library(Rsubread)

bam.files <- list.files("BAMoutput")
bam.files <- bam.files[grep("Control", bam.files, invert=T)]
bam.files <- bam.files[grep("vcf|summary", bam.files, invert=T)]
bam.files <- bam.files[grep("TRIMMED", bam.files)]
gff <- "miRBase/hsa.v2.gff3" 

bam.files <- paste0("BAMoutput/",bam.files)

fCounts <- featureCounts(files=bam.files,
                         isGTFAnnotationFile=T,
                         GTF.attrType='Name',
                         annot.ext=gff,
                         GTF.featureType="miRNA",
                         allowMultiOverlap=T,
                         countMultiMappingReads=T
)


