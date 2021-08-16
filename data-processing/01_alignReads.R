##############
## align fastqs to ref ##
## get summaryStats ##
##############

library(Rsubread)

fq.dir <- "trimmed" #source directory
fq.files <- list.files(fq.dir, pattern="_TRIMMED.fastq.gz")
dest.dir<-"BAMoutput/"
dir.create(dest.dir)

for(f in fq.files){
   cat("File: ", f,"\n")
   f.out.name <- paste0(dest.dir, gsub(".fastq.gz", ".bam", f))
   try(align(index="indices/genome.fa",  #points to location
         readfile1=file.path(fq.dir, f),
         nBestLocations=10, 
         unique=F, 
         indels=0,
         nsubreads=50,
         TH1=2,
         maxMismatches=4,
         type=1, 
         annot.ext="miRBase/hsa.gff3", #points to location
         isGTF=T,
         GTF.featureType="miRNA",
         output_file=f.out.name))
}

