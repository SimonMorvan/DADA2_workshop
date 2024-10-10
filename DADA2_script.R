#### Libraries ####
library(dada2); packageVersion("dada2")
library(ggplot2)
library(reshape2)
library(tictoc)

#### Get your raw material ####
tic()
### Locate where the reads are and set your directory there

# On cluster
path <- "~/projects/rrg-yergeaue/projects/DADA2_Workshop/data/raw_reads/subset_2000/"
setwd(path)

# On personal computer
path <- "data/raw_reads/"


list.files(path)


### Store location of R1 fastq files in the Forward_reads
### Store location of R2 fastq files in the Reverse_reads

Forward_reads <- sort(list.files(path, pattern="_R1.fastq")) #149 
Reverse_reads <- sort(list.files(path, pattern="_R2.fastq")) #149
### Some samples do not have both the R1 and R2 files


### Shorten sample names 
sample.names <- sapply(strsplit(Forward_reads, "---"), function(x) {
  sub("_.*$", "", x[2])  # Remove everything after (and including) the underscore
})
sample.names 


### Now indicate the path of the files that we will analyse
Forward_reads <- file.path(path, Forward_reads)
Reverse_reads <- file.path(path, Reverse_reads)


#### Quality profile ####

plotQualityProfile(Forward_reads[1])
plotQualityProfile(Reverse_reads[1])
#
#plotQualityProfile(Forward_reads, aggregate = TRUE)
#plotQualityProfile(Reverse_reads, aggregate = TRUE)


#### Filter and trim ####

# First create a new path for the filtered reads we are going to generate
# 
# On cluster
filt_path <- file.path("~/projects/rrg-yergeaue/projects/DADA2_Workshop/data/filtered_reads") 

# On local computer
filt_path <- file.path("data/filtered_reads") 
 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

### Name the location of the filtered files with the sample name


# Amplified region : 16S V4-V5
# Primers : Forward_515F-Y: GTGYCAGCMGCCGCGGTAA 19 bases
#           Reverse_926R: CCGYCAATTYMTTTRAGTTT  20 bases
# Amplicon length ~405


out <- filterAndTrim(Forward_reads, filtFs, Reverse_reads, filtRs, 
                     truncQ=6,
                     truncLen = c(280,250),
                     trimLeft=c(19,20),
                     maxEE=c(2,2))
                     #multithread=TRUE)

# Check output
out

#Reads are filtered as pairs, so the number of R1 reads that pass the filter is the same as the number of R2 reads.

plotQualityProfile(filtFs[1])
plotQualityProfile(filtRs[1])
#


#### Learn error rates ####
errF <- learnErrors(filtFs)
errR <- learnErrors(filtRs)
#multithread=TRUE

# See what errors look like
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


#### Dereplication ####
derepFs <- derepFastq(filtFs)
names(derepFs) <- sample.names
derepFs[1]

derepRs <- derepFastq(filtRs)
names(derepRs) <- sample.names
derepRs[1]


##### Dada2 ASVs ####
dadaFs <- dada(derepFs, 
               err = errF, 
               #multithread=TRUE,
               pool=TRUE)


dadaRs <- dada(derepRs, 
               err=errR,
               #multithread=TRUE,
               pool=TRUE)

dadaFs[[1]]
dadaRs[[1]]

#save(dadaRs, file="data/dadaRs.rdata")
#save(dadaFs, file="data/dadaFs.rdata")


##### Merging R1 and R2s #####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      minOverlap = 12, 
                      maxMismatch = 0)

max(mergers[[1]]$nmatch) # Largest overlap 
min(mergers[[1]]$nmatch) # Smallest overlap

#### Abundance table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
seqtab[,1]
seqtab_df <- as.data.frame(seqtab)

hist(nchar(getSequences(seqtab)),xlab="Size", ylab="Frequency", main = "ASVs length",ylim=c(0,400)) 


#### Remove chimeras #####
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "pooled", 
                                    #multithread = TRUE,
                                    verbose = TRUE) 

#### Presence / absence of ASVs #####
seqtab.nochim.bin <- ifelse(seqtab.nochim>0,1,0) 
seqtab.nochim.bin.df <- as.data.frame(t(seqtab.nochim.bin))

#### Inspect effect of all filtering steps #####
getN <- function(x) sum(getUniques(x))
track <- data.frame(Input=as.numeric(out[,1]), # input
                    Filtered=as.numeric(out[,2]), # filtered
                    "Filt//In"=as.numeric(round(((out[,2]/out[,1])*100),2)),# % (Filtered / Input)
                    Merge = as.numeric(sapply(mergers, getN)), # Merged 
                    "Mer//In"=as.numeric(round(((sapply(mergers, getN)/out[,1])*100),2)),# % (Merged / Input)
                    Nonchim = as.numeric(rowSums(seqtab.nochim)),# Non-chimeric                       
                    "Nonchim//In"=as.numeric(round(((rowSums(seqtab.nochim)/out[,1])*100),2)),# % (Non-chimeric / Input)
                    ASV = as.numeric(rowSums(seqtab.nochim.bin))) # Number of ASVs per sample 
rownames(track) <- sample.names # Row names
head(track)

#### Plot to see all filtering steps #####
gtrack<- track[,c(1,2,4,6)]
gtrack$ID <- rownames(gtrack)

lgtrack <- melt(gtrack, id.vars="ID")


bar_track <- ggplot(lgtrack ,aes(x=ID, y=as.numeric(value), fill=variable)) +
  geom_bar(stat="identity", position = "identity") + 
  theme_classic() + # Theme
  theme(axis.ticks.length=unit(0.3,"cm"), # Ticks size
        axis.text.x = element_text(angle=45,hjust=1),# Changes the x labels orientation 
        axis.title = element_text(size=18),
        axis.text.y = element_text(size=15),
        title = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.background = element_rect(color="black"))+ 
  scale_x_discrete(name ="Sample ID", limits=rownames(track))+ # Changes x-axis title & sorts the x label names
  scale_y_continuous(name="Abundance",expand = expansion(mult = c(0, 0.05)))+ #Changes y-axis titl
  scale_fill_manual(name="Reads kept at each step",values=c("#003049", "#d62828", "#f77f00", "#fcbf49"))+
  ggtitle("Tracking table")# Main title
bar_track 

##### Taxonomy assignment

# If on cluster : 
taxa <- assignTaxonomy(seqtab.nochim, "~/projects/rrg-yergeaue/projects/DADA2_Workshop/data/raw_reads/Silva_taxonomy/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/projects/rrg-yergeaue/projects/DADA2_Workshop/data/raw_reads/Silva_taxonomy/silva_species_assignment_v138.1.fa.gz")

# If on personal computer 
taxa <- assignTaxonomy(seqtab.nochim, "data/silva_taxonomy/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "data/silva_taxonomy/silva_species_assignment_v138.1.fa.gz")
  
taxa <- as.data.frame(taxa)

table(taxa$Genus)
toc()

library(phyloseq)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               #sample_data(samdf), 
               tax_table(taxa))

palette <- c("#355070","#6d597a","#b56576","#e56b6f","#eaac8b",
             "#ff9f1c","#ffbf69","#ffffff","#cbf3f0","#2ec4b6",
             "#ff99c8","#fcf6bd","#d0f4de","#a9def9","#e4c1f9","red")

plot_bar(ps,fill="Family")+
geom_bar(stat = "identity")+
  theme_bw()+
  theme(title = element_text(size=22),
        axis.text.x = element_text(size=5,angle=45,hjust=1),# Changes the x labels orientation 
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=12,angle=90,hjust=0.5),
        legend.title = element_text(size=18),
        legend.text = element_text(size=8,face = "italic"),
        legend.background = element_rect(color="black"))+
  scale_fill_manual(values=palette)+
  scale_y_continuous(name="Abundance",expand = expansion(mult = c(0, 0.05)))+ #Changes y-axis titl
  labs(fill="Taxa")

