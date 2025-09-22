###PREPPING YOUR R ENVIRONMENT###

##Load up SummarySE function##

#To run this put cursor in front of summarySE below and click Run above. This is a function created by someone else that we borrow, never seen a citation for it.#

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##Install packages##

#Some packages require BiocManager to be installed before being able to actually install the package you want, think of it as like a non-default conda channel#
install.packages("BiocManager")

#phyloseq is the main package you will use to do most of the analysis scripts. Like QIIME2 it is a wrapper of a lot of other tools. You will also get vegan (statistics), ggplot2(graphics), and plyr(data manipulation) by installing this#
BiocManager::install("phyloseq")

#Install other useful packages#
install.packages(c("tidyverse", "FSA", "dplyr", "reshape", "rcompanion", "devtools", "ggpubr"))
BiocManager::install("DESeq2")
BiocManager::install("metagenomeSeq")
BiocManager::install("decontam")
devtools::install_github("yanlinlin82/ggvenn")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
BiocManager::install("TreeSummarizedExperiment")
BiocManager::install("phangorn")


##Load packages##
library(phyloseq)
library(FSA)
library(vegan)
library(plyr)
library(reshape)
library(DESeq2)
library(metagenomeSeq)
library(ggvenn)
library(rcompanion)
library(car)
library(ggpubr)
library(ggtext)
library(decontam)
library(pairwiseAdonis)
library(tidyverse)
library(phangorn)

####Uploading data into R####

#Upload files from dada2 analysis from Chris Anderson
taxa <- as.matrix(read.delim("ASVs_taxonomy.tsv", row.names = 1))
asv_tab <- read.delim("ASVs_counts.tsv", row.names = 1)
samdf <- read.delim("16S_Chicken_Cecal_Vacc_SE_Chall_metadata.txt", row.names = 1)

#Select the same rownames as there are in the sample dataframe
seqtab.nochim <- 
  select(asv_tab, row.names(samdf))

#Make it a matrix
seqtab.nochim <- as.matrix(seqtab.nochim)

#Create a phyloseq object
post_dada2_data <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=TRUE), 
                            sample_data(samdf), 
                            tax_table(taxa))
post_dada2_data
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2066 taxa and 81 samples ]
# sample_data() Sample Data:       [ 81 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 2066 taxa by 6 taxonomic ranks ]

####Filtering data bDataframe manipulation and summarizing sequence####

##Change taxa names if needed##

#See what present names are#
colnames(tax_table(post_dada2_data))

#Switch them to classic format#
colnames(tax_table(post_dada2_data))=c("Kingdom","Phylum","Class","Order","Family","Genus")

#Check that it worked#
colnames(tax_table(post_dada2_data))

#Check number of taxa and samples#
ntaxa(post_dada2_data) #2066#
nsamples(post_dada2_data) #81#

#Get number of sequences
as.data.frame(taxa_sums(post_dada2_data)) %>%
  dplyr::rename("Sequences" = "taxa_sums(post_dada2_data)") %>%
  dplyr::summarize(sum(Sequences))
#2620765

#Remove taxonomy exclusive to other project
post_dada2_data_wo_zeros <- filter_taxa(post_dada2_data, function(x) sum(x) >0, TRUE)
ntaxa(post_dada2_data_wo_zeros) #804#

#Get number of sequences
as.data.frame(taxa_sums(post_dada2_data_wo_zeros)) %>%
  dplyr::rename("Sequences" = "taxa_sums(post_dada2_data_wo_zeros)") %>%
  dplyr::summarize(sum(Sequences))
#2620765

#Remove ASVs with no annotation
post_dada2_data_wo_zeros_unknowns <- subset_taxa(post_dada2_data_wo_zeros, !is.na(Kingdom))
ntaxa(post_dada2_data_wo_zeros_unknowns) #764#

#Get number of sequences
as.data.frame(taxa_sums(post_dada2_data_wo_zeros_unknowns)) %>%
  dplyr::rename("Sequences" = "taxa_sums(post_dada2_data_wo_zeros_unknowns)") %>%
  dplyr::summarize(sum(Sequences))
#2618790

#Remove ASVs with Archaea designation
post_dada2_data_wo_zeros_unknowns_archaea <- subset_taxa(post_dada2_data_wo_zeros_unknowns, Kingdom != "Archaea")
ntaxa(post_dada2_data_wo_zeros_unknowns_archaea) #763#

#Get number of sequences
as.data.frame(taxa_sums(post_dada2_data_wo_zeros_unknowns_archaea)) %>%
  dplyr::rename("Sequences" = "taxa_sums(post_dada2_data_wo_zeros_unknowns_archaea)") %>%
  dplyr::summarize(sum(Sequences))
#2618786

####Decontam####


##Take a look at library size
# Put sample_data into a ggplot-friendly data.frame (go to matrix first to remove sample_data object label)
lib_size_df <- as.data.frame(as.matrix(sample_data(post_dada2_data_wo_zeros_unknowns_archaea)))

#Add the number of sequences
lib_size_df$LibrarySize <- sample_sums(post_dada2_data_wo_zeros_unknowns_archaea)

#Put rows in order by library size
lib_size_df <- lib_size_df[order(lib_size_df$LibrarySize),]

#Add a column designating the Rank of increasing library size from 1 to x
lib_size_df$Index <- seq(nrow(lib_size_df))

#Plot the library size
ggplot(data=lib_size_df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + 
  geom_point()

####Identify Contaminants - Frequency####

#Get number of samples
nsamples(post_dada2_data_wo_zeros_unknowns_archaea)
#81/81

#Remove control samples without DNA quantification
post_dada2_data_wo_zeros_unknowns_archaea_controls <- 
  subset_samples(post_dada2_data_wo_zeros_unknowns_archaea, Sample_or_Control=="True Sample")

#Check have removed 2 control samples
nsamples(post_dada2_data_wo_zeros_unknowns_archaea_controls) #Should be 6#
#79/79

#Use DNA quantification to check for contaminants i.e. frequency method
contamdf.freq <- isContaminant(post_dada2_data_wo_zeros_unknowns_archaea_controls, method="frequency", conc="quant_reading")

#Check to see what the resulting anlaysis looks like
head(contamdf.freq)
#the p column contains the probability that was used to classify contaminants
#Default threshold is 0.1
#The contaminant column helps determine if the ASV was a contaminant or not in TRUE/FALSE designation

#Get a summary of the ASVs identified as contaminants
table(contamdf.freq$contaminant)
# FALSE  TRUE 
# 758     5 

#Get the ranks of the ASV frequencies and display the first 6 
head(which(contamdf.freq$contaminant))
#335 452 533 549 616
#Just 5 out of the 763 ASVs were classified as contaminants, 

#Compare a non-contaminant (2) to a contaminant (335)
plot_frequency(post_dada2_data_wo_zeros_unknowns_archaea_controls, taxa_names(post_dada2_data_wo_zeros_unknowns_archaea_controls)[c(2,335)], conc="quant_reading") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")
#Black line is model of a noncontaminant sequence i.e. independent of DNA concentration
#Red line shows model of contaminant i.e. frequency is inversely proportional to input DNA concentration

#Remove contamnant reads
post_dada2_data_wo_zeros_unknowns_archaea_contam <- prune_taxa(!contamdf.freq$contaminant, post_dada2_data_wo_zeros_unknowns_archaea)

#Check it worked
ntaxa(post_dada2_data_wo_zeros_unknowns_archaea_contam)
#758/758

#Get number of sequences
as.data.frame(taxa_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam)) %>%
  dplyr::rename("Sequences" = "taxa_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam)") %>%
  dplyr::summarize(sum(Sequences))
#2618064


####Analysis of Filtered Data####

##Remove controls samples before statistical analysis
#Remove control samples from analysis
post_dada2_data_wo_zeros_unknowns_archaea_contam_controls <- subset_samples(post_dada2_data_wo_zeros_unknowns_archaea_contam, Sample_or_Control=="True Sample")
nsamples(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls)
#79/79

#Get number of sequences
as.data.frame(taxa_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls)) %>%
  dplyr::rename("Sequences" = "taxa_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls)") %>%
  dplyr::summarize(sum(Sequences))
#2567747

#Check number of ASVs
ntaxa(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls)
#758

#Remove taxa that are no longer present with removal of controls
post_dada2_data_wo_zeros_unknowns_archaea_contam_controls <- 
  filter_taxa(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls, function(x) sum(x) >0, TRUE)
ntaxa(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls)
#720


##Remove samples with less than 10,000 reads

#Save a variable of the number of sequences you want to filter by
MINREADS=10000

#Filter for samples with less than 10,000 reads
post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold <- 
  prune_samples(sample_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls) >= MINREADS, post_dada2_data_wo_zeros_unknowns_archaea_contam_controls)

#Check number of samples
nsamples(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold)
#78/78

#Get number of sequences
as.data.frame(taxa_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold)) %>%
  dplyr::rename("Sequences" = "taxa_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold)") %>%
  dplyr::summarize(sum(Sequences))
#2557791

#Check number of ASVs
ntaxa(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold)
#720

#Remove taxa that are no longer present with removal of below threshold sample(s)
post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold <- 
  filter_taxa(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold, function(x) sum(x) >0, TRUE)
ntaxa(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold)
#717



#Save the preanalysis phyloseq in a shorter name for conciseness moving forward
data <- post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold


##Filter low abundance ASVs

#Remove rare ASVs by removing those that occur less than 11 times
Fdata <- filter_taxa(data, function(x) sum(x) >10, TRUE)
ntaxa(Fdata) #532#

#Print out taxonomy table for future reference
write.csv(tax_table(Fdata), 'Fdata_taxa.csv')

####Track the changes in ASVs per sample over time - Supplementary Table 1####

#Combine the first sets of filters that just removed ASVs with cbind
post_dada2_sequences_per_sample_track <- 
  cbind(sample_sums(post_dada2_data), 
        sample_sums(post_dada2_data_wo_zeros),
        sample_sums(post_dada2_data_wo_zeros_unknowns),
        sample_sums(post_dada2_data_wo_zeros_unknowns_archaea),
        sample_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam)) %>%
  data.frame()

#Note: Need to use data frames and thus merge after removing samples

#Add in ASVs per sample post control sample removal
post_dada2_sequences_per_sample_track <- 
  merge(post_dada2_sequences_per_sample_track,
        sample_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls), 
        by=0,
        all.x=TRUE)

#Add in ASVs per sample post removal of samples below sequence threshold
post_dada2_sequences_per_sample_track <- 
  merge(post_dada2_sequences_per_sample_track,
        sample_sums(post_dada2_data_wo_zeros_unknowns_archaea_contam_controls_threshold), 
        by.x="Row.names", by.y = 0,
        all.x=TRUE)

#Add in ASVs per sample post removal of low abundance ASVs
post_dada2_sequences_per_sample_track <- 
  merge(post_dada2_sequences_per_sample_track,
        sample_sums(Fdata), 
        by.x="Row.names", by.y = 0,
        all.x=TRUE)

#Add column names
colnames(post_dada2_sequences_per_sample_track) <-
  c("Sample.name", "post_dada2", "no_zero_asvs", "no_unknown_asvs", 
    "no_archaea_asvs", "no_contaminant_asvs", "no_control_samples", "no_threshold_samples", "no_low_abd_asvs")



##Summarize sequences by ASVs (filtered and unfiltered)## 

#Unfiltered version first#

#export table of sequence sums#
data_seqs_per_ASV <- as.data.frame(taxa_sums(data))

#Change colnames to Sequences#
colnames(data_seqs_per_ASV) <- c("Sequences")
sum(data_seqs_per_ASV$Sequences) #2557791#

#Use cbind to get taxonomy added in#
data_seqs_per_ASV = cbind(as(data_seqs_per_ASV, "data.frame"), as(tax_table(data)[rownames(data_seqs_per_ASV), ], "matrix"))

#export it as a csv, and tell it that the the rownames should be named OTUID#
write.csv(data.frame("OTUID" =rownames(data_seqs_per_ASV), data_seqs_per_ASV) , "seqs_per_ASV.csv", row.names=FALSE)

#Filtered version second#

Fdata_seqs_per_ASV <- as.data.frame(taxa_sums(Fdata))
colnames(Fdata_seqs_per_ASV) <- c("Filtered_Sequences")
sum(Fdata_seqs_per_ASV$Filtered_Sequences) #2556945#
Fdata_seqs_per_ASV = cbind(as(Fdata_seqs_per_ASV, "data.frame"), as(tax_table(Fdata)[rownames(Fdata_seqs_per_ASV), ], "matrix"))
write.csv(data.frame("OTUID" =rownames(Fdata_seqs_per_ASV), Fdata_seqs_per_ASV) , "filtered_seqs_per_ASV.csv", row.names=FALSE)


####data Version - Alpha & Beta Analysis####

####Determine Various Taxonomic Level Numbers - data phyloseq####

datataxa <- as.data.frame(tax_table(data))

datataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>%nrow #13
datataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>% filter(!is.na(Phylum), Phylum!="uncultured") %>%nrow #12

datataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>%nrow #25
datataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>% 
  filter(!is.na(Class), !grepl(pattern = "uncultured", x = Class)) %>% 
  filter(!grepl(pattern = "unidentified", x = Class))%>% 
  filter(!grepl(pattern = "metagenome", x = Class)) %>% 
  filter(as.character(Phylum) != as.character(Class))%>%
  nrow #23

datataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>%nrow #49
datataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>% 
  filter(!is.na(Order), !grepl(pattern = "uncultured", x = Order)) %>% 
  filter(!grepl(pattern = "unidentified", x = Order))%>% 
  filter(!grepl(pattern = "metagenome", x = Order))%>%
  filter(!grepl(pattern = "Unknown", x = Order))%>%
  filter(as.character(Class) != as.character(Order))%>%
  nrow #45

datataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>%nrow #80
datataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>% 
  filter(!is.na(Family), !grepl(pattern = "uncultured", x = Family)) %>% 
  filter(!grepl(pattern = "unidentified", x = Family))%>% 
  filter(!grepl(pattern = "metagenome", x = Family))%>%
  filter(!grepl(pattern = "Unknown", x = Family))%>%
  filter(as.character(Order) != as.character(Family))%>%
  nrow #69

datataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>%nrow #152
datataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>% 
  filter(!is.na(Genus), !grepl(pattern = "uncultured", x = Genus)) %>% 
  filter(!grepl(pattern = "unidentified", x = Genus))%>% 
  filter(!grepl(pattern = "metagenome", x = Genus))%>%
  filter(!grepl(pattern = "Unknown", x = Genus))%>%
  filter(as.character(Family) != as.character(Genus))%>%
  nrow #119

#Extract out and copy otu/taxonomy tables for analysis in Excel
Uniq_data_taxa_genus_level <- 
  distinct(datataxa, Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE) #152

##Add three columns
#last_non_na_taxa_name = Name of last non-NA value in a row for taxonomy designation
#last_non_na_taxa_category = Column name aka taxonomic category of last_non_na_taxa_name
#name_col = combination of previous columns for visualization
Uniq_data_taxa_genus_level <-
  Uniq_data_taxa_genus_level %>% 
  dplyr::mutate(last_non_na_taxa_name = coalesce(Genus, Family, Order, Class, Phylum, Kingdom)) %>% 
  rowwise() %>% 
  dplyr::mutate(last_non_na_taxa_category = rev(names(Uniq_data_taxa_genus_level)[!is.na(c_across(Kingdom:Genus))])[1]) %>%
  dplyr::mutate(name_col = paste0(last_non_na_taxa_category, ":", last_non_na_taxa_name))



#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g_uncultured, with different preceding taxonomy strings
Uniq_data_taxa_genus_level$Taxonomy <- paste(Uniq_data_taxa_genus_level$Kingdom,
                                             Uniq_data_taxa_genus_level$Phylum,
                                             Uniq_data_taxa_genus_level$Class,
                                             Uniq_data_taxa_genus_level$Order,
                                             Uniq_data_taxa_genus_level$Family,
                                             Uniq_data_taxa_genus_level$Genus,
                                             sep=":")

#Fuse Family and Genus columns
Uniq_data_taxa_genus_level<- unite(Uniq_data_taxa_genus_level, Family_Genus, Family:Genus, sep=';', remove=FALSE)

####Determine Various Taxonomic Level Numbers - Fdata phyloseq####

Fdatataxa <- as.data.frame(tax_table(Fdata))

Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>%nrow #5
Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>% filter(!is.na(Phylum), Phylum!="uncultured") %>%nrow #4

Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>%nrow #8
Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>% 
  filter(!is.na(Class), !grepl(pattern = "uncultured", x = Class)) %>% 
  filter(!grepl(pattern = "unidentified", x = Class))%>% 
  filter(!grepl(pattern = "metagenome", x = Class)) %>% 
  filter(as.character(Phylum) != as.character(Class))%>%
  nrow #6

Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>%nrow #23
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>% 
  filter(!is.na(Order), !grepl(pattern = "uncultured", x = Order)) %>% 
  filter(!grepl(pattern = "unidentified", x = Order))%>% 
  filter(!grepl(pattern = "metagenome", x = Order))%>%
  filter(!grepl(pattern = "Unknown", x = Order))%>%
  filter(as.character(Class) != as.character(Order))%>%
  nrow #19

Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>%nrow #40
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>% 
  filter(!is.na(Family), !grepl(pattern = "uncultured", x = Family)) %>% 
  filter(!grepl(pattern = "unidentified", x = Family))%>% 
  filter(!grepl(pattern = "metagenome", x = Family))%>%
  filter(!grepl(pattern = "Unknown", x = Family))%>%
  filter(as.character(Order) != as.character(Family))%>%
  nrow #33

Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>%nrow #105
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>% 
  filter(!is.na(Genus), !grepl(pattern = "uncultured", x = Genus)) %>% 
  filter(!grepl(pattern = "unidentified", x = Genus))%>% 
  filter(!grepl(pattern = "metagenome", x = Genus))%>%
  filter(!grepl(pattern = "Unknown", x = Genus))%>%
  filter(as.character(Family) != as.character(Genus))%>%
  nrow #85

#Extract out and copy otu/taxonomy tables for analysis in Excel
Uniq_Fdata_taxa_genus_level <- 
  distinct(Fdatataxa, Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE) #120

##Add three columns
#last_non_na_taxa_name = Name of last non-NA value in a row for taxonomy designation
#last_non_na_taxa_category = Column name aka taxonomic category of last_non_na_taxa_name
#name_col = combination of previous columns for visualization
Uniq_Fdata_taxa_genus_level <-
  Uniq_Fdata_taxa_genus_level %>% 
  dplyr::mutate(last_non_na_taxa_name = coalesce(Genus, Family, Order, Class, Phylum, Kingdom)) %>% 
  rowwise() %>% 
  dplyr::mutate(last_non_na_taxa_category = rev(names(Uniq_Fdata_taxa_genus_level)[!is.na(c_across(Kingdom:Genus))])[1]) %>%
  dplyr::mutate(name_col = paste0(last_non_na_taxa_category, ":", last_non_na_taxa_name))



#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g_uncultured, with different preceding taxonomy strings
Uniq_Fdata_taxa_genus_level$Taxonomy <- paste(Uniq_Fdata_taxa_genus_level$Kingdom,
                                              Uniq_Fdata_taxa_genus_level$Phylum,
                                              Uniq_Fdata_taxa_genus_level$Class,
                                              Uniq_Fdata_taxa_genus_level$Order,
                                              Uniq_Fdata_taxa_genus_level$Family,
                                              Uniq_Fdata_taxa_genus_level$Genus,
                                              sep=":")

#Fuse Family and Genus columns
Uniq_Fdata_taxa_genus_level<- unite(Uniq_Fdata_taxa_genus_level, Family_Genus, Family:Genus, sep=';', remove=FALSE)



####Combining Taxonomic Level Dataframe - Supplementary Table 2####

#Add a column to designate it as removed in Fdata or not
Uniq_data_taxa_genus_level$Fdata <-
  Uniq_data_taxa_genus_level$Taxonomy %in% Uniq_Fdata_taxa_genus_level$Taxonomy

#Get distribution
table(Uniq_data_taxa_genus_level$Fdata)
# FALSE  TRUE 
# 47   105 

#Send to Excel 
clipr::write_clip(Uniq_data_taxa_genus_level)


####Alpha Diversity####

##Calculate alpha diversity##
alphadiv = estimate_richness(data, split = TRUE)
write.csv(alphadiv, file='alphadiv.csv')
alphadiv <- read.csv("alphadiv.csv", row.names = 1)

##cbind in the metadata from the filtered phyloseq
alphadiv_metadata = cbind(as(alphadiv, "data.frame"), as(sample_data(data)[rownames(alphadiv), ], "data.frame"))

####Alpha Diversity with Metadata - Supplementary Table 3 - Used to Create Figure 4####
#Send to Excel
clipr::write_clip(alphadiv_metadata)

##Test for normality##
shapiro.test(alphadiv_metadata$Shannon)
# W = 0.96872, p-value = 0.05193

shapiro.test(alphadiv_metadata$Observed)
#W = 0.98617, p-value = 0.5632

shapiro.test(alphadiv_metadata$Fisher)
#W = 0.98956, p-value = 0.7808

shapiro.test(alphadiv_metadata$Simpson)
#W = 0.88052, p-value = 2.574e-06

##Correlation Tests##
cor.test(alphadiv_metadata$Shannon, alphadiv_metadata$Observed, method = "spearman", exact = FALSE)
#S = 44261, p-value = 5.485e-05, rho = 0.4402978   

cor.test(alphadiv_metadata$Shannon, alphadiv_metadata$Fisher, method = "spearman", exact = FALSE)
#S = 31012, p-value = 3.6e-09, rho = 0.6078352   

cor.test(alphadiv_metadata$Shannon, alphadiv_metadata$Simpson, method = "spearman", exact = FALSE)
#S = 4718, p-value < 2.2e-16, rho = 0.9403381   

##Basic sumamry of Shannon statistics##
mean(alphadiv_metadata$Shannon)#3.763385
sd(alphadiv_metadata$Shannon) #0.2638508

#Create a summary of Shannon Diversity by Vaccination Status and dpi - Reference Table
alphadiv_metadata %>%
  group_by(VS_DPC) %>%
  summarize(mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon), count=n()) %>%
  data.frame() %>%
  clipr::write_clip()
# VS_DPC     mean        sd   median       IQR count
# 1   M_D0 3.542406 0.2730713 3.478568 0.3399132    12
# 2  M_D14 3.744577 0.2647807 3.727120 0.3185737    13
# 3   M_D7 3.759635 0.2973733 3.892001 0.3774272    13
# 4   V_D0 3.735634 0.1985803 3.770352 0.2903547    11
# 5  V_D14 3.982450 0.1817325 3.965613 0.1258484    15
# 6   V_D7 3.760835 0.1885009 3.779070 0.2139858    14



####Tests For Total And Pairwise Significance####

###Kruskal-Wallis####
kruskal.test(Shannon ~ VS_DPC, data = alphadiv_metadata)
# Kruskal-Wallis chi-squared = 19.03, df = 5, p-value = 0.001897

##Dunn:
Dunn <- FSA::dunnTest(Shannon ~ VS_DPC, data = alphadiv_metadata, method="bh") 
Dunn
# Comparison                Z   P.unadj     P.adj
# 1   M_D0 - M_D14 -1.77225008 7.635305e-02 0.1636136866
# 2    M_D0 - M_D7 -2.01816038 4.357456e-02 0.1089364028
# 3   M_D14 - M_D7 -0.25098115 8.018287e-01 0.8591021627
# 4    M_D0 - V_D0 -1.36473548 1.723362e-01 0.2872270058
# 5   M_D14 - V_D0  0.34123560 7.329262e-01 0.9994448481
# 6    M_D7 - V_D0  0.58153176 5.608821e-01 0.8413231814
# 7   M_D0 - V_D14 -4.23104927 2.326037e-05 0.0003489056
# 8  M_D14 - V_D14 -2.45218030 1.419935e-02 0.0532475598
# 9   M_D7 - V_D14 -2.19239013 2.835135e-02 0.0850540451
# 10  V_D0 - V_D14 -2.69299333 7.081369e-03 0.0531102638
# 11   M_D0 - V_D7 -1.75474423 7.930307e-02 0.1486932516
# 12  M_D14 - V_D7  0.04973239 9.603356e-01 0.9603356485
# 13   M_D7 - V_D7  0.30531908 7.601231e-01 0.9501539363
# 14   V_D0 - V_D7 -0.29942048 7.646192e-01 0.8822529637
# 15  V_D14 - V_D7  2.55203404 1.070961e-02 0.0535480271
#6/15 are unneeded

Dunn = Dunn$res

####Dunn Shannon Diversity - Supplementary Table 4####
#Send down to Excel
clipr::write_clip(Dunn)

#Get the Compact Letter Display for VS_DPC designation
Dunn_VS_DPC_cld <- rcompanion::cldList(comparison = Dunn$Comparison,
                                       p.value    = Dunn$P.adj,
                                       threshold  = 0.05)
Dunn_VS_DPC_cld

# Group Letter MonoLetter
# 1   M_D      a         a 
# 2 M_D14     ab         ab
# 3  M_D7     ab         ab
# 4   V_D     ab         ab
# 5 V_D14      b          b
# 6  V_D7     ab         ab

#Original Order: "a", "ab", "ab", "ab", "b", "ab"
#Actual Order: "a", "ab", "ab", "ab", "ab", "b"

#Send down to Excel - Reference Table
clipr::write_clip(FDunn_VS_DPC_cld)

####Beta Diversity Analysis Prep####

##Determining best transformation##

#Transform data using square root
SR_data <- transform_sample_counts(data, function(x){x^(1/2)})


#Alternative:
#https://github.com/joey711/phyloseq/issues/299
#https://github.com/joey711/phyloseq/issues/229

#Install and load DESeq2 so you can use their rlog based transfomation
library("DESeq2")

#Load alternative gm_mean function to handle zeros based upon github issue shown above#
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
data_dds = phyloseq_to_deseq2(data, ~1)

#Calculate geometric means prior to estimate size factors#
#https://github.com/joey711/phyloseq/issues/445
geoMeans = apply(counts(data_dds), 1, gm_mean)
data_dds = estimateSizeFactors(data_dds, geoMeans = geoMeans)

#Conduct DESEQ2 test#
data_dds = DESeq(data_dds, fitType="local")

#Make a copy of data so you can have a designated DESeq transformed copy
data
DS_data = data
DS_data

#Switch the asv table with the DESeq2 transformed data
otu_table(DS_data) <- otu_table(getVarianceStabilizedData(data_dds), taxa_are_rows = TRUE)

#Check to see if your basic phyloseq information was not changed
DS_data

#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
#https://github.com/joey711/phyloseq/issues/445

#Make any negative values equal to 0 since the more negative they are the more likely they were to be zero over very small and unlikely to affect your results

#Make a copy of DS_data to manipulate
ZDS_data <- DS_data
ZDS_data

#extract out asv table
DESeq2_otu_table <- as.data.frame(otu_table(ZDS_data))

#Change all negatives to zero
DESeq2_otu_table[DESeq2_otu_table < 0.0] <- 0.0

#Switch out the asv table in phyloseq object
otu_table(ZDS_data) <-otu_table(DESeq2_otu_table, taxa_are_rows = TRUE)

#Check to make sure basic phyloseq info did not change
ZDS_data

#Show how amount of positive numbers changed throughout the transformation 
z <- otu_table(data)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE      TRUE 
# 0.7610593 0.2389407 

z <- otu_table(DS_data)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# 0             1 

z <- otu_table(ZDS_data)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# 0             1 

#Standardize with CSS:
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

#Install and load
library(metagenomeSeq)

#Convert file types to use in metagenomeSeq
MGS <- phyloseq_to_metagenomeSeq(data) 

#Perform normalization following: https://bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf
CSSp <- cumNormStatFast(MGS)
CSSp
MGS <- cumNorm(MGS, p = CSSp)

#Store post-transformation abundance table
normmybiom <- MRcounts(MGS, norm = T)

#Create copy of filtered data to be modified
MGS_data = data
MGS_data

#Switch out old ASV table for transformed one
otu_table(MGS_data) <-otu_table(normmybiom, taxa_are_rows = TRUE)
MGS_data

#Plot each of the different methods to see the differences#
plot_ordination(SR_data, ordinate(SR_data, "PCoA", "bray"), color = "VS_DPC") 
plot_ordination(ZDS_data, ordinate(ZDS_data, "PCoA", "bray"), color = "VS_DPC")
plot_ordination(MGS_data, ordinate(MGS_data, "PCoA", "bray"), color = "VS_DPC")

#Check your library size differences by dividing the sequences sum of the sample with the most sequences by the one with the least sequences typically want it to be <10x difference

max(sample_sums(data))/min(sample_sums(data))#6.531277
max(sample_sums(SR_data))/min(sample_sums(SR_data)) #2.800138
max(sample_sums(ZDS_data))/min(sample_sums(ZDS_data)) #1.084707
max(sample_sums(MGS_data))/min(sample_sums(MGS_data)) #3.385987

#Check to see if they are statistically normal, if >0.05 then it is likely normal
shapiro.test(sample_sums(data))
#W = 0.96745, p-value = 0.04329

shapiro.test(sample_sums(SR_data)) 
#W = 0.98914, p-value = 0.7546

shapiro.test(sample_sums(ZDS_data)) 
#W = 0.98498, p-value = 0.4914

shapiro.test(sample_sums(MGS_data)) 
#W = 0.94904, p-value = 0.003526

#Visual representation of the sample sums
hist(sample_sums(data))
hist(sample_sums(SR_data))
hist(sample_sums(ZDS_data))
hist(sample_sums(MGS_data))

####One-Way PERMANOVA Analysis - Bray-Curtis - Supplementary Table 5####

#Run total PERMANOVA on  data based on variable of choice#
#Make a distance matrix#
#Note that we must use the phyloseq::distance here instead of just distance 
#because other libraries we are using also have the distance function. 
#This allows us to tell R to use the phyloseq version. 
ZDS_data_Bray<-phyloseq::distance(ZDS_data, "bray")

#convert metadata of the phyloseq to a dataframe#
ZDS_data_sampledf<-data.frame(sample_data(ZDS_data))

#run a PERMANVOA using vegan's adonis function#
vegan::adonis2(ZDS_data_Bray~VS_DPC, data=ZDS_data_sampledf, permutations=9999)
#           Df SumOfSqs     R2      F Pr(>F)    
# Model     5  0.11846 0.50586 14.741  1e-04 ***
#   Residual 72  0.11571 0.49414                  
# Total    77  0.23417 1.00000 

##Pairwise PERMANOVA
ZDS_data_pairwise_adonis_vs_dpc_BC <-
  pairwise.adonis(ZDS_data_Bray, ZDS_data_sampledf$VS_DPC, perm = 9999, p.adjust.m="BH")
ZDS_data_pairwise_adonis_vs_dpc_BC

# pairs Df   SumsOfSqs   F.Model        R2 p.value p.adjusted sig
# 1   M_D0 vs M_D14  1 0.019502131 11.920316 0.3413576   1e-04      1e-04 ***
#   2    M_D0 vs M_D7  1 0.012128593  8.296121 0.2650846   1e-04      1e-04 ***
#   3    M_D0 vs V_D0  1 0.029155558 25.045685 0.5439312   1e-04      1e-04 ***
#   4   M_D0 vs V_D14  1 0.040356421 27.961487 0.5279589   1e-04      1e-04 ***
#   5    M_D0 vs V_D7  1 0.030962482 21.882282 0.4769223   1e-04      1e-04 ***
#   6   M_D14 vs M_D7  1 0.005679625  2.947923 0.1093933   1e-04      1e-04 ***
#   7   M_D14 vs V_D0  1 0.028846868 17.124148 0.4376874   1e-04      1e-04 ***
#   8  M_D14 vs V_D14  1 0.027888554 14.890161 0.3641502   1e-04      1e-04 ***
#   9   M_D14 vs V_D7  1 0.027355370 14.683959 0.3700225   1e-04      1e-04 ***
#   10   M_D7 vs V_D0  1 0.025159430 16.744178 0.4321728   1e-04      1e-04 ***
#   11  M_D7 vs V_D14  1 0.033247037 19.341409 0.4265727   1e-04      1e-04 ***
#   12   M_D7 vs V_D7  1 0.025813142 15.159349 0.3774799   1e-04      1e-04 ***
#   13  V_D0 vs V_D14  1 0.026228926 17.725347 0.4248100   1e-04      1e-04 ***
#   14   V_D0 vs V_D7  1 0.014709705 10.132289 0.3058131   1e-04      1e-04 ***
#   15  V_D14 vs V_D7  1 0.008367177  5.017323 0.1567065   1e-04      1e-04 ***

#Send to Excel
clipr::write_clip(ZDS_data_pairwise_adonis_vs_dpc_BC)

####PCoA Analysis - Bray-Curtis - Axis 1 & 2 - Figure 5A####

#Sample Type by Vaccination Status by dpi
ZDS_data_VS_DPC_PCoA_BC <- 
  plot_ordination(ZDS_data, ordinate(ZDS_data, "PCoA", "bray"), 
                  shape = "Vaccination_Status", color = "Days_Post_Challenge")+ 
  scale_color_manual(name="Days Post Inoculation", 
                     breaks = c("D0", "D7", "D14"),
                     labels=c("0", "7", "14"),
                     values=c("purple", "red", "orange"))+
  scale_shape_discrete("Vaccination Status",
                       labels=c("Mock","BBS 1134")) +
  geom_point(size=7)+
  ggtitle("Bray-Curtis") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=15, face = "bold")) +
  guides(color = guide_legend(override.aes = list(shape = c(15, 15, 15))))+
  theme(axis.text=element_text(size=14, face = "bold", color = "black"), axis.title=element_text(size=15, face = "bold"), 
        legend.text=element_text(size=12, face = "bold", color = "black"), legend.key.size = unit(0.7,"cm"), 
        legend.title=element_text(size=13, face = "bold")) 
ZDS_data_VS_DPC_PCoA_BC

ggsave("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Manuscripts/Bearson_Vacc_vs_DSC_Enteritidis_Challenge/Figures/Figure 5A - Bray-Curtis PCoA VS by DPI.tiff",
       ZDS_data_VS_DPC_PCoA_BC,
       device = "tiff", width = 10, height = 6)

####PCoA Analysis - Bray-Curtis - Axis 1 & 2 - Supplementary Figure 1A####

ZDS_data_VS_DPC_PCoA_BC_Axes_1_3 <- 
  plot_ordination(ZDS_data, ordinate(ZDS_data, "PCoA", "bray"), 
                  shape = "Vaccination_Status", color = "Days_Post_Challenge", axes = c(1,3))+ 
  scale_color_manual(name="Days Post Inoculation", 
                     breaks = c("D0", "D7", "D14"),
                     labels=c("0", "7", "14"),
                     values=c("purple", "red", "orange"))+
  scale_shape_discrete("Vaccination Status",
                       labels=c("Mock","BBS 1134")) +
  geom_point(size=7)+
  ggtitle("Bray-Curtis") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=15, face = "bold")) +
  guides(color = guide_legend(override.aes = list(shape = c(15, 15, 15))))+
  theme(axis.text=element_text(size=14, face = "bold", color = "black"), axis.title=element_text(size=15, face = "bold"), 
        legend.text=element_text(size=12, face = "bold", color = "black"), legend.key.size = unit(0.7,"cm"), 
        legend.title=element_text(size=13, face = "bold")) 
ZDS_data_VS_DPC_PCoA_BC_Axes_1_3

ggsave("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Manuscripts/Bearson_Vacc_vs_DSC_Enteritidis_Challenge/Figures/Supplementary Figure 1A - Bray-Curtis PCoA VS by DPI - Axes 1 and 3.tiff",
       ZDS_data_VS_DPC_PCoA_BC_Axes_1_3,
       device = "tiff", width = 10, height = 6)

####data Phylogentic Tree - RESTART HERE####

#Upload the fasta sequences
seqtab <- readDNAStringSet("ASVs.fa")

#Create a dataframe with ASV number as name and sequence as the DNA sequence
seqtab_df <-
  data.frame(name = names(seqtab), sequence =paste(seqtab))

#Filter for ASVs in the phyloseq
data_seqtab_df <-
  filter(seqtab_df, name %in% rownames(otu_table(data)))

#Create a sequences list
data_seqs <- data_seqtab_df$sequence

#Add names based upon ASV values
names(data_seqs) <- data_seqtab_df$name # This propagates to the tip labels of the tree

#Make an alignment of the sequences
data_alignment <- AlignSeqs(DNAStringSet(data_seqs), anchor=NA)

#Transform alignment to a phyDat format
data_phang.align <- phyDat(as(data_alignment, "matrix"), type="DNA")

#Compute pairwise distance
data_dm <- dist.ml(data_phang.align)

#Construct a neighbor-joining tree
data_treeNJ <- NJ(data_dm) # Note, tip order != sequence order

#Compute the likelhood of a tree given a squence alighment and model
data_fit = pml(data_treeNJ, data=data_phang.align)

## negative edges length changed to 0!
#
data_fitGTR <- update(data_fit, k=4, inv=0.2)

#Fit a GTR+G+I (generalized time-reversible with Gamma rate variation) maximum likelihood tree using the NJT as a starting point
data_fitGTR <- optim.pml(data_fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                         rearrangement = "stochastic", control = pml.control(trace = 0))

#Detach phagorn package
detach("package:phangorn", unload=TRUE)

#Make a copy of the main phyloseq
ZDS_data_tree <- ZDS_data

#Create a tree phyloseq object
data_tree <- phy_tree(data_fitGTR$tree)

#Change the tip labels to match copied phyloseq
data_tree$tip.label <- data_seqtab_df[[1]][match(data_tree$tip.label, data_seqtab_df[[2]])]

#Add to the copied phyolseq
ZDS_data_tree <- phyloseq(otu_table(ZDS_data),
                          phy_tree(data_tree),
                          sample_data(ZDS_data),
                          tax_table(ZDS_data))

#Check to see if it is rooted
is.rooted(phy_tree(ZDS_data_tree)) 
#FALSE

#Create a midpoint rooted tree
phy_tree(ZDS_data_tree)<-phangorn::midpoint(phy_tree(ZDS_data_tree))

#Or
#phy_tree(ZDS_data_tree) <- root(phy_tree(ZDS_data_tree), sample(taxa_names(ZDS_data_tree), 1), resolve.root = TRUE)

#Chick if it is rooted
is.rooted(phy_tree(ZDS_data_tree)) 
##TRUE

####PERMANOVA Analysis - Unweighted UniFrac####

#Run total PERMANOVA on  data based on variable of choice#
#Make a distance matrix#
#Note that we must use the phyloseq::distance here instead of just distance 
#because other libraries we are using also have the distance function. 
#This allows us to tell R to use the phyloseq version. 
ZDS_data_UW_UniFrac<-phyloseq::distance(ZDS_data_tree, "unifrac")

#convert metadata of the phyloseq to a dataframe#
ZDS_data_sampledf<-data.frame(sample_data(ZDS_data_tree))

#run a PERMANVOA using vegan's adonis function#
vegan::adonis2(ZDS_data_UW_UniFrac~VS_DPC, data=ZDS_data_sampledf, permutations=9999)
#Not working - Potential Reasons
#https://github.com/joey711/phyloseq/issues/936
#https://github.com/joey711/phyloseq/issues/299
#https://github.com/joey711/phyloseq/issues/255


####PERMANOVA Analysis - Weighted UniFrac - Supplementary Table 5B####

#Run total PERMANOVA on  data based on variable of choice#
#Make a distance matrix#
#Note that we must use the phyloseq::distance here instead of just distance 
#because other libraries we are using also have the distance function. 
#This allows us to tell R to use the phyloseq version. 
ZDS_data_W_UniFrac<-phyloseq::distance(ZDS_data_tree, "wunifrac")

#convert metadata of the phyloseq to a dataframe#
ZDS_data_sampledf<-data.frame(sample_data(ZDS_data_tree))

#run a PERMANVOA using vegan's adonis function#
vegan::adonis2(ZDS_data_W_UniFrac~VS_DPC, data=ZDS_data_sampledf, permutations=9999)
#           Df SumOfSqs     R2      F Pr(>F)    
# Model     5 0.0088956 0.393 9.3232  1e-04 ***
#   Residual 72 0.0137396 0.607                  
# Total    77 0.0226352 1.000 

##Pairwise PERMANOVA
ZDS_data_pairwise_adonis_vs_dpc_W_UniFrac <-
  pairwise.adonis(ZDS_data_W_UniFrac, ZDS_data_sampledf$VS_DPC, perm = 9999, p.adjust.m="BH")
ZDS_data_pairwise_adonis_vs_dpc_W_UniFrac

clipr::write_clip(ZDS_data_pairwise_adonis_vs_dpc_W_UniFrac)

####PCoA Analysis - Weighted Unifrac - Axes 1 & 2 - Figure 5B####

#Sample Type by Vaccination Status by dpi
ZDS_data_VS_DPC_PCoA_W_UniFrac <- 
  plot_ordination(ZDS_data_tree, ordinate(ZDS_data_tree, "PCoA", "wunifrac"), 
                  shape = "Vaccination_Status", color = "Days_Post_Challenge")+ 
  scale_color_manual(name="Days Post Inoculation", 
                     breaks = c("D0", "D7", "D14"),
                     labels=c("0", "7", "14"),
                     values=c("purple", "red", "orange"))+
  scale_shape_discrete("Vaccination Status",
                       labels=c("Mock","BBS 1134")) +
  geom_point(size=7)+
  ggtitle("Weighted UniFrac") +
  theme(plot.title = element_text(hjust = 0.5, size=15)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=15, face = "bold")) +
  guides(color = guide_legend(override.aes = list(shape = c(15, 15, 15))))+
  theme(axis.text=element_text(size=14, face = "bold", color = "black"), axis.title=element_text(size=15, face = "bold"), 
        legend.text=element_text(size=12, face = "bold", color = "black"), legend.key.size = unit(0.7,"cm"), 
        legend.title=element_text(size=13, face = "bold")) 
ZDS_data_VS_DPC_PCoA_W_UniFrac

tiff('ZDS data - Vaccination Status by Days Post Inoculation PCoA - Weighted Unifrac.tiff', 
     units="in", width=10, height=6, res=300)
ZDS_data_VS_DPC_PCoA_W_UniFrac
dev.off()

ggsave("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Manuscripts/Bearson_Vacc_vs_DSC_Enteritidis_Challenge/Figures/Figure 5B - Weighted UniFrac PCoA VS by DPI.tiff",
       ZDS_data_VS_DPC_PCoA_W_UniFrac,
       device = "tiff", width = 10, height = 6)

####PCoA Analysis - Weighted Unifrac - Axes 1 & 3 - Supplementary Figure 1B####
ZDS_data_VS_DPC_PCoA_W_UniFrac_Axes_1_3 <- 
  plot_ordination(ZDS_data_tree, ordinate(ZDS_data_tree, "PCoA", "wunifrac"), 
                  shape = "Vaccination_Status", color = "Days_Post_Challenge", axes = c(1,3))+ 
  scale_color_manual(name="Days Post Inoculation", 
                     breaks = c("D0", "D7", "D14"),
                     labels=c("0", "7", "14"),
                     values=c("purple", "red", "orange"))+
  scale_shape_discrete("Vaccination Status",
                       labels=c("Mock","BBS 1134")) +
  geom_point(size=7)+
  ggtitle("Weighted UniFrac") +
  theme(plot.title = element_text(hjust = 0.5, size=15)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=15, face = "bold")) +
  guides(color = guide_legend(override.aes = list(shape = c(15, 15, 15))))+
  theme(axis.text=element_text(size=14, face = "bold", color = "black"), axis.title=element_text(size=15, face = "bold"), 
        legend.text=element_text(size=12, face = "bold", color = "black"), legend.key.size = unit(0.7,"cm"), 
        legend.title=element_text(size=13, face = "bold")) 
ZDS_data_VS_DPC_PCoA_W_UniFrac_Axes_1_3

tiff('ZDS data - Vaccination Status by Days Post Inoculation PCoA - Weighted Unifrac - Axes 1_3.tiff', 
     units="in", width=10, height=6, res=300)
ZDS_data_VS_DPC_PCoA_W_UniFrac_Axes_1_3
dev.off()

ggsave("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Manuscripts/Bearson_Vacc_vs_DSC_Enteritidis_Challenge/Figures/Supplementary Figure 1B - Weighted UniFrac PCoA VS by DPI - Axes 1 and 3.tiff",
       ANCOMBC2_dot_Plot_VS_DPC_RA,
       device = "tiff", width = 10, height = 6)

####Joint Beta Diversity Graphics - Figure 5 and Supplementary Figure 1####

#Make a joint figure of 50 ns simulations results 
Figure_Beta_Diversity_PCoAs <- ggarrange(ZDS_data_VS_DPC_PCoA_BC,
                                         ZDS_data_VS_DPC_PCoA_W_UniFrac,
                                         labels = c("A", "B"),
                                         ncol = 2, nrow = 1,
                                         common.legend = TRUE,
                                         legend = "right")
Figure_Beta_Diversity_PCoAs

ggsave("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Manuscripts/Bearson_Vacc_vs_DSC_Enteritidis_Challenge/Figures/Figure 5 - PCoAs VS by DPI.tiff",
       Figure_Beta_Diversity_PCoAs,
       device = "tiff", width = 10, height = 6, bg="white")

#Make a joint figure of 50 ns simulations results 
Figure_Beta_Diversity_PCoAs_Axes_1_3 <- ggarrange(ZDS_data_VS_DPC_PCoA_BC_Axes_1_3,
                                                  ZDS_data_VS_DPC_PCoA_W_UniFrac_Axes_1_3,
                                                  labels = c("A", "B"),
                                                  ncol = 2, nrow = 1,
                                                  common.legend = TRUE,
                                                  legend = "right")
Figure_Beta_Diversity_PCoAs_Axes_1_3


ggsave("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Manuscripts/Bearson_Vacc_vs_DSC_Enteritidis_Challenge/Figures/Supplementary Figure 1 - PCoAs VS by DPI - Axes 1 and 3.tiff",
       Figure_Beta_Diversity_PCoAs_Axes_1_3,
       device = "tiff", width = 10, height = 6)

####Taxonomic Level Summary Statistics####

#Get abundance in relative percentage#
phy <- transform_sample_counts(Fdata, function(x) 100*x/sum(x))

##Create table with all Phyla##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Phylum', NArm = FALSE)
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)
#remove unnessary columns
dat <- subset(dat, select=c(Abundance, Phylum))
#Summarize based upon target parameter
Phylum_Summary <- summarySE(data=dat, measurevar="Abundance", groupvars="Phylum", na.rm = TRUE)

##Create table with all Genera##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus', NArm = FALSE)
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine all taxonomy columns#
dat$Taxonomy <- paste(dat$Kingdom,dat$Phylum,dat$Class,dat$Order,dat$Family,dat$Genus,sep=":")
# convert Genus to a character vector from a factor because R
dat$Taxonomy <- as.character(dat$Taxonomy)
#remove unnessary columns
dat <- subset(dat, select=c(Abundance, Taxonomy))
#Summarize based upon target parameter
Genus_Summary <- summarySE(data=dat, measurevar="Abundance", groupvars="Taxonomy", na.rm = TRUE)

#Add taxonomic information relevant to differential abundance anlysis
Genus_Summary <- 
  merge(dplyr::select(Uniq_data_taxa_genus_level, Taxonomy, Family_Genus, name_col), Genus_Summary)

#Split Taxonomy into all of the taxonomic levels
Genus_Summary <-
  Genus_Summary %>%
  separate(Taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=":", remove=FALSE)

#Add an Uknown Genus row
Genus_Summary_ANCOMBC_merging <-
  rbind(Genus_Summary, c("Bacteria:Mixed:Mixed:Mixed:Mixed:Unknown", "Bacteria", "Mixed", "Mixed", "Mixed", "Mixed", 
                         "Unknown", "Mixed_Unknown", "Mixed:Mixed", NA, NA, NA, NA, NA))

####Determine Fdata Means####

##Create tabled based upon means of all Genera in Fdata##
glom <- tax_glom(phy, taxrank = 'Genus', NArm = FALSE)
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine all taxonomy columns#
dat<- unite(dat, Taxonomy, Kingdom:Genus, sep=';', remove = FALSE)
# convert to a character vector from a factor because R
dat$Taxonomy <- as.character(dat$Taxonomy)
#Combine Family and Genus columns#
dat<- unite(dat, Family_Genus, Family:Genus, sep=';', remove = FALSE)
# convert to a character vector from a factor because R
dat$Family_Genus <- as.character(dat$Family_Genus)

# group dataframe by Taxonomy, calculate mean rel. abundance
Fdata_RA_means_Taxonomy <- plyr::ddply(dat, ~Taxonomy, function(x) c(mean=mean(x$Abundance)))
#105

# find Genus whose rel. abund. is less than 1%
Fdata_RA_Other_LTET_1_Taxonomy <- Fdata_RA_means_Taxonomy[Fdata_RA_means_Taxonomy$mean <= 1,]$Taxonomy
#85

# group dataframe by Family_Genus, calculate mean rel. abundance
Fdata_RA_means_Family_Genus <- plyr::ddply(dat, ~Family_Genus, function(x) c(mean=mean(x$Abundance)))
#99 - 7 taxonomies were NA for both Family and Genus and are combined into NA/NA

# find Genus whose rel. abund. is less than 1%
Fdata_RA_Other_LTET_1_Family_Genus <- Fdata_RA_means_Family_Genus[Fdata_RA_means_Family_Genus$mean <= 1,]$Family_Genus
#81

# group dataframe by Family_Genus, calculate mean rel. abundance
Fdata_RA_means_Genus <- plyr::ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
#86 - 20 taxonomies were NA for Genus and are combined into NA

#Change NA to Unknown per ANCOMBC2
Fdata_RA_means_Genus$Genus[is.na(Fdata_RA_means_Genus$Genus)] = "Unknown"


# find Genus whose rel. abund. is less than 1%
Fdata_RA_Other_LTET_1_Genus <- Fdata_RA_means_Genus[Fdata_RA_means_Genus$mean <= 1,]$Genus
#72

####Differentially Abundant Genera###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

library(ANCOMBC)
library(tidyverse)
library(caret)
library(DT)


####Fdata - V514_D0 vs M514_D0 - Genus - Supplementary Table 7####

#Filter for D0 Samples
Fdata_D0 <- 
  subset_samples(Fdata, Days_Post_Challenge=="D0")


#Check number of taxa and samples#
ntaxa(Fdata_D0) #532#
nsamples(Fdata_D0) #23#

#Remove taxonomy exclusive to other samples
Fdata_D0 <- filter_taxa(Fdata_D0, function(x) sum(x) >0, TRUE)
ntaxa(Fdata_D0) 
#397

#Make a TSE from phyloseq
tse_D0 = mia::makeTreeSummarizedExperimentFromPhyloseq(Fdata_D0)

#Run ANCOMBC2
set.seed(123)
D0_output = ancombc2(data = tse_D0, tax_level = "Genus",
                     fix_formula = "Vaccination_Status", rand_formula = NULL,
                     p_adj_method = "holm", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                     group = "Vaccination_Status", struc_zero = TRUE, neg_lb = FALSE,
                     alpha = 0.05, n_cl = 2, verbose = TRUE)

##ANCOMBC2 Test Results

#Extract results
D0_output_res = D0_output$res

#Merge with Genus summary
D0_output_res_genus_summary <- merge(Genus_Summary_ANCOMBC_merging, D0_output_res,
                                     by.x="Genus", by.y="taxon")

#Write out the results
clipr::write_clip(D0_output_res)

####Fdata - V514_D7 vs M514_D7 - Genus - Supplementary Table 7####

#Filter for D7 Samples
Fdata_D7 <- 
  subset_samples(Fdata, Days_Post_Challenge=="D7")


#Check number of taxa and samples#
ntaxa(Fdata_D7) #532#
nsamples(Fdata_D7) #27#

#Remove taxonomy exclusive to other samples
Fdata_D7 <- filter_taxa(Fdata_D7, function(x) sum(x) >0, TRUE)
ntaxa(Fdata_D7) 
#436

#Make a TSE from phyloseq
tse_D7 = mia::makeTreeSummarizedExperimentFromPhyloseq(Fdata_D7)

#Run ANCOMBC2
set.seed(123)
D7_output = ancombc2(data = tse_D7, tax_level = "Genus",
                     fix_formula = "Vaccination_Status", rand_formula = NULL,
                     p_adj_method = "holm", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                     group = "Vaccination_Status", struc_zero = TRUE, neg_lb = FALSE,
                     alpha = 0.05, n_cl = 2, verbose = TRUE)

##ANCOMBC2 Test Results

#Extract results
D7_output_res = D7_output$res

#Merge with Genus summary
D7_output_res_genus_summary <- merge(Genus_Summary_ANCOMBC_merging, D7_output_res,
                                     by.x="Genus", by.y="taxon")

#Write out the results
clipr::write_clip(D7_output_res)

####Fdata - V514_D14 vs M514_D14 - Genus - Supplementary Table 7####

#Filter for D14 Samples
Fdata_D14 <- 
  subset_samples(Fdata, Days_Post_Challenge=="D14")


#Check number of taxa and samples#
ntaxa(Fdata_D14) #532#
nsamples(Fdata_D14) #28#

#Remove taxonomy exclusive to other samples
Fdata_D14 <- filter_taxa(Fdata_D14, function(x) sum(x) >0, TRUE)
ntaxa(Fdata_D14) 
#468

#Make a TSE from phyloseq
tse_D14 = mia::makeTreeSummarizedExperimentFromPhyloseq(Fdata_D14)

#Run ANCOMBC2
set.seed(123)
D14_output = ancombc2(data = tse_D14, tax_level = "Genus",
                      fix_formula = "Vaccination_Status", rand_formula = NULL,
                      p_adj_method = "holm", pseudo_sens = TRUE,
                      prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                      group = "Vaccination_Status", struc_zero = TRUE, neg_lb = FALSE,
                      alpha = 0.05, n_cl = 2, verbose = TRUE)

##ANCOMBC2 Test Results

#Extract results
D14_output_res = D14_output$res

#Merge with Genus summary
D14_output_res_genus_summary <- merge(Genus_Summary_ANCOMBC_merging, D14_output_res,
                                      by.x="Genus", by.y="taxon")

#Write out the results
clipr::write_clip(D14_output_res)


####General Taxonomic Bar Plots Prep####

#Note: These graphics and statistics will be paired with differential abundance analysis,
#Thus they will be run on filtered data (Fdata)

##Sample by Vaccination Status by DPC

##Create table ready for making stacked bar graph for Genuses <1%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus', NArm = FALSE)
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- plyr::ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(VS_DPC, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("VS_DPC","Genus"), na.rm = FALSE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- plyr::ddply(dat, ~VS_DPC, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(VS_DPC, Abundance, Genus))
#combine with original table
Genus_Fdata_VS_DPC <- rbind(dat, Abundance)

##Genus Sample ID

##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus', NArm = FALSE)
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- plyr::ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Sample, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Sample","Genus"), na.rm = FALSE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- plyr::ddply(dat, ~Sample, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Sample, Abundance, Genus))
#combine with original table
Genus_Fdata_Sample <- rbind(dat, Abundance)


####ANCOMBC2 Bar Plots Prep####
#All graphics are exploratory and not included in manuscript
##D0

#Get list of Taxa Category: Taxa Name taxas that were significantly different between Mock and Vaccinated samples
ANCOMBC2_D0_sig_Family_Genus <- 
  filter(D0_output_res_genus_summary, 
         diff_robust_Vaccination_StatusVaccinated==TRUE)$Family_Genus

#Keep only differentially abundant taxa from the unfiltered Fdata by VS/DPC subcategory relative percentages
ANCOMBC2_D0_sig_Full_Taxonomy_Fdata_VS_DPC <- 
  Full_Taxonomy_Fdata_VS_DPC %>%
  filter(Family_Genus %in% ANCOMBC2_D0_sig_Family_Genus) %>%
  filter(VS_DPC %in% c("M_D0", "V_D0"))

##D7

#Get list of Taxa Category: Taxa Name taxas that were significantly different between Mock and Vaccinated samples
ANCOMBC2_D7_sig_Family_Genus <- 
  filter(D7_output_res_genus_summary, 
         diff_Vaccination_StatusVaccinated==TRUE & passed_ss_Vaccination_StatusVaccinated == TRUE)$Family_Genus

#Keep only differentially abundant taxa from the unfiltered Fdata by VS/DPC subcategory relative percentages
ANCOMBC2_D7_sig_Full_Taxonomy_Fdata_VS_DPC <- 
  Full_Taxonomy_Fdata_VS_DPC %>%
  filter(Family_Genus %in% ANCOMBC2_D7_sig_Family_Genus) %>%
  filter(VS_DPC %in% c("M_D7", "V_D7"))

##D14

#Get list of Taxa Category: Taxa Name taxas that were significantly different between Mock and Vaccinated samples
ANCOMBC2_D14_sig_Family_Genus <- 
  filter(D14_output_res_genus_summary, 
         diff_Vaccination_StatusVaccinated==TRUE & passed_ss_Vaccination_StatusVaccinated == TRUE)$Family_Genus

#Keep only differentially abundant taxa from the unfiltered Fdata by VS/DPC subcategory relative percentages
ANCOMBC2_D14_sig_Full_Taxonomy_Fdata_VS_DPC <- 
  Full_Taxonomy_Fdata_VS_DPC %>%
  filter(Family_Genus %in% ANCOMBC2_D14_sig_Family_Genus) %>%
  filter(VS_DPC %in% c("M_D14", "V_D14"))

##All

#Keep all differentially abundant taxa from the unfiltered Fdata by VS/DPC subcategory relative percentages
ANCOMBC2_all_sig_Full_Taxonomy_Fdata_VS_DPC <- 
  Full_Taxonomy_Fdata_VS_DPC %>%
  filter(Family_Genus %in% c(ANCOMBC2_D0_sig_Family_Genus, 
                             ANCOMBC2_D7_sig_Family_Genus,
                             ANCOMBC2_D14_sig_Family_Genus))


####Determine colors for Family_Genus Bar Charts####
#Greater than 1% across all VS/DPC subcategories and ANCOMBC2 taxa

#Determine the length of list of all Family_Genus taxas greater than 1 % across all VS_DPC categories and
#All Family_Genus taxas found to be significant by ANCOMBC2
ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_names <- 
  sort(unique(c(Genus_Fdata_VS_DPC$Genus, 
                ANCOMBC2_all_sig_Full_Taxonomy_Fdata_VS_DPC$Family_Genus)))
ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_names
#29

#polychrome from pals has 36 values, just use the first 29
ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_colors <- unname(pals::polychrome(n=29))

dput(ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_colors)
c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", 
  "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", 
  "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D", 
  "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", 
  "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C")

#Use the list of serotypes to give each color a name
ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_colors <- 
  setNames(ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_colors, ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_names)
ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_colors
[Eubacterium] coprostanoligenes group;Incertae Sedis_234                               Bacteroidaceae;Bacteroides                              Barnesiellaceae;Barnesiella 
"#5A5156"                                                "#E4E1E3"                                                "#F6222E" 
Eggerthellaceae_2;Asaccharobacter_2                             Eggerthellaceae_2;CHKCI002_2             Erysipelatoclostridiaceae;Massiliomicrobiota 
"#FE00FA"                                                "#16FF32"                                                "#3283FE" 
Erysipelatoclostridiaceae;Thomasclavelia                           Erysipelotrichaceae;Holdemania                         Erysipelotrichaceae;Turicibacter 
"#FEAF16"                                                "#B00068"                                                "#1CFFCE" 
Incertae Sedis_119;Incertae Sedis_217          Lachnospiraceae;[Ruminococcus] gauvreauii group                          Lachnospiraceae;Anaerobutyricum 
"#90AD1C"                                                "#2ED9FF"                                                "#DEA0FD" 
Lachnospiraceae;CHKCI001                           Lachnospiraceae;Eisenbergiella                           Lachnospiraceae;Frisingicoccus 
"#AA0DFE"                                                "#F8A19F"                                                "#325A9B" 
Lachnospiraceae;GCA-900066575                        Lachnospiraceae;Lachnoclostridium                       Lachnospiraceae;Mediterraneibacter 
"#C4451C"                                                "#1C8356"                                                "#85660D" 
Lachnospiraceae;NA                               Lachnospiraceae;Sellimonas                               Lachnospiraceae;Tyzzerella 
"#B10DA1"                                                "#FBE426"                                                "#1CBE4F" 
Lactobacillaceae;Ligilactobacillus                                      Lactobacillaceae;NA                                 Monoglobaceae;Monoglobus 
"#FA0087"                                                "#FC1CBF"                                                "#F7E1A0" 
Oscillospiraceae;NA                                        Other Prokaryotes                                 Peptostreptococcaceae;NA 
"#C075A6"                                                "#782AB6"                                                "#AAF400" 
Ruminococcaceae;Faecalibacterium                                       Ruminococcaceae;NA 
"#BDCDFF"                                                "#822E1C" 

#Adjust the colors pending what is happening on the graphs
ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_colors['Other Prokaryotes'] = 'rosybrown'
ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_colors['Lactobacillaceae;Ligilactobacillus'] = '#FA0087'


####General Taxonomic Bar plot####

Genus_Fdata_VS_DPC_Family_genus_Order <- sort(unique(Genus_Fdata_VS_DPC$Genus))

# Find the index of the element to move (e.g., "b")
index_to_move <- which(Genus_Fdata_VS_DPC_Family_genus_Order == "Other Prokaryotes")  # Returns 2

# Remove the element from its current position
Genus_Fdata_VS_DPC_Family_genus_Order <- Genus_Fdata_VS_DPC_Family_genus_Order[-index_to_move]

# Move the element to the end using append()
Genus_Fdata_VS_DPC_Family_genus_Order <- append(Genus_Fdata_VS_DPC_Family_genus_Order, "Other Prokaryotes")

#Genera greater than 1% across all VS/DPC subcategories Only
spatial_plot_Genus_Fdata_VS_DPC_FP <- 
  ggplot(data=Genus_Fdata_VS_DPC, aes(x=VS_DPC, y=Abundance, 
                                      fill=factor(Genus, levels=Genus_Fdata_VS_DPC_Family_genus_Order))) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_colors)+
  scale_x_discrete("Mock Vaccination by Days Post Challenge",
                   limits=c("M_D0", "M_D7", "M_D14", "V_D0", "V_D7", "V_D14"),
                   labels=c("Mock <br> 0 DPC", "Mock <br> 7 DPC", "Mock <br> 14 DPC",
                            "Vaccinated <br> 0 DPC", "Vaccinated <br> 7 DPC", "Vaccinated <br> 14 DPC"))+
  ylab("Relative Percentage") +
  #theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_markdown(size=13), axis.title=element_text(size=15), 
        legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), 
        legend.title=element_text(size=13))
spatial_plot_Genus_Fdata_VS_DPC_FP

tiff('Genera by Vaccination Status by DPC Bar Plots.tiff', units="in", width=12, height=6, res=300)
spatial_plot_Genus_Fdata_VS_DPC_FP
dev.off()

####Greater than 1% across all VS/DPC subcategories and ANCOMBC2 Relative Abundance Stacked Bar Graph - Figure 6 Prep####
##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus', NArm = FALSE)
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- plyr::ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
#Remove ANCOMBC2 DAGs
Other <- setdiff(Other, ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_names)
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(VS_DPC, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("VS_DPC","Genus"), na.rm = FALSE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- plyr::ddply(dat, ~VS_DPC, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(VS_DPC, Abundance, Genus))
#combine with original table
Genus_Fdata_VS_DPC_ANCOMBC2 <- rbind(dat, Abundance)


#Make a wide version for reference
wide_Genus_Fdata_VS_DPC_ANCOMBC2 <- pivot_wider(Genus_Fdata_VS_DPC_ANCOMBC2, 
                                                names_from = VS_DPC, values_from = Abundance)

#Get the alphabetical genus order
Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order <- sort(unique(Genus_Fdata_VS_DPC_ANCOMBC2$Genus))

# Find the index of the element to move (e.g., "b")
index_to_move <- which(Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order == "Other Prokaryotes")  # Returns 2

# Remove the element from its current position
Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order <- Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order[-index_to_move]

# Move the element to the end using append()
Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order <- append(Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order, "Other Prokaryotes")

#Add a * to the end of every genus found by ANCOMBC2 to be differentially abundant at any dpi
Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order_V2 <- c()

for(genus in Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order){
  if (genus %in% unique(ANCOMBC2_all_sig_Full_Taxonomy_Fdata_VS_DPC$Family_Genus)){
    genus_V2 <- paste0(genus,"*") 
  }
  else{
    genus_V2 <- genus
  }
  Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order_V2 <-
    append(Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order_V2, genus_V2)
}

#Determine the length of list of all Family_Genus taxas greater than 1 % across all VS_DPC categories and
#All Family_Genus taxas found to be significant by ANCOMBC2
ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_names <- 
  sort(unique(c(Genus_Fdata_VS_DPC$Genus, 
                ANCOMBC2_all_sig_Full_Taxonomy_Fdata_VS_DPC$Family_Genus)))
ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_names
#29

#Create alternate version with just Genus bolded
dput(ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_names)
Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order_V2_bold_Genus <- c("[Eubacterium] coprostanoligenes group;**Incertae Sedis_234***", 
                                                                  "Bacteroidaceae;**Bacteroides***", "Barnesiellaceae;**Barnesiella**", 
                                                                  "Eggerthellaceae_2;**Asaccharobacter_2***", "Eggerthellaceae_2;**CHKCI002_2***", 
                                                                  "Erysipelatoclostridiaceae;**Massiliomicrobiota**", "Erysipelatoclostridiaceae;**Thomasclavelia**", 
                                                                  "Erysipelotrichaceae;**Holdemania***", "Erysipelotrichaceae;**Turicibacter***", 
                                                                  "Incertae Sedis_119;**Incertae Sedis_217**", "Lachnospiraceae;**[Ruminococcus] gauvreauii group***", 
                                                                  "Lachnospiraceae;**Anaerobutyricum***", "Lachnospiraceae;**CHKCI001**", 
                                                                  "Lachnospiraceae;**Eisenbergiella**", "Lachnospiraceae;**Frisingicoccus***", 
                                                                  "Lachnospiraceae;**GCA-900066575***", "Lachnospiraceae;**Lachnoclostridium***", 
                                                                  "Lachnospiraceae;**Mediterraneibacter**", "Lachnospiraceae;**NA**", "Lachnospiraceae;**Sellimonas**", 
                                                                  "Lachnospiraceae;**Tyzzerella***", "Lactobacillaceae;**Ligilactobacillus**", 
                                                                  "Lactobacillaceae;**NA**", "Monoglobaceae;**Monoglobus***", "Oscillospiraceae;**NA**", 
                                                                  "Peptostreptococcaceae;**NA**", "Ruminococcaceae;**Faecalibacterium***", 
                                                                  "Ruminococcaceae;**NA**", "Other Prokaryotes")

####Greater than 1% across all VS/DPC subcategories and ANCOMBC2 Relative Abundance Stacked Bar Graph - Figure 6####
spatial_plot_Genus_Fdata_VS_DPC_ANCOMBC2_FP <- 
  ggplot(data=Genus_Fdata_VS_DPC_ANCOMBC2, aes(x=VS_DPC, y=Abundance, 
                                               fill=factor(Genus, levels=Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order))) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=ANCOMBC2_and_GT_1_Genus_Fdata_VS_DPC_colors,
                    labels=Genus_Fdata_VS_DPC_ANCOMBC2_Family_genus_Order_V2)+
  scale_x_discrete("Vaccination Status by Days Post Inoculation",
                   limits=c("M_D0", "M_D7", "M_D14", "V_D0", "V_D7", "V_D14"),
                   labels=c("Mock <br> 0 dpi", "Mock <br> 7 dpi", "Mock <br> 14 dpi",
                            "BBS 1134 <br> 0 dpi", "BBS 1134 <br> 7 dpi", "BBS 1134 <br> 14 dpi"))+
  ylab("Relative Percentage") +
  #theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() +
  theme(axis.text=element_markdown(size=14, face = "bold", color = "black"), axis.title=element_markdown(size=15, face = "bold"), 
        legend.text=element_markdown(size=12, face = "bold", color = "black"), legend.key.size = unit(0.7,"cm"), 
        legend.title=element_markdown(size=13, face = "bold")) 

spatial_plot_Genus_Fdata_VS_DPC_ANCOMBC2_FP

tiff('Genera (GT 1 and DAGs) by Vaccination Status by dpi Bar Plots.tiff', units="in", width=14, height=6, res=300)
spatial_plot_Genus_Fdata_VS_DPC_ANCOMBC2_FP
dev.off()

ggsave("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Manuscripts/Bearson_Vacc_vs_DSC_Enteritidis_Challenge/Figures/Figure 6 - Genera Relative Percentage VS by DPI.tiff",
       spatial_plot_Genus_Fdata_VS_DPC_ANCOMBC2_FP,
       device = "tiff", width = 16, height = 6)

####Relative Percentages Table for all Genera - Supplementary Table 6####

# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus', NArm = FALSE)
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Family_Genus, Family:Genus, sep=';')
# convert Family_Genus to a character vector from a factor because R
dat$Family_Genus <- as.character(dat$Family_Genus)
#remove unnessary columns
dat <- subset(dat, select=c(VS_DPC, Abundance, Family_Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("VS_DPC","Family_Genus"), na.rm = FALSE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(VS_DPC, Abundance, Family_Genus))
#Pivot wider by VS_DPC
wide_Family_Genus_Fdata_VS_DPC <- pivot_wider(dat, names_from = VS_DPC, values_from = Abundance)

#Add the means
wide_Family_Genus_Fdata_VS_DPC <-
  merge(wide_Family_Genus_Fdata_VS_DPC, Fdata_RA_means_Family_Genus)

#Add a column to designate it as Other Prokaryote or not
wide_Family_Genus_Fdata_VS_DPC$Other_Prokaryote <-
  ifelse(wide_Family_Genus_Fdata_VS_DPC$Family_Genus %in% Fdata_RA_Other_LTET_1_Family_Genus, TRUE, FALSE)

#Sanity Check
table(wide_Family_Genus_Fdata_VS_DPC$Other_Prokaryote)
# FALSE  TRUE 
# 18    81 
# TRUE value matches length of Fdata_RA_Other_LTET_1_Family_Genus


#Add column to relative abundance table denoting ANCOMBC2 identified genera
wide_Family_Genus_Fdata_VS_DPC$ANCOMBC2 <-
  ifelse(wide_Family_Genus_Fdata_VS_DPC$Family_Genus %in% unique(ANCOMBC2_all_sig_Full_Taxonomy_Fdata_VS_DPC$Family_Genus), "TRUE", "FALSE")

#Sanity Check
table(wide_Family_Genus_Fdata_VS_DPC$ANCOMBC2)
# FALSE  TRUE 
# 85    14 

#Sent to Excel to create Supplementary Table 6
clipr::write_clip(wide_Family_Genus_Fdata_VS_DPC)

#Create a version focused on bar graph visualized genera
wide_Family_Genus_Fdata_VS_DPC_bar_graph <- 
  wide_Family_Genus_Fdata_VS_DPC %>%
  filter(Other_Prokaryote == FALSE)

####ANCOMBC2 LFC Dot Plot - Figure 7 General Prep####

##Get the significant ANCOMBC Results

#DO
D0_output_sig_res_genus_summary <-
  D0_output_res_genus_summary %>%
  filter(diff_Vaccination_StatusVaccinated == TRUE &
           passed_ss_Vaccination_StatusVaccinated == TRUE) %>%
  select(name_col, Family_Genus, lfc_Vaccination_StatusVaccinated)
#7

#Add a column for dpi
D0_output_sig_res_genus_summary$DPI <- "D0"

#D7
D7_output_sig_res_genus_summary <-
  D7_output_res_genus_summary %>%
  filter(diff_Vaccination_StatusVaccinated == TRUE &
           passed_ss_Vaccination_StatusVaccinated == TRUE) %>%
  select(name_col, Family_Genus, lfc_Vaccination_StatusVaccinated)
#5

#Add a column for dpi
D7_output_sig_res_genus_summary$DPI <- "D7"


#D14
D14_output_sig_res_genus_summary <-
  D14_output_res_genus_summary %>%
  filter(diff_Vaccination_StatusVaccinated == TRUE &
           passed_ss_Vaccination_StatusVaccinated == TRUE) %>%
  select(name_col, Family_Genus, lfc_Vaccination_StatusVaccinated)
#5

#Add a column for dpi
D14_output_sig_res_genus_summary$DPI <- "D14"
#4

#Row bind all dpis
sig_res_genus_summary <-
  rbind(D0_output_sig_res_genus_summary, 
        D7_output_sig_res_genus_summary,
        D14_output_sig_res_genus_summary)

#Create a Vaccination Status column
sig_res_genus_summary$Vaccination_Status <-
  ifelse(sig_res_genus_summary$lfc_Vaccination_StatusVaccinated >0, "Vaccinated", "Mock")


#Pivot the greater than 1% and ANCOMBC2 genera relative abundance data frame to long format
Genus_Fdata_VS_DPC_ANCOMBC2_V2 <-
  Genus_Fdata_VS_DPC_ANCOMBC2 %>%
  separate(col = VS_DPC, into = c("Vaccination_Status_Short", "DPI"))

#Create a vaccination status column that matches sig_res_genus_summary
Genus_Fdata_VS_DPC_ANCOMBC2_V2$Vaccination_Status <-
  ifelse(Genus_Fdata_VS_DPC_ANCOMBC2_V2$Vaccination_Status_Short == "V", "Vaccinated", "Mock")


####ANCOMBC2 LFC Dot Plot - Figure 7 Relative Abundance Size Prep####
#Pivot the greater than 1% and ANCOMBC2 genera relative abundance data frame to long format
Genus_Fdata_VS_DPC_ANCOMBC2_V2 <-
  Genus_Fdata_VS_DPC_ANCOMBC2 %>%
  separate(col = VS_DPC, into = c("Vaccination_Status_Short", "DPI"))

#Create a vaccination status column that matches sig_res_genus_summary
Genus_Fdata_VS_DPC_ANCOMBC2_V2$Vaccination_Status <-
  ifelse(Genus_Fdata_VS_DPC_ANCOMBC2_V2$Vaccination_Status_Short == "V", "Vaccinated", "Mock")

#Select only needed columns
Genus_Fdata_VS_DPC_ANCOMBC2_V2 <-
  select(Genus_Fdata_VS_DPC_ANCOMBC2_V2, -Vaccination_Status_Short)

#Change column names
colnames(Genus_Fdata_VS_DPC_ANCOMBC2_V2) <- c("DPI", "Abundance", "Family_Genus", "Vaccination_Status")

#Merge with sig_res_genus_summary
sig_res_genus_summary_RA <-
  merge(sig_res_genus_summary, 
        Genus_Fdata_VS_DPC_ANCOMBC2_V2)

#Put rows in order by relative percentage size
Fdata_RA_means_Family_Genus_ANCOMBC2 <- filter(Fdata_RA_means_Family_Genus,
                                               Family_Genus %in% unique(sig_res_genus_summary_RA$Family_Genus))

Fdata_RA_means_Family_Genus_ANCOMBC2 <- Fdata_RA_means_Family_Genus_ANCOMBC2[order(Fdata_RA_means_Family_Genus_ANCOMBC2$mean),]

sig_res_genus_summary_RA$Family_Genus <-factor(sig_res_genus_summary_RA$Family_Genus, 
                                               levels=Fdata_RA_means_Family_Genus_ANCOMBC2$Family_Genus)

####ANCOMBC2 LFC Dot Plot - Figure 7####
ANCOMBC2_dot_Plot_VS_DPC_RA <-
  ggplot(sig_res_genus_summary_RA, aes(lfc_Vaccination_StatusVaccinated, Family_Genus)) +
  geom_point(aes(color=DPI, shape = Vaccination_Status, size=Abundance)) +
  scale_shape_discrete("Vaccination Status",
                       label=c("Mock", "BBS 1134")) +
  scale_color_manual("Days Post Inoculation",
                     limits=c("D0", "D7", "D14"),
                     labels=c("0", "7", "14"),
                     values=c("purple", "red", "orange")) +
  scale_size_continuous("Group Abundance (%)", range= c(2, 6), breaks = c(1, 5, 10)) +
  xlab("ANCOMBC2 Log2 Fold Change") +
  ylab("Family; Genus") +
  geom_vline(xintercept=0) +
  theme_bw()+
  guides(color = guide_legend(order = 2, override.aes = list(shape = 15, size=5)))+
  guides(size = guide_legend(order=3, override.aes = list(shape = 15)))+
  guides(shape = guide_legend(order =1, override.aes = list(size=5)))+
  theme(axis.text=element_text(size=14, face = "bold", color = "black"), axis.title=element_text(size=15, face = "bold"), 
        legend.text=element_text(size=12, face = "bold", color = "black"), legend.key.size = unit(0.7,"cm"), 
        legend.title=element_text(size=13, face = "bold")) 
ANCOMBC2_dot_Plot_VS_DPC_RA

ggsave("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Manuscripts/Bearson_Vacc_vs_DSC_Enteritidis_Challenge/Figures/Figure 7 - ANCOMBC2 Log2 Fold Change Dot Plot.tiff",
       ANCOMBC2_dot_Plot_VS_DPC_RA,
       device = "tiff", width = 10, height = 6)
