rm(list = ls())

#### PRELIMINARIES ############################################################################################# 

#*Uploads the needed libraries --------------------------------------------------------------------------------

#require(DESeq2)

require(ggplot2)

require(data.table)

require(usedist)

require(DT)

require(R2HTML)

#*Set the number of significant digits for the output --------------------------
sig_dig = 4

#*Set the working directory ---------------------------------------------------------------------------------

#get the default wd
default_wd <- getwd()

#Set the directory where all the saved outputs will be stored
#setwd("/media/davide/Dropbox/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Anne_Claire_Project/Analysis/Downstream_analysis/Finding_correlations/")# <--- insert here the path to the working (output) directory
setwd("~/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Stefania_Carobbio_project/First_dataset/Analysis_old/Downstream_analysis/Finding_correlations/")# <--- insert here the path to the working (output) directory

new_wd <- getwd()

#### DATA UPLOAD ###############################################################################################################

annotation_table <- read.csv("~/Dropbox (Cambridge University)/Bioinformatics_core/Useful_bash_scripts_and_references/Homo_sapiens.GRCh38.95_gene_annotation_table.txt", sep = "\t")
#annotation_table_entrez <- read.csv("~/Dropbox (Cambridge University)/Bioinformatics_core/Useful_bash_scripts_and_references/", sep = "\t")
#annotation_table_entrez <- read.csv("/media/davide/Dropbox/Dropbox (Cambridge University)/Bioinformatics_core/Useful_bash_scripts_and_references/Mus_musculus.NCBI_ENTREZ_gene_annotation_table.txt", sep = "\t")
#annotation_table <- read.csv("/media/davide/Dropbox/Dropbox (Cambridge University)/Bioinformatics_core/Useful_bash_scripts_and_references/Homo_sapiens.GRCh38.95_gene_annotation_table.txt", sep = "\t")

# inputs the list of the count files
#input_files <- list.files(path = "/media/davide/Dropbox/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Stefania_Carobbio_project/First_dataset/Analysis_old/Counts/Pluripotent/", pattern = "*.tab", full.names = TRUE) #<--- insert here the path to the working directory
input_files <- list.files(path = "~/Dropbox (Cambridge University)/Bioinformatics_core/Projects/Stefania_Carobbio_project/First_dataset/Analysis_old/Counts/Pluripotent/", pattern = "*.tab", full.names = TRUE) #<--- insert here the path to the working directory

#Reads the count files (input_files)
#create a list; each element of a list (named ''sample'' in the following) is a count table
counts_tables <- lapply(input_files, fread, header = FALSE, stringsAsFactors = FALSE, col.names=c("GeneID","Counts"))

#takes the names of the genes from the counts_tables
genes_names <- counts_tables[[1]]$GeneID
genes_number = length(genes_names)

### BUILDS THE COUNTS TABLES DATAFRAME ###############################################################################################################
# Each column of the counts tables dataframe contain a replicate of the various samples; Each row contain a gene;
# The raw counts for each replicate and for each gene are reported

#take the names of each element (sample) of the list from the input files 
samples_names <- substr(noquote(unlist(lapply(basename(input_files), tools::file_path_sans_ext))) ,1,5)

#assign the names to the elements of the counts_tables list, composed by the counts tables; NOTE: each replicate has its ID
names(counts_tables) <- samples_names

#Creates a single data frame with all the samples as columns, for reporting and clustering purposes -- see heatmaps below
counts_tables_dataframe <- sapply(counts_tables, '[[', 2)

#assigns the genes names to the rows of the counts_tables_dataframe
rownames(counts_tables_dataframe) <- genes_names

#Cuts off the last 5 rows of the counts_tables_dataframe; the columns are taken from feature counts outputs and the last 5 rows don't contain any gene
counts_tables_dataframe<- head(counts_tables_dataframe, n=nrow(counts_tables_dataframe) -5 )

### BUILDS THE COUNTS TABLES DATAFRAME AVERAGE ###############################################################################################################
# It is similar to the counts_tables_dataframe but the columns contain the averages over the replicates

# extracts the root of the replicates names, obtaining one name for each sample.
samples_names_unique <- unique(substr(samples_names,1,3))

#initializes the matrix
counts_tables_dataframe_average <- matrix(ncol=length(samples_names_unique), nrow=nrow(counts_tables_dataframe))

# Builds the matrix by columns
for(i in 1:length(samples_names_unique)){
  submatrix<- counts_tables_dataframe[, grep(colnames(counts_tables_dataframe), pattern = samples_names_unique[i])]
  dummy <- (apply(submatrix,1,mean))
  counts_tables_dataframe_average[,i]<- dummy
}
colnames(counts_tables_dataframe_average) <- samples_names_unique
rownames(counts_tables_dataframe_average) <- head(genes_names,n=length(genes_names)-5 )

# FILTERS THE COUNTS TABLES DATAFRAME AVERAGE ###############################################################################################################
# Filters for low counts

CPM <- apply(counts_tables_dataframe_average,2, function(x){x*10^6/sum(x)})

count_zeroes<- function(x){length(which(x<=1.5))}
numberofzeroes <- apply(CPM, 1, count_zeroes)
CPM_filtered <- CPM[which(numberofzeroes<=2),]
CPM_filtered_sorted<- CPM_filtered[order(CPM_filtered[,1], decreasing = TRUE),]

# COMPUTES THE DISTANCES MATRIX ###############################################################################################################

# Standardises, according to different criteria, the CPM_filteres matrix tandardises the  
#
standardiser <- function(x) (x-mean(x))/sd(x)
#
CPM_filtered_standardised_by_columns <- (apply(CPM_filtered, 2, standardiser))
#
CPM_filtered_standardised_by_rows <- t(apply(CPM_filtered, 1, standardiser))
#
CPM_filtered_logged_standardised_by_rows <- t(apply(log(CPM_filtered), 1, standardiser))
#
CPM_filtered_standardised_by_columns_rows <- t(apply(CPM_filtered_standardised_by_columns, 1, standardiser))

# Computes the euclidean distance
#
distances_logged_standardised_by_rows<- as.matrix(dist((CPM_filtered_logged_standardised_by_rows)))
#write.csv(distances, file="distances_matrix.csv")

# Computes the ''absolute distance'' with a customised formula; 
# in this case both the perfectly correlated or anticorrelated trajectories
# will result in a distance =0 . Very ''different'' trajectories will have 
# distances close to 1
#
CPM_filtered_logged_standardised_by_rows <- t(CPM_filtered_logged_standardised_by_rows) # note: the rows become columns here !
correlation_matrix <- cor(CPM_filtered_logged_standardised_by_rows) # computes the covariance between the columns; 
absolute_distance <- 1-abs(correlation_matrix)                      # being the columns expressed as Z scores (standardised) the covariance is the correlation coefficient

#abs_euclid <- function(x,y) 1-(abs(cor(x,y, method="pearson")))
#aaaa<- dist_make(distances_logged_standardised_by_rows, abs_euclid) 
#
# Uploads the distance matrix computed as above: the coputation is intensive and the result has been saved for convenience
#absolute_distance<- data.frame(fread("./absolute_distance.csv"), row.names = 1)

# FINDS THE INFLUENCERS ###############################################################################################################
# Finds the genes that have more other genes at a closer distance

# Computes the influencers according to the euclidean distance
#
influencer<- function(x) length(which(x<=1|x>=5.2))
#
influencers_1- apply(distances_logged_standardised_by_rows, 2,influencer) 
#
hist(ttt, breaks=200)

# Computes the influencers according to the absolute distance
#
influencer<- function(x) length(which(x<=0.05))
#
influencer_2 <- apply(absolute_distance, 2, influencer)
influencer_4 <- apply(absolute_distance, 2, influencer)

write.table(names(influencer_2[which(influencer_2>11000)]), quote=FALSE, row.names = FALSE )

write.csv(influencer_2, "influencers_list")


#hist(distances_logged_standardised_by_rows[,"PPARA"], ylim=c(0, 300), breaks = 100)

hist(absolute_distance[,"PPARA"], breaks = 500)

hist(distances_logged_standardised_by_rows[,"PPARD"], ylim=c(0, 300), breaks = 100)

hist(absolute_distance[,"PPARD"], breaks=500)

hist(absolute_distance[,"SLC38A1"],breaks = 500 )

hist(absolute_distance[,"ZNF490"], )

hist(absolute_distance[,"TXNDC5"], )

write.table(names(absolute_distance[,"PPARA"])[1:100], row.names = F, quote = F)

#write.table(names(, row.names = F, quote = F)



plot(CPM_filtered_standardised_by_columns["ZNF490",], type="o", col="red")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["ZFP28",], type="o", col="blue")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["CCDC61",], type="o", col="green")
#
plot(CPM_filtered_standardised_by_columns["SOX2",], type="o", col="red")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["SNX16",], type="o", col="blue")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["LAMB1",], type="o", col="blue")
#
plot(CPM_filtered_standardised_by_columns["PPARA",], type="o", col="red")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["THAP2",], type="o", col="blue")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["P4HA1",], type="o", col="green")
#
#
plot(CPM_filtered_standardised_by_columns["PPARG",], type="o", col="red")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["GRASP",], type="o", col="blue")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["CLUH",], type="o", col="green")
#
plot(CPM_filtered_standardised_by_columns["ZIC1",], type="o", col="red")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["RIN3",], type="o", col="blue")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["BMP7",], type="o", col="green")

write.table(names(ttt[which(ttt>1000&ttt<1100)]), quote=FALSE, row.names = FALSE )

plot(CPM_filtered["PPARA",], type="o", col="red", ylim=c(0,80))
par(new=TRUE)
plot(CPM_filtered["UPF3A",], type="o", col="blue", ylim=c(0,80))

plot(CPM_filtered["PAX3",], type="o", col="red", ylim=c(0,80))
par(new=TRUE)
plot(CPM_filtered["PDZD4",], type="o", col="blue", ylim=c(0,80))

plot(CPM_filtered_standardised_by_columns["SOX2",], type="o", col="red", ylim=c(0,80))
par(new=TRUE)
plot(CPM_filtered["SNX16",], type="o", col="blue", ylim=c(0,80))



plot(counts_tables_dataframe_average["PPARA",], type="o", col="red", ylim=c(100, 1000))
par(new=TRUE)
plot(counts_tables_dataframe_average["UPF3A",], type="o", col="blue", ylim=c(100, 1000))


max(ttt)



genes_list <- c("UCP1", "PDGRFRA", "DIO2", "EBF2", "PRDM16", "PGC1", "ADRP", "FATP4", "PPARG", "PPARG2", "TBOX", "PAX3", "PAX7", "OCT4", "PAT2","P2RX5","PREF1")
genes_list_1<-c("PPARA", "ZIC1", "SOX2", "DLK1","PEL1", "PAX3")

genes_list[which(genes_list %in% names(influencer_2))]

genes_list_1[which(genes_list_1 %in% names(influencer_2))]


distances_for_genes_of_interest <- t(absolute_distance[c("PPARA", "PPARG", "PAX3", "PAX7", "ZIC1", "SOX2", "DLK1"),])

plot(CPM_filtered_standardised_by_columns["PPARA",], type="o", col="red")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["THAP2",], type="o", col="blue")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["TRIQK",], type="o", col="green")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["LEPR",], type="o", col="yellow")
par(new=TRUE)
plot(CPM_filtered_standardised_by_columns["CAMK2A",], type="o", col="purple")

standardised_rows_of_interest <- data.frame(cbind(CPM_filtered_standardised_by_rows["PPARA",], CPM_filtered_standardised_by_rows["THAP2",], CPM_filtered_standardised_by_rows["TRIQK",], CPM_filtered_standardised_by_rows["LEPR",], CPM_filtered_standardised_by_rows["CAMK2A",]))
colnames(standardised_rows_of_interest) <-c("PPARA", "THAP2", "TRIQK", "LEPR", "CAMK2A")

x<- seq(1:8)

ggplot(standardised_rows_of_interest, aes(x)) +                    # basic graphical object
  geom_line(aes(y=standardised_rows_of_interest$PPARA), colour="red") +  # first layer
  geom_line(aes(y=standardised_rows_of_interest$THAP2), colour="black")+  # second layer
  geom_line(aes(y=standardised_rows_of_interest$TRIQK), colour="purple")+  # second layer
  geom_line(aes(y=standardised_rows_of_interest$LEPR), colour="yellow")+  # second layer
  geom_line(aes(y=standardised_rows_of_interest$CAMK2A), colour="blue") +
  scale_color_discrete(name = "Genes", labels=colnames(standardised_rows_of_interest)) # second layer
  
  
  


