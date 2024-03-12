################### START of PE DADA2 protocol for V3V4 2016 Shrimp samples ###################
UPDATE 2020-03-03: Once more, I used dada2 to use its cleaning filters so I took the cutadapt-treated samples (no primers nor adapters), this time using the 2016 V3V4 samples.
cd /home/dada2/
  mv Dada2_tutorial Dada2_tutorial_ok
mkdir -p 2016_shrimp_V3V4/01_import_files/02_paired
cd 2016_shrimp_V3V4/
  cp -R /home/luigui/a3/Shrimp_2015-2018/V3-V4_comparison/03_cutadapt_trim16S/cutadapt-2016_V3V4_* 01_import_files/02_paired
# Use the same name as before
rename 's/cutadapt-2016_//' 01_import_files/02_paired/*
  
  ### 1.- Loading dada and input file information ###
  # R --no-save # if requiered to initialize R with savind disabled
  library(dada2)
PE.path <- "/home/dada2/2016_shrimp_V3V4/01_import_files/02_paired" #Set the path containing the input single-end files (joined). Full paths are better to avoid issues
setwd("/home/dada2/2016_shrimp_V3V4") # set the working directory manually
getwd() # Check current directory
# [1] "/home/dada2/Dada2_tutorial"
PE.files_full <- list.files(PE.path, pattern="fastq.gz", full.names = TRUE) #Get full paths and filenames for the input data folder
PE.files_full_R1 <- sort(list.files(PE.path, pattern="_R1.fastq.gz", full.names = TRUE))
PE.files_full_R2 <- sort(list.files(PE.path, pattern="_R2.fastq.gz", full.names = TRUE))
# Next, extract the sample names (this will require specific adjustment to match the sample naming convention of each project.
library("stringr") # In order to use the string replacement function, we load the stringr library
PE.sam_name <- str_replace(basename(PE.files_full_R1), "_R1.fastq.gz", "")
### 2.- Filter reads ###
test <- c(1,7,13,19) # we can select a subset
pdf("PE_Tutorial-01-plotQualityProfile.pdf") # Due to some missing fonts, we get an error if we try to run it directly, so we output it to a pdf file
plotQualityProfile(c(rbind(PE.files_full_R1[test],PE.files_full_R2[test])))
dev.off()
PE.fil_name_R1 <- file.path("02_filtered/02_paired", paste0(PE.sam_name, "_filt_R1.fastq.gz")) # Paste with no separator the new sufix to all file names
PE.fil_name_R2 <- file.path("02_filtered/02_paired", paste0(PE.sam_name, "_filt_R2.fastq.gz"))
names(PE.fil_name_R1) <- PE.sam_name # we can now rename the rows with the sample names
names(PE.fil_name_R2) <- PE.sam_name

### Optional.- Evaluate read lengths distribution ###
library("Biostrings") #we now load the Biostring library for string manipulation
# The next code is for creating an example of the lenght density plots
pdf("PE_Tutorial-02-read_length_R1.pdf") # Create output file
par(mfrow=c(2,2)) # plot 2x2 graphs per page
test <- c(1,7,13,19) # Define which samples to use (subset)
colors <- c("black","coral","turquoise","cornflowerblue") # Define color name vector 
for(i in 1:4){ # For each subset item
  plot(stats::density(fastq.seqlengths(PE.files_full_R1[test[i]])),main=PE.sam_name[test[i]]) # Create the density graph
  polygon(stats::density(fastq.seqlengths(PE.files_full_R1[test[i]])),col=colors[i],border=colors[i]) # and color it accordingly
}
dev.off()
pdf("PE_Tutorial-02-read_length_R2.pdf") # Create output file
par(mfrow=c(2,2)) # plot 2x2 graphs per page
for(i in 1:4){ # For each subset item
  plot(stats::density(fastq.seqlengths(PE.files_full_R2[test[i]])),main=PE.sam_name[test[i]]) # Create the density graph
  polygon(stats::density(fastq.seqlengths(PE.files_full_R2[test[i]])),col=colors[i],border=colors[i]) # and color it accordingly
}
par(mfrow=c(1,1)) # Reset 1x1 graph per page
dev.off()
# We can also get the same info for all sets
# The next is an for plotting length distribution of all samples
groups <- PE.sam_name # create a vector for coloring the groups
groups[grep("_HE", groups)] <- "black" #assign a color to each year
groups[grep("_HL", groups)] <- "coral"
groups[grep("_IE", groups)] <- "aquamarine"
groups[grep("_IL", groups)] <- "cornflowerblue"

pdf("PE_Tutorial-03-read_length-all_R1.pdf")
par(mfrow=c(4,4)) # Plot 4x4 graphs per page
all_lengths <- 0
for(i in 1:length(PE.sam_name)){ # For each file
  lengths <- fastq.seqlengths(PE.files_full_R1[i]) # exctract the actual lengths of all sequences
  plot(stats::density(lengths),main=PE.sam_name[i]) # now create density graphs
  polygon(stats::density(fastq.seqlengths(PE.files_full_R1[i])),col=as.character(groups[i]),border=as.character(groups[i])) # and use the corresponding color per group
  all_lengths <- c(all_lengths,lengths) # create a large cummulative vector with all lengths
}
par(mfrow=c(1,1)) # Reset 1x1 graph per page
plot(stats::density(all_lengths),main=paste("All samples-mean",round(mean(all_lengths),2))) # Finally, add a single graph summarizing all samples
polygon(stats::density(all_lengths),col=1,border=1)
dev.off()
pdf("PE_Tutorial-03-read_length-all_R2.pdf")
par(mfrow=c(4,4)) # Plot 4x4 graphs per page
all_lengths <- 0
for(i in 1:length(PE.sam_name)){ # For each file
  lengths <- fastq.seqlengths(PE.files_full_R2[i]) # exctract the actual lengths of all sequences
  plot(stats::density(lengths),main=PE.sam_name[i]) # now create density graphs
  polygon(stats::density(fastq.seqlengths(PE.files_full_R2[i])),col=as.character(groups[i]),border=as.character(groups[i])) # and use the corresponding color per group
  all_lengths <- c(all_lengths,lengths) # create a large cummulative vector with all lengths
}
par(mfrow=c(1,1)) # Reset 1x1 graph per page
plot(stats::density(all_lengths),main=paste("All samples-mean",round(mean(all_lengths),2))) # Finally, add a single graph summarizing all samples
polygon(stats::density(all_lengths),col=1,border=1)
dev.off()

### Optional.- Example: Evaluate the total reads passing variable maxEE filters ###
library("R.utils")
series <- seq(0.5,10,.5) # create a float vector to test Expected Errors
MaxEE.compare_R1 <- data.frame(matrix(NA,nrow=length(series),ncol=4));names(MaxEE.compare_R1) <- PE.sam_name[test]; rownames(MaxEE.compare_R1) <- series #Initialize empty table
sZ1 <-countLines(paste0(PE.path,"/",PE.sam_name[test[1]],"_R1.fastq.gz"))[1]/4 # Get the total lines from each of the files
sZ2 <-countLines(paste0(PE.path,"/",PE.sam_name[test[2]],"_R1.fastq.gz"))[1]/4
sZ3 <-countLines(paste0(PE.path,"/",PE.sam_name[test[3]],"_R1.fastq.gz"))[1]/4
sZ4 <-countLines(paste0(PE.path,"/",PE.sam_name[test[4]],"_R1.fastq.gz"))[1]/4
j=0
for(i in series){
  j=j+1
  PE.dadaclean <- filterAndTrim(PE.files_full_R1[test], PE.fil_name_R1[test], truncQ=2, truncLen = 0, trimLeft = 0, trimRight = 0, maxLen = Inf, minLen = 20, maxN = 0, minQ = 0, rm.phix=TRUE, compress=TRUE, multithread=10, maxEE=i)
  MaxEE.compare_R1[j,1] <- countLines(paste0("02_filtered/02_paired/",PE.sam_name[test[1]],"_filt_R1.fastq.gz"))[1]/4
  MaxEE.compare_R1[j,2] <- countLines(paste0("02_filtered/02_paired/",PE.sam_name[test[2]],"_filt_R1.fastq.gz"))[1]/4
  MaxEE.compare_R1[j,3] <- countLines(paste0("02_filtered/02_paired/",PE.sam_name[test[3]],"_filt_R1.fastq.gz"))[1]/4
  MaxEE.compare_R1[j,4] <- countLines(paste0("02_filtered/02_paired/",PE.sam_name[test[4]],"_filt_R1.fastq.gz"))[1]/4
}
# Plot sequence totals
pdf("PE_Tutorial-04-Adjust_MaxEE_R1.pdf")
plot(MaxEE.compare_R1[,1],type="o",col=1,pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example",ylim=c(0,max(MaxEE.compare_R1)),xaxt="n")
axis(1,at=1:length(series),labels=series)
lines(MaxEE.compare_R1[,2],type="o",col="coral",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
lines(MaxEE.compare_R1[,3],type="o",col="turquoise",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
lines(MaxEE.compare_R1[,4],type="o",col="cornflowerblue",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
abline(v=8,col="gray", lty=2)
legend("bottomright", pch=16, lty=1, col=c("black","coral","turquoise","cornflowerblue"), legend=names(MaxEE.compare_R1),bg="white")
# Plot %
plot(MaxEE.compare_R1[,1]*100/sZ1,type="o",col=1,pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter (%)", main="Effect of Max Expected Error filtering example",ylim=c(0,100),xaxt="n")
axis(1,at=1:length(series),labels=series)
lines(MaxEE.compare_R1[,2]*100/sZ2,type="o",col="coral",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
lines(MaxEE.compare_R1[,3]*100/sZ3,type="o",col="turquoise",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
lines(MaxEE.compare_R1[,4]*100/sZ4,type="o",col="cornflowerblue",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
legend("bottomright", pch=16, lty=1, col=c("black","coral","turquoise","cornflowerblue"), legend=names(MaxEE.compare_R1),bg="white")
abline(v=8,col="gray", lty=2)
dev.off()
write.table(MaxEE.compare_R1, file="PE_Tutorial-04-EE_compare_R1.tsv", quote=F, sep="\t",col.names=NA) 

# Now, repeat for R2
MaxEE.compare_R2 <- data.frame(matrix(NA,nrow=length(series),ncol=4));names(MaxEE.compare_R2) <- PE.sam_name[test]; rownames(MaxEE.compare_R2) <- series #Initialize empty table
sZ1 <-countLines(paste0(PE.path,"/",PE.sam_name[test[1]],"_R2.fastq.gz"))[1]/4 # Get the total lines from each of the files
sZ2 <-countLines(paste0(PE.path,"/",PE.sam_name[test[2]],"_R2.fastq.gz"))[1]/4
sZ3 <-countLines(paste0(PE.path,"/",PE.sam_name[test[3]],"_R2.fastq.gz"))[1]/4
sZ4 <-countLines(paste0(PE.path,"/",PE.sam_name[test[4]],"_R2.fastq.gz"))[1]/4
j=0
for(i in series){
  j=j+1
  PE.dadaclean <- filterAndTrim(PE.files_full_R2[test], PE.fil_name_R2[test], truncQ=2, truncLen = 0, trimLeft = 0, trimRight = 0, maxLen = Inf, minLen = 20, maxN = 0, minQ = 0, rm.phix=TRUE, compress=TRUE, multithread=10, maxEE=i)
  MaxEE.compare_R2[j,1] <- countLines(paste0("02_filtered/02_paired/",PE.sam_name[test[1]],"_filt_R2.fastq.gz"))[1]/4
  MaxEE.compare_R2[j,2] <- countLines(paste0("02_filtered/02_paired/",PE.sam_name[test[2]],"_filt_R2.fastq.gz"))[1]/4
  MaxEE.compare_R2[j,3] <- countLines(paste0("02_filtered/02_paired/",PE.sam_name[test[3]],"_filt_R2.fastq.gz"))[1]/4
  MaxEE.compare_R2[j,4] <- countLines(paste0("02_filtered/02_paired/",PE.sam_name[test[4]],"_filt_R2.fastq.gz"))[1]/4
}
# Plot sequence totals
pdf("PE_Tutorial-04-Adjust_MaxEE_R2.pdf")
plot(MaxEE.compare_R2[,1],type="o",col=1,pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example",ylim=c(0,max(MaxEE.compare_R2)),xaxt="n")
axis(1,at=1:length(series),labels=series)
lines(MaxEE.compare_R2[,2],type="o",col="coral",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
lines(MaxEE.compare_R2[,3],type="o",col="turquoise",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
lines(MaxEE.compare_R2[,4],type="o",col="cornflowerblue",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
abline(v=8,col="gray", lty=2)
legend("bottomright", pch=16, lty=1, col=c("black","coral","turquoise","cornflowerblue"), legend=names(MaxEE.compare_R2),bg="white")
# Plot %
plot(MaxEE.compare_R2[,1]*100/sZ1,type="o",col=1,pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter (%)", main="Effect of Max Expected Error filtering example",ylim=c(0,100),xaxt="n")
axis(1,at=1:length(series),labels=series)
lines(MaxEE.compare_R2[,2]*100/sZ2,type="o",col="coral",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
lines(MaxEE.compare_R2[,3]*100/sZ3,type="o",col="turquoise",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
lines(MaxEE.compare_R2[,4]*100/sZ4,type="o",col="cornflowerblue",pch=16, xlab="MaxEE cutoff", ylab="Sequences remaining after filter", main="Effect of Max Expected Error filtering example")
legend("bottomright", pch=16, lty=1, col=c("black","coral","turquoise","cornflowerblue"), legend=names(MaxEE.compare_R2),bg="white")
abline(v=8,col="gray", lty=2)
dev.off()
write.table(MaxEE.compare_R2, file="PE_Tutorial-04-EE_compare_R2.tsv", quote=F, sep="\t",col.names=NA) 
### End of Optional procedures ###
# Error will be tailored to emulate that in the 150 nt sets (V3 and V4) which is 3 and 3.5 errors in all their span (0.02 and 0.02333 respectively for R1 and R2). This is equivalent to 5 and 5.83333 max errors in 250 nt sets. Evidently, these should not be linear but it is just an approximate approach.
PE.dadaclean <- filterAndTrim(PE.files_full_R1, PE.fil_name_R1, PE.files_full_R2, PE.fil_name_R2, trimLeft = 0, trimRight = 0, minLen = c(150,150), maxLen = c(240,240), truncLen=0, maxN=0, maxEE=c(5,5.8), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
write.table(PE.dadaclean, file="PE_Tutorial-05-Dada_clean.tsv", quote=F, sep="\t",col.names=NA) # Table with the totals before and after cleaning

### 3.- Create error models and run Dada ###
PE.derep_R1 <- derepFastq(PE.fil_name_R1, verbose=TRUE)
PE.derep_R2 <- derepFastq(PE.fil_name_R2, verbose=TRUE) # same thing for R2

# Now run the error model creation separately (optionally, nbases can be increased from 1e8 default to a large number such as 1e20 [e.g. nbases=1e20] to force it to check use sequence for model calculation if fastq files have more sequences than the default, which takes more time but is expected to create a more adjusted)
PE.err_R1 <- learnErrors(PE.derep_R1, multithread=10, nbases=1e20, MAX_CONSIST=10,randomize=TRUE, OMEGA_C=0, verbose=TRUE)
# 124711175 total bases in 952553 reads from 24 samples will be used for learning the error rates.
PE.err_R2 <- learnErrors(PE.derep_R2, multithread=10, nbases=1e20, MAX_CONSIST=10,randomize=TRUE, OMEGA_C=0, verbose=TRUE)
# 124711175 total bases in 952553 reads from 24 samples will be used for learning the error rates.

# Now plot the quality models fit (regression) to actual errors
pdf("PE_Tutorial-06-Plot_error_model_R1.pdf")
plotErrors(PE.err_R1, nominalQ=TRUE)
dev.off()
pdf("PE_Tutorial-06-Plot_error_model_R2.pdf")
plotErrors(PE.err_R2, nominalQ=TRUE)
dev.off()

# Now run dada for ASV calculation (The actual dada algorithm is based on the estimated error)
PE.dada_R1 <- dada(PE.fil_name_R1, err=PE.err_R1, multithread=10, pool=FALSE)
# Sample 1 - 35432 reads in 15613 unique sequences.
# Sample 2 - 20881 reads in 10068 unique sequences.
# Sample 3 - 73171 reads in 44898 unique sequences.
# Sample 4 - 53010 reads in 29671 unique sequences.
# Sample 5 - 20465 reads in 12018 unique sequences.

PE.dada_R2 <- dada(PE.fil_name_R2, err=PE.err_R2, multithread=10, pool=FALSE)
# Sample 1 - 35432 reads in 18342 unique sequences.
# Sample 2 - 20881 reads in 10000 unique sequences.
# Sample 3 - 73171 reads in 42719 unique sequences.
# Sample 4 - 53010 reads in 32265 unique sequences.
# Sample 5 - 20465 reads in 13541 unique sequences.

save.image(file = "my_work_space_paired_ends_clean.RData") # Save point to stop for now

### 4.- Merge corrected pairs ###
PE.merged <- mergePairs(PE.dada_R1, PE.derep_R1, PE.dada_R2, PE.derep_R2, maxMismatch = 0, minOverlap = 30, verbose=TRUE)
# 10698 paired-reads (in 66 unique pairings) successfully merged out of 32165 (in 844 pairings) input.
# 1956 paired-reads (in 56 unique pairings) successfully merged out of 19158 (in 910 pairings) input.
# 34695 paired-reads (in 507 unique pairings) successfully merged out of 65485 (in 7883 pairings) input.
# 24460 paired-reads (in 336 unique pairings) successfully merged out of 48571 (in 3582 pairings) input.
# 8368 paired-reads (in 170 unique pairings) successfully merged out of 17823 (in 1670 pairings) input.


### 5.- Create ASVs table ###
# #IMPORTANT: Contrary to PE processing, SE uses the dada ASVs for table reconstruction (no merge intermediary-step is required)
PE.seqtab <- makeSequenceTable(PE.merged) # Sample names as rows and denoised sequences as columns
# Table is ordered by total ASV frequence (across all samples) and has no singletons
dim(PE.seqtab) # The totals table includes no singletons
# [1]   24 10310
table(nchar(getSequences(PE.seqtab)))

#OPTIONAL: subset by length:
PE.seqtab_369_431 <- PE.seqtab[,nchar(colnames(PE.seqtab)) %in% 369:431]
write.table(table(nchar(getSequences(PE.seqtab_369_431))), "PE_Tutorial-07_lenghts.tsv", sep="\t", row.names=FALSE)

We can create a table with each ASV number and its representative sequence
PE.table_tsv_output <- PE.seqtab_369_431
PE.table_tsv_output[PE.table_tsv_output==1]=0 # Don't consider those values that have a single observation per sample, make them 0 (sample singletons)
PE.table_tsv_output <- PE.table_tsv_output[,colSums(PE.table_tsv_output)>1] # filter singleton ASVs across the table
# Export sequences as in fasta format
uniquesToFasta(PE.table_tsv_output, fout='PE_Tutorial-07_ASV.fasta', ids=paste("ASV_",1:ncol(PE.table_tsv_output), sep=""))
ok_tab = PE.table_tsv_output
write.table(cbind("ASVs"=1:nrow(t(PE.table_tsv_output)),"rep_seq"=rownames(t(PE.table_tsv_output))), file="PE_Tutorial-07_ASV_to_seqs-default.tsv", quote=F, sep="\t",row.names=FALSE)
# Now replace the rep_seq with an incremental ASV number
PE.table_tsv_output <- t(PE.table_tsv_output)
rownames(PE.table_tsv_output) <- paste0("ASV_",1:nrow(PE.table_tsv_output))
# and print the output ASV table
write.table(PE.table_tsv_output, file="PE_Tutorial-07_ASV_Table-default.tsv", quote=F, sep="\t",col.names=NA)

### 6.- Chimera checking ###
PE.seqtab.nochim <- removeBimeraDenovo(PE.seqtab_369_431, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 862 bimeras out of 2002 input sequences.
# We can create a new table with each ASV number and its representative sequence
PE.table_tsv_output <- PE.seqtab.nochim
PE.table_tsv_output[PE.table_tsv_output==1]=0 # Don't consider those values that have a single observation per sample, make them 0 (sample singletons)
PE.table_tsv_output <- PE.table_tsv_output[,colSums(PE.table_tsv_output)>1] # filter singleton ASVs across the table
# Export sequences as in fasta format
uniquesToFasta(PE.table_tsv_output, fout='PE_Tutorial-08_ASV.fasta', ids=paste("ASV_",1:ncol(PE.table_tsv_output), sep=""))
nochim=PE.table_tsv_output
write.table(cbind("ASVs"=1:nrow(t(PE.table_tsv_output)),"rep_seq"=rownames(t(PE.table_tsv_output))), file="PE_Tutorial-08_ASV_to_seqs-nochim.tsv", quote=F, sep="\t",row.names=FALSE)
# Now replace the rep_seq with an incremental ASV number
PE.table_tsv_output <- t(PE.table_tsv_output)
rownames(PE.table_tsv_output) <- paste0("ASV_",1:nrow(PE.table_tsv_output))
# and print the output ASV table
write.table(PE.table_tsv_output, file="PE_Tutorial-08_ASV_Table-nochim.tsv", quote=F, sep="\t",col.names=NA)

# We can evaluate the total table dimensions
dim(nochim)
# [1]   24 3994
# We can calclate the % of reads passing the chimera filters
sum(nochim)*100/sum(ok_tab)
# [1] 89.87097
table(nchar(getSequences(nochim)))
# 369  370  372  373  375  381  384  386  387  388  390  391  392  393  394  395 
#   18    7    8    1    1    1    2    3   87    4    3    1    1    5    6    2 
#  396  398  399  400  401  402  403  404  405  406  407  408  409  410  411  412 
#    1    2    6   37   97  493  281  183  128   40   69   11   14    6    3   10 
#  414  415  417  418  419  420  421  422  423  424  425  426  427  428  429  430 
#   14    5   12    9   13   49   53  197   20  272    8   94 1097  593   24    2 
#  431 
#    1 

write.table(table(nchar(getSequences(PE.seqtab.nochim))), "PE_Tutorial-08_nochim_lenghts.tsv", row.names=FALSE, sep="\t")

### 7.- Track reads lost per step ###
# getUniques(PE.derep$V3_HE1) # we can use the getUniques to get the total items per different sequence and we can sum it to calculate total abundance per sample
# By using this, we can create a function to automate this for all samples in a set:
getN <- function(x) sum(getUniques(x)) # Where getUniques gets non-repeated sequences from a dada2 object or merger object (joined reads)
sequences_lost <- cbind(PE.dadaclean, sapply(PE.derep_R1, getN), sapply(PE.dada_R1, getN), sapply(PE.dada_R2, getN), rowSums(ok_tab), rowSums(nochim))
colnames(sequences_lost) <- c("Raw", "Qual_filter", "Derep", "ASVs R1", "ASVs R2", "Merged", "nonchim")
rownames(sequences_lost) <- PE.sam_name
write.table(sequences_lost, "PE_Tutorial-09_Seqs_lost_in_ASVs_processing.tsv", col.names=NA, sep="\t")
# Create a quick assesment of sequences lost throughout the process
pdf("PE_Tutorial-10_preview_reads_passing_ASV_processing.pdf")
matplot(t(sequences_lost[,-5]),type='l',xaxt='n', main="Sequences remaining after each step - R1", xlab="Step", ylab=" Sequences remaining")
axis(1,at=1:ncol(sequences_lost[,-5]),labels=colnames(sequences_lost[,-5]))
# Now R2
matplot(t(sequences_lost[,-4]),type='l',xaxt='n', main="Sequences remaining after each step - R2", xlab="Step", ylab=" Sequences remaining")
axis(1,at=1:ncol(sequences_lost[,-4]),labels=colnames(sequences_lost[,-4]))
# And same thing for the percentage remaining
matplot(t(sequences_lost[,-5]/sequences_lost[,1]*100),type='l',xaxt='n', main="Sequences remaining after each step  - R1 (%)", xlab="Step", ylab=" Percentage of Sequences remaining")
axis(1,at=1:ncol(sequences_lost[,-5]),labels=colnames(sequences_lost[,-5]))
# Now R2
matplot(t(sequences_lost[,-4]/sequences_lost[,1]*100),type='l',xaxt='n', main="Sequences remaining after each step  - R2 (%)", xlab="Step", ylab=" Percentage of Sequences remaining")
axis(1,at=1:ncol(sequences_lost[,-4]),labels=colnames(sequences_lost[,-4]))
dev.off()

# Save work so far
save.image(file = "PE_clean_tutorial.RData")

### 8.- Import ASVs into QIIME 2 (tables and sequence sets) ###
# QIIME-enabled output according to https://forum.qiime2.org/t/exporting-dada2-r-package-results-to-work-with-qiime2/2573/5:
# Import fasta into qiime artifact qza
conda activate qiime2-2019.10
qiime tools import --input-path PE_Tutorial-08_ASV.fasta --type 'FeatureData[Sequence]' --output-path PE.Tutorial-11_ASV_rep_seq.qza
# append missing header to the table for import
cat <(echo -n "#OTU Table") PE_Tutorial-08_ASV_Table-nochim.tsv > temp.txt
# convert to biom
biom convert -i temp.txt -o temp.biom --table-type="OTU table" --to-hdf5
# and create table-type qza
qiime tools import --input-path temp.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path PE.Tutorial-11_ASV_table.qza

################### END of PE DADA2 protocol for V3V4 2016 Shrimp samples ###################