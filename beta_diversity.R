
#Load the packages
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")


library(tidyverse)
library(vegan)
library(qiime2R)


##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza
#
# To get these files you need to scp them from the cluster:
#
# first on  your laptop cd to the directory where you want to save them.
# Then use this code for our example dataset today:
# mkdir core-metrics-results/
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/* .
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/core-metrics-results/* core-metrics-results/.
##############################################


###Set your working directory
#setwd ("~/Documents/ANSC516/ANSC516-repo/ANSC516/data/moving-pictures/")

list.files()

if(!dir.exists("output"))
  dir.create("output")

#How to load a file into R
metadata2 <- read.delim("MetadataS2.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadata2[1,]
metadata2[,1]
# When subsetting, the first number is the row and after the comma is the column
metadata2 <- metadata2[-1,]

#Now the qiime2R method
metadata<-read_q2metadata("MetadataS2.txt")
str(metadata)
levels(metadata$`treatment`)
colnames(metadata)[2] <- "treatment"
colnames(metadata)[6] <- "donor"
str(metadata)

row.names(metadata) <- metadata[,1]
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)

bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")
j_PCoA<-read_qza("core-metrics-results/jaccard_pcoa_results.qza")
wUF <- read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")
uwUF <- read_qza("core-metrics-results/unweighted_unifrac_pcoa_results.qza")

body_colors <- c("Black", "Blue", "Green", "Gray", "Red", "Yellow")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Now we are going to make an ordination plot
BC_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")

ggplot(bc_meta, aes(x=PC1, y=PC2, color=treatment)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #theme_q2r() +
  xlab("PC1 (68.49%)") +
  ylab("PC2 (7.91%)") +
  scale_color_manual(values=c("Blue", "Black", "Green", "Gray", "Red","Yellow"), name = "treatment")

# Now we are going to make our code a little more re-usable
my_column <- "treatment"
my_column2 <- "donor"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #theme_q2r() +
  facet_grid(~treatment) +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/BC-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "treatment"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= treatment)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
  #scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,"-subject.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

Jaccard

J_PCoA<-read_qza("core-metrics-results/jaccard_pcoa_results.qza")

J_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Now we are going to make an ordination plot
ggplot(J_meta, aes(x=PC1, y=PC2, color=treatment)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #theme_q2r() +
  xlab("PC1 (63.07%)") +
  ylab("PC2 (5.88%)") +
  scale_color_manual(values=c("Blue", "Black", "Green", "Gray", "Red","Yellow"), name = "treatment")

# Now we are going to make our code a little more re-usable
my_column <- "treatment"
#my_column <- "DietTreatment"

ggplot(J_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #theme_q2r() +
  facet_grid(~treatment) +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/Jd-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),J_meta,mean)
colnames(centroids)[1] <- "treatment"

ggplot(J_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/Jd-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(J_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= treatment)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
#scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/Jd-ellipse_", my_column,"-subject.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches














## SAME thing but with weighted UniFrac

Wuni_PCoA<-read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)
colnames(centroids)[1] <- "treatment"

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "treatment")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= donor), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "treatment")
ggsave(paste0("output/Wuni-ellipse_", my_column,"-treatment.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Unweighted Unifrac


UWuni_PCoA<-read_qza("core-metrics-results/unweighted_unifrac_pcoa_results.qza")

UWuni_meta <- UWuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),UWuni_meta,mean)
colnames(centroids)[1] <- "treatment"

ggplot(UWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*UWuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*UWuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "treatment")
ggsave(paste0("output/UWuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(UWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= donor), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*UWuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*UWuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "treatment")
ggsave(paste0("output/UWuni-ellipse_", my_column,"-treatment.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#Bray curtis

bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "treatment"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "treatment")
ggsave(paste0("output/BC-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#USE to divide IBD and Healthy.
ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= donor), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "treatment")
ggsave(paste0("output/BC-ellipse_", my_column,"-treatment.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Bray curtis by donor


bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")

bc2_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc2_meta,mean)
colnames(centroids)[1] <- "donor"

ggplot(bc2_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "donor")
ggsave(paste0("output/BC2-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(bc2_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= treatment,), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "donor")
ggsave(paste0("output/BC2-ellipse_", my_column,"-donor.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

# Jaccard

j_PCoA<-read_qza("core-metrics-results/jaccard_pcoa_results.qza")

j_meta <- j_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),j_meta,mean)
colnames(centroids)[1] <- "treatment"

ggplot(j_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "treatment")
ggsave(paste0("output/j-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(j_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= donor), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  #theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "treatment")
ggsave(paste0("output/j-ellipse_", my_column,"-treatment.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

bc_dist_mat<-read_qza("core-metrics-results/bray_curtis_distance_matrix.qza")
bc_dm <- as.matrix(bc_dist_mat$data) 
rownames(bc_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(bc_dm),metadata$SampleID),]
rownames(bc_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out1 <- adonis2(bc_dm ~ donor, data = metadata_sub)

j_dist_mat<-read_qza("core-metrics-results/jaccard_distance_matrix.qza")
j_dm <- as.matrix(j_dist_mat$data) 
rownames(j_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(j_dm),metadata$SampleID),]
rownames(j_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out2 <- adonis2(j_dm ~ donor, data = metadata_sub)

uw_dist_mat<-read_qza("core-metrics-results/unweighted_unifrac_distance_matrix.qza")
uw_dm <- as.matrix(uw_dist_mat$data) 
rownames(uw_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(uw_dm),metadata$SampleID),]
rownames(uw_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out4 <- adonis2(uw_dm ~ donor, data = metadata_sub)

w_dist_mat<-read_qza("core-metrics-results/weighted_unifrac_distance_matrix.qza")
w_dm <- as.matrix(w_dist_mat$data) 
rownames(w_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(w_dm),metadata$SampleID),]
rownames(w_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out3 <- adonis2(w_dm ~ donor, data = metadata_sub)




write.table(PERMANOVA_out,"output/Body.site_Adonis_overall.csv",sep=",", row.names = TRUE) 

######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

body.site_Pair <- pairwise.adonis2(bc_dm ~ body.site, data = metadata_sub)
write.table(body.site_Pair,"output/Body.site_Adonis_pairwise.csv",sep=",", row.names = TRUE) 

