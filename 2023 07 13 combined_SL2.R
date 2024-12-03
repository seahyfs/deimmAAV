# Setting up --------------------------------------------------------------
pacman::p_load(DESeq2, limma, tidyverse, readxl, janitor, ggplot2, ComplexHeatmap, circlize, ggpmisc, ggExtra, Cairo, gridExtra, patchwork, cowplot)
rm(list=ls())

#loading dataset 1
working_directory <- "C:/Users/seahyfs/Desktop/20230624_SL2_Cas_only_analysis/SL2 combined"
setwd(working_directory)

countsfile1 <- "counts_no_input.tsv"
countsfile2 <- "counts2.tsv"

cts1 <- read.csv(countsfile1, sep = "\t", header=TRUE, row.names=1) %>% as.matrix(row.names="gene_id")
cts2 <- read.csv(countsfile2, sep = "\t", header=TRUE, row.names=1) %>% as.matrix(row.names="gene_id")

#generate actual cts1 dataset (to do because cts1 still has data from different lanes)

pair_columns <- c("Neg1", "Neg2", "Neg3", "Neg4", "Neg5", "Neg6", "Neg7", "Neg8", "Cas9.1", "Cas9.2",
                  "MTP009.1", "MTP022.1", "MTP030.1", "MTP031.1", "MTP039.1", "MTP043.1", "MTP045.1", "MTP053.1",
                  "MTP054.1", "MTP056.1", "MTP063.1", "MTP065.1", "MTP085.1", "MTP100.1", "MTP103.1", "MTP109.1",
                  "MTP115.1", "MTP135.1", "MTP145.1", "MTP146.1", "MTP147.1", "MTP148.1", "MTP154.1", "MTP156.1", 
                  "MTP009.2", "MTP022.2", "MTP030.2", "MTP031.2", "MTP039.2", "MTP043.2", "MTP045.2", "MTP053.2",
                  "MTP054.2", "MTP056.2", "MTP063.2", "MTP065.2", "MTP085.2", "MTP100.2", "MTP103.2", "MTP109.2",
                  "MTP115.2", "MTP135.2", "MTP145.2", "MTP146.2", "MTP147.2", "MTP148.2", "MTP154.2", "MTP156.2",
                  "OP1",  "SP1", "OP2", "SP2")

cts1 <- cts1 %>% data.frame()

for (col in pair_columns) {
  col1 <- paste0(col, "_L1")
  col2 <- paste0(col, "_L3")
  
  cts1 <- cts1 %>%
    mutate(!!col := !!sym(col1) + !!sym(col2)) %>%
    select(-one_of(col1, col2))
}

#duplicate negative samples (cts1)
ctsadd <- cts1[, grepl("Neg", colnames(cts1))]

oldnames <-  colnames(ctsadd)
newnames <- NULL

for (i in oldnames){
  newname <- paste0(i,"A")
  newnames <- c(newnames, newname)
}

colnames(ctsadd) <- newnames
cts1a <- cbind(cts1, ctsadd)

#duplicate negative samples (cts2)
ctsadd <- cts2[, grepl("Neg", colnames(cts2))]

oldnames <-  colnames(ctsadd)
newnames <- NULL

for (i in oldnames){
  newname <- paste0(i,"A")
  newnames <- c(newnames, newname)
}

colnames(ctsadd) <- newnames
cts2a <- cbind(cts2, ctsadd)

#rename columns too for cts2
cts2oldnames <- colnames(cts2a)
cts2newnames <- NULL

for (i in cts2oldnames){
  newname <- paste0("2",i)
  cts2newnames <- c(cts2newnames, newname)
}

colnames(cts2a) <- cts2newnames

cts2a <- data.frame(cts2a)
cts <- cbind(cts1a, cts2a)

write.csv(cts,'combinedcounts.csv',row.names=TRUE) 

#load colfile
colfile <- "20230713_comblabels.csv"
coldata <- read.csv(colfile, header=TRUE, row.names=1) %>% 
  select(Sample, Category, Batch) %>%  
  mutate_at(vars(Sample, Category, Batch), list(factor))

#check if rownames are there, check if names match
if (!all(rownames(coldata) %in% colnames(cts))){
  stop("There is a problem with your input data!")
}else if (all(rownames(coldata) != colnames(cts))) {
  stop("There is a problem with your input data!")
}else{print("You may proceed")}

all(rownames(coldata) %in% colnames(cts))#run DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Sample)

#pre-filtering (remove those with less than 10 counts total == in the end it's removing 6 rows)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$Sample <- relevel(dds$Sample, ref = "Neg") ##change this to change what the data is normalised to!
dds <- DESeq(dds)

#check for batch effects
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup="Sample")

pcaData <- plotPCA(vsd, intgroup=c("Sample", "Batch", "Category"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Category)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#subset data that is relevant
ddsdata<- as.data.frame(dds@rowRanges@elementMetadata@listData)
pvals <- ddsdata %>% select_if(grepl("WaldPvalue_Sample_", names(.)))
l2fc <- ddsdata %>% select_if(grepl("\\bSample_", names(.)))

#write as csv (just to save)
write.csv(pvals,'pvalsall.csv',row.names=TRUE)
write.csv(l2fc,'l2fc.csv',row.names=TRUE)

# Functions ---------------------------------------------------------------
#function to convert to neg log10
convertToNegativeLog10 <- function(data) {
  # Identify numeric columns
  numeric_cols <- sapply(data, is.numeric)
  
  # Loop through each column and convert numeric values
  for (col in names(data)[numeric_cols]) {
    data[[col]] <- ifelse(data[[col]] > 0, -log10(data[[col]]), 0)
  }
  
  return(data)
}

#function to return nothing for every nth
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

# Data manipulation -----------------------------------------------------------
df <- pvals %>% 
  convertToNegativeLog10() %>% 
  select_all(~gsub("WaldPvalue_Sample_", "", .)) %>% 
  select_all(~gsub("_vs_Neg", "", .)) 

## Combining pvals and l2fc - where l2fc is negative
int_values <- l2fc[, !(names(l2fc) %in% "...1")]

# Create a mask for values less than or equal to 0
mask <- int_values <= 0

# Apply the mask to the numeric values of array1
masked_values <- df[, !(names(df) %in% "...1")]

for (i in 1:nrow(masked_values)) {
  for (j in 1:ncol(masked_values)) {
    if (mask[i, j]) {
      masked_values[i, j] <- 0
    }
  }
}

df <- masked_values %>% tibble::rownames_to_column("epitope")

#separate name into protein and location
df <- df %>% separate_wider_delim(epitope, delim="|", names = c("protein", "loc"))

## rename samples
colnames(df) <- c("protein", "loc", "Cas",  "FECD#009",  "FECD#022", 
                  "FECD#030", "FECD#031", "FECD#039", "FECD#043", "FECD#045", "FECD#053",
                  "FECD#054", "FECD#056", "FECD#063", "FECD#065", "FECD#085", "FECD#100", "FECD#103", "FECD#109",
                  "FECD#115", "FECD#135", "FECD#145", "FECD#146", "FECD#147", "FECD#148", "FECD#154", "FECD#156", 
                  "Neg1","Neg2", "Neg3", "Neg4", "Neg5", "Neg6", "Neg7", "Neg8", 
                  "HD", "HDD", "SG", "SGD")


# Sample orders -----------------------------------------------------------
#Sample order for AAV
sample_order1 <- c("Neg1","Neg2", "Neg3", "Neg4", "Neg5", "Neg6", "Neg7", "Neg8", "Cas")

sample_order2 <- c(
  "FECD#009",  "FECD#022", 
  "FECD#030", "FECD#031", "FECD#039", "FECD#043", "FECD#045", "FECD#053",
  "FECD#054", "FECD#056", "FECD#063", "FECD#065", "FECD#085", "FECD#100", "FECD#103", "FECD#109",
  "FECD#115", "FECD#135", "FECD#145", "FECD#146", "FECD#147", "FECD#148", "FECD#154", "FECD#156", 
  "HD", "SG")

sample_order <- c(sample_order1, sample_order2)

#Sample order for Cas - typically only the Cas is different?
sample_orderc1 <- c("Neg1","Neg2", "Neg3", "Neg4", "Neg5", "Neg6", "Neg7", "Neg8")#, "hAlb")

sample_orderc2 <- c("FECD#009",  "FECD#022", 
                    "FECD#030", "FECD#031", "FECD#039", "FECD#043", "FECD#045", "FECD#053",
                    "FECD#054", "FECD#056", "FECD#063", "FECD#065", "FECD#085", "FECD#100", "FECD#103", "FECD#109",
                    "FECD#115", "FECD#135", "FECD#145", "FECD#146", "FECD#147", "FECD#148", "FECD#154", "FECD#156", 
                    "HD", "SG")

#capsid order
caporder = c("AAV_DJ", "AAV1", "AAV2", "AAV3", "AAV4", "AAV5", "AAV6", "AAV7", "AAV8", "AAV9", "AAV10", "AAV11", "AAV12", "AAV13" )

# Extracting relevant subset of data --------------------------------------
#extract AAVs only
AAV <-  df %>% filter(grepl('AAV', protein))
serotypes <- c(1:13, "_DJ")

#extracting for each AAV serotype
AAVsero <- list()
AAVserononeg <- list()
AAVtemp <- NULL
AAVmat <- NULL
AAVserocontrols <- list()
AAVserolowAb <- list()
AAVserohighAb <- list()

for (i in serotypes){
  sero <- paste0("AAV",i,"_capsid", "")
  AAVsero[[i]] <- AAV %>% filter(grepl(sero, protein))
}

for (i in serotypes){
  AAVtemp <- AAVsero[[i]]
  AAVmat<- AAVtemp %>% select (-protein ) %>% select(-loc) %>% as.matrix()
  rownames(AAVmat) <- AAVtemp$loc
  AAVmat <- AAVmat[str_order(AAVtemp$loc, numeric = TRUE),] 
  AAVsero[[i]] <- AAVmat[, sample_order]
  AAVserononeg[[i]] <- AAVmat[, sample_order2] 
  AAVserocontrols[[i]] <- AAVmat[, sample_order1]
#  AAVserolowAb[[i]] <- AAVmat[, lowAb] 
#  AAVserohighAb[[i]] <- AAVmat[, highAb]
}


#extract for each sample, for different serotypes
AAVtemp2 <- NULL
AAVbysample <- list()
AAVmat2 <- NULL

for (i in sample_order){
  AAVtemp2 <- AAV %>% select(protein, loc, i) %>% 
    pivot_wider(names_from = protein, values_from = i) %>% 
    select_all(~gsub("_capsid", "", .)) 
  AAVmat2 <- AAVtemp2 %>% select(-loc) %>% as.matrix()
  rownames(AAVmat2) <- AAVtemp2$loc
  AAVmat2 <- AAVmat2[str_order(AAVtemp2$loc, numeric = TRUE),]
  AAVbysample[[i]] <- AAVmat2[, caporder] 
}

#extract Cas only
Cas_types <- c("spCas9", "saCas9", "CasRX")
Cas <- list()
Casneg <- list()
Casnoneg <- list()
Castemp <- NULL
Casmat <- NULL

for (i in Cas_types){
  Castemp <-  df %>% 
    filter(grepl(i, protein)) %>% 
    as.data.frame()
  Castemp <- Castemp %>% 
    select (- protein)%>% 
    arrange(match(loc,str_sort(Castemp$loc, numeric = TRUE)))
  Casmat <- Castemp %>% select(-loc) %>% as.matrix()
  rownames(Casmat) <- Castemp$loc
  row.names(Casmat)[row.names(Casmat) =="160-2099"] <- "160-209"
  Cas[[i]] <- Casmat[, sample_order]
  Casnoneg[[i]] <- Casmat[, sample_orderc2] 
  Casneg[[i]] <- Casmat[, sample_orderc1]
}

#to include plots in a different folder
wdir<- getwd()
newwdir <- paste0(wdir, "/plots")
setwd(newwdir)

# Heatmaps ----------------------------------------------------------------
#Heatmap conditions
col_fun1 = colorRamp2(c(4, 0 ), c("#8F8F8F", "white"))
col_fun2 = colorRamp2(c(4, 0 ), c("blue", "white"))
col_fun3 = colorRamp2(c(8, 0 ), c("blue", "white"))
col_fun1(seq(0, 4))
col_fun2(seq(0, 4))
col_fun3(seq(0, 8))

dir1 <- paste0(newwdir, "/AAV by serotype_20240129")
setwd(dir1)

#for loop to generate AAV heatmaps (for individual AAVs, including all samples)
for (i in serotypes){
  #defined column order
  ht_list <-draw( Heatmap(AAVserocontrols[[i]], 
                          col = col_fun1, 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE, 
                          row_names_side = "left",
                          heatmap_legend_param = list(title = "-log10(p-value)"),  
                          column_names_gp = grid::gpar(fontfamily = "sans", fontsize = 10),
                          row_names_gp = grid::gpar(fontfamily = "sans", fontsize = 8)) +
                    Heatmap(AAVserononeg[[i]], 
                            col = col_fun2, 
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE, 
                            row_names_side = "left", 
                            heatmap_legend_param = list(title = "-log10(p-value)"),
                            column_names_gp = grid::gpar(fontfamily = "sans", fontsize = 10),
                            row_names_gp = grid::gpar(fontfamily = "sans", fontsize = 8)), 
                  column_title = paste0("AAV",i), column_title_gp = gpar(fontfamily = "sans", fontsize = 22), 
                  ht_gap = unit(0, "cm"))
  
  
  pdf(paste0("AAV",i,"_capsid", " ", "deforder.pdf"))
  print(ht_list)
  dev.off() 
  
  #still split into 2 different segments, each individually undergoing clustering
  ht_list2 <-draw( Heatmap(AAVserocontrols[[i]], 
                           col = col_fun1, 
                           cluster_rows = FALSE, 
                           cluster_columns = TRUE, 
                           column_dend_reorder = TRUE,
                           row_names_side = "left",
                           heatmap_legend_param = list(title = "-log10(p-value)"),  
                           column_names_gp = grid::gpar(fontsize = 10),
                           row_names_gp = grid::gpar(fontsize = 8)) +
                     Heatmap(AAVserononeg[[i]], 
                             col = col_fun2, 
                             cluster_rows = FALSE, 
                             cluster_columns = TRUE, 
                             column_dend_reorder = TRUE,
                             row_names_side = "left", 
                             heatmap_legend_param = list(title = "-log10(p-value)"),
                             column_names_gp = grid::gpar(fontsize = 10),
                             row_names_gp = grid::gpar(fontsize = 8)), 
                   column_title = paste0("AAV",i,"_capsid"), column_title_gp = gpar(fontsize = 22), 
                   ht_gap = unit(0, "cm"))
  
  
  pdf(paste0("AAV",i,"_capsid", " ", "clustordersep.pdf"))
  print(ht_list2)
  dev.off() 
  
  #clustered columns
  clust <- Heatmap(AAVsero[[i]], 
                   col = col_fun2, 
                   cluster_rows = FALSE, 
                   cluster_columns = TRUE, 
                   column_dend_reorder = TRUE,
                   row_names_side = "left",
                   heatmap_legend_param = list(title = "-log10(p-value)"),
                   column_names_gp = grid::gpar(fontsize = 10),
                   row_names_gp = grid::gpar(fontsize = 8), 
                   column_title = paste0("AAV",i,"_capsid"), column_title_gp = gpar(fontsize = 22))
  
  pdf(paste0("AAV",i,"_capsid", " ", "clustorder.pdf"))
  print(clust)
  dev.off()
}

#clustering epitopes
dir0 <- paste0(newwdir, "/AAV by serotype clustered epitopes")
setwd(dir0)

Nterm <- c("0-49","10-59","25-74","35-84","50-99","60-109","75-124","85-134","100-149","110-159","125-174","135-184","150-199","160-209","175-224", "185-234", "200-249", "210-259")
Cterm <- c("385-434", "400-449", "410-459", "425-474","435-484","450-499","460-509","475-524","485-534","500-549","510-559","525-574","535-584","550-599","560-609","575-624","585-634", "600-649", "610-659", "625-674")

#for loop to generate AAV heatmaps (for individual AAVs, including all samples)
for (i in serotypes){
  
  # Set stylings for row names and make our selected rows unique
  Nterm_idx <- which(rownames(AAVserononeg[[i]]) %in% Nterm)
  Cterm_idx <- which(rownames(AAVserononeg[[i]]) %in% Cterm)
  fontcolors <- rep('black', nrow(AAVserononeg[[i]]))
  fontcolors[Nterm_idx] <- '#70AD47'
  fontcolors[Cterm_idx] <- '#ED7D31'
  
  # Create text annotation object for displaying row names
  rowAnno <- rowAnnotation(rows = anno_text(rownames(AAVserononeg[[i]]), gp = gpar(col = fontcolors, fontsize = 8)))

  #defined column order
  ht_list <-draw( Heatmap(AAVserononeg[[i]], 
                            col = col_fun2, 
                            cluster_rows = TRUE, 
                            cluster_columns = FALSE, 
                            show_row_names = FALSE,
                            #right_annotation = rowAnno,
                            heatmap_legend_param = list(title = "-log10(p-value)"),
                            column_names_gp = grid::gpar(fontsize = 10)) +
                    Heatmap(AAVserocontrols[[i]], 
                            col = col_fun1, 
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE, 
                            show_row_names = FALSE,
                            right_annotation = rowAnno,
                            heatmap_legend_param = list(title = "-log10(p-value)"),  
                            column_names_gp = grid::gpar(fontsize = 10),
                            row_names_gp = grid::gpar(fontsize = 8)), 
                  column_title = paste0("AAV",i,"_capsid"), column_title_gp = gpar(fontsize = 22), 
                  ht_gap = unit(0, "cm"))
  
  
  pdf(paste0("AAV",i,"_capsid", " ", "deforder+epi.pdf"))
  print(ht_list)
  dev.off() 
  
  #still split into 2 different segments, each individually undergoing clustering
  ht_list2 <-draw(Heatmap(AAVserononeg[[i]], 
                             col = col_fun2, 
                             cluster_rows = TRUE, 
                             cluster_columns = TRUE, 
                             column_dend_reorder = TRUE,
                            show_row_names = FALSE,
                             heatmap_legend_param = list(title = "-log10(p-value)"),
                             column_names_gp = grid::gpar(fontsize = 10),
                             row_names_gp = grid::gpar(fontsize = 8)) +
                      Heatmap(AAVserocontrols[[i]], 
                              col = col_fun1, 
                              cluster_rows = FALSE, 
                              cluster_columns = TRUE, 
                              show_row_names = FALSE,
                              right_annotation = rowAnno,
                              heatmap_legend_param = list(title = "-log10(p-value)"),  
                              column_names_gp = grid::gpar(fontsize = 10),
                              row_names_gp = grid::gpar(fontsize = 8)), 
                   column_title = paste0("AAV",i,"_capsid"), column_title_gp = gpar(fontsize = 22), 
                   ht_gap = unit(0, "cm"))
  
  
  pdf(paste0("AAV",i,"_capsid", " ", "clustordersep+epi.pdf"))
  print(ht_list2)
  dev.off() 
}

dir2 <- paste0(newwdir, "/AAV by sample_20240129")
setwd(dir2)

#for loop to generate AAV heatmaps (for individual samples, all AAVs)
for (i in sample_order){
  #defined column order
  sample1 <-Heatmap(AAVbysample[[i]], 
                    col = col_fun2, 
                    cluster_rows = FALSE, 
                    cluster_columns = FALSE, 
                    row_names_side = "left",
                    heatmap_legend_param = list(title = "-log10(p-value)"),  
                    column_names_gp = grid::gpar(fontfamily = "sans", fontsize = 10),
                    row_names_gp = grid::gpar(fontfamily = "sans", fontsize = 8),  
                    column_title = paste0(i), column_title_gp = gpar(fontfamily = "sans", fontsize = 18))
  
  pdf(paste0(i,"_for_all_AAVs.pdf"), width = 6, height = 8)
  print(sample1)
  dev.off() 
}

dir3 <- paste0(newwdir, "/Cas by protein")
setwd(dir3)

#for loop to generate Cas heatmaps
for (i in Cas_types){
  #defined column order
  ht_list <-draw( Heatmap(Casneg[[i]], 
                          col = col_fun1, 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE, 
                          row_names_side = "left",
                          heatmap_legend_param = list(title = "-log10(p-value)"),  
                          column_names_gp = grid::gpar(fontsize = 10),
                          row_names_gp = grid::gpar(fontsize = 8)) +
                    Heatmap(Casnoneg[[i]], 
                            col = col_fun2, 
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE, 
                            row_names_side = "left", 
                            heatmap_legend_param = list(title = "-log10(p-value)"),
                            column_names_gp = grid::gpar(fontsize = 10),
                            row_names_gp = grid::gpar(fontsize = 6)), 
                  column_title = i, column_title_gp = gpar(fontsize = 22), 
                  ht_gap = unit(0, "cm"))
  
  
  pdf(paste0(i, " ", "deforder.pdf"), width = 8, height = 11)
  print(ht_list)
  dev.off() 
  
  #still split into 2 different segments, each individually undergoing clustering
  ht_list2 <-draw( Heatmap(Casneg[[i]], 
                           col = col_fun1, 
                           cluster_rows = FALSE, 
                           cluster_columns = TRUE, 
                           row_names_side = "left",
                           heatmap_legend_param = list(title = "-log10(p-value)"),  
                           column_names_gp = grid::gpar(fontsize = 10),
                           row_names_gp = grid::gpar(fontsize = 8)) +
                     Heatmap(Casnoneg[[i]], 
                             col = col_fun2, 
                             cluster_rows = FALSE, 
                             cluster_columns = TRUE, 
                             row_names_side = "left", 
                             heatmap_legend_param = list(title = "-log10(p-value)"),
                             column_names_gp = grid::gpar(fontsize = 10),
                             row_names_gp = grid::gpar(fontsize = 6)), 
                   column_title = i, column_title_gp = gpar(fontsize = 22), 
                   ht_gap = unit(0, "cm"))
  
  
  pdf(paste0(i, " ", "clustordersep.pdf"), width = 8, height = 11)
  print(ht_list2)
  dev.off() 
  
  
  #clustered columns
  clust <- Heatmap(Cas[[i]], 
                   col = col_fun2, 
                   cluster_rows = FALSE, 
                   cluster_columns = TRUE, 
                   column_dend_reorder = TRUE,
                   row_names_side = "left",
                   heatmap_legend_param = list(title = "-log10(p-value)"),
                   column_names_gp = grid::gpar(fontsize = 10),
                   row_names_gp = grid::gpar(fontsize = 6), 
                   column_title = i, column_title_gp = gpar(fontsize = 22))
  
  pdf(paste0(i, " ", "clustorder.pdf"), width = 8, height = 11)
  print(clust)
  dev.off()
  
}

dir4 <- paste0(newwdir, "/AAV serum categories")
setwd(dir4)

#for serum samples in different categories based on Amine's data
i = "_DJ"

lowhigh <- draw(Heatmap(AAVserocontrols[[i]], 
                         col = col_fun1, 
                         cluster_rows = FALSE, 
                         cluster_columns = TRUE, 
                         column_dend_reorder = TRUE,
                         row_names_side = "left",
                         column_title = "Negative controls",
                         heatmap_legend_param = list(title = "-log10(p-value)"),  
                         column_names_gp = grid::gpar(fontsize = 10),
                         row_names_gp = grid::gpar(fontsize = 8)) +
                   Heatmap(AAVserolowAb[[i]], 
                           col = col_fun2, 
                           cluster_rows = FALSE, 
                           cluster_columns = TRUE, 
                           column_dend_reorder = TRUE,
                           row_names_side = "left",
                           column_title = "Titer < 1:100",
                           heatmap_legend_param = list(title = "-log10(p-value)"),  
                           column_names_gp = grid::gpar(fontsize = 10),
                           row_names_gp = grid::gpar(fontsize = 8)) +
                   Heatmap(AAVserohighAb[[i]], 
                           col = col_fun2, 
                           cluster_rows = FALSE, 
                           cluster_columns = TRUE, 
                           column_dend_reorder = TRUE,
                           row_names_side = "left", 
                           column_title = "Titer ≥ 1:100",
                           show_heatmap_legend = FALSE,
                           column_names_gp = grid::gpar(fontsize = 10),
                           row_names_gp = grid::gpar(fontsize = 8)), 
                 column_title = paste0("AAV",i,"_capsid"), 
                 column_title_gp = gpar(fontsize = 22), 
                 ht_gap = unit(0, "cm"))

cairo_pdf(paste0("AAV",i,"_capsid", " ", "nAbgroups_withnegs.pdf"))
print(lowhigh)
dev.off() 

lowhigh2 <-draw( Heatmap(AAVserolowAb[[i]], 
                          col = col_fun2, 
                          cluster_rows = FALSE, 
                          cluster_columns = TRUE, 
                          column_dend_reorder = TRUE,
                          row_names_side = "left",
                          column_title = "Titer < 1:100",
                          heatmap_legend_param = list(title = "-log10(p-value)"),  
                          column_names_gp = grid::gpar(fontsize = 10),
                          row_names_gp = grid::gpar(fontsize = 8)) +
                  Heatmap(AAVserohighAb[[i]], 
                          col = col_fun2, 
                          cluster_rows = FALSE, 
                          cluster_columns = TRUE, 
                          column_dend_reorder = TRUE,
                          row_names_side = "left", 
                          column_title = "Titer ≥ 1:100",
                          show_heatmap_legend = FALSE,
                          column_names_gp = grid::gpar(fontsize = 10),
                          row_names_gp = grid::gpar(fontsize = 8)), 
                column_title = paste0("AAV",i,"_capsid"), 
                column_title_gp = gpar(fontsize = 22), 
                ht_gap = unit(0, "cm"))



cairo_pdf(paste0("AAV",i,"_capsid", " ", "nAbgroups_nonegs.pdf"))
print(lowhigh2)
dev.off() 

# Comparing AAV1 and AAV6 -------------------------------------------------
#extract AAVs only
AAV <-  df %>% filter(grepl('AAV', protein))

col_order <- c("Neg1-AAV1", "Neg1-AAV6",   "Neg2-AAV1",   "Neg2-AAV6",   "Neg3-AAV1",   "Neg3-AAV6",   "Neg4-AAV1",  "Neg4-AAV6",
               "Neg5-AAV1",   "Neg5-AAV6",   "Neg6-AAV1",   "Neg6-AAV6",  "Neg7-AAV1",   "Neg7-AAV6",   "Neg8-AAV1",   "Neg8-AAV6", 
                "Cas-AAV1",    "Cas-AAV6",
               "MTP009-AAV1", "MTP009-AAV6", "MTP022-AAV1", "MTP022-AAV6", "MTP030-AAV1",
               "MTP030-AAV6", "MTP031-AAV1", "MTP031-AAV6", "MTP039-AAV1", "MTP039-AAV6", "MTP043-AAV1", "MTP043-AAV6", "MTP045-AAV1", "MTP045-AAV6", "MTP053-AAV1", "MTP053-AAV6",
               "MTP054-AAV1", "MTP054-AAV6", "MTP056-AAV1", "MTP056-AAV6", "MTP063-AAV1", "MTP063-AAV6", "MTP065-AAV1", "MTP065-AAV6", "MTP085-AAV1", "MTP085-AAV6", "MTP100-AAV1",
               "MTP100-AAV6", "MTP103-AAV1", "MTP103-AAV6", "MTP109-AAV1", "MTP109-AAV6", "MTP115-AAV1", "MTP115-AAV6", "MTP135-AAV1", "MTP135-AAV6", "MTP145-AAV1", "MTP145-AAV6",
               "MTP146-AAV1", "MTP146-AAV6", "MTP147-AAV1", "MTP147-AAV6", "MTP148-AAV1", "MTP148-AAV6" ,"MTP154-AAV1", "MTP154-AAV6", "MTP156-AAV1", "MTP156-AAV6",
               "OP-AAV1",     "OP-AAV6", "OPD-AAV1",     "OPD-AAV6", "SP-AAV1"   ,  "SP-AAV6", "SPD-AAV1"   ,  "SPD-AAV6" )


#stuff of interest
AAV1 <- AAV %>% 
  filter(grepl("AAV1_capsid", protein)) %>% 
  mutate(protein = gsub("_capsid", "", protein)) %>% 
  pivot_longer(cols = !c(protein, loc),
               names_to = "Sample",
               values_to = "p-values")


AAV6 <- AAV %>% 
  filter(grepl("AAV6_capsid", protein)) %>% 
  mutate(protein = gsub("_capsid", "", protein)) %>% 
  pivot_longer(cols = !c(protein, loc),
               names_to = "Sample",
               values_to = "p-values")

AAVcomb <- bind_rows(AAV1, AAV6)
AAVcomb$Sample_AAV <- str_c(AAVcomb$Sample, '-', AAVcomb$protein)
AAVtemp <- AAVcomb %>% 
  select(-c(protein, Sample)) %>% 
  pivot_wider(names_from = Sample_AAV, values_from = 'p-values') 

AAVtemp <- AAVtemp%>% 
  arrange(match(loc,str_sort(AAVtemp$loc, numeric = TRUE)))  %>% 
  select(loc, col_order)

AAVfinal <- AAVtemp %>% select(-loc) %>% select(col_order) %>% as.matrix()
rownames(AAVfinal) <- AAVtemp$loc   

#Heatmap conditions
col_fun2 = colorRamp2(c(4, 0 ), c("blue", "white"))
col_fun2(seq(0, 4))
col_fun3 = colorRamp2(c(4, 0, -4 ), c("green", "white", "red"))
col_fun3(seq(-4, 4))

myRows <- c("85-134", "100-149", "110-159", "125-174", #L129F
            "375-424", "385-434", "400-449", "410-459", #E418D
            "485-534", "500-549", "510-559", "525-574", #E531K
            "535-584", "550-599", "560-609", "575-624", #F584L
            "585-634", #A598V
            "600−649", "610-659", "625-674", "635-684") #N642H


# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(AAVfinal) %in% myRows)
fontcolors <- rep('black', nrow(AAVfinal))
fontcolors[row_idx] <- '#8F8F8F'

# Create text annotation object for displaying row names
rowAnno <- rowAnnotation(rows = anno_text(rownames(AAVfinal), gp = gpar(col = fontcolors, fontsize = 8)))

#clean up column names
columnnames <- c("Neg1", "Neg1-AAV6",   "Neg2",   "Neg2-AAV6",   "Neg3",   "Neg3-AAV6",   "Neg4",  "Neg4-AAV6",
                 "Neg5",   "Neg5-AAV6",   "Neg6",   "Neg6-AAV6",  "Neg7",   "Neg7-AAV6",   "Neg8",   "Neg8-AAV6", 
                 "Cas",    "Cas-AAV6",
                  "MTP009", "MTP009-AAV6", "MTP022", "MTP022-AAV6", "MTP030",
                 "MTP030-AAV6", "MTP031", "MTP031-AAV6", "MTP039", "MTP039-AAV6", "MTP043", "MTP043-AAV6", "MTP045", "MTP045-AAV6", "MTP053", "MTP053-AAV6",
                 "MTP054", "MTP054-AAV6", "MTP056", "MTP056-AAV6", "MTP063", "MTP063-AAV6", "MTP065", "MTP065-AAV6", "MTP085", "MTP085-AAV6", "MTP100",
                 "MTP100-AAV6", "MTP103", "MTP103-AAV6", "MTP109", "MTP109-AAV6", "MTP115", "MTP115-AAV6", "MTP135", "MTP135-AAV6", "MTP145", "MTP145-AAV6",
                 "MTP146", "MTP146-AAV6", "MTP147", "MTP147-AAV6", "MTP148", "MTP148-AAV6" ,"MTP154", "MTP154-AAV6", "MTP156", "MTP156-AAV6",
                 "OP",     "OP-AAV6", "OPD",     "OPD-AAV6", "SP"   ,  "SP-AAV6", "SPD"   ,  "SPD-AAV6")

columnnamesF <- every_nth(columnnames, 2, inverse=TRUE)
columnAnno <- HeatmapAnnotation(text = anno_text(columnnamesF, gp = gpar(fontsize = 8)))

#creating a column_split
column_split=NULL
for(i in 1: 41){
  val = 2*i
  column_split[(val-1):val] = paste0(i)
}
column_split <- as.numeric(column_split)

AAVnum <- rep(c("AAV1", "AAV6"), times=41)
col_vector <-(c("AAV1" ="blue", "AAV6" = "grey"))
column_ha = HeatmapAnnotation(AAV=AAVnum, col=list(AAV=col_vector))

#comparison heatmap 1
sample<-Heatmap(AAVfinal, 
                col = col_fun2, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                show_row_names = FALSE,
                show_column_names = FALSE,
                column_split = column_split,
                left_annotation = rowAnno,
                top_annotation = column_ha,
                bottom_annotation = columnAnno,
                column_gap = unit(0, "mm"),
                border = TRUE,
                heatmap_legend_param = list(title = "-log10(p-value)"),  
                column_title = "Comparing AAV1 and AAV6", column_title_gp = gpar(fontsize = 18)) 

pdf(paste0("ComparingAAV1_6.pdf"), width = 12, height = 8)
sample
dev.off() 


# l2fc plots --------------------------------------------------------------
setwd('..')

df <- l2fc %>% 
  select_all(~gsub("Sample_", "", .)) %>% 
  select_all(~gsub("_vs_Neg", "", .)) %>% 
  tibble::rownames_to_column("epitope")

#separate name into protein and location
df <- df %>% separate_wider_delim(epitope, delim="|", names = c("protein", "loc"))

# Extracting relevant subset of data --------------------------------------
#extract AAVs only
AAV <-  df %>% filter(grepl('AAV', protein))
serotypes <- c(1:13, "_DJ")

#extracting for each AAV serotype
AAVsero <- list()
AAVserononeg <- list()
AAVtemp <- NULL
AAVmat <- NULL
AAVserocontrols <- list()

for (i in serotypes){
  sero <- paste0("AAV",i,"_capsid", "")
  AAVsero[[i]] <- AAV %>% filter(grepl(sero, protein))
}

for (i in serotypes){
  AAVtemp <- AAVsero[[i]]
  AAVmat<- AAVtemp %>% select (-protein ) %>% select(-loc) %>% as.matrix()
  rownames(AAVmat) <- AAVtemp$loc
  AAVmat <- AAVmat[str_order(AAVtemp$loc, numeric = TRUE),] 
  AAVsero[[i]] <- AAVmat[, sample_order]
  AAVserononeg[[i]] <- AAVmat[, sample_order2] 
  AAVserocontrols[[i]] <- AAVmat[, sample_order1]
}

#extract for each sample, for different serotypes
AAVtemp2 <- NULL
AAVbysample <- list()
AAVmat2 <- NULL

for (i in sample_order){
  AAVtemp2 <- AAV %>% select(protein, loc, i) %>% 
    pivot_wider(names_from = protein, values_from = i) %>% 
    select_all(~gsub("_capsid", "", .)) 
  AAVmat2 <- AAVtemp2 %>% select(-loc) %>% as.matrix()
  rownames(AAVmat2) <- AAVtemp2$loc
  AAVmat2 <- AAVmat2[str_order(AAVtemp2$loc, numeric = TRUE),]
  AAVbysample[[i]] <- AAVmat2[, caporder] 
}

#extract Cas only
Cas_types <- c("spCas9", "saCas9", "CasRX")
Cas <- list()
Casneg <- list()
Casnoneg <- list()
Castemp <- NULL
Casmat <- NULL

for (i in Cas_types){
  Castemp <-  df %>% 
    filter(grepl(i, protein)) %>% 
    as.data.frame()
  Castemp <- Castemp %>% 
    select (- protein)%>% 
    arrange(match(loc,str_sort(Castemp$loc, numeric = TRUE)))
  Casmat <- Castemp %>% select(-loc) %>% as.matrix()
  rownames(Casmat) <- Castemp$loc
  row.names(Casmat)[row.names(Casmat) =="160-2099"] <- "160-209"
  Cas[[i]] <- Casmat[, sample_order]
  Casnoneg[[i]] <- Casmat[, sample_orderc2] 
  Casneg[[i]] <- Casmat[, sample_orderc1]
}

#to include plots in a different folder
wdir<- getwd()
newwdir <- paste0(wdir, "/plots")
setwd(newwdir)

# Heatmaps ----------------------------------------------------------------
#Heatmap conditions
col_fun1 = colorRamp2(c(1.5, 0, -1.5 ), c("#8F8F8F", "white", "white"))
col_fun2 = colorRamp2(c(1.5, 0, -1.5 ), c("blue", "white", "white"))
col_fun3 = colorRamp2(c(8, 0 ), c("blue", "white"))
col_fun1(seq(-1.5, 1.5))
col_fun2(seq(-1.5, 1.5))
col_fun3(seq(0, 8))

dir1 <- paste0(newwdir, "/l2fc AAV by serotype")
setwd(dir1)

#for loop to generate AAV heatmaps (for individual AAVs, including all samples)
for (i in serotypes){
  #defined column order
  ht_list <-draw( Heatmap(AAVserocontrols[[i]], 
                          col = col_fun1, 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE, 
                          row_names_side = "left",
                          heatmap_legend_param = list(title = "l2fc"),  
                          column_names_gp = grid::gpar(fontsize = 10),
                          row_names_gp = grid::gpar(fontsize = 8)) +
                    Heatmap(AAVserononeg[[i]], 
                            col = col_fun2, 
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE, 
                            row_names_side = "left", 
                            heatmap_legend_param = list(title = "l2fc"),
                            column_names_gp = grid::gpar(fontsize = 10),
                            row_names_gp = grid::gpar(fontsize = 8)), 
                  column_title = paste0("AAV",i,"_capsid"), column_title_gp = gpar(fontsize = 22), 
                  ht_gap = unit(0, "cm"))
  
  
  pdf(paste0("AAV",i,"_capsid", " ", "deforder.pdf"))
  print(ht_list)
  dev.off() 
  
  #still split into 2 different segments, each individually undergoing clustering
  ht_list2 <-draw( Heatmap(AAVserocontrols[[i]], 
                           col = col_fun1, 
                           cluster_rows = FALSE, 
                           cluster_columns = TRUE, 
                           column_dend_reorder = TRUE,
                           row_names_side = "left",
                           heatmap_legend_param = list(title = "l2fc"),  
                           column_names_gp = grid::gpar(fontsize = 10),
                           row_names_gp = grid::gpar(fontsize = 8)) +
                     Heatmap(AAVserononeg[[i]], 
                             col = col_fun2, 
                             cluster_rows = FALSE, 
                             cluster_columns = TRUE, 
                             column_dend_reorder = TRUE,
                             row_names_side = "left", 
                             heatmap_legend_param = list(title = "l2fc"),
                             column_names_gp = grid::gpar(fontsize = 10),
                             row_names_gp = grid::gpar(fontsize = 8)), 
                   column_title = paste0("AAV",i,"_capsid"), column_title_gp = gpar(fontsize = 22), 
                   ht_gap = unit(0, "cm"))
  
  
  pdf(paste0("AAV",i,"_capsid", " ", "clustordersep.pdf"))
  print(ht_list2)
  dev.off() 
  
  #clustered columns
  clust <- Heatmap(AAVsero[[i]], 
                   col = col_fun2, 
                   cluster_rows = FALSE, 
                   cluster_columns = TRUE, 
                   column_dend_reorder = TRUE,
                   row_names_side = "left",
                   heatmap_legend_param = list(title = "l2fc"),
                   column_names_gp = grid::gpar(fontsize = 10),
                   row_names_gp = grid::gpar(fontsize = 8), 
                   column_title = paste0("AAV",i,"_capsid"), column_title_gp = gpar(fontsize = 22))
  
  pdf(paste0("AAV",i,"_capsid", " ", "clustorder.pdf"))
  print(clust)
  dev.off()
}

dir2 <- paste0(newwdir, "/l2fc AAV by sample")
setwd(dir2)

#for loop to generate AAV heatmaps (for individual samples, all AAVs)
for (i in sample_order){
  #defined column order
  sample1 <-Heatmap(AAVbysample[[i]], 
                    col = col_fun2, 
                    cluster_rows = FALSE, 
                    cluster_columns = FALSE, 
                    row_names_side = "left",
                    heatmap_legend_param = list(title = "l2fc"),  
                    column_names_gp = grid::gpar(fontsize = 10),
                    row_names_gp = grid::gpar(fontsize = 8),  
                    column_title = paste0("AAV serotypes for ", i), column_title_gp = gpar(fontsize = 18))
  
  
  pdf(paste0(i,"_for_all_AAVs.pdf"), width = 6, height = 8)
  print(sample1)
  dev.off() 
}

dir3 <- paste0(newwdir, "/l2fc Cas by protein")
setwd(dir3)

#for loop to generate Cas heatmaps
for (i in Cas_types){
  #defined column order
  ht_list <-draw( Heatmap(Casneg[[i]], 
                          col = col_fun1, 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE, 
                          row_names_side = "left",
                          heatmap_legend_param = list(title = "l2fc"),  
                          column_names_gp = grid::gpar(fontsize = 10),
                          row_names_gp = grid::gpar(fontsize = 8)) +
                    Heatmap(Casnoneg[[i]], 
                            col = col_fun2, 
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE, 
                            row_names_side = "left", 
                            heatmap_legend_param = list(title = "l2fc"),
                            column_names_gp = grid::gpar(fontsize = 10),
                            row_names_gp = grid::gpar(fontsize = 6)), 
                  column_title = i, column_title_gp = gpar(fontsize = 22), 
                  ht_gap = unit(0, "cm"))
  
  
  pdf(paste0(i, " ", "deforder.pdf"), width = 8, height = 11)
  print(ht_list)
  dev.off() 
  
  #still split into 2 different segments, each individually undergoing clustering
  ht_list2 <-draw( Heatmap(Casneg[[i]], 
                           col = col_fun1, 
                           cluster_rows = FALSE, 
                           cluster_columns = TRUE, 
                           row_names_side = "left",
                           heatmap_legend_param = list(title = "l2fc"),  
                           column_names_gp = grid::gpar(fontsize = 10),
                           row_names_gp = grid::gpar(fontsize = 8)) +
                     Heatmap(Casnoneg[[i]], 
                             col = col_fun2, 
                             cluster_rows = FALSE, 
                             cluster_columns = TRUE, 
                             row_names_side = "left", 
                             heatmap_legend_param = list(title = "l2fc"),
                             column_names_gp = grid::gpar(fontsize = 10),
                             row_names_gp = grid::gpar(fontsize = 6)), 
                   column_title = i, column_title_gp = gpar(fontsize = 22), 
                   ht_gap = unit(0, "cm"))
  
  
  pdf(paste0(i, " ", "clustordersep.pdf"), width = 8, height = 11)
  print(ht_list2)
  dev.off() 
  
  
  #clustered columns
  clust <- Heatmap(Cas[[i]], 
                   col = col_fun2, 
                   cluster_rows = FALSE, 
                   cluster_columns = TRUE, 
                   column_dend_reorder = TRUE,
                   row_names_side = "left",
                   heatmap_legend_param = list(title = "l2fc"),
                   column_names_gp = grid::gpar(fontsize = 10),
                   row_names_gp = grid::gpar(fontsize = 6), 
                   column_title = i, column_title_gp = gpar(fontsize = 22))
  
  pdf(paste0(i, " ", "clustorder.pdf"), width = 8, height = 11)
  print(clust)
  dev.off()
  
}



setwd('..')


#average mxlp for all serotypes
OP <- AAVbysample$OP %>% 
  as.data.frame()
OP <- tibble::rownames_to_column(OP, "epitopes")
OP <- separate(OP, epitopes, into = c("start", "end"), sep = "-")
#OP <- column_to_rownames(OP, var = "start")

ser1 <- colnames(OP)
ser2 <- ser1[3:16]
ser <- ser2[-5] #remove AAV4, because the CTERM is exactly the same as the second last fragment. Run that separately without the CTERM fragment
cterm <- c(688, 687, 686, 687, 675, 687, 688, 689, 687, 689, 684, 693, 684)
CTERM <- data.frame(serotype = ser, start = cterm)
CTERM <- CTERM %>% mutate(end = start + 49)

fill <- 1:742
mxlpall <- data.frame(epi=fill)

for (i in ser){
  mxlp <- data.frame(epi = fill)
  counts <- data.frame(epi = fill)

  CTERMstart <- CTERM %>% filter(serotype == i) %>% select(start) %>% as.numeric()
  CTERMend <- CTERM %>% filter(serotype == i) %>% select(end) %>% as.numeric()
  
  OPtemp <- OP %>% 
    select(c(start, end, i)) %>% 
    mutate(start = ifelse(row_number()==57, CTERMstart, start)) %>% 
    mutate(end = ifelse(row_number()==57, CTERMend, end))
  
  epistarts <- OPtemp$start
  
  for (j in epistarts){
    j2 <- as.numeric(j)
    end <- j2 + 49
    seq <- (j2+1):(j2+49) 
    value <- OPtemp %>% 
      filter(start == j2) %>% 
      pull(!!i) %>% 
      as.numeric()
    
    names <- c("epi", as.character(j))
    mxlptemp <- data.frame(matrix(nrow = 0, ncol = length(names)))
    counttemp <- data.frame(matrix(nrow = 0, ncol = length(names)))
    colnames(mxlptemp) <- names
    colnames(counttemp) <- names
    
    for (x in 1:742){
      rowtoadd <- ifelse(x %in% seq, value, 0) 
      colname <- as.character(j)
      mxlptemp <- add_row(mxlptemp, epi = x, !!colname := rowtoadd)
      
      counttoadd <- ifelse(x %in% seq, 1, 0)
      colname <- as.character(j)
      counttemp <- add_row(counttemp, epi = x, !!colname := counttoadd)
    }
    
    mxlp <- merge(mxlp, mxlptemp, by = "epi")
    counts <- merge(counts, counttemp, by = "epi")
    
  }
  mxlp <- mxlp %>% 
    mutate(sum = rowSums(across(!epi))) %>% 
    select(c(epi, sum))
  
  counts <- counts %>% 
    mutate(div= rowSums(across(!epi))) %>% 
    select(c(epi, div))
  
  comb <- merge(mxlp, counts) %>% 
    mutate(ave = sum/div) %>% 
    select(c(epi, ave)) 
  
  names(comb) <- c("epi", paste0(i))
  mxlpall <- merge(mxlpall, comb)
}

#AAV4
i = "AAV4"
mxlp <- data.frame(epi = fill)
counts <- data.frame(epi = fill)

OPtemp <- OP %>% 
  select(c(start, end, i)) %>% 
  drop_na(end)

epistarts <- OPtemp$start

for (j in epistarts){
  j2 <- as.numeric(j)
  end <- j2 + 49
  seq <- (j2+1):(j2+49) 
  value <- OPtemp %>% 
    filter(start == j2) %>% 
    select(i) %>% 
    as.numeric()
  
  names <- c("epi", as.character(j))
  mxlptemp <- data.frame(matrix(nrow = 0, ncol = length(names)))
  counttemp <- data.frame(matrix(nrow = 0, ncol = length(names)))
  colnames(mxlptemp) <- names
  colnames(counttemp) <- names
  
  for (x in 1:742){
    rowtoadd <- ifelse(x %in% seq, value, 0) 
    colname <- as.character(j)
    mxlptemp <- add_row(mxlptemp, epi = x, !!colname := rowtoadd)
    
    counttoadd <- ifelse(x %in% seq, 1, 0)
    colname <- as.character(j)
    counttemp <- add_row(counttemp, epi = x, !!colname := counttoadd)
  }
  
  mxlp <- merge(mxlp, mxlptemp, by = "epi")
  counts <- merge(counts, counttemp, by = "epi")
  
}
mxlp <- mxlp %>% 
  mutate(sum = rowSums(across(!epi))) %>% 
  select(c(epi, sum))

counts <- counts %>% 
  mutate(div= rowSums(across(!epi))) %>% 
  select(c(epi, div))

comb <- merge(mxlp, counts) %>% 
  mutate(ave = sum/div) %>% 
  select(c(epi, ave)) 

names(comb) <- c("epi", paste0(i))
mxlpall <- merge(mxlpall, comb)

#export as excel file
write.csv(mxlpall,'mlxp_for_all_sero_epi.csv',row.names=TRUE)

mxlp2 <- mxlpall %>% 
  select(-epi)
min(mxlp2, na.rm=T)
max(mxlp2, na.rm=T)


#repeat for l2fc
l2fc <- l2fc %>% 
  select_all(~gsub("Sample_", "", .)) %>% 
  select_all(~gsub("_vs_Neg", "", .)) %>% 
  tibble::rownames_to_column("epitope") %>% 
  separate_wider_delim(epitope, delim="|", names = c("protein", "loc"))


# Extracting relevant subset of data --------------------------------------
#extract AAVs only
AAVl2fc <-  l2fc %>% 
  filter(grepl('AAV', protein))

OP <- AAVl2fc %>%
  select(protein, loc, OP) %>% 
  as.data.frame()
OP4 <- separate(OP, loc, into = c("start", "end"), sep = "-") %>% 
  pivot_wider(names_from = protein, values_from = OP) %>% 
  select_all(~gsub("_capsid", "", .)) 

OP <- OP4 %>% 
  select(c(start, end, ser))
OP4 <- OP4 %>% 
  select(c(start, end, AAV4))
#OP <- column_to_rownames(OP, var = "start")

cterm <- c(688, 687, 686, 687, 675, 687, 688, 689, 687, 689, 684, 693, 684)
CTERM <- data.frame(serotype = ser, start = cterm)
CTERM <- CTERM %>% mutate(end = start + 49)

fill <- 1:742
mxlpall <- data.frame(epi=fill)

for (i in ser){
  mxlp <- data.frame(epi = fill)
  counts <- data.frame(epi = fill)
  
  CTERMstart <- CTERM %>% filter(serotype == i) %>% select(start) %>% as.numeric()
  CTERMend <- CTERM %>% filter(serotype == i) %>% select(end) %>% as.numeric()
  
  OPtemp <- OP %>% 
    select(c(start, end, i)) %>% 
    mutate(start = ifelse(row_number()==57, CTERMstart, start)) %>% 
    mutate(end = ifelse(row_number()==57, CTERMend, end))
  
  epistarts <- OPtemp$start
  
  for (j in epistarts){
    j2 <- as.numeric(j)
    end <- j2 + 49
    seq <- (j2+1):(j2+49) 
    value <- OPtemp %>% 
      filter(start == j2) %>% 
      pull(!!i) %>% 
      as.numeric()
    
    names <- c("epi", as.character(j))
    mxlptemp <- data.frame(matrix(nrow = 0, ncol = length(names)))
    counttemp <- data.frame(matrix(nrow = 0, ncol = length(names)))
    colnames(mxlptemp) <- names
    colnames(counttemp) <- names
    
    for (x in 1:742){
      rowtoadd <- ifelse(x %in% seq, value, 0) 
      colname <- as.character(j)
      mxlptemp <- add_row(mxlptemp, epi = x, !!colname := rowtoadd)
      
      counttoadd <- ifelse(x %in% seq, 1, 0)
      colname <- as.character(j)
      counttemp <- add_row(counttemp, epi = x, !!colname := counttoadd)
    }
    
    mxlp <- merge(mxlp, mxlptemp, by = "epi")
    counts <- merge(counts, counttemp, by = "epi")
    
  }
  mxlp <- mxlp %>% 
    mutate(sum = rowSums(across(!epi))) %>% 
    select(c(epi, sum))
  
  counts <- counts %>% 
    mutate(div= rowSums(across(!epi))) %>% 
    select(c(epi, div))
  
  comb <- merge(mxlp, counts) %>% 
    mutate(ave = sum/div) %>% 
    select(c(epi, ave)) 
  
  names(comb) <- c("epi", paste0(i))
  mxlpall <- merge(mxlpall, comb)
}

#AAV4
i = "AAV4"
mxlp <- data.frame(epi = fill)
counts <- data.frame(epi = fill)

OPtemp <- OP4 %>% 
  drop_na(end)

epistarts <- OPtemp$start

for (j in epistarts){
  j2 <- as.numeric(j)
  end <- j2 + 49
  seq <- (j2+1):(j2+49) 
  value <- OPtemp %>% 
    filter(start == j2) %>% 
    select(i) %>% 
    as.numeric()
  
  names <- c("epi", as.character(j))
  mxlptemp <- data.frame(matrix(nrow = 0, ncol = length(names)))
  counttemp <- data.frame(matrix(nrow = 0, ncol = length(names)))
  colnames(mxlptemp) <- names
  colnames(counttemp) <- names
  
  for (x in 1:742){
    rowtoadd <- ifelse(x %in% seq, value, 0) 
    colname <- as.character(j)
    mxlptemp <- add_row(mxlptemp, epi = x, !!colname := rowtoadd)
    
    counttoadd <- ifelse(x %in% seq, 1, 0)
    colname <- as.character(j)
    counttemp <- add_row(counttemp, epi = x, !!colname := counttoadd)
  }
  
  mxlp <- merge(mxlp, mxlptemp, by = "epi")
  counts <- merge(counts, counttemp, by = "epi")
  
}
mxlp <- mxlp %>% 
  mutate(sum = rowSums(across(!epi))) %>% 
  select(c(epi, sum))

counts <- counts %>% 
  mutate(div= rowSums(across(!epi))) %>% 
  select(c(epi, div))

comb <- merge(mxlp, counts) %>% 
  mutate(ave = sum/div) %>% 
  select(c(epi, ave)) 

names(comb) <- c("epi", paste0(i))
mxlpall <- merge(mxlpall, comb)

#export as excel file
write.csv(mxlpall,'l2fc_for_all_sero_epi.csv',row.names=TRUE)

mxlp2 <- mxlpall %>% 
  select(-epi)
min(mxlp2, na.rm=T)
max(mxlp2, na.rm=T)

## analysing OP and SP and which are enriched for
dir1 <- paste0(newwdir, "/OPSP proportion of fragments enriched")
setwd(dir1)

Nterm <- c("0-49","10-59","25-74","35-84","50-99","60-109","75-124","85-134","100-149","110-159","125-174","135-184","150-199","160-209","175-224", "185-234", "200-249", "210-259")
Cterm <- c("385-434", "400-449", "410-459", "425-474","435-484","450-499","460-509","475-524","485-534","500-549","510-559","525-574","535-584","550-599","560-609","575-624","585-634", "600-649", "610-659", "625-674")

for (i in serotypes){
  AAVtemp <- AAVsero[[i]]
  AAVOPSP<- AAVtemp %>% 
    as.data.frame() %>% 
    select(c(OP, SP)) %>% 
    mutate(cat = case_when(row.names(.) %in% Nterm ~ "N-term",
                           row.names(.) %in% Cterm ~ "C-term",
                           TRUE ~ "other")) %>% 
    mutate(OP = case_when(OP > 0.0 ~ 1,
                               OP == 0.0 ~ 0)) %>% 
    mutate(SP = case_when(SP > 0.0 ~ 1,
                               SP == 0.0 ~ 0)) 
  OPSPplot <- AAVOPSP %>%
    pivot_longer(!cat, names_to = "Sample", values_to = "count") 
  
  OPSPplot$cat <- factor(OPSPplot$cat, levels = c("N-term", "C-term", "other"))
  
  data_ratio <- OPSPplot %>%
    group_by(cat, Sample) %>%
    summarise(ratio = sum(count == 1) / n())
  
  sample1 <- ggplot(data_ratio, aes(x = cat, y = ratio, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = NULL, y = "Ratio (enriched/all)", title = paste0("AAV",i," - Proportion of fragments enriched")) +
    scale_fill_manual(values = c("OP" = "blue", "SP" = "darkgrey")) + 
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal()
  
  pdf(paste0("AAV",i,"_bargraph.pdf"), width = 6, height = 6)
  print(sample1)
  dev.off() 
}

#create them as combined plot instead:
plot_list <- list()
legend_list <- list()


# Loop through your data
for (i in serotypes) {
  AAVtemp <- AAVsero[[i]]
  AAVOPSP <- AAVtemp %>%
    as.data.frame() %>%
    select(c(OP, SP)) %>%
    mutate(
      cat = case_when(
        row.names(.) %in% Nterm ~ "N-term",
        row.names(.) %in% Cterm ~ "C-term",
        TRUE ~ "other"
      ),
      OP = case_when(OP > 0.0 ~ 1, OP == 0.0 ~ 0),
      SP = case_when(SP > 0.0 ~ 1, SP == 0.0 ~ 0)
    )
  
  OPSPplot <- AAVOPSP %>%
    pivot_longer(!cat, names_to = "Sample", values_to = "count")
  
  OPSPplot$cat <- factor(OPSPplot$cat, levels = c("N-term", "C-term", "other"))
  
  data_ratio <- OPSPplot %>%
    group_by(cat, Sample) %>%
    summarise(ratio = sum(count == 1) / n())
  
  sample_plot <- ggplot(data_ratio, aes(x = cat, y = ratio, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      x = NULL,
      y = "Ratio (enriched/all)",
      title = paste0("AAV", i)
    ) +
    scale_fill_manual(values = c("OP" = "blue", "SP" = "darkgrey")) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() 
  
  # Combine the current plot with the accumulated plots
    plot_list[[i]] <- sample_plot
}

# Create a custom legend for one of the plots (here, we use the first plot)
legend_plot <- get_legend(plot_list[[1]])

# Remove all legends from individual plots
for (i in seq_along(plot_list)) {
  plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
}

# Arrange all plots in a grid
grid <- wrap_plots(plotlist = plot_list, ncol = 3)  # You can adjust ncol as needed

# Combine the grid of plots with the custom legend
final_plot <- grid + legend_plot

# Save the combined plot as a PDF
pdf("combined_plots.pdf", width = 8, height = 12)
final_plot
dev.off()

### Subsetting OP and SP data
subset <- df %>% 
  filter(grepl('AAV', protein)) %>% 
  select(c(protein, loc, HD, SG))
write.csv(subset,'mlxpHDSG.csv',row.names=TRUE)

# Identify the rows to average
rows_to_average <- c("425-474", "435-484", "450-499", "460-509", "475-524", "485-534", "500-549", "510-559", "525-574", "535-584", "550-599", "560-609", "575-624", "585-634", "600-649")

HDsub <- AAVbysample$HD 
average_row <- HDsub[rows_to_average, ]
average_values <- colMeans(average_row[,], na.rm = TRUE)
new_row <- data.frame(t(average_values))
rownames(new_row) <- "Cterm_ave"
HDsub <- rbind(HDsub, new_row)
write.csv(HDsub,'mlxpHD.csv',row.names=TRUE)

SGsub <- AAVbysample$SG 
average_row <- SGsub[rows_to_average, ]
average_values <- colMeans(average_row[,], na.rm = TRUE)
new_row <- data.frame(t(average_values))
rownames(new_row) <- "Cterm_ave"
SGsub <- rbind(SGsub, new_row)
write.csv(SGsub,'mlxpSG.csv',row.names=TRUE)

HDsmall <- HDsub %>%
  select(c(protein, Cterm_ave)) %>% 
  rename(HDCterm_ave = Cterm_ave)

SGsmall <- SGsub %>%
  select(c(protein, Cterm_ave)) %>% 
  rename(SGCterm_ave = Cterm_ave)

summary <- merge(HDsmall, SGsmall, by = "protein", all=TRUE)
write.csv(summary,'mlxpHDSGCterm_ave.csv',row.names=TRUE)


