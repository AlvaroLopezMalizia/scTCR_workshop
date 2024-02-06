
#Disclaimer----
#this script is only meant for teaching purposes
#it is not an exhaustive analysis of a dataset
#the parameters, thresholds and procedures here should not be directly used in a research task without serious consideration
#in a real situation, every decision taken needs to be thoroughly considered,
#using deep knowledge about how the data was collected, and what questions we may want to answer
#about the organization of the script, this is merely the way I work and
#it's surely plagued with things that can be optimized and/or personalized to your needs and likings

#INTRODUCTION TO THE WORKSHOP----
#On the top-right corner of the script windows, to the right of the "Source" dropdown menu,
#you will find the outline button, press it, it will unfold a clickable index to different parts of the script
#Always try to keep your scripts organized in sections and subsections
#you can create a section by starting a line with a hashtag and ending it with 4 -
#this is a section----
#you can create subsections by adding more hashtags
##this is nested section----
###this is another nested section----

#I tend to name objects starting with temp_ 
#for global environment objects that I will sooner or later want to delete to clean up.
#I use the starting string for_ to identify objects created within different levels of a loop
#because it helps me testing the loop level by level and then easily removing the loop temporal objects
#I sometimes use the starting string fun_ to name objects created within a function
#so that I can test the function part by part and then easily remove those objects
#I may use a starting string named after the name of the section when I need that object
#to remain in the working environment for the following sections

#I chose not to use Rmarkdown for this workshop because it can get buggy and or very slow
#and, in reality, in our lab we rarely use it

#01) SCRIPT SETUP----
#set the working directory to the location of this script automatically
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#check that the working directory is indeed the location of this script
getwd()
#if something went wrong, manually set the directory to the location of this script
#using setwd() function, for example:
#setwd('path_to_directory/TCR_bundle')

library(Seurat) #our rnaseq base pipeline
library(dplyr) # %>% pipes for data management 
library(tidyr) #more tools for data management
library(clustree) #explore nested clustering
library(patchwork) #arrange plots into one thing
library(ggplot2) #standard plotting library
library(DescTools) #Gini index
library(divo) #morisita horn index
library(ggrepel) #repel labels on ggplot
library(scGate) #auto annotation
library(UCell) #auto annotation dependency
library(stringdist) #string manipualtion
library(pheatmap) #simple heatmap library (complexHeatmaps is the one to go if you need more complex things)

#Set random number generator seed for reproducibility
set.seed(123)

#02) LOAD DATA----

#load RNAseq counts, cell names and gene names from cellranger output
temp_rna <- Read10X(data.dir = "scTCR_workshop_data/filtered_gene_bc_matrices")

#"project" will set the sample name of this seurat object
#avoid using strings that start with a number, because that will make downstream analysis more complicated in many situations
temp_seurat <- CreateSeuratObject(counts = temp_rna, project = "s01",
                              min.cells = 5, min.features = 300)

#when working with multiple samples and combining their data,
#the same cell barcodes can be found in different samples more often than expected
#seurat is supposed to properly manage this issue, but I prefer to solve it manually
#by pasting the sample name stored in orig.ident to the cell barcodes

colnames(temp_seurat) <- paste0(temp_seurat$orig.ident, ':', colnames(temp_seurat))

##add MT content----
#the ^ indicates that the genes need to start with that pattern
#sometimes the "-" character is not allowed in feature names
#and it gets replaced with "."
#in that case, the proper search pattern is going to be "^MT\\."
#"\\" is needed then because "." is a special character that has other meaning otherwise
temp_seurat[["percent.mt"]] <- PercentageFeatureSet(temp_seurat, pattern = "^MT-")

#notice how temp_seurat[["percent.mt"]] results in a dataframe instead of in a vector
class(temp_seurat[["percent.mt"]])
dim(temp_seurat[["percent.mt"]])

#then check that we got information that makes sense
plot(density(log10(temp_seurat[["percent.mt"]][[1]])))

#we can imagine this MT content distribution as the combination of two distributions:
#one coming from the natural MT content of the cells, 
#and one coming from the MT content derived from mithochondrial rupture. 
#This is expected in scRNAseq

##TCR data----
#load tcr data from the csv output
temp_tcr <- read.delim('scTCR_workshop_data/vdj_v1_hs_nsclc_t_filtered_contig_annotations.csv', sep = ',')

#Select only full length and productive contigs
#in older cellranger versions, the TCR data is not filtered by full length or being productive
temp_tcr <- temp_tcr[grep("[Tt]rue",temp_tcr$full_length),]
temp_tcr <- temp_tcr[grep("[Tt]rue",temp_tcr$productive),]

#Select rows with TRA or TRB chains
temp_tcr <- temp_tcr[which((temp_tcr$chain=="TRA")|(temp_tcr$chain=="TRB")),]

temp_tcr=temp_tcr[order(temp_tcr$chain),]

#Notice that each droplet barcode will be repeated several times 
#depending on the TCR information found in it. 
#This needs to be processed in order to make it compatible with our scRNAseq, 
#where each barcode is present only once

#we can include the V gene information in a string along with the sequence information
#and we can include the V gene information when defining the clonotypes
temp_tcr <- temp_tcr %>% 
  mutate('v_cdr3_nt'=paste0(v_gene,"_",cdr3_nt),
         'full_cdr3_aa' = paste0(chain,":", v_cdr3_nt),
         'full_cdr3_nt' = paste0(chain,":", v_cdr3_nt),
         'b_cdr3_aa' = ifelse(chain == 'TRB', cdr3, NA),
         'b_cdr3_nt' = ifelse(chain == 'TRB', cdr3_nt, NA),
         'b_v_cdr3_nt' = ifelse(chain == 'TRB', v_cdr3_nt, NA)) %>% 
  group_by(barcode) %>% 
  mutate(
    'nTRA' = sum(chain == 'TRA'),
    'nTRB' = sum(chain == 'TRB')) %>% 
  mutate('tcr_data' = 
    ifelse(nTRA > 2 | nTRB > 1, 'tcr_doublet', 'tcr_singlet')
  ) %>% 
  mutate_at(vars('full_cdr3_aa', 'full_cdr3_nt', 'chain'), ~ paste0(., collapse = ';')) %>% 
  ungroup() %>% 
  dplyr::select(barcode, v_gene, full_cdr3_aa, full_cdr3_nt, b_cdr3_aa, b_cdr3_nt, b_v_cdr3_nt, nTRA, nTRB, tcr_data) %>% 
  dplyr::filter(!(is.na(b_cdr3_aa) & nTRB > 0)) %>% 
  distinct(barcode, .keep_all = TRUE)

#explore the information that we generated
#we have defined as doublets the droplets that contain more than 1 TRB chain, why?

##add tcr to seurat----
#we need to generate the same type of barcode that includes the sample name 's01'
#this will prevent mixing TCR information between samples when there is a repeated barcode
temp_tcr$barcode[1]
colnames(temp_seurat)[1]
temp_tcr$cell_id <- paste0('s01', ':', temp_tcr$barcode)

#I like adding a $cell_id column to the seurat metadata
temp_seurat$cell_id <- colnames(temp_seurat)

#whenever I edit the seurat metadata, I prefer to use a temporary object in case I mess up
temp_metadata <- left_join(temp_seurat@meta.data, temp_tcr, by = 'cell_id')
temp_metadata$tcr_data[is.na(temp_metadata$tcr_data)] <- 'tcr_absence'
rownames(temp_metadata)[1:5]
rownames(temp_metadata) <- temp_metadata$cell_id

#it's important to make sure that we haven't altered the metadata row order or names
#and that the metadata is a data.frame and not a tibble or something else
#because otherwise the internal connectivity of the seurat object may be broken

all(rownames(temp_metadata) == colnames(temp_seurat))
temp_seurat@meta.data <- temp_metadata

##remove doublets----
#part of the quality control involves removing technical doublets
#it's common practice to consider droplets with more than one TRB chain as a doublet
temp_seurat <- subset(temp_seurat, cells = temp_seurat$cell_id[temp_seurat$tcr_data != 'tcr_doublet'])
dim(temp_seurat)

##TCR metadata----
#after we have removed the droplets with more than one TRB chain,
#we can use the V gene + TRB chain to define the clonotypes
#why not using both TRA and TRB chains to define the clonotype?
#in part because it's easier to use only the TRB chain

#what is a clonotype?
#using 10x definition:
#Clonotype.
#A set of adaptive immune cells that are the clonal progeny of a fully recombined, 
#unmutated common ancestor.

names(temp_seurat@meta.data)
temp_tcrcount <- temp_seurat@meta.data %>% 
  group_by(orig.ident) %>% 
  mutate('n_droplets_sample' = n()) %>% 
  ungroup() %>% 
  dplyr::filter(!is.na(b_v_cdr3_nt)) %>% 
  group_by(b_v_cdr3_nt) %>% 
  mutate('b_v_cdr3_nt_count' = n()) %>% 
  ungroup() %>% 
  distinct(b_v_cdr3_nt, .keep_all = TRUE) %>% 
  mutate('rank' = row_number(desc(b_v_cdr3_nt_count))) %>% 
  mutate('rank' = sprintf('%04d', rank)) %>% 
  mutate('b_nt_id' =  paste0('b_nt_', rank, '_', round(b_v_cdr3_nt_count, 0))) %>% 
  dplyr::select(b_v_cdr3_nt, b_v_cdr3_nt_count, b_v_cdr3_nt_count, b_nt_id) %>% 
  as.data.frame()

#explore the information generated, notice how we've put a name to each clonotype
#that can be easily ordered, contains its ranking in terms of expansion relative
#to all other clonotypes, and also contains the number of droplets detected
#with that clonotype

#explore tcr abundance
table(temp_tcrcount$b_v_cdr3_nt_count)

#this is the normal thing to see, with a lot of clonotypes with just one clone, then few to none
#clonotypes observed as we increase the abundance, and finally a small peak of several highly abundant clonotypes

#to make observation easier, we can categorize the clonotypes according to their abundance
temp_tcrcount$b_cdr3_nt_cat <- cut(x = temp_tcrcount$b_v_cdr3_nt_count, c(0, 1, 5, 10, 20, Inf))
table(temp_tcrcount$b_cdr3_nt_cat)

#next we should add this information to the metadata of the seurat object, so that
#we can have everything neatly organized

temp_metadata <- left_join(temp_seurat@meta.data, temp_tcrcount, by = 'b_v_cdr3_nt')
rownames(temp_metadata)[1:5] #still lost rowname info
rownames(temp_metadata) <- temp_metadata$cell_id #this is one of the reasons I like keeping a $cell_id column
all(rownames(temp_metadata) == rownames(temp_seurat@meta.data))

#why not Seurat::AddMetaData?
#it has limited functionalities for my taste. I'd just rather treat the metadata as any other dataframe
#and put extra attention on not messing up the cell names or order
#I try not to directly modify the seurat metadata until I am sure that the modifications are correct
#so always be very careful not to be messing up the order of the cells,
#and only then I replace it. Although this takes extra memory

table(temp_metadata$b_cdr3_nt_cat, useNA = 'ifany')
table(temp_metadata$b_cdr3_nt_cat, temp_metadata$b_v_cdr3_nt_count, useNA = 'ifany')
temp_seurat@meta.data <- temp_metadata

#Why didn't I perform a low quality cells control?
#because we have a very heterogeneous mix of cell types
#each cell type can have a different threshold for what is considered low quality
#Ideally, I'd split the dataset into major celltypes and perform quality control separately for each on of them

#how many cells with RNA and now TCR data do we have?
table(temp_seurat$tcr_data)

#03) PROCESSING----
#If memory is sufficient, I like keeping a backup of the seurat object, 
#then clear the working environment
#and create new temp objects to work on
#if memory is insufficient, you can always store the object on a temporary file on-disk
#and load it as needed
p02_seurat <- temp_seurat
rm(list=ls(pattern = 'temp_'))
gc() #and garbage collect to make sure that objects are properly unloaded from memory

#clustering
temp_seurat <- NormalizeData(p02_seurat)
temp_seurat <- FindVariableFeatures(temp_seurat, nfeatures = 2000)
temp_varfeats <- VariableFeatures(temp_seurat)
temp_keepindex <- grepl("^trav|^trac|^traj|^trbv|^trbc|^trbj|^tradv", tolower(temp_varfeats))
temp_varfeats <- temp_varfeats[!temp_keepindex]
VariableFeatures(temp_seurat) <- temp_varfeats
temp_seurat <- ScaleData(temp_seurat)
temp_seurat <- RunPCA(temp_seurat, npcs = 30)

temp_seurat <- FindNeighbors(temp_seurat, reduction = "pca", dims = 1:30)
temp_seurat <- FindClusters(temp_seurat, resolution = seq(0, 1, 0.1), random.seed	= 1234)

names(temp_seurat@meta.data)

#seurat default cluster names have two issues:
#they start with a number, and eventually we will have those as column names, which is ill advised
#the cluster numbers will be treated as characters for plotting order, 
#generating ordering such as 1, 10, 11, 2, 3; instead of 1, 2, 3, 10, 11

#rename clustering variable
temp_clusterindex <- grep('RNA_snn_res.', names(temp_seurat@meta.data))
temp_clusternames <- names(temp_seurat@meta.data)[temp_clusterindex]

#fix clusternames
for(temp_res in temp_clusternames){
  
  for_newnames <- sprintf("%02s", as.character(temp_seurat@meta.data[[temp_res]]))
  for_newnames <- paste0('c', for_newnames)
  temp_seurat@meta.data[[temp_res]] <- for_newnames
  
}
rm(list=ls(pattern = '^for_')) #the clearer you keep the env, the better

names(temp_seurat@meta.data)
table(temp_seurat@meta.data[[temp_res]])

#umap
temp_seurat <- RunUMAP(temp_seurat, reduction = "pca", dims = 1:30, reduction.name = "umap")

#the heaviest processing part is mostly over
#if everything worked just fine, we can replace the backup copy of the object and clear up the working environ

p03_seurat <- temp_seurat
rm(list=ls(pattern = '^temp_'))
rm(p02_seurat)
gc()

#04) EXPLORATION----
#again my working copy of the object
#this is memory inefficient, but it's just a way to be able to rollback in case I screw up something

temp_seurat <- p03_seurat

##plots clustree----

#create a clustree
temp_clustree <-   clustree(temp_seurat, 
                            prefix="RNA_snn_res.",
                            highlight_core = TRUE, node_size=5, 
                            edge_arrow = TRUE,edge_arrow_ends="first", 
                            edge_width = 0.5)


#create dimplots of each resolution
temp_clusteres <- grepl('^RNA_snn_res.', names(temp_seurat@meta.data))
temp_clusteres <- names(temp_seurat@meta.data)[temp_clusteres]

temp_dimplots <- lapply(temp_clusteres, function(one_clusteres){
  
  DimPlot(object = temp_seurat,
          reduction = 'umap',
          group.by = one_clusteres,
          label = TRUE) +
    ggtitle(one_clusteres)
  
})

#using patchwork we can create a single organized plot
temp_plot <- wrap_plots(temp_dimplots, nrow = 2)
temp_plot <- wrap_plots(temp_clustree, temp_plot, ncol = 2, widths = c(2, length(temp_dimplots)))
temp_plot <- temp_plot

#print clustree with umaps of every resolution
png("scTCR_workshop_clustree_1.png", width = 5+10*length(temp_dimplots), height = 30, units = 'cm', res = 320)
plot(temp_plot)
dev.off()

#we can't say much right now about this clustree
#we need to consider this along with other variables to see which resolution is better suited to our needs

##variable plots-----
temp_genes <- c('CD2','CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B',
                'MKI67', 'TOP2A',
                'IGHD','CD38', 'CD19', 'MS4A1',
                'ITGA2B', 'GP1BA', 'ITGB3', 'PECAM1',
                'NCAM1', 'NKG7', 'FCGR3A', 'KLRC1',
                'CD14', 'CD86', 'CD14', 'LYZ',
                'CD68', 'CD163', 'MRC1', 'ITGAX',
                'CD80','CD83','CD1C', 'CD83',
                'ITGAM', 'CCR3', 'IL3RA', 'IL5RA',
                "CR1", "ITGAM", "FUT4", "FCGR3B",
                "ENPP3", "ANPEP", "MS4A3", "CD69",
                "IL7R", "KIT", "IL2RA", "THY1",
                "SDC1", 'XBP1', "CD27", "IRF4",
                'MALAT1', 'nCount_RNA', 'nFeature_RNA', 'percent.mt')

temp_genes %in% rownames(temp_seurat)

temp_gene_subtitles <- c('T cells', 'T cells', 'T cells', 'T cells' , 'T cells' , 'T cells',
                         'Proliferation', 'Proliferation',
                         'B cells', 'B cells', 'B cells', 'B cells',
                         'platelets', 'platelets', 'platelets', 'platelets',
                         'NK phenotype', 'NK phenotype', 'NK phenotype', 'NK phenotype',
                         'Monocytes', 'Monocytes', 'Monocytes', 'Monocytes',
                         'Macrophages', 'Macrophages', 'Macrophages', 'Macrophages',
                         'Dendritic', 'Dendritic', 'Dendritic', 'Dendritic',
                         'Eosinophils', 'Eosinophils', 'Eosinophils', 'Eosinophils',
                         'Neutrophils', 'Neutrophils', 'Neutrophils', 'Neutrophils',
                         'Basophils', 'Basophils', 'Basophils', 'Basophils',
                         'innate like', 'innate like', 'innate like', 'innate like',
                         'plasma', 'plasma', 'plasma', 'plasma',
                         'QC', 'QC', 'QC', 'QC')

length(temp_gene_subtitles) == length(temp_genes)

names(temp_gene_subtitles) <- temp_genes

#feature plots
#You can pass a list to features to FeaturePlot
#but I want to add a subtitle to each of the plots
#so I plot each feature separately and then combine them with patchwork
temp_feat_list <- list()
for(for_gene in temp_genes){
  
  temp_feat_list[[for_gene]] <- FeaturePlot(object = temp_seurat,
                                                 features = for_gene,
                                                 min.cutoff = 'q5', #the cutoffs sometimes make it easier to see imporant things
                                                 max.cutoff = 'q95',
                                                 order = T) + ggtitle(label = for_gene, 
                                                                      subtitle = temp_gene_subtitles[for_gene])
  
  
}

temp_feat_plot <- wrap_plots(temp_feat_list, ncol = 4)

png("scTCR_workshop_features_1.png", width = 60, height = length(temp_feat_list) * 2.5, units = 'cm', res = 120)
plot(temp_feat_plot)
dev.off()

rm(list=ls(pattern = "^for_"))

##plots by TCR----
#The TCR themselves are also helpful at understanding our initial UMAP
#so we should plot them too

table(temp_seurat$b_cdr3_nt_cat, useNA = 'ifany')

#additionally I'd like to individualize the top 4 most expanded clonotypes
temp_selection <- sort(unique(temp_seurat$b_nt_id), decreasing = FALSE)[1:5]

temp_seurat$b_nt_id_top <- temp_seurat$b_nt_id
temp_seurat$b_nt_id_top[!(temp_seurat$b_nt_id_top %in% temp_selection)] <- NA

temp_selection <- sort(unique(temp_seurat$b_nt_id), decreasing = FALSE)[6:11]

temp_seurat$b_nt_id_top2 <- temp_seurat$b_nt_id
temp_seurat$b_nt_id_top2[!(temp_seurat$b_nt_id_top2 %in% temp_selection)] <- NA

##auto annotation----
#This auto annotation tool will let us perform a quick broad anotation on our dataset
#I haven't experimented much with it, but it doesn't seem to be too accurate to subdivide each immune compartment
#into specific cell types, but it does a very nice job at separating broad immune types

temp_sc_models <- readRDS('scTCR_workshop_data/scGate_models.RDS')

temp_seurat <- scGate(temp_seurat, 
                      model = temp_sc_models$human$TME_HiRes, 
                      reduction = "pca", ncores = 4)

##extra dimplots----
temp_tcrplots <- DimPlot(temp_seurat,
                         group.by = c('b_cdr3_nt_cat', 'b_nt_id_top2', 'b_nt_id_top', 'scGate_multi'),
                         pt.size = 1,
                         order = TRUE,
                         alpha = 0.7,
                         ncol = 2,
                         label = TRUE)

png("scTCR_workshop_umaps_1.png", width = 40, height = length(temp_tcrplots) * 7.5, units = 'cm', res = 120)
plot(temp_tcrplots)
dev.off()

#it's clear that we have TCR data way beyond Tcells
#like the cycling cluster has an obvious explanation, they are all cycling cells, including Tcells
#The NK-phenotype cluster can be a mix of actually NK cells, NK-T cells and Gamma_delta T cells
#since this dataset doesn't have the TRG and TRD genes, spotting the gama delta can be a bit complicated
#then we have multiple TCR clonotypes spread all over other celltypes
#what are the possible explanations for this?

#05) ASSESMENT----

#Closing the TCR debate, we can talk a bit about the quality of our clustering

table(temp_seurat$scGate_multi, temp_seurat$RNA_snn_res.1)
#we can see some potential mixing of CD4 and CD8 cells, sometimes also with B cells. This isn't necessarily true
#scGate could be mislabeling some cells, but if there is discrepancy between the two classifications, then we should
#pay more attention to what is going on

table(temp_seurat$b_nt_id_top, temp_seurat$RNA_snn_res.1)
#It's not ideal to have one of the most expanded clonotypes of the CD8T cells present in the CD4T cells

table(temp_seurat$b_nt_id_top2, temp_seurat$RNA_snn_res.1)

#we are getting a bit of discrepancies
#I am kind of worried about the apparent B cells and NK cells getting mixed with Tcells
#one clonotype with a degree of sharing between cluster 01 of tcd4 and 04 of tcd8 may
#be indicating some kind of problem with the clustering
#I see low quality cells in the CD4 area that I would like to get clustered apart
#so I'll subset the seurat object to Tcells, NK cells and other clusters containing rather abundant TCR data

##quality control thresholds with celltypes----
#why didn't I just subset it all to cells with TCR data?
#because we will be skewing Tcell types abudances

table(temp_seurat$RNA_snn_res.1, temp_seurat$tcr_data) / rowSums(table(temp_seurat$RNA_snn_res.1, temp_seurat$tcr_data))

#changing topic, earlier I sad that we need to be careful with harsh quality control thresholds
#now lets check the scatterplots of qc measures by cell type

png("scTCR_workshop_qc_celltype.png", width = 20, height = 15, units = 'cm', res = 120)
ggplot(data = temp_seurat@meta.data, aes(x = nCount_RNA, y = percent.mt, color = scGate_multi)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()
dev.off()

#lets save progress and clear space
p05_seurat <- temp_seurat
rm(list=ls(pattern = 'temp_'))
rm(p03_seurat)
gc()
temp_seurat <- p05_seurat

#06) SUBSET SEURAT----
#before that I want to save some of the information and keep it along the next stage, to have some kind of comparison
names(temp_seurat@meta.data)

#probably a very inefficient and prone to errors way of selecting desired columns
temp_seurat@meta.data <- temp_seurat@meta.data[,c("cell_id", "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "barcode", 
                                                  "v_gene", "full_cdr3_aa", "full_cdr3_nt", "b_cdr3_aa", "b_cdr3_nt", "b_v_cdr3_nt",
                                                  "nTRA", "nTRB", "tcr_data", "b_v_cdr3_nt_count", "b_nt_id", "b_cdr3_nt_cat", "b_nt_id_top", "b_nt_id_top2",
                                                  "RNA_snn_res.1", "scGate_multi")]
rownames(temp_seurat@meta.data)[1]
names(temp_seurat@meta.data)

#we can change the name of the variables that we will re-compute and would otherwise get overwritten, but we want to keep
temp_metadata <- temp_seurat@meta.data %>% 
  rename('full_RNA_snn_res.1' = 'RNA_snn_res.1',
         'full_scGate_multi' = 'scGate_multi')

names(temp_metadata)
rownames(temp_metadata)[1]
class(temp_metadata) #when using dplyr sometimes the output is a tibble, which is not always compatible with functions designed for seurat objects

temp_seurat@meta.data <- temp_metadata[temp_seurat@meta.data$cell_id,] #this is a good way of making sure the new metadata is in the same order as the old one

all(colnames(temp_seurat) == temp_seurat@meta.data$cell_id)

##subset seurat and recluster----
#im keeping cluster 19 because it was very close to the Tcells cluster, and c14 for it's relationship with the cluster containing the cycling Tcells
temp_cells <- temp_seurat$cell_id[temp_seurat$full_RNA_snn_res.1 %in% c('c01', 'c04', 'c09', 'c10', 'c14', 'c16','c17', 'c19')]
temp_seurat <- subset(temp_seurat, cells = temp_cells)
dim(temp_seurat)
names(temp_seurat@meta.data)

#clustering
temp_seurat <- NormalizeData(temp_seurat) 
temp_seurat <- FindVariableFeatures(temp_seurat, nfeatures = 2000) 
#now that our dataset is different, the variable features will be different
#and we want our new variable features to be more focused on the variability of our t-cell enriched dataset
temp_varfeats <- VariableFeatures(temp_seurat)
temp_keepindex <- grepl("^trav|^trac|^traj|^trbv|^trbc|^trbj|^tradv", tolower(temp_varfeats))
temp_varfeats <- temp_varfeats[!temp_keepindex]
VariableFeatures(temp_seurat) <- temp_varfeats
temp_seurat <- ScaleData(temp_seurat)
temp_seurat <- RunPCA(temp_seurat, npcs = 30)

ElbowPlot(temp_seurat, ndims = 30)
#30 dims may be a bit of an overkill

#umap
temp_seurat <- RunUMAP(temp_seurat, reduction = "pca", dims = 1:20, reduction.name = "umap")
temp_20dims <- DimPlot(temp_seurat,
                       group.by = c('full_RNA_snn_res.1', 'full_scGate_multi', 'b_cdr3_nt_cat', 'b_v_cdr3_nt_count', 'b_nt_id_top2', 'b_nt_id_top'),
                       ncol = 6)

temp_seurat <- RunUMAP(temp_seurat, reduction = "pca", dims = 1:30, reduction.name = "umap")
temp_30dims <- DimPlot(temp_seurat,
                       group.by = c('full_RNA_snn_res.1', 'full_scGate_multi', 'b_cdr3_nt_cat', 'b_v_cdr3_nt_count', 'b_nt_id_top2', 'b_nt_id_top'),
                       ncol = 6)

temp_dimchoice_plot <- wrap_plots(temp_20dims, temp_30dims, ncol = 1) + plot_annotation(subtitle = '20 dims on first row, 30 dims on second row')

png("scTCR_workshop_dimchoice.png", width = 70, height = 20, units = 'cm', res = 120)
temp_dimchoice_plot
dev.off()

#there isn't much difference, and people usually like more dims, so we can keep 30

temp_seurat <- FindNeighbors(temp_seurat, reduction = "pca", dims = 1:30)
temp_seurat <- FindClusters(temp_seurat, resolution = seq(0, 1, 0.1), random.seed	= 1234)

names(temp_seurat@meta.data)

#Again we will rename the seurat clusters

#rename clustering variable
temp_clusterindex <- grep('^RNA_snn_res.', names(temp_seurat@meta.data)) #notice the ^ to make sure I don't match the full_RNA_snn_res.1
temp_clusternames <- names(temp_seurat@meta.data)[temp_clusterindex]

#fix clusternames
for(temp_res in temp_clusternames){
  
  temp_for_newnames <- sprintf("%02s", as.character(temp_seurat@meta.data[[temp_res]]))
  temp_for_newnames <- paste0('c', temp_for_newnames)
  temp_seurat@meta.data[[temp_res]] <- temp_for_newnames
  
}
rm(list=ls(pattern = 'temp_for')) #the clearer you keep the env, the better

names(temp_seurat@meta.data)
table(temp_seurat@meta.data[[temp_res]])


##plot clustree----
#create a clustree
names(temp_seurat@meta.data)

temp_clustree <-   clustree(temp_seurat, 
                            prefix="RNA_snn_res.",
                            highlight_core = TRUE, node_size=5, 
                            edge_arrow = TRUE,edge_arrow_ends="first", 
                            edge_width = 0.5)


#create dimplots of each resolution
temp_clusteres <- grepl('^RNA_snn_res.', names(temp_seurat@meta.data))
temp_clusteres <- names(temp_seurat@meta.data)[temp_clusteres]

temp_dimplots <- lapply(temp_clusteres, function(one_clusteres){
  
  DimPlot(object = temp_seurat,
          reduction = 'umap',
          group.by = one_clusteres,
          label = TRUE) +
    ggtitle(one_clusteres)
  
})

#using patchwork we can create a single organized plot
temp_plot <- wrap_plots(temp_dimplots, nrow = 2)
temp_plot <- wrap_plots(temp_clustree, temp_plot, ncol = 2, widths = c(2, length(temp_dimplots)))
temp_plot <- temp_plot

#print clustree
png("scTCR_workshop_clustree_2.png", width = 5+10*length(temp_dimplots), height = 30, units = 'cm', res = 320)
plot(temp_plot)
dev.off()

#the clustree looks pretty stable on the bottom half, but there is one thing that we should check:
#between resolutions 0.8 and 0.9 there is one cluster that dissapears
#the clustree indicates that this cluster 9 of resoltuion 0.8 is getting merged into cluster 3 of resolution 0.9
#we can see the same thing in the umap
#we can't really say a lot more with this information so we need more information
#we can start by looking at the same variables as before

##variable plots-----
temp_genes <- c('CD2','CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B',
                'MKI67', 'TOP2A',
                'IGHD','CD38', 'CD19', 'MS4A1', 'TRBC1', 'TRBC2',
                'ITGA2B', 'GP1BA', 'ITGB3', 'PECAM1',
                'NCAM1', 'NKG7', 'FCGR3A', 'KLRC1',
                'CD14', 'CD86', 'CD14', 'LYZ',
                'CD68', 'CD163', 'MRC1', 'ITGAX',
                'CD80','CD83','CD1C', 'CD83',
                'ITGAM', 'CCR3', 'IL3RA', 'IL5RA',
                "CR1", "ITGAM", "FUT4", "FCGR3B",
                "ENPP3", "ANPEP", "MS4A3", "CD69",
                "IL7R", "KIT", "IL2RA", "THY1",
                "SDC1", 'XBP1', "CD27", "IRF4",
                'MALAT1', 'nCount_RNA', 'nFeature_RNA', 'percent.mt')

temp_genes %in% rownames(temp_seurat)

temp_gene_subtitles <- c('T cells', 'T cells', 'T cells', 'T cells' , 'T cells' , 'T cells',
                         'Proliferation', 'Proliferation',
                         'B cells', 'B cells', 'B cells', 'B cells', 'B cells', 'B cells',
                         'platelets', 'platelets', 'platelets', 'platelets',
                         'NK phenotype', 'NK phenotype', 'NK phenotype', 'NK phenotype',
                         'Monocytes', 'Monocytes', 'Monocytes', 'Monocytes',
                         'Macrophages', 'Macrophages', 'Macrophages', 'Macrophages',
                         'Dendritic', 'Dendritic', 'Dendritic', 'Dendritic',
                         'Eosinophils', 'Eosinophils', 'Eosinophils', 'Eosinophils',
                         'Neutrophils', 'Neutrophils', 'Neutrophils', 'Neutrophils',
                         'Basophils', 'Basophils', 'Basophils', 'Basophils',
                         'innate like', 'innate like', 'innate like', 'innate like',
                         'plasma', 'plasma', 'plasma', 'plasma',
                         'QC', 'QC', 'QC', 'QC')

length(temp_gene_subtitles) == length(temp_genes)

names(temp_gene_subtitles) <- temp_genes

#feature plots
temp_feat_list <- list()
for(temp_for_gene in temp_genes){
  
  temp_feat_list[[temp_for_gene]] <- FeaturePlot(object = temp_seurat,
                                                 features = temp_for_gene,
                                                 min.cutoff = 'q5', #the cutoffs sometimes make it easier to see imporant things
                                                 max.cutoff = 'q95',
                                                 order = T) + ggtitle(label = temp_for_gene, 
                                                                      subtitle = temp_gene_subtitles[temp_for_gene])
  
  
}

temp_feat_plot <- wrap_plots(temp_feat_list, ncol = 4)

png("scTCR_workshop_features_2.png", width = 60, height = length(temp_feat_list) * 2.5, units = 'cm', res = 420)
plot(temp_feat_plot)
dev.off()

##auto annotation----
temp_sc_models <- readRDS('scTCR_workshop_data/scGate_models.RDS')

temp_seurat <- scGate(temp_seurat, 
                      model = temp_sc_models$human$TME_HiRes, 
                      reduction = "pca", ncores = 4)

#is our new auto annotation concordant with the previous one?

table(temp_seurat$scGate_multi, temp_seurat$full_scGate_multi)

#there is evidently some confusion between CD8 and CD4 cells
#so when plotting the new classification, we might want to also use the old one to check how things changed

##plots by TCR----
table(temp_seurat$b_cdr3_nt_cat, useNA = 'ifany')

#now that we have subset the dataset and changed the total counts of different clonotypes
#should we use the original counts or should we re-count them?
#debatable
#for simplicity, I am going to be using the things we calculated before

##extra dimplots----
#we can look at what the data looks like in our new reduced dataset
temp_tcrplots <- DimPlot(temp_seurat,
                         group.by = c('b_cdr3_nt_cat', 'b_nt_id_top2', 'b_nt_id_top', 'b_v_cdr3_nt_count', 'scGate_multi', 'full_scGate_multi', 'full_RNA_snn_res.1', 'RNA_snn_res.0.8'),
                         pt.size = 1,
                         order = TRUE,
                         alpha = 0.7,
                         ncol = 2)

png("scTCR_workshop_umaps_2.png", width = 40, height = length(temp_tcrplots) * 7.5, units = 'cm', res = 420)
plot(temp_tcrplots)
dev.off()

table(temp_seurat$scGate_multi, temp_seurat$RNA_snn_res.1)
table(temp_seurat$b_nt_id_top, temp_seurat$RNA_snn_res.1)
table(temp_seurat$b_nt_id_top2, temp_seurat$RNA_snn_res.1) 

#checking what happened to the cycling group
table(temp_seurat$RNA_snn_res.0.3, temp_seurat$full_RNA_snn_res.1)


#what other things changed?
table(temp_seurat$RNA_snn_res.1, temp_seurat$full_RNA_snn_res.1)

#what about the auto annotation?
table(temp_seurat$scGate_multi, temp_seurat$full_scGate_multi)

#what about the cd4/cd8 cluster separation?
table(temp_seurat$RNA_snn_res.1, temp_seurat$scGate_multi)
table(temp_seurat$RNA_snn_res.1, temp_seurat$full_scGate_multi)

#the cd4/cd8 new classification seems to be a bit more coherent with the new clustering
#compared to the previous cd4/cd8 classification compared to the previous clustering
#there is a small group of cells that get classified as cd8, but we don't seem to manage to separate them
#from the cd4 cells
#cluster 10 that were the NK phenotype cells got split, but we still didn't get a separated cluster with the TCR data
#Some of the cells of proliferating cluster 16 now got grouped with the plasma cells, so maybe we cleaned a bit the proliferation cluster
#but we still see some macro/mono markers, maybe they are biological doublets
#we are getting cluster 6 that is probably the low quality cells getting clustered together
#we got a new cluster, cluster 8, with cells coming from old clusters 1, 4, 10 and 14
#even though they have TCR data, they also have high counts, high features and B cell markers
#they couldn't be identified by the auto annotation, and they seem to be a clearly separated cluster
#they may be doublets or other kind of biological or tehcnical artifact
#So we are most likely removing noise from our clusters of interest

#Now we will plot some more genes to help us choose a resolution for the Tcells
#but the resolution that first splits the low quality cells is 0.6. We can't really go lower than this

##T cells variable plots-----
temp_genes <- c('CD2','CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B',
                'MKI67', 'TOP2A',
                'CCR7','LEF1', 'SELL', 'TCF7',
                'TIGIT', 'PDCD1', 'CTLA4', 'LAG3',
                'GZMA', 'GZMH', 'PRF1', 'NKG7',
                'CCL3', 'CCL4', 'CCL5', 'IFNG',
                'FOXP3', 'IL2RA', 'CD40LG', 'IL17A')

temp_genes %in% rownames(temp_seurat)

temp_gene_subtitles <- c('T cells', 'T cells', 'T cells', 'T cells' , 'T cells' , 'T cells',
                         'Proliferation', 'Proliferation',
                         'Naive', 'Naive', 'Naive', 'Naive',
                         'terminal', 'terminal', 'terminal', 'terminal',
                         'cytolitic', 'cytolitic', 'cytolitic', 'cytolitic',
                         'inflamatory', 'inflamatory', 'inflamatory', 'inflamatory',
                         'Treg', 'Treg', 'Tconv', 'Th17')

length(temp_gene_subtitles) == length(temp_genes)

names(temp_gene_subtitles) <- temp_genes

#feature plots
temp_feat_list <- list()
for(temp_for_gene in temp_genes){
  
  temp_feat_list[[temp_for_gene]] <- FeaturePlot(object = temp_seurat,
                                                 features = temp_for_gene,
                                                 min.cutoff = 'q5', #the cutoffs sometimes make it easier to see imporant things
                                                 max.cutoff = 'q95',
                                                 order = T) + ggtitle(label = temp_for_gene, 
                                                                      subtitle = temp_gene_subtitles[temp_for_gene])
  
  
}

temp_feat_plot <- wrap_plots(temp_feat_list, ncol = 4)

png("scTCR_workshop_features_3.png", width = 60, height = length(temp_feat_list) * 2.5, units = 'cm', res = 420)
plot(temp_feat_plot)
dev.off()

#vemos que el cluster 9 tiene un perfil marcadamente terminal comparado con el 2
#por lo que me quedo con la reoslucion 07. La 08 no parece aportar mucho

#07)CLUSTER ANNOTATION----
##cluster c06----
temp_seurat$cluster <- temp_seurat$RNA_snn_res.0.7
temp_seurat$cluster[temp_seurat$cluster == 'c08'] <- 'TCD4_reg'

temp_c06 <- FindMarkers(object = temp_seurat,
                        ident.1 = 'c06',
                        ident.2 = 'c00',
                        group.by = 'RNA_snn_res.0.7',
                        logfc.threshold = 1,
                        min.pct = 0.7)

temp_seurat$cluster[temp_seurat$cluster == 'c06'] <- 'LowQuality'

##cluster c11----
#evidently cycling

temp_seurat$cluster[temp_seurat$cluster == 'c11'] <- 'Cycling'

##cluster c08----
#definitely Tregs


##cluster c09----
#cluster9 seems to be radically different from 00 and 02, so lets check what is making it special
#but we split it form cluster 2, so we may want to know if this split makes any sense

temp_c09 <- FindMarkers(object = temp_seurat,
                        ident.1 = 'c09',
                        ident.2 = 'c02',
                        group.by = 'RNA_snn_res.0.7',
                        logfc.threshold = 0.5,
                        min.pct = 0.5)

temp_c09_all <- FindMarkers(object = temp_seurat,
                        ident.1 = 'c09',
                        ident.2 = c('c02', 'c00', 'c12'),
                        group.by = 'RNA_snn_res.0.7',
                        logfc.threshold = 0.7,
                        min.pct = 0.5)

temp_seurat$cluster[temp_seurat$cluster == 'c09'] <- 'TCD4_TIGIT'

##cluster c02----
#now let's try to figure out what cl00 and 02 are
temp_c00 <- FindMarkers(object = temp_seurat,
                        ident.1 = 'c00',
                        ident.2 = 'c02',
                        group.by = 'RNA_snn_res.0.7',
                        logfc.threshold = 0.5,
                        min.pct = 0.5)

#cluster 2 is evidently more activated
temp_seurat$cluster[temp_seurat$cluster == 'c02'] <- 'TCD4_memory'

#its rare to find real naives in the tumor
temp_seurat$cluster[temp_seurat$cluster == 'c00'] <- 'TCD4_early_activation'

##cluster c12----
#it's hard to think of an useful direct comparison, so maybe against all the other tcd4 clusters
temp_c12 <- FindMarkers(object = temp_seurat,
                        ident.1 = 'c12',
                        ident.2 = c('c00', 'c02', 'c09'),
                        group.by = 'RNA_snn_res.0.7',
                        logfc.threshold = 1,
                        min.pct = 0.5)

#it comes out as the interferon cluster
temp_seurat$cluster[temp_seurat$cluster == 'c12'] <- 'TCD4_interferon'

##cluster c10----

temp_c10 <- FindMarkers(object = temp_seurat,
                        ident.1 = 'c10',
                        ident.2 = c('c01', 'c05'),
                        group.by = 'RNA_snn_res.0.7',
                        logfc.threshold = 1,
                        min.pct = 0.5)

temp_seurat$cluster[temp_seurat$cluster == 'c10'] <- 'TCD8_CD39'

##clusters c01 and c05----

temp_c05 <- FindMarkers(object = temp_seurat,
                        ident.1 = 'c05',
                        ident.2 = 'c01',
                        group.by = 'RNA_snn_res.0.7',
                        logfc.threshold = 0.5,
                        min.pct = 0.5)

#we could call them c17 and GZMK TCD8 clusters

temp_seurat$cluster[temp_seurat$cluster == 'c05'] <- 'TCD8_c17'
temp_seurat$cluster[temp_seurat$cluster == 'c01'] <- 'TCD8_GZMK'

##clusters c03 and 13----
#cluster 03 may still be of interest because it is a mix of NK and Tcells

temp_seurat$cluster[temp_seurat$cluster == 'c03'] <- 'NK/NKT_mix'

#we will rejoin c13 to c03 in NK/NKT_mix
temp_seurat$cluster[temp_seurat$cluster == 'c13'] <- 'NK/NKT_mix'

##cluster c04----
#cluster 04 are plasma cells according to scgate

temp_seurat$cluster[temp_seurat$cluster == 'c04'] <- 'Plasma_cells'

##cluster c07----
#cluster 07 are probably doublets
table(temp_seurat$cluster, temp_seurat$full_scGate_multi)
table(temp_seurat$cluster, temp_seurat$full_RNA_snn_res.1)
table(temp_seurat$cluster, temp_seurat$b_cdr3_nt_cat, useNA = 'ifany')
temp_seurat$cluster[temp_seurat$cluster == 'c07'] <- 'doublets'


#lets check our classification
table(temp_seurat$cluster, useNA = 'ifany')

#lets save progress and clear space
p07_seurat <- temp_seurat
rm(list = ls(pattern = 'temp_'))
rm(p05_seurat)
gc()

temp_seurat <- p07_seurat

#saving RDS object so that we can start the analysis from this point
saveRDS(temp_seurat, 'scTCR_workshop_data/seurat_object_annotated.RDS')

#08) TCR DIVERSITY----
#reny entropy----

f_renyi_entropy <- function(freq_vector, q_value) {
  
  if(q_value == 1) {
    
    fun_freq_vector <- freq_vector[freq_vector != 0 & !is.na(freq_vector)]
    fun_total <- sum(fun_freq_vector, na.rm = TRUE)
    fun_props <- fun_freq_vector / fun_total
    fun_shannon <- -sum(fun_props * log(fun_props))
    
    return(fun_shannon)
    
    
  } else {
    
    fun_freq_vector <- freq_vector[freq_vector != 0 & !is.na(freq_vector)]
    fun_total <- sum(fun_freq_vector, na.rm = TRUE)
    fun_props <- fun_freq_vector / fun_total
    fun_props_exp <- fun_props ^ q_value
    fun_props_exp_sum <- sum(fun_props_exp, na.rm = TRUE)
    fun_prop_exp_sum_log <- log(fun_props_exp_sum)
    fun_prop_exp_sum_weighted <- fun_prop_exp_sum_log / (1 - q_value)
    
    return(fun_prop_exp_sum_weighted)
    
  }
}

#clone count dataframe----
temp_clone_by_cluster <- temp_seurat@meta.data %>% 
  dplyr::filter(!is.na(b_nt_id)) %>% 
  group_by(cluster, b_nt_id) %>% 
  mutate('b_v_cdr3_nt_count_cluster' = n()) %>% 
  ungroup() %>% 
  distinct(b_nt_id, cluster, b_v_cdr3_nt_count_cluster) %>%
  ungroup() %>% 
  pivot_wider(id_cols = b_nt_id,
              names_from = cluster,
              values_from = b_v_cdr3_nt_count_cluster,
              values_fill = NA) %>% 
  as.data.frame()

names(temp_clone_by_cluster)

#diversity calculation----
#now we can calculate the different entropies
#for plotting purposes, I am going to get the exp(entropy) and then
#log transform the variables in ggplot, so that we get y and x axis with the raw numbers but in log scale

temp_diversity <- data.frame(
  'richness' = exp(sapply(temp_clone_by_cluster[,-1], f_renyi_entropy, q_value = 0)),
  'shannon' = sapply(temp_clone_by_cluster[,-1], f_renyi_entropy, q_value = 1),
  'reciprocal_simpson' = exp(sapply(temp_clone_by_cluster[,-1], f_renyi_entropy, q_value = 2)),
  'renyi_h5' = sapply(temp_clone_by_cluster[,-1], f_renyi_entropy, q_value = 5),
  'cluster' = names(temp_clone_by_cluster[,-1])
)

temp_diversity$complement_pielou <- round(1 - temp_diversity$shannon / log(temp_diversity$richness), 4)
temp_diversity$complement_h2_h0 <- round(1 - log(temp_diversity$reciprocal_simpson) / log(temp_diversity$richness), 4)
temp_diversity$gini <- sapply(temp_clone_by_cluster[, -1], function(fun_one_cluster){
  
  fun_one_cluster <- fun_one_cluster[fun_one_cluster != 0]
  fun_one_gini <- Gini(fun_one_cluster, na.rm = TRUE)
  return(fun_one_gini)
  
})

temp_divplots <- list()

temp_divplots$shannon <- ggplot(data = temp_diversity, aes(x = richness, y = shannon, label = cluster)) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_point() +
  geom_label_repel() +
  scale_x_continuous(trans = 'log', labels = function(x){round(x, 0)}) +
  ggtitle('shannon (diversity)')

temp_divplots$simpson <- ggplot(data = temp_diversity, aes(x = richness, y = reciprocal_simpson, label = cluster)) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_point() +
  geom_label_repel() +
  scale_y_continuous(trans = 'log', labels = function(x){round(x, 0)}) +
  scale_x_continuous(trans = 'log', labels = function(x){round(x, 0)}) +
  ggtitle('inverse_simpson (diversity)')

temp_divplots$renyi_h5 <- ggplot(data = temp_diversity, aes(x = richness, y = renyi_h5, label = cluster)) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_point() +
  geom_label_repel() +
  scale_x_continuous(trans = 'log', labels = function(x){round(x, 0)}) +
  ggtitle('renyi_h5 (diversity)')

temp_divplots$rec_pielou <- ggplot(data = temp_diversity, aes(x = richness, y = complement_pielou, label = cluster)) +
  geom_point() +
  geom_label_repel() +
  scale_y_continuous(trans = 'sqrt', labels = function(x){round(x, 2)}) +
  scale_x_continuous(trans = 'log', labels = function(x){round(x, 0)}) +
  ggtitle('pielou complementary (inequality)')

temp_divplots$complement_h2_h0 <- ggplot(data = temp_diversity, aes(x = richness, y = complement_h2_h0, label = cluster)) +
  geom_point() +
  geom_label_repel() +
  scale_y_continuous(trans = 'sqrt', labels = function(x){round(x, 2)}) +
  scale_x_continuous(trans = 'log', labels = function(x){round(x, 0)}) +
  ggtitle('J2 complementary (inequality)')

temp_divplots$gini <- ggplot(data = temp_diversity, aes(x = richness, y = gini, label = cluster)) +
  geom_point() +
  geom_label_repel() +
  scale_x_continuous(trans = 'log', labels = function(x){round(x, 0)}) +
  ggtitle('gini (inequality)')

temp_divplots <- wrap_plots(temp_divplots, ncol = 3)

png("scTCR_workshop_tcr_diversity.png", width = length(temp_divplots) * 3 * 3, height = length(temp_divplots) * 3 * 2, units = 'cm', res = 420)
plot(temp_divplots)
dev.off()

#we can see that the difference between shannon and inverse simpson is mostly a matter of scale
#although some differences can happen in terms of order
cor(temp_diversity$shannon, temp_diversity$reciprocal_simpson, method = 'spearman')
cor(temp_diversity$shannon, temp_diversity$renyi_h5, method = 'spearman')

#and what about the different inequality indeces?
cor(temp_diversity$complement_pielou, temp_diversity$complement_h2_h0, method = 'spearman')
cor(temp_diversity$gini, temp_diversity$complement_pielou, method = 'spearman')
cor(temp_diversity$gini, temp_diversity$complement_h2_h0, method = 'spearman')

#why am I using spearman correlation?

#there are many things to discuss here, like what is the potential influence of the size of the cluster
#as in number of cells in the cluster
#and how that could be affecting the diversity and evenness indexes

#09) TCR OVERLAP----
#we may also be interested in studying the associations between the clonotypes of different clusters
#one thing we can do, is to construct a table that, for each pair of clusters, counts how many clonotypes are shared
#we can start that from the counts by species and cluster

#temp_freq_mat a matrix with species in rows, groups in columns and number of individuals by species and group in cells,
#rownames with species, no other variables
#species with no abundance as 0, not as NA

temp_freq_mat <- temp_clone_by_cluster
rownames(temp_freq_mat) <- temp_freq_mat$b_nt_id
temp_freq_mat <- as.matrix(temp_freq_mat[,-1])

f_species_overlap <- function(fun_freq_mat){
  
  fun_freq_mat[is.na(fun_freq_mat)] <- 0
  
  fun_df_split <- split(fun_freq_mat, f = rownames(fun_freq_mat))
  
  fun_individuals_list <- lapply(fun_df_split, function(one_fun_split){
    
    fun_res <- matrix(one_fun_split, ncol = length(one_fun_split), nrow = length(one_fun_split))
    colnames(fun_res) <- colnames(fun_freq_mat)
    rownames(fun_res) <- colnames(fun_freq_mat)
    return(fun_res)
    
  })
  
  fun_species_list <- lapply(fun_individuals_list, function(one_fun_split){
    
    fun_res <- matrix(as.numeric(one_fun_split != 0), ncol = sqrt(length(one_fun_split)))
    colnames(fun_res) <- colnames(fun_freq_mat)
    rownames(fun_res) <- colnames(fun_freq_mat)
    return(fun_res)
    
  })
  
  fun_species_t <- lapply(fun_species_list, function(one_fun_split){
    
    t(one_fun_split)
    
  })
  
  #get total counts by species
  fun_individuals_list <- Map(function(x, y) {x * y}, fun_individuals_list, fun_species_t)
  fun_species_list <- Map(function(x, y) {x * y}, fun_species_list, fun_species_t)
  
  fun_individuals <- Reduce('+', fun_individuals_list)
  fun_species <- Reduce('+', fun_species_list)

  
  fun_result <- list('species_overlap' = fun_species,
                     'individuals_overlap' = fun_individuals)
  
  return(fun_result)
  
}

temp_species_overlap <- f_species_overlap(temp_freq_mat)

##Morisita-Horn index----
#Now, these are very raw numbers, we may want a measure of overlap that is more summarized and can
#be compared between clusters with very different number of individuals and species
#one of the things we use a lot is the morisita index

temp_mhdata <- table(temp_seurat$b_nt_id, temp_seurat@meta.data$cluster)
temp_dimnames <- dimnames(temp_mhdata)
temp_mhdata <- as.matrix.data.frame(temp_mhdata)
dimnames(temp_mhdata) <- temp_dimnames

temp_mh <- mh(temp_mhdata, PlugIn = TRUE)$PlugIn

#since morisita is a measure of similarity, we could do some clustering to improve visualization

#get order of clusters so that we can plot them nicely
temp_order <- colnames(temp_mh)[hclust(as.dist(temp_mh), method="complete")$order]

#merge the whole data into a long format dataframe

temp_plot_data <- list()

temp_plot_data$mh <- as.data.frame(temp_mh + t(temp_mh)) %>%
  mutate('cluster' = colnames(.)) %>%
  pivot_longer(-cluster, names_to = "cluster_y", values_to = "mh_overlap") %>% 
  mutate('cluster_pair' = paste0(cluster, cluster_y)) %>% 
  distinct(cluster_pair, .keep_all = TRUE) %>% 
  dplyr::filter(cluster != cluster_y) %>% 
  mutate('mh_overlap' = round(mh_overlap, 3))

temp_plot_data$species <- data.frame(temp_species_overlap$species_overlap) %>% 
  mutate('cluster' = colnames(.)) %>%
  pivot_longer(-cluster, names_to = "cluster_y", values_to = "species_overlap") %>% 
  mutate('cluster_pair' = paste0(cluster, cluster_y)) %>% 
  distinct(cluster_pair, .keep_all = TRUE) %>% 
  dplyr::filter(cluster != cluster_y) %>% 
  dplyr::select(cluster_pair, species_overlap)

temp_plot_data$individuals <- data.frame(temp_species_overlap$individuals_overlap) %>% 
  mutate('cluster' = colnames(.)) %>%
  pivot_longer(-cluster, names_to = "cluster_y", values_to = "individuals_overlap") %>% 
  mutate('cluster_pair' = paste0(cluster, cluster_y)) %>% 
  distinct(cluster_pair, .keep_all = TRUE) %>% 
  dplyr::filter(cluster != cluster_y) %>% 
  dplyr::select(cluster_pair, individuals_overlap)

sapply(temp_plot_data, dim)

temp_plot_data <- merge(temp_plot_data$mh, temp_plot_data$species, by = 'cluster_pair') %>% 
  merge(., temp_plot_data$individuals, by = 'cluster_pair')

temp_plot_data <- temp_plot_data %>% 
  pivot_longer(-c(cluster, cluster_pair, cluster_y),
               names_to = 'index_type',
               values_to = 'index_value') %>% 
  group_by(index_type) %>% 
  mutate('index_value_scaled' = scale(index_value)[,1]) %>% 
  ungroup()

temp_plot_data <- temp_plot_data[temp_plot_data$index_value != 0,]

#set the order of the factors
temp_plot_data$cluster <- factor(temp_plot_data$cluster, levels = temp_order)
temp_plot_data$cluster_y <- factor(temp_plot_data$cluster_y, levels = temp_order)

unique(temp_plot_data$index_type)
temp_plot_data$index_type <- factor(temp_plot_data$index_type, levels = c('individuals_overlap', 'species_overlap', 'mh_overlap'))

png('scTCR_workshop_tcr_overlap.png', width = 60, height = 20, units = 'cm', res = 120)
ggplot(data = temp_plot_data, aes(x = cluster, y = cluster_y, size = index_value_scaled, label = index_value)) +
  geom_point() +
  geom_text(vjust = -0.5) +
  facet_wrap(~ index_type) +
  scale_size_continuous(range = c(3, 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        axis.title.y = element_blank(),
        legend.position = 'none')
dev.off()

#explore most abundant species
temp_hm <- temp_freq_mat
temp_hm[is.na(temp_hm)] <- 0 
temp_hm <- temp_hm[rowSums(temp_hm) > 3,] #keep species with at least 3 individuals
temp_hm <- temp_hm[, colSums(temp_hm) > 0] #keep clusters with some information left

dev.off()
png('scTCR_workshop_tcr_heatmap.png', width = 20, height = 20, units = 'cm', res = 420)
print(pheatmap::pheatmap(temp_hm, scale = 'column', display_numbers = temp_hm))
dev.off()

#explore cycling vs cd39
temp_clone_by_cluster %>% 
  dplyr::select(b_nt_id, TCD8_GZMK, TCD4_early_activation) %>% 
  dplyr::filter(complete.cases(.))

temp_seurat$temp_variable <- temp_seurat$b_nt_id
temp_seurat$temp_variable[!(temp_seurat$b_nt_id %in% c('b_nt_0003_16', 'b_nt_0070_2', 'b_nt_0106_2', 'b_nt_0145_2'))] <- NA
temp_seurat$temp_scgate <- temp_seurat$scGate_multi
temp_seurat$temp_scgate[is.na(temp_seurat$temp_variable)] <- NA

png('scTCR_workshop_tcr_umap_clarification.png', width = 30, height = 10, units = 'cm', res = 420)
DimPlot(temp_seurat, 
        group.by = c('temp_variable', 'temp_scgate'),
        order = TRUE,
        pt.size = 1.5)
dev.off()

table(temp_seurat$temp_variable, temp_seurat$temp_scgate)

#10) TCR matching----
temp_vdjdb <- readRDS('scTCR_workshop_data/vdjdb.rds')

# Explore the database
head(temp_vdjdb)
colnames(temp_vdjdb)

table(temp_vdjdb$Pathology)

#filtered vdjdb to common clonotypes with vdjdb score higher than 0
temp_vdj_filter <- temp_vdjdb[temp_vdjdb$vdjdb.score > 0 ,]
temp_vdj_filter <- temp_vdj_filter[temp_vdj_filter$cdr3 %in% temp_seurat$b_cdr3_aa,]

#we can add this to our metadata
temp_vdj_filter <-  temp_vdj_filter %>% 
  dplyr::rename('b_cdr3_aa' = 'cdr3')
temp_metadata <- left_join(temp_seurat@meta.data, temp_vdj_filter, by = 'b_cdr3_aa')
rownames(temp_metadata) <- temp_metadata$cell_id

all(temp_metadata$cell_id == temp_seurat$cell_id)

temp_seurat@meta.data <- temp_metadata
names(temp_seurat@meta.data)

#Check relevant information about these clonotypes
temp_detail <- temp_seurat@meta.data %>% 
  dplyr::select(b_cdr3_aa, b_nt_id, scGate_multi, cluster, antigen.species) %>% 
  dplyr::filter(!is.na(antigen.species))
View(temp_detail)

#calculate levenshtein distances to vdjdb sequences with vdjdb score > 0
#and for the most expanded clonotypes
temp_vdj_filter <- temp_vdjdb[temp_vdjdb$vdjdb.score > 0 ,]
temp_sequences <- na.omit(unique(temp_seurat$b_cdr3_aa[temp_seurat$b_v_cdr3_nt_count > 5]))

table(temp_vdj_filter$antigen.species)

temp_stringdist <- stringdistmatrix(a = temp_sequences,
                                    b = unique(temp_vdj_filter$cdr3),
                                    useNames = TRUE,
                                    method = "lv")

#For all of our expanded clonotypes, we can know what is the smallest distance to a sequence in the vdjdb
#and keep only the information when the minimum distance is 2 or less
#that 2 is arbitrary, what is a relevant distance needs to be thorougly discussed in each case

temp_stringdist_long <- data.frame(temp_stringdist) %>% 
  mutate('b_cdr3_aa' = rownames(.)) %>%
  pivot_longer(-b_cdr3_aa, names_to = "cdr3", values_to = "levenshtein_dist") %>% 
  dplyr::filter(levenshtein_dist <= 2) %>%
  left_join(., temp_vdj_filter, by = 'cdr3')

temp_metadata <- temp_seurat@meta.data %>% 
  dplyr::select(v_gene, full_cdr3_nt, b_cdr3_aa, b_nt_id, cluster) %>% 
  distinct(b_cdr3_aa, .keep_all = TRUE)

temp_detail <- left_join(temp_stringdist_long, temp_metadata, by = 'b_cdr3_aa') %>% 
  dplyr::select(b_nt_id,  b_cdr3_aa, cdr3, antigen.species,levenshtein_dist)

View(temp_detail)

temp_tabledata <- temp_seurat@meta.data[temp_seurat$b_nt_id %in% c('b_nt_0005_11', 'b_nt_0006_11', 'b_nt_0008_10', 'b_nt_0010_7', 'b_nt_0012_6'),]
table(temp_tabledata$b_nt_id, temp_tabledata$cluster)

#what about the expanded clonotypes with larger differences?
temp_tabledata <- temp_seurat@meta.data[!(temp_seurat$b_nt_id %in% c('b_nt_0005_11', 'b_nt_0006_11', 'b_nt_0008_10', 'b_nt_0010_7', 'b_nt_0012_6')),]
temp_tabledata <- temp_tabledata[temp_tabledata$b_v_cdr3_nt_count > 5,]
table(temp_tabledata$b_nt_id, temp_tabledata$cluster)

