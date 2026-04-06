# Functions and code that will run to format the input files to be ready for the 
# make_circos.CMVD.rmd

rm(list = ls())

################################################################################
#                                 functions                                    #
################################################################################



#' Formats data for gene input
#' @param data data.frame containing biofilter output

#data = genes_gr_c
format_genes <- function(data){
  
  data <- as.data.frame(data)
  data = data[order(data$P),]
  data2 <- data[,c(2,3,3,16,24)] # label = search gene (not gene)
  data2$POS.1 <- data2$POS.1 + 1
  names(data2)[1:5] <- c("chr", "start", "end", "label","group") 
  data2$chr <- paste0("chr", data2$chr)
  data2$start <- as.numeric(data2$start)
  data3 = unique(data2)
  data3 = data3[!duplicated(data3$label), ]
  
  data.bed <- data3[data3$chr != "chrNA", ]
  
  return(data.bed)
}



#' Formats data for heterozygosity track
#' @param data data.frame containing meta analysis output
format_het <- function(data){
  
  data <- as.data.frame(data)
  
  data$end <- data$BP + 1
  
  data.het <- data[, c("CHR", "BP", "end", "I")]
  
  names(data.het) <- c("chr", "start", "end", "value1")
  
  data.het$chr <- paste0("chr", data.het$chr)
  
  return(data.het)
}

#' Formats data for Manhattan plot
#' @param data data.frame containing meta analysis output
#' @param annotation data.frame containing annotations
#' @param p.threshold p-value to filter out

#data = mega.data2 # for testing
#annotations = anno.data2 # for testing
format_manhattan <- function(data, annotations, 
                             p.threshold = 0.001,
                             col = "black",
                             anno.col = "#009E73"){
  

  data <- as.data.frame(data)
  annotations <- as.data.frame(annotations)
  
  data$end <- data$BP + 1
  data$value1 <- -log10(data$P)
  data$value2 <- ifelse(data$SNP %in% annotations$ID, 
                        anno.col, col)

  data.fmt <- data[data$value1 > -log10(p.threshold) | data$value2 == 1, 
                   c("CHR", "BP", "end", "value1", "value2")]

  
  names(data.fmt)[1:2] <- c("chr", "start")
  
  data.fmt$chr <- paste0("chr", data.fmt$chr)

  return(data.fmt)  
}

################################################################################
#                                 file paths                                   #
################################################################################

# inputs 
pheno.path <-file.path("/path/CMVD_ALL.global_reserve_merged.meta")
path_AFR <- file.path("/path/global_reserve_maf-0.01.regenie.gz")
path_EUR <- file.path("/path/global_reserve_maf-0.01.regenie.gz")
# Correct meta for heterozygosity
meta.path_AFR <- file.path("/path/CMVD_AFR_MF.global_reserve_merged.meta")
meta.path_EUR <- file.path("/path/CMVD_EUR_MF.global_reserve_merged.meta")
anno.path <- file.path("/path/sumstats_all_replication_genes_CAD.csv")
twas.path <- file.path("/path/collated_TWAS_results_withancestrypops_nopcutoff_GTEXpos.csv")

# outputs
genes.out <- file.path("/path/Circos/processed", "gene_track_data_gr")
het.out_AFR <- file.path("/path/Circos/processed", "heterozygosity_data_AFR_5e-4_gr")
het.out_EUR <- file.path("/path/Circos/processed", "heterozygosity_data_EUR_5e-4_gr")
manhattan.out <- file.path("/path/Circos/processed", "manhattan_data_gr")
twas.out <- file.path("/path/Circos/processed", "twas_data_gr")

################################################################################
#                               loading zone                                   #
################################################################################

#mega.data <- data.table::fread(mega.path)
mega.data2 <- data.table::fread(pheno.path)
colnames(mega.data2) <- gsub("#", "", colnames(mega.data2))

# Load the data ancestry stratefied data for the GWAS annotation
data_AFR <- data.table::fread(path_AFR)
colnames(data_AFR) <- gsub("#", "", colnames(data_AFR))
data_EUR <- data.table::fread(path_EUR)
colnames(data_EUR) <- gsub("#", "", colnames(data_EUR))

# Load the data meta-analyzed for SEX
meta.data_AFR <- data.table::fread(meta.path_AFR)
meta.data_EUR <- data.table::fread(meta.path_EUR)

#gene.data <- data.table::fread(gene.path)

anno.data2 <- data.table::fread(anno.path)
colnames(anno.data2) <- gsub("#", "", colnames(anno.data2))
#anno.data <- readxl::read_excel(anno.path)

twas.data <- data.table::fread(twas.path)
#rsid <- data.table::fread(rsid.path)



################################################################################
#                            format for gene track                             #
################################################################################

# Extract the global reserve data
genes_gr = anno.data2[grep("global_reserve", anno.data2$PHENO),]


# Check how many AFR GWAS and EUR GWAS are in the meta-analysis


genes_gr_AFR = genes_gr[genes_gr$COHORT == "CMVD_AFR", ]
genes_gr_AFR$group <- "AFR"
genes_gr_EUR = genes_gr[genes_gr$COHORT == "CMVD_EUR", ]
genes_gr_EUR$group <- "EUR"

# Threshold genes to be labeled that are above study-wide significance 
genes_gr_c = rbind(genes_gr_AFR, genes_gr_EUR)
genes_gr_c$group <- ifelse((genes_gr_c$COHORT == "CMVD_AFR" & genes_gr_c$P <= 0.0002), "AFR", 
                           ifelse((genes_gr_c$COHORT == "CMVD_EUR" & genes_gr_c$P <= 0.0002), "EUR", "NS"))


genes_new = 
mult = read.table("multiple_loci.txt")
afr = read.table("afr_loci.txt")
eur = read.table("eur_loci.txt")
# Subset the data columns to have only
# chr, start, end, label, group
gene.bed <- format_genes(genes_gr_c) # 370

# Add additional important genes
ADK = c("chr10", "74731160", "74731161", "ADK", "Important")
CAPN2 = c("chr1", "223825378", "223825379", "CAPN2", "Important")
ERC1 = c("chr12", "1447550", "1447551", "ERC1", "Important")
AADAC = c("chr3", "151823495", "151823496", "AADAC", "Important")
PLA2G5 = c("chr1", "20041185", "20041186", "PLA2G5", "Important")

# Create final formatted dataset
gene.bed = rbind(gene.bed, ADK, CAPN2, ERC1, AADAC, PLA2G5)
gene.bed$start = as.numeric(gene.bed$start)
gene.bed$end = as.numeric(gene.bed$end)


save(gene.bed, file = genes.out)

################################################################################
#                      format for manhattan plot track                         #
################################################################################

mega.bed_gr <- format_manhattan(data = mega.data2,
                             annotations = anno.data2,
                             p.threshold = 1) # changed from 0.01 to 1 to include all SNPs
save(mega.bed_gr, file = manhattan.out)

################################################################################
#                      format for heterozygosity track                         #
################################################################################

meta.data_AFR = meta.data_AFR[meta.data_AFR$P<=5e-4,] # 2771
meta.data_EUR = meta.data_EUR[meta.data_EUR$P<=5e-4,] # 901

het.data_AFR <- format_het(meta.data_AFR)
het.data_EUR <- format_het(meta.data_EUR)

save(het.data_AFR, file = het.out_AFR)
save(het.data_EUR, file = het.out_EUR)

# Uniquify
length(unique(meta.data_AFR$SNP, meta.data_EUR$SNP))

# Annotation and SNP file
gene_snp = genes_gr[,c(1,23)]

gene_snp_AFR = inner_join(gene_snp, meta.data_AFR, by = c("ID"="SNP"))
gene_snp_EUR = inner_join(gene_snp, meta.data_EUR, by = c("ID"="SNP"))

length(unique(gene_snp_AFR$Gene, gene_snp_EUR$Gene))


################################################################################
#                              format for twas                              #
################################################################################
# Extract the global reserve data
twas_gr = twas.data[grep("global_reserve", twas.data$phenotype),]
# p-value threshold
#twas_gr = twas_gr[twas_gr$pvalue<=0.05,]


# Remove TUBGCP5 because it's crazy
twas_gr3 = twas_gr[!grepl("TUBGCP5", gene_name),]

twas_gr4 = twas_gr3[,c("gene_name", "tissue_weights", "log10_pvalue", "gene_chr", "gene_start")]
# Flip the orientation for the plot
twas_gr4$end_position = twas_gr4$gene_start + 1
# Make pretty colors
twas_gr4$tissue = ifelse(twas_gr4$tissue_weights=="MAGNET_weights_genetype_filtered_allelefix_rsid", "#D81B60", 
                         ifelse(twas_gr4$tissue_weights=="mashr_Heart_Left_Ventricle", "#9C9B65",
                                ifelse(twas_gr4$tissue_weights=="mashr_Heart_Atrial_Appendage", "#FFC107", "")))


twas_gr5 = twas_gr4[,c(4,5,6,3,7)]
colnames(twas_gr5) <- c("chr", "start", "end", "value1", "value2")
twas_gr5$value1 <- twas_gr5$value1 * (-1)


save(twas_gr5, file = twas.out)
