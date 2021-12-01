# breast cancer demo
BRCA_Phenotype <- read.csv('TCGA-BRCA.GDC_phenotype.tsv', header = T, 
                           row.names=1,sep="\t", check.names=F, quote = "")


which(BRCA_Phenotype$sample_type_id.samples == 2)
BRCA_type <- data.frame(ID = rownames(BRCA_Phenotype),
                        category = BRCA_Phenotype$sample_type_id.samples,
                        type = BRCA_Phenotype$sample_type.samples,
                        stringsAsFactors = F)
BRCA_miRNA <- read.csv('TCGA-BRCA.mirna.tsv', header = T,
                       row.names=1,sep="\t", check.names=F, quote = "")

BRCA_mRNA <- read.csv('TCGA-BRCA.htseq_counts.tsv', header = T,
                      row.names=1,sep="\t", check.names=F, quote = "")

rownames(BRCA_mRNA)

process2 <- function(x){
  y <- strsplit(x, split = "\\.")[[1]][1]
  y
}

ensemble_ID <- sapply(rownames(BRCA_mRNA), process2)

rownames(BRCA_mRNA) <- ensemble_ID



library("biomaRt")
listMarts()


## ----ID transformation entrez---------------------------------------------------------------------------------------------------------
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

SYMBOL <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                filters="ensembl_gene_id",values = rownames(BRCA_mRNA),
                mart = ensembl)

entrez <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
                filters="ensembl_gene_id",values = rownames(BRCA_mRNA),
                mart = ensembl)
entrez1 <- entrez[!is.na(entrez$entrezgene),]

##########################3
hsa <- read.csv('hsa.csv', stringsAsFactors = F)

miRNA <- unique(hsa$miRNA)
mRNA <- unique(hsa$Target.Gene..Entrez.Gene.ID.)

relation <- data.frame(miRNA = hsa$miRNA, mRNA = hsa$Target.Gene, gene_id = hsa$Target.Gene..Entrez.Gene.ID., stringsAsFactors = F)
#delete duplicated rows
library(plyr)
relation1 <- ddply(relation,.(miRNA,gene_id),nrow)



##############################differential expression
##############################differential expression
library(limma)
library(affy)

id_type1 <- BRCA_type[BRCA_type$ID %in% colnames(BRCA_miRNA),]

id_type2 <- id_type1[match(colnames(BRCA_miRNA),id_type1$ID),]
#id_type2 <- L_type[colnames(L_miRNA) %in% L_type$ID,]
table(id_type2$category)
strain <- id_type2$type
#strain1 <- c()
for (i in 1:length(strain)){
  if(strain[i] == 'Metastatic'){
    strain[i] <- 'Primary Tumor'
  }
}

design <- model.matrix(~0+factor(strain))
colnames(design) <- c("Tumor","Normal")
design
contrast.matrix<-makeContrasts('Tumor-Normal',levels = design)

fit <- lmFit(BRCA_miRNA, design)

fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
options(digits=2)

#top1<-topTable(fit2,adjust="BH")
top1<-topTable(fit, n=1000*1000, adjust="BH")



de_mirna <- top1[top1$adj.P.Val < 0.01,]


################################mRNA differential expression
id_type3 <- BRCA_type[BRCA_type$ID %in% colnames(BRCA_mRNA),]

id_type4 <- id_type3[match(colnames(BRCA_mRNA),id_type3$ID),]
#id_type2 <- L_type[colnames(L_miRNA) %in% L_type$ID,]
table(id_type4$category)
strain <- id_type4$type
#strain1 <- c()
for (i in 1:length(strain)){
  if(strain[i] == 'Metastatic'){
    strain[i] <- 'Primary Tumor'
  }
}

design <- model.matrix(~0+factor(strain))
colnames(design) <- c("Tumor","Normal")
design


contrast.matrix<-makeContrasts('Tumor-Normal',levels = design)



#top1<-topTable(fit2,adjust="BH")
#top1<-topTable(fit, n=1000*1000, adjust="BH")
fit <- lmFit(BRCA_mRNA, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
options(digits=2)
top1_mRNA<-topTable(fit, n=1000*1000, adjust="BH")

de_mRNA <- top1_mRNA[top1_mRNA$adj.P.Val < 0.01,]

##########################################################
#####match miRNA
miRNA_p <- read.csv('miRNA_P.csv', stringsAsFactors = F)
miRNA_p1 <- miRNA_p[miRNA_p$pre.miRNA %in% rownames(BRCA_miRNA),]


rownames(de_mRNA)
entrez2 <- entrez1[entrez1$ensembl_gene_id %in% rownames(de_mRNA),]
miRNA <- BRCA_miRNA

miRNA_gene_num1 <- c()
de_miRNA_gene_num <- c()
pvalue_miRNA_BRCA <- c()
miRNA_delta_BRCA <- c()
miRNA_name <- c()

#entrez_p1 <- read.csv('entrez_p1.csv', header = T, stringsAsFactors = F)
for (i in 1:length(rownames(top1))){
  cat('i:',i,'\n')
  len1 = length(which(miRNA_p1$pre.miRNA == rownames(top1)[i]))
  if(len1 != 0){
    miRNA_name <- c(miRNA_name, rownames(top1)[i])
    id1 = which(miRNA_p1$pre.miRNA == rownames(top1)[i])
    m_miRNA <- as.character(miRNA_p1[id1,]$mature.miRNA)
    gene <- c()
    for(j in 1:length(m_miRNA)){
      id_gene <- which(relation1$miRNA == m_miRNA[j])
      gene_name <- relation1[id_gene,]$gene_id
      gene <- c(gene, gene_name)
    }
    
    overlap_gene <- intersect(gene, entrez2$entrezgene)
    overlap_gene_num <- length(overlap_gene)
    all_gene_num <- length(entrez1$ensembl_gene_id)
    de_gene_num <- length(entrez2$ensembl_gene_id)
    miRNA_gene_num <- length(gene)
    # 
    miRNA_gene_num1 <- c(miRNA_gene_num1, miRNA_gene_num)
    de_miRNA_gene_num <- c(de_miRNA_gene_num, overlap_gene_num)
    
    
    #enrichment analysis
    pvalue_miRNA <- phyper(overlap_gene_num-1,miRNA_gene_num, all_gene_num-miRNA_gene_num, 
                           de_gene_num, lower.tail=FALSE)
    
    
    pvalue_miRNA_BRCA <- c(pvalue_miRNA_BRCA, pvalue_miRNA)
    
    
    
  }
}

top_miRNA <- data.frame(name = miRNA_name, miRNA_gene_num = miRNA_gene_num1,
                        de_miRNA_gene_num = de_miRNA_gene_num,
                        pvalue_miRNA = pvalue_miRNA_BRCA, stringsAsFactors = F)

top_miRNA_order <- top_miRNA[order(top_miRNA$pvalue_miRNA),]

top2 <- top1
top2$name <- rownames(top2)
library(SPIA)
combine <- "fisher"


#pcombFDR = p.adjust(pG, "fdr")

top_miRNA_all <- merge(top2, top_miRNA_order, by = 'name')
pG <- combfunc(top_miRNA_all$P.Value, top_miRNA_all$pvalue_miRNA, combine)
top_miRNA_all$pG <- combfunc(top_miRNA_all$P.Value, top_miRNA_all$pvalue_miRNA, combine)
top_miRNA_all$pG_fdr <- p.adjust(pG, "fdr")
top_miRNA_all_order <- top_miRNA_all[order(top_miRNA_all$pG_fdr),]

#score <- top_miRNA_all_order$logFC * (-log10(top_miRNA_all_order$pG))
#top_miRNA_all_order$score <- score

write.csv(top_miRNA_all_order,'top_miRNA_all_order_BRCA_1.csv')
