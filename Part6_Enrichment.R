#########################
#
# Part 6: Enrichment Analyses
#
#########################



#####################
# Part 6.1 : GO enrichment
#####################

R

library(data.table)
library(dplyr)
library(topGO) 

# Global genes with GO annotation
genes=readMappings("FromJose2023Aug11_FINAL_gff3_extracted_GO_FINAL_2024June12.txt", sep = "\t", IDsep = ",")

# Candidate gene names
gene=read.table("SELECTED_MERGED.500kbpwindow.genesGO_2024July22.txt",header=F)


#------------ BP: biological process --------------

  # Prep 
  geneNames=names(genes)
  outlier_gene = unique(as.vector(gene$V1))
  geneList=factor(as.integer(geneNames %in% outlier_gene))
  names(geneList)=geneNames

  # Make a database: GOdata
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot=annFUN.gene2GO, gene2GO = genes)

  # Run stats test 
  restRes=runTest(GOdata, algorithm="classic", statistic="fisher")

  # Final output table with all info 
  sigterm=restRes@geneData[['Significant']];sigterm
  
  if (sigterm<1){gene_table=data.frame()}else{

    all_sigterms=restRes@geneData[['SigTerms']]
    gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=length(restRes@score),numChar=999999) #topNodes how many to display (all here)

    pvalues=GenTable(GOdata,Fisher.p=restRes,topNodes=length(restRes@score))$Fisher.p
    fdr=p.adjust(pvalues,method="BH")

    gene_table$adjust.p=round(fdr,6) #round the number 

	}

write.table(gene_table,file="Candidate_Haploblock_Gene_GO_enrich_forBP_all.txt", sep="\t", quote=F, row.names=F, col.names=T)



#------------ MF: Molecular Function --------------
rm(gene_table)

# Make a database: GOdata
  GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot=annFUN.gene2GO, gene2GO = genes)

  # Run stats test 
  restRes=runTest(GOdata,algorithm="classic",statistic="fisher")

  # Final output table with all info 
  sigterm=restRes@geneData[['Significant']];sigterm
  
  if (sigterm<1){gene_table=data.frame()} else{
    gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=length(restRes@score),numChar=999999) #topNodes how many to display (all here)
    pvalues=GenTable(GOdata,Fisher.p=restRes,topNodes=length(restRes@score))$Fisher.p
    fdr=p.adjust(pvalues,method="BH")
    gene_table$adjust.p=round(fdr,6) #round the number 
}


write.table(gene_table,file="Candidate_Haploblock_Gene_GO_enrich_forMF_all.txt", sep="\t", quote=F, row.names=F, col.names=T)


#------------ CC: cellular component --------------

rm(gene_table)
# Make a database: GOdata
  GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, annot=annFUN.gene2GO, gene2GO = genes)

  # Run stats test 
  restRes=runTest(GOdata,algorithm="classic",statistic="fisher")

  # Final output table with all info 
  sigterm=restRes@geneData[['Significant']];sigterm
 
  if (sigterm<1){gene_table=data.frame()} else{
    gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=length(restRes@score),numChar=999999) #topNodes how many to display (all here)
    pvalues=GenTable(GOdata,Fisher.p=restRes,topNodes=length(restRes@score))$Fisher.p
    fdr=p.adjust(pvalues,method="BH")
    gene_table$adjust.p=round(fdr,6) #round the number 
}

write.table(gene_table,file="Candidate_Haploblock_Gene_GO_enrich_forCC_all.txt", sep="\t", quote=F, row.names=F, col.names=T)












#####################
# Part 6.2 : KEGG enrichment
#####################

#Annotate at KAAS:https://www.genome.jp/kegg/kaas/

R

library("clusterProfiler")

# Global genes with KEGG annotation
all_keggs <- read.delim("gene_K_KO_multirow_2024June14.txt", header=T, sep="\t")
head(all_keggs)

# Candidate gene names
gene <- read.table("SELECTED_MERGED.500kbpwindow.genes_forKEGG_2024July22.txt",header=F)
colnames(gene)[1] <- "gene_ID"
dim(gene)
# 51 selected genes
# But some might not match to KEGG file, becuase they did not map to a KO path


#------ Run clusterProfiler:Enricher 
term2gene <- data.frame(KO=all_keggs$KO_Path,Gene=all_keggs$gene_ID)
term2name=data.frame(KO=all_keggs$KO_Path,description=all_keggs$KO_Path_desp)

gene <- gene$gene_ID

table <- enricher(gene,TERM2GENE = term2gene,TERM2NAME =term2name,
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH",
                 maxGSSize = length(gene))


kegg_table <- as.data.frame(table)
dim(kegg_table)
kegg_table



#END