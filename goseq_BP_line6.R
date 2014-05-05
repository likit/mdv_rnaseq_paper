library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(ggplot2)
library(biomaRt)
library(GO.db)

degenes.table<-read.table(
  'line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.txt',
  stringsAsFactors=F, sep="\t", header=T)

annots<-select(org.Gg.eg.db, keys=degenes.table$geneID,
               columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

annotated.degenes<-merge(degenes.table, annots,
                         by.x="geneID", by.y="ENSEMBL")

# remove duplicated Entrez ID
uniq.annotated.degenes<-annotated.degenes[
  !duplicated(annotated.degenes$geneID),]

# remove gene with no Entrez ID
uniq.annotated.degenes<-uniq.annotated.degenes[
  !is.na(uniq.annotated.degenes$ENTREZID),]

mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
allgenes<-getBM(attributes='ensembl_gene_id', mart=mart)
allgenes<-allgenes$ensembl_gene_id

#write.table(uniq.annotated.degenes$geneID,
#            'lin6u_vs_i.uniq.annot.degenes.txt', sep='\t', quote=F)

gene.vector<-as.integer(allgenes%in%uniq.annotated.degenes$geneID)
names(gene.vector)<-allgenes

pwf=nullp(gene.vector, 'galGal4', 'ensGene')

cat("BP pathway analysis..\n")
# KEGG Pathway analysis
BP = goseq(pwf, "galGal4", "ensGene", test.cats="GO:BP")

# Adjust P-value using BH method
BP$padjust = p.adjust(BP$over_represented_pvalue, method="BH")

# Get pathway names for significant patways
BP_SIG = BP[BP$padjust<0.05,]
terms = stack(lapply(mget(BP[BP$padjust<0.05,]$category, GOTERM), Term))
BP_SIG$terms = terms$values
xx = as.list(org.Gg.egPATH2EG)
xx = xx[!is.na(xx)] # remove BP IDs that do not match any gene

cat("Writing genes to files..\n")
# Write genes in each pathway to separate files
get_genes_BP = function(cat, data, prefix)
{
  m = match(xx[[cat]], data$ENTREZID)
  mm = m[!is.na(m)]
  d = data.frame(cat, data[mm,]$geneID, data[mm,]$ENTREZID)
  
  filename = paste(prefix, cat, sep="_")
  write.table(d, filename, sep="\t", row.names=F, col.names=F, quote=F)
  return(d)
}
df = lapply(BP_SIG$category, get_genes_BP,
            uniq.annotated.degenes, "line6_goseq_BP_genes")

# Get number of genes in each pathway
#BP_SIG$ngenes = sapply(df, nrow) # get a number of genes from each category

# Writing pathway information to a file
cat("Writing pathways' info to a file...\n")
write.table(BP_SIG, 'line6u_vs_i.degenes.BP.txt',
            sep='\t', row.names=F, quote=F)
