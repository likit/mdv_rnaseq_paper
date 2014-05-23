library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(ggplot2)
library(biomaRt)

degenes.table<-read.table('combined_deu_gene_ensembl.txt',
                          stringsAsFactors=F, sep="\t", header=F)
colnames(degenes.table)<-c("geneID", "ENSEMBL")
mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
allgenes<-getBM(attributes=c('ensembl_gene_id', 'start_position',
                             'end_position'), mart=mart)
allgenes$length<-allgenes$end_position - allgenes$start_position

gene.vector<-as.integer(allgenes$ensembl_gene_id%in%degenes.table$ENSEMBL)
names(gene.vector)<-allgenes$ensembl_gene_id
gene.vector.length<-allgenes$length
names(gene.vector.length)<-allgenes$ensembl_gene_id

pwf=nullp(gene.vector, bias.data=gene.vector.length)

# KEGG Pathway analysis
KEGG = goseq(pwf, "galGal4", "ensGene", test.cats="KEGG")

# Adjust P-value using BH method
KEGG$padjust = p.adjust(KEGG$over_represented_pvalue, method="BH")

# Get pathway names for significant patways
KEGG_SIG = KEGG[KEGG$padjust<0.1,]
pathway = stack(mget(KEGG[KEGG$padjust<0.1,]$category, KEGGPATHID2NAME))
KEGG_SIG$pathway = pathway$values
xx = as.list(org.Gg.egPATH2EG)
xx = xx[!is.na(xx)] # remove KEGG IDs that do not match any gene

annots<-select(org.Gg.eg.db, keys=degenes.table$ENSEMBL,
               columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

annotated.degenes<-merge(degenes.table, annots,
                                    by.x="ENSEMBL", by.y="ENSEMBL")

# remove duplicated Entrez ID
#uniq.annotated.degenes<-annotated.degenes.table.uniq[
#  !duplicated(annotated.degenes.table.uniq$geneID),]

# remove gene with no Entrez ID
#uniq.annotated.degenes<-uniq.annotated.degenes[
#  !is.na(uniq.annotated.degenes$ENTREZID),]

# Write genes in each pathway to separate files

get_genes_kegg = function(cat, data, prefix)
{
  m = match(xx[[cat]], data$ENTREZID)
  mm = m[!is.na(m)]
  d = data.frame(cat, data[mm, ]$geneID, data[mm, ]$ENSEMBL,
                 data[mm, ]$ENTREZID, data[mm, ]$SYMBOL)

  filename = paste(prefix, cat, sep="_")
  write.table(d, filename, sep="\t", row.names=F, col.names=F, quote=F)
  return(d)
}
df = lapply(KEGG_SIG$category, get_genes_kegg, annotated.degenes, "deu_goseq_KEGG_genes")

#mapply(write.table, x=df, file="sample.goseq.txt",
#        MoreArgs=list(row.names=FALSE, sep="\t", quote=FALSE,
#        col.names=FALSE, append=TRUE))

# Writing pathway information to a file
write.table(KEGG_SIG, 'deu.combined.KEGG.txt', sep='\t', row.names=F, quote=F)
