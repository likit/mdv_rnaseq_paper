library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(biomaRt)

degenes.table<-read.table('results/line7u_vs_i.degenes.fdr.05.fa.tophits.txt',
                          stringsAsFactors=F, sep="\t", header=T)

# remove gene with no Ensembl ID
degenes.table<-degenes.table[!is.na(degenes.table$geneID),]

degenes.table<-degenes.table[order(degenes.table$geneID,
                          degenes.table$bits, decreasing=T),]

degenes.table.uniq<-degenes.table[!duplicated(degenes.table$geneID),]

mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
allgenes<-getBM(attributes=c('ensembl_gene_id', 'start_position',
                             'end_position'), mart=mart)
allgenes$length<-allgenes$end_position - allgenes$start_position

gene.vector<-as.integer(allgenes$ensembl_gene_id%in%degenes.table.uniq$geneID)
names(gene.vector)<-allgenes$ensembl_gene_id
gene.vector.length<-allgenes$length
names(gene.vector.length)<-allgenes$ensembl_gene_id

pwf=nullp(gene.vector, bias.data=gene.vector.length)

# KEGG Pathway analysis
KEGG = goseq(pwf, "galGal4", "ensGene", test.cats="KEGG")

# Adjust P-value using BH method
KEGG$padjust = p.adjust(KEGG$over_represented_pvalue, method="BH")

# Get pathway names for significant patways
KEGG_SIG = KEGG[KEGG$padjust<0.5,]
pathway = stack(mget(KEGG[KEGG$padjust<0.5,]$category, KEGGPATHID2NAME))
KEGG_SIG$pathway = pathway$values
xx = as.list(org.Gg.egPATH2EG)
xx = xx[!is.na(xx)] # remove KEGG IDs that do not match any gene

degenes.table.uniq$id<-rownames(degenes.table.uniq)
annots<-select(org.Gg.eg.db, keys=degenes.table.uniq$geneID,
               columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

annotated.degenes.table.uniq<-merge(degenes.table.uniq, annots,
                         by.x="geneID", by.y="ENSEMBL")

# remove duplicated Entrez ID
uniq.annotated.degenes<-annotated.degenes.table.uniq[
                      !duplicated(annotated.degenes.table.uniq$geneID),]

# remove gene with no Entrez ID
uniq.annotated.degenes<-uniq.annotated.degenes[
                      !is.na(uniq.annotated.degenes$ENTREZID),]

# Write genes in each pathway to separate files
get_genes_kegg = function(cat, data, prefix)
{
    m = match(xx[[cat]], data$ENTREZID)
    mm = m[!is.na(m)]
    d = data.frame(cat, data[mm, ]$id,
                   data[mm,]$geneID, data[mm,]$ENTREZID, data[mm,]$SYMBOL)

    filename = paste(prefix, cat, sep="_")
    write.table(d, filename, sep="\t", row.names=F, col.names=F, quote=F)
    return(d)
}
df = lapply(KEGG_SIG$category,
            get_genes_kegg, uniq.annotated.degenes,
            "results/line7_goseq_KEGG_genes")

# Get number of genes in each pathway
#KEGG_SIG$ngenes = sapply(df, nrow) # get a number of genes from each category

# Writing pathway information to a file
write.table(KEGG_SIG, 'results/line7u_vs_i.degenes.KEGG.txt',
            sep='\t', row.names=F, quote=F)
