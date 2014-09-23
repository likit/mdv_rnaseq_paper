library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(biomaRt)
library(GO.db)

degenes.table<-read.table('results/line6u_vs_i.fold_change.txt',
                          stringsAsFactors=F, sep="\t", header=F)
colnames(degenes.table)<-c('id', 'ENSEMBL', 'FC')

degenes.table.uniq<-degenes.table[!duplicated(degenes.table$ENSEMBL),]

mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
allgenes<-getBM(attributes=c('ensembl_gene_id', 'start_position',
                             'end_position'), mart=mart)
allgenes$length<-allgenes$end_position - allgenes$start_position

run_goseq<-function(gene.list) {
  
  gene.vector<-as.integer(allgenes$ensembl_gene_id%in%gene.list$ENSEMBL)
  names(gene.vector)<-allgenes$ensembl_gene_id
  gene.vector.length<-allgenes$length
  names(gene.vector.length)<-allgenes$ensembl_gene_id

  pwf=nullp(gene.vector, bias.data=gene.vector.length)

  ontology = goseq(pwf, "galGal4", "ensGene", test.cats="GO:BP")

  # Adjust P-value using BH method
  ontology$padjust = p.adjust(ontology$over_represented_pvalue, method="BH")

  # Get pathway names for significant patways
  ontology_sig = ontology[ontology$padjust<0.05,]
  terms = stack(lapply(mget(ontology[ontology$padjust<0.05,]$category,
                            GOTERM), Term))
  ontology_sig$terms = terms$values
  ontology_sig$termids = terms$ind
  return(ontology_sig)
}

enriched<-run_goseq(degenes.table.uniq[which(degenes.table.uniq$FC>0.0),])
write.table(enriched, 'results/line6u_vs_i.up_genes.BP.txt',
            sep='\t', row.names=F, quote=F)

enriched<-run_goseq(degenes.table.uniq[which(degenes.table.uniq$FC<0.0),])
write.table(enriched, 'results/line6u_vs_i.down_genes.BP.txt',
            sep='\t', row.names=F, quote=F)

degenes.table<-read.table('results/line7u_vs_i.fold_change.txt',
                          stringsAsFactors=F, sep="\t", header=F)
colnames(degenes.table)<-c('id', 'ENSEMBL', 'FC')
degenes.table.uniq<-degenes.table[!duplicated(degenes.table$ENSEMBL),]

enriched<-run_goseq(degenes.table.uniq[which(degenes.table.uniq$FC>0.0),])
write.table(enriched, 'results/line7u_vs_i.up_genes.BP.txt',
            sep='\t', row.names=F, quote=F)

enriched<-run_goseq(degenes.table.uniq[which(degenes.table.uniq$FC<0.0),])
write.table(enriched, 'results/line7u_vs_i.down_genes.BP.txt',
            sep='\t', row.names=F, quote=F)

# xx = as.list(org.Gg.egGO2ALLEGS)
# xx = xx[!is.na(xx)] # remove GO IDs that do not match any gene
# #degenes = temp[temp$V2>0,]$V3
# 
# annots<-select(org.Gg.eg.db, keys=degenes.table.uniq$ENSEMBL,
#                columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
# 
# annotated.degenes.table.uniq<-merge(annots, degenes.table.uniq,
#                                     by.x="ENSEMBL", by.y="ENSEMBL")

# # remove duplicated Entrez ID
# uniq.annotated.degenes<-annotated.degenes.table.uniq[
#   !duplicated(annotated.degenes.table.uniq$ENTREZID),]
# 
# # remove gene with no Entrez ID
# uniq.annotated.degenes<-uniq.annotated.degenes[
#   !is.na(uniq.annotated.degenes$ENTREZID),]

# Write genes in each pathway to separate files
# get_genes_list = function(cat, data, prefix)
# {
#   m = match(xx[[cat]], data$ENTREZID)
#   mm = m[!is.na(m)]
#   d = data.frame(cat, data[mm, ]$id, data[mm,]$ENSEMBL, 
#                  data[mm, ]$FC, data[mm,]$ENTREZID, data[mm,]$SYMBOL)
#   #d = data.frame(data[mm,]$SYMBOL)
#   filename = paste(prefix, cat, sep="_")
#   write.table(d, filename, sep="\t", row.names=F, col.names=F, quote=F)
#   return(d)
# }
# df = lapply('GO:0007155', get_genes_list,
#             uniq.annotated.degenes[which(uniq.annotated.degenes$FC<0.0), ],
#             "line7_down_goseq_BP_genes")
