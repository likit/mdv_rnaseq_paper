library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(biomaRt)

degenes.table<-read.table('line6u_vs_i.fold_change.txt',
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
  terms = stack(lapply(mget(ontology[ontology$padjust<0.05,]$category, GOTERM), Term))
  ontology_sig$terms = terms$values
  return(ontology_sig)
}

enriched<-run_goseq(degenes.table.uniq[which(degenes.table.uniq$FC>0.0),])
write.table(enriched, 'line6u_vs_i.up_genes.BP.txt', sep='\t', row.names=F, quote=F)

enriched<-run_goseq(degenes.table.uniq[which(degenes.table.uniq$FC<0.0),])
write.table(enriched, 'line6u_vs_i.down_genes.BP.txt', sep='\t', row.names=F, quote=F)


degenes.table<-read.table('line7u_vs_i.fold_change.txt',
                          stringsAsFactors=F, sep="\t", header=F)
colnames(degenes.table)<-c('id', 'ENSEMBL', 'FC')
degenes.table.uniq<-degenes.table[!duplicated(degenes.table$ENSEMBL),]

enriched<-run_goseq(degenes.table.uniq[which(degenes.table.uniq$FC>0.0),])
write.table(enriched, 'line7u_vs_i.up_genes.BP.txt', sep='\t', row.names=F, quote=F)

enriched<-run_goseq(degenes.table.uniq[which(degenes.table.uniq$FC<0.0),])
write.table(enriched, 'line7u_vs_i.down_genes.BP.txt', sep='\t', row.names=F, quote=F)