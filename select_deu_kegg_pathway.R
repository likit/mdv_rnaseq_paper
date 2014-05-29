library(GO.db)
library(KEGG.db)
library(org.Gg.eg.db)

degenes.table<-read.table('select_deu_gene_ensembl.txt', stringsAsFactors=F, sep="\t", header=F)
colnames(degenes.table)<-c("geneID", "ENSEMBL")

degenes.table.path <- select(org.Gg.eg.db, keys=as.vector(degenes.table$ENSEMBL),
                           keytype="ENSEMBL", columns=c("ENSEMBL", "SYMBOL", "ENTREZID", "PATH"))

degenes.table.path<-degenes.table.path[!is.na(degenes.table.path$PATH),]
pathway = stack(mget(degenes.table.path$PATH, KEGGPATHID2NAME))
degenes.table.path$pathway = pathway$values

write.table(degenes.table.path, 'select_deu_kegg_pathway.txt', sep='\t', quote=F, row.names=F, col.names=F)
