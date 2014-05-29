library(GO.db)

degenes.table.go <- select(org.Gg.eg.db, keys=as.vector(degenes.table$ENSEMBL),
                           keytype="ENSEMBL", columns=c("ENSEMBL", "SYMBOL", "GO"))

degenes.table.go<-degenes.table.go[!is.na(degenes.table.go$ONTOLOGY),]
degenes.table.go.bp<-degenes.table.go[which(degenes.table.go$ONTOLOGY=="BP"),]
degenes.table.go.mf<-degenes.table.go[which(degenes.table.go$ONTOLOGY=="MF"),]

degenes.table.go.bp.terms<-select(GO.db, keytype="GOID", keys=as.vector(degenes.table.go.bp$GO),
                                  columns=c("DEFINITION", "GOID", "TERM"))
degenes.table.go.mf.terms<-select(GO.db, keytype="GOID", keys=as.vector(degenes.table.go.mf$GO),
                                  columns=c("DEFINITION", "GOID", "TERM"))
write.table(degenes.table.go, 'select_deu_ontology.txt', sep='\t', quote=F, row.names=F, col.names=F)
write.table(degenes.table.go.bp.terms, 'select_deu_BP_terms.txt', sep='\t', quote=F, row.names=F, col.names=F)
write.table(degenes.table.go.mf.terms, 'select_deu_MF_terms.txt', sep='\t', quote=F, row.names=F, col.names=F)

head(degenes.table.go.mf.terms)
