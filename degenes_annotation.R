retreive_annot<-function(input, output) {
  degenes.table<-read.table(input, stringsAsFactors=F, sep="\t", header=T)
  degenes.table <- data.frame(rownames(degenes.table), degenes.table$geneID)
  colnames(degenes.table) <- c("geneID", "ENSEMBL")
  annots<-select(org.Gg.eg.db, keys=as.vector(degenes.table$ENSEMBL),
               columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

  annotated.degenes<-merge(degenes.table, annots,
                           by.x="ENSEMBL", by.y="ENSEMBL",)
  head(annotated.degenes)
  # remove duplicated Entrez ID
  uniq.annotated.degenes<-annotated.degenes[
    !duplicated(annotated.degenes$geneID),]

  # remove gene with no Entrez ID
  uniq.annotated.degenes<-uniq.annotated.degenes[
    !is.na(uniq.annotated.degenes$ENTREZID),]

  write.table(uniq.annotated.degenes, output, sep='\t', quote=F)
}

retreive_annot('line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.txt',
               'line6u_vs_i.uniq.annot.degenes.txt')
retreive_annot('line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.txt',
               'line7u_vs_i.uniq.annot.degenes.txt')