library('org.Gg.eg.db')

retrieve_annot<-function(input, output) {
  degenes.table<-read.table(input, stringsAsFactors=F, sep="\t", header=T, row.names=NULL)

  # remove gene with no Ensembl ID
  degenes.table<-degenes.table[!is.na(degenes.table$geneID),]
  
  degenes.table <- degenes.table[order(degenes.table$geneID,
                                       degenes.table$bits, decreasing=T),]

  # remove duplicated Ensembl ID
  degenes.table<-degenes.table[!duplicated(degenes.table$geneID),]
  
  degenes.table <- degenes.table[order(degenes.table$row.names,
                                       degenes.table$bits, decreasing=T),]
  # remove duplicated gene ID
  degenes.table<-degenes.table[!duplicated(degenes.table$row.names),]
  
  degenes.table <- data.frame(degenes.table$row.names, degenes.table$geneID)
  colnames(degenes.table) <- c("geneID", "ENSEMBL")  
  
  annots<-select(org.Gg.eg.db, keys=as.vector(degenes.table$ENSEMBL),
               columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

  annotated.degenes<-merge(degenes.table, annots,
                           by.x="ENSEMBL", by.y="ENSEMBL",)
  write.table(annotated.degenes, output, sep='\t', quote=F)
}

retrieve_annot('results/line6u_vs_i.degenes.fdr.05.fa.tophits.txt',
               'results/line6u_vs_i.uniq.annot.degenes.txt')
retrieve_annot('results/line7u_vs_i.degenes.fdr.05.fa.tophits.txt',
               'results/line7u_vs_i.uniq.annot.degenes.txt')

retrieve_annot('results/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt.se'
               ,'results/line6i_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt.se')

retrieve_annot('results/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt.a3ss'
               ,'results/line6i_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt.a3ss')

retrieve_annot('results/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt.a5ss'
               ,'results/line6i_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt.a5ss')
