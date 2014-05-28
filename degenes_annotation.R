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

retrieve_annot('line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.txt',
               'line6u_vs_i.uniq.annot.degenes.txt')
retrieve_annot('line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.txt',
               'line7u_vs_i.uniq.annot.degenes.txt')

retrieve_annot('miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')

retrieve_annot('miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')

#retrieve_annot('miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt'
#               ,'miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
