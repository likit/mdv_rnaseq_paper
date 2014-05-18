retrieve_annot<-function(input, output) {
  degenes.table<-read.table(input, stringsAsFactors=F, sep="\t", header=T, row.names=NULL)
  degenes.table <- data.frame(degenes.table$row.names, degenes.table$geneID)
  colnames(degenes.table) <- c("geneID", "ENSEMBL")
  annots<-select(org.Gg.eg.db, keys=as.vector(degenes.table$ENSEMBL),
               columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

  annotated.degenes<-merge(degenes.table, annots,
                           by.x="ENSEMBL", by.y="ENSEMBL",)
  # remove gene with no Entrez ID
  annotated.degenes<-annotated.degenes[
    !is.na(annotated.degenes$ENTREZID),]
  
  # remove duplicated Entrez ID
  uniq.annotated.degenes<-annotated.degenes[
    !duplicated(annotated.degenes$ENTREZID),]
  
  write.table(uniq.annotated.degenes, output, sep='\t', quote=F)
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

retrieve_annot('miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
retrieve_annot('miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt'
               ,'miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt.fa.gallus.uniq.annot.miso.txt')
