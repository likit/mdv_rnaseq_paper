library('pathview')
library(org.Gg.eg.db)

line6.expr<-read.table('line6u_vs_i.cuffref.degenes.fdr.05')
line6.tophits<-read.table(
  'line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.txt',
  stringsAsFactors=F, sep="\t", header=T)
line6.expr<-merge(line6.tophits, line6.expr, by=0)
line6.expr<-line6.expr[order(line6.expr$PostFC),]
line6.expr.uniq<-line6.expr[!duplicated(line6.expr$geneID),]
line6.data <- data.frame(line6.expr.uniq$PostFC,
                         row.names=line6.expr.uniq$geneID)
#line 7
line7.expr<-read.table('line7u_vs_i.cuffref.degenes.fdr.05')
line7.tophits<-read.table(
  'line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.txt',
  stringsAsFactors=F, sep="\t", header=T)
line7.expr<-merge(line7.tophits, line7.expr, by=0)
line7.expr<-line7.expr[order(line7.expr$PostFC),]
line7.expr.uniq<-line7.expr[!duplicated(line7.expr$geneID),]
line7.data <- data.frame(line7.expr.uniq$PostFC,
                         row.names=line7.expr.uniq$geneID)

plot_pathview <- function(pathway.id) {
line6.pathway.file<-paste("line6_goseq_KEGG_genes", pathway.id, sep="_")
line7.pathway.file<-paste("line7_goseq_KEGG_genes", pathway.id, sep="_")
line6.kegg<-read.table(line6.pathway.file, sep='\t', header=F)
colnames(line6.kegg)<-c("pathway", "geneID", "ENTREZID")
line6.kegg<-merge(line6.expr.uniq, line6.kegg,
                  by.x="geneID",by.y="geneID", all.y=T)

line7.kegg<-read.table(line7.pathway.file, sep='\t', header=F)
colnames(line7.kegg)<-c("pathway", "geneID", "ENTREZID")
line7.kegg<-merge(line7.expr.uniq, line7.kegg,
                  by.x="geneID",by.y="geneID", all.y=T)

line6.kegg.data<-data.frame(line6.kegg$PostFC,
                            row.names=line6.kegg$geneID)
line7.kegg.data<-data.frame(line7.kegg$PostFC,
                            row.names=line7.kegg$geneID)
line67.kegg.data<-merge(line6.kegg.data, line7.kegg.data, by=0, all=T)
line67.kegg.data[is.na(line67.kegg.data)]<-0
line67.kegg.data$line6.kegg.PostFC<-(-1)*log2(line67.kegg.data$line6.kegg.PostFC)
line67.kegg.data$line7.kegg.PostFC<-(-1)*log2(line67.kegg.data$line7.kegg.PostFC)
line67.kegg.data[line67.kegg.data==Inf]<- 0.0
line67.kegg.data<-data.frame(line67.kegg.data$line6.kegg.PostFC,
                             line67.kegg.data$line7.kegg.PostFC,
                             row.names=line67.kegg.data$Row.names)
pv.out<-pathview(gene.data=line67.kegg.data[,1:2], pathway.id=pathway.id,
        species="gga", out.suffix="line67", multi.state=TRUE,
        gene.idtype="ENSEMBL",
        limit=list(gene=2,cpd=1))
}

pathway.id = "04630"
plot_pathview(pathway.id)
