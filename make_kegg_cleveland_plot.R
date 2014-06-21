library(ggplot2)
line7 <- read.table('line7u_vs_i.degenes.KEGG.txt',
                    sep="\t", header=T)
line6 <- read.table('line6u_vs_i.degenes.KEGG.txt',
                    sep="\t", header=T)
line7$sample = "Line 7"
line6$sample = "Line 6"

line6 <- line6[line6$padjust<0.1,]
line7 <- line7[line7$padjust<0.1,]

line67 <- rbind(line6, line7)
line67$log10padjust = (-1)*log10(line67$padjust)
ggplot(line67, aes(x=log10padjust, y=reorder(pathway, log10padjust))) +
  geom_point(aes(size=numDEInCat, colour=sample)) +
  facet_grid(sample ~ ., scales="free_y", space="free_y") +
  geom_segment(aes(yend=pathway), xend=0,
                 colour="grey60",
                 linetype="dashed") +
  theme_bw() +
  scale_colour_brewer(palette="Set1", limits=c("Line 6", "Line 7")) +
  theme(panel.grid.major.y = element_blank()) +
  scale_size_area(max_size=12, name="Number of genes") +
  labs(list(title="Enriched KEGG Pathways",
              x=expression(paste(-log[10](adj.pvalue))),
            y="Pathway")) +
  guides(colour=FALSE)

ggsave("line67_KEGG_cleveland_2.pdf")
