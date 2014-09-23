annotate:

	Rscript $(protocol)/degenes_annotation.R
	Rscript $(protocol)/degenes_annotation_human.R

run-kegg:

	# Rscript $(protocol)/goseq_line6.R
	# Rscript $(protocol)/goseq_line7.R
	Rscript $(protocol)/make_kegg_cleveland_plot.R

run-bp-ontology:

	Rscript $(protocol)/goseq_BP_up_down_both_lines.R
