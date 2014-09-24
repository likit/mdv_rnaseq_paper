###############################################
#####     Prepare reference genome       ######
###############################################

init:

	grep '>' gal4selected.fa | sort | sed s/\>// > ref.sort.list
	python $(protocol)/sort-fasta.py gal4selected.fa \
		ref.sort.list gal4selected.sorted.fa

	java -jar /opt/software/picardTools/1.113/CreateSequenceDictionary.jar \
		R=gal4selected.sorted.fa O=gal4selected.sorted.dict

	samtools faidx gal4selected.sorted.fa

###############################################
####    RSEM calculate expression      ########
###############################################

rsem-prepare-ref:

	cd gimme; \
		cat gimme.bed.fa | python $(protocol)/fasta-to-gene-list.py > gimme.bed.txt; \
		qsub -v "list=gimme.bed.txt,input=gimme.bed.fa,sample=models-rsem" \
		$(protocol)/rsem_prepare_reference.sh

rsem-calc-expression:

	cd gimme; \
	qsub -v "input_read=../reads/line6u.se.fq,\
		sample_name=line6u-single-rsem,index=models-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd gimme; \
	qsub -v "input_read=../reads/line6i.se.fq,\
		sample_name=line6i-single-rsem,index=models-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd gimme; \
	qsub -v "input_read=../reads/line7u.se.fq,\
		sample_name=line7u-single-rsem,index=models-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd gimme; \
	qsub -v "input_read=../reads/line7i.se.fq,\
		sample_name=line7i-single-rsem,index=models-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd gimme; \
	qsub -v "input_read1=../reads/line6u.pe.1,input_read2=../reads/line6u.pe.2,\
		sample_name=line6u-paired-rsem,index=models-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh

	cd gimme; \
	qsub -v "input_read1=../reads/line6i.pe.1,input_read2=../reads/line6i.pe.2,\
		sample_name=line6i-paired-rsem,index=models-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh

	cd gimme; \
	qsub -v "input_read1=../reads/line7u.pe.1,input_read2=../reads/line7u.pe.2,\
		sample_name=line7u-paired-rsem,index=models-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh

	cd gimme; \
	qsub -v "input_read1=../reads/line7i.pe.1,input_read2=../reads/line7i.pe.2,\
		sample_name=line7i-paired-rsem,index=models-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh

ebseq-line6:

	cd gimme; \
	rsem-generate-data-matrix line6u-single-rsem.genes.results \
		line6u-paired-rsem.genes.results line6i-single-rsem.genes.results \
		line6i-paired-rsem.genes.results > line6u_vs_i.gene.counts.matrix; \
	rsem-run-ebseq line6u_vs_i.gene.counts.matrix 2,2 line6u_vs_i.degenes; \
	rsem-control-fdr line6u_vs_i.degenes 0.05 line6u_vs_i.degenes.fdr.05

ebseq-line7:

	cd gimme; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix; \
	rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes; \
	rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

rsem-output-to-fasta:

	cd gimme; \
	python $(protocol)/rsem-output-to-fasta.py \
		line6u_vs_i.degenes.fdr.05 gimme.bed.fa > \
		line6u_vs_i.degenes.fdr.05.fa

	cd gimme; \
	python $(protocol)/rsem-output-to-fasta.py \
		line7u_vs_i.degenes.fdr.05 gimme.bed.fa > \
		line7u_vs_i.degenes.fdr.05.fa

get-longest-sequences:

	cd gimme; \
	python $(protocol)/gene-rep.py line6u_vs_i.degenes.fdr.05.fa > \
		line6u_vs_i.degenes.fdr.05.fa.longest

	cd gimme; \
	python $(protocol)/gene-rep.py line7u_vs_i.degenes.fdr.05.fa > \
		line7u_vs_i.degenes.fdr.05.fa.longest

create-gallus-blastdb:

	wget ftp://ftp.ensembl.org/pub/release-73/fasta/gallus_gallus/pep/Gallus_gallus.Galgal4.73.pep.all.fa.gz
	gunzip Gallus_gallus.Galgal4.73.pep.all.fa.gz
	formatdb -i Gallus_gallus.Galgal4.73.pep.all.fa -o T -p T -n Gallus_prot

run-blast-gallus:

	cd gimme; \
	qsub -v "db=$(projdir)/Gallus_prot,\
		input=line6u_vs_i.degenes.fdr.05.fa.longest,program=blastx,\
		output=line6u_vs_i.degenes.fdr.05.fa.longest.gallus.xml" \
		$(protocol)/blast.sh

	cd gimme; \
	qsub -v "db=$(projdir)/Gallus_prot,\
		input=line7u_vs_i.degenes.fdr.05.fa.longest,program=blastx,\
		output=line7u_vs_i.degenes.fdr.05.fa.longest.gallus.xml" \
		$(protocol)/blast.sh

get-tophits:

	cd gimme; \
	python $(protocol)/get_top_hits.py line6u_vs_i.degenes.fdr.05.fa.longest.gallus.xml \
		> line6u_vs_i.degenes.fdr.05.fa.tophits.txt

	cd gimme; \
	python $(protocol)/get_top_hits.py line7u_vs_i.degenes.fdr.05.fa.longest.gallus.xml \
		> line7u_vs_i.degenes.fdr.05.fa.tophits.txt
