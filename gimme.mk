# Requires Gimme
# Please set gimmetools to gimme/src/utils
#
combine-transcripts:

	cd local_assembly; \
		for d in chr*_[0-9][0-9]; \
			do python $(gimmetools)/rename_fasta.py $$d/transcripts.fa local_$$d >> local_merged.fa; \
	done

	cd assembly; \
		for d in global_[0-9][0-9]; \
			do python $(gimmetools)/rename_fasta.py $$d/transcripts.fa global_$$d >> global_merged.fa; \
	done

clean-transcripts:

	cd local_assembly; seqclean local_merged.fa -c 8
	cd assembly; seqclean global_merged.fa -c 8

remove-redundant-seq:

	if [ ! -d gimme ]; then mkdir gimme; fi
	cat local_assembly/local_merged.fa.clean \
		assembly/global_merged.fa.clean >> gimme/all.fa.clean
	cd gimme; \
		qsub -v input="all.fa.clean",output="all.fa.clean.nr",c="1.0" \
		$(protocol)/cdhit_job.sh

	cd assembly; qsub -v "input=global_merged.fa.clean,\
		output=global_merged.fa.clean.nr,c=1.0" \
		$(protocol)/cdhit_job.sh

align-transcripts:

	cd gimme; \
		python $(protocol)/split-fa.py all.fa.clean.nr 10000 subsets
	cd gimme; \
		for f in subsets*.fa; do \
			qsub -v input="$$f" $(protocol)/blat_job.sh; \
		done
		# wait till all the jobs are done before proceeding

sort-alignments:

	cd gimme; \
		cat subsets*.fa.psl > all.fa.clean.nr.psl; \
		sort -k 10 all.fa.clean.nr.psl > all.fa.clean.nr.psl.sorted; \
		pslReps -nohead -singleHit all.fa.clean.nr.psl.sorted all.fa.clean.nr.psl.best info; \
		rm subsets*.fa.psl; \
		rm subsets*.fa

build-merged-gene-models:

	cd gimme; \
	python $(gimmedir)/src/utils/gff2bed.py ../tophat/merged_cuff_ref/merged.gtf > \
		../tophat/merged_cuff_ref/merged.bed; \
	qsub -v "output=gimme.bed,input1=all.fa.clean.nr.psl.best,\
		input2=../tophat/merged_cuff_ref/merged.bed,\
		gimme_dir=$(gimmedir)/src/,ref=../gal4selected.fa" $(protocol)/run_gimme.sh

merged-models-to-transcripts:

	cd gimme; python $(gimmedir)/src/utils/get_transcript_seq.py \
		gimme.bed ../gal4selected.fa > gimme.bed.fa

annotate-merged-genes:

	cd gimme; python $(protocol)/gene-rep-merged.py gimme.bed.fa > gimme.longest.fa

	cd gimme; \
		qsub -v db="Gallus_prot",input="gimme.longest.fa",program="blastx",output="gimme-gga.xml" \
			$(protocol)/blast.sh

	cd gimme; \
		qsub -v db="Human_prot",input="gimme.longest.fa",program="blastx",output="gimme-hsa.xml" \
			$(protocol)/blast.sh

rsem-prepare-reference-merged-models:

	cd gimme; \
		python $(protocol)/get_top_hits.py gimme-gga.xml > gimme-gga-tophits.txt; \
		python $(protocol)/get_best_ensembl_hits_combined.py gimme-gga-tophits.txt \
		gimme.bed.fa > gimme-gga.fa

	cd gimme; \
		cat gimme-gga.fa | python $(protocol)/fasta-to-gene-list.py > gimme-gga-list.txt

	cd gimme; \
	qsub -v "knownIsoforms=gimme-gga-list.txt,input=gimme-gga.fa,\
		output=gimme-gga-rsem" $(protocol)/rsem_prepare_reference.sh

rsem-calc-expression-merged-models:

	cd gimme; \
	qsub -v "input_read=../reads/line7u.se.fq,\
		sample_name=line7u-single-rsem,index=gimme-gga-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd gimme; \
	qsub -v "input_read=../reads/line7i.se.fq,\
		sample_name=line7i-single-rsem,index=gimme-gga-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd gimme; \
	qsub -v "input_read1=../reads/line7u.pe.1,\
		input_read2=../reads/line7u.pe.2,sample_name=line7u-paired-rsem,\
		index=gimme-gga-rsem" $(protocol)/rsem_calculate_expr_paired.sh

	cd gimme; \
	qsub -v "input_read1=../reads/line7i.pe.1,\
		input_read2=../reads/line7i.pe.2,sample_name=line7i-paired-rsem,\
		index=gimme-gga-rsem" $(protocol)/rsem_calculate_expr_paired.sh

ebseq-gimme:

	cd gimme; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix

	cd gimme; \
	rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes

	cd gimme; \
	rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-tophits-degenes:

	cd gimme; \
		python $(protocol)/tophits-to-degenes-gimme.py \
		line7u_vs_i.degenes.fdr.05 gimme-gga-tophits.txt \
		> line7u_vs_i.degenes.fdr.05.gga.tophits

	cd gimme; \
		python $(protocol)/get_top_hits.py gimme-hsa.xml \
		> gimme-hsa-tophits.txt

	cd gimme; \
		python $(protocol)/tophits-to-degenes-gimme.py \
		line7u_vs_i.degenes.fdr.05 gimme-hsa-tophits.txt \
		> line7u_vs_i.degenes.fdr.05.hsa.tophits
