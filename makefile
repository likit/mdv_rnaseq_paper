###############################################
##### De novo assembly with Velvet+Oases ######
###############################################

run-quality-trim-pe:

	qsub -v left="reads/line7u.pe.1,right=reads/line7u.pe.2" protocol/quality_trim_pe_job.sh
	qsub -v left="reads/line7i.pe.1,right=reads/line7i.pe.2" protocol/quality_trim_pe_job.sh
	qsub -v left="reads/line6u.pe.1,right=reads/line6u.pe.2" protocol/quality_trim_pe_job.sh
	qsub -v left="reads/line6i.pe.1,right=reads/line6i.pe.2" protocol/quality_trim_pe_job.sh

run-quality-trim-se:

	for r in reads/*.se.fq; do qsub -v input="$$r" protocol/quality_trim_se_job.sh; done

interleave-reads:

	cd assembly; ~/velvet_1.2.03/shuffleSequences_fastq.pl pe.1.fastq pe.2.fastq paired.fastq

run-velveth:

	cd assembly; qsub -v pe_input="paired.fastq",se_input="single.fastq" ~/rnaseq-comp-protocol/velveth_job.sh

run-velvetg:

	cd assembly; qsub ~/rnaseq-comp-protocol/velvetg_job.sh

run-oases:

	cd assembly; qsub ~/rnaseq-comp-protocol/oases_job.sh

###############################################
###### Local assembly with Velvet+Oases #######
###############################################

local-assembly: run-tophat-pe run-tophat-se extract-reads
run-tophat-pe:

	cd tophat; qsub -v outdir="line6u_pe",index="gal4selected",left="../reads/line6u.1_trim1.fastq",right="../reads/line6u.1_trim2.fastq",unpaired="../reads/line6u.1_trim_unpaired.fastq" ~/mdv-protocol/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line6i_pe",index="gal4selected",left="../reads/line6i.1_trim1.fastq",right="../reads/line6i.1_trim2.fastq",unpaired="../reads/line6i.1_trim_unpaired.fastq" ~/mdv-protocol/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7u_pe",index="gal4selected",left="../reads/line7u.1_trim1.fastq",right="../reads/line7u.1_trim2.fastq",unpaired="../reads/line7u.1_trim_unpaired.fastq" ~/mdv-protocol/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7i_pe",index="gal4selected",left="../reads/line7i.1_trim1.fastq",right="../reads/line7i.1_trim2.fastq",unpaired="../reads/line7i.1_trim_unpaired.fastq" ~/mdv-protocol/tophat_pe_job.sh 

run-tophat-se:

	cd tophat; qsub -v outdir="line6u_se",index="gal4selected",input="../reads/line6u.fq_trim.fastq" ~/mdv-protocol/tophat_se_job.sh
	cd tophat; qsub -v outdir="line6i_se",index="gal4selected",input="../reads/line6i.fq_trim.fastq" ~/mdv-protocol/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7u_se",index="gal4selected",input="../reads/line7u.fq_trim.fastq" ~/mdv-protocol/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7i_se",index="gal4selected",input="../reads/line7i.fq_trim.fastq" ~/mdv-protocol/tophat_se_job.sh

extract-reads:

	cd tophat; for dir in line??_?e; \
		do ~/mdv-protocol/extract_reads.sh $$dir/accepted_hits.bam ~/mdv-protocol/chromosomes.txt; \
	done

merge-bams:

	cd tophat; \
	for chr in $(cat ~/mdv-protocol/chromosomes.txt); do printf "merging %s..\n" "$chr";  \
		samtools merge -n merged/"$chr".bam \
		line6u_pe/"$chr".bam line6u_se/"$chr".bam \
		line6i_pe/"$chr".bam line6i_se/"$chr".bam \
		line7u_pe/"$chr".bam line7u_se/"$chr".bam \
		line7i_pe/"$chr".bam line7i_se/"$chr".bam; \
	done

run-velveth-local:

	cd tophat/merged; \
	for f in *.bam; \
		do qsub -v outdir=$$(basename "$$f" .bam),input="$$f" ~/mdv-protocol/velveth_local_job.sh; \
	done

run-velvetg-local:

	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" ~/mdv-protocol/velvetg_local_job.sh; \
	done

run-oases-local:

	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" ~/mdv-protocol/oases_local_job.sh; \
	done

combine-transcripts:

	cd tophat/merged; \
		for d in chr*_[0-9][0-9]; \
			do python ~/gimme/src/utils/rename_fasta.py $$d/transcripts.fa local_$$d >> local_merged.fa; \
	done

	cd assembly; \
		for d in global_[0-9][0-9]; \
			do python ~/gimme/src/utils/rename_fasta.py $$d/transcripts.fa global_$$d >> global_merged.fa; \
	done

clean-transcripts:

	cd tophat/merged; ~/seqclean-x86_64/seqclean local_merged.fa
	cd assembly; ~/seqclean-x86_64/seqclean global_merged.fa

remove-redundant-seq:

	cat tophat/merged/local_merged.fa.clean assembly/global_merged.fa.clean >> all.fa.clean
	qsub -v input="all.fa.clean",output="all.fa.clean.nr",c="1.0" ~/mdv-protocol/cdhit_job.sh

	cd assembly; qsub -v input="global_merged.fa.clean",output="global_merged.fa.clean.nr",c="1.0" ~/mdv-protocol/cdhit_job.sh

align-transcripts:

	python ~/mdv-protocol/split-fa.py all.fa.clean.nr
	for f in subsets*.fa; do \
		qsub -v input="$$f" ~/mdv-protocol/blat_job.sh; \
	done
	cat subsets*.fa.psl > all.fa.clean.nr.psl
	sort -k 10 all.fa.clean.nr.psl > all.fa.clean.nr.psl.sorted
	pslReps -nohead -singleHit all.fa.clean.nr.psl.sorted all.fa.clean.nr.psl.best info
	rm subsets*.fa.psl
	rm subsets*.fa

run-cufflinks:

	cd tophat; for d in line??_?e; do qsub -v outdir="$$d",input="$$d/accepted_hits.bam" \
		~/mdv-protocol/cufflinks_job.sh; echo $$d; done

run-cuffmerge-ref:

	cd tophat; cuffmerge -g ../Gallus_UCSC_ensembl_73.gtf.removed -o merged_cuff_ref -s gal4selected.fa -p 4 ~/mdv-protocol/merge_list.txt

build-gene-models-with-cufflinks-ref:

	python ~/gimme/src/utils/gff2bed.py tophat/merged_cuff_ref/merged.gtf > tophat/merged_cuff_ref/merged.bed
	qsub -v output="asm_cuff_ref_models.bed",input1="all.fa.clean.nr.psl.best",input2="tophat/merged_cuff_ref/merged.bed",ref="tophat/gal4selected.fa" ~/mdv-protocol/run_gimme2.sh

models-to-transcripts:

	python ~/gimme/src/utils/get_transcript_seq.py asm_cuff_ref_models.bed tophat/gal4selected.fa | sed 1d > asm_cuff_ref_models.bed.fa

rsem-gimme-models-ref:

	cat asm_cuff_ref_models.bed.fa | python ~/mdv-protocol/fasta-to-gene-list.py > asm_cuff_ref_models.txt
	qsub -v list="asm_cuff_ref_models.txt",input="asm_cuff_ref_models.bed.fa",sample="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_prepare_reference.sh

rsem-calc-gimme-models-ref:

	qsub -v input_read="reads/line6u.se.fq",sample_name="line6u-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line6i.se.fq",sample_name="line6i-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7u.se.fq",sample_name="line7u-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7i.se.fq",sample_name="line7i-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line6u.pe.1",input_read2="reads/line6u.pe.2",sample_name="line6u-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line6i.pe.1",input_read2="reads/line6i.pe.2",sample_name="line6i-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7u.pe.1",input_read2="reads/line7u.pe.2",sample_name="line7u-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7i.pe.1",input_read2="reads/line7i.pe.2",sample_name="line7i-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh

rsem-calc-gimme-models-ref-rspd:

	qsub -v input_read="reads/line6u.se.fq",sample_name="line6u-single-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line6i.se.fq",sample_name="line6i-single-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line6u.pe.1",input_read2="reads/line6u.pe.2",sample_name="line6u-paired-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line6i.pe.1",input_read2="reads/line6i.pe.2",sample_name="line6i-paired-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh

	#qsub -v input_read="reads/line7u.se.fq",sample_name="line7u-single-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" \
	#	~/mdv-protocol/rsem_calculate_expr_single.sh
	#qsub -v input_read="reads/line7i.se.fq",sample_name="line7i-single-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" \
	#	~/mdv-protocol/rsem_calculate_expr_single.sh

	#qsub -v input_read1="reads/line7u.pe.1",input_read2="reads/line7u.pe.2",sample_name="line7u-paired-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh
	#qsub -v input_read1="reads/line7i.pe.1",input_read2="reads/line7i.pe.2",sample_name="line7i-paired-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh

ebseq-line6-models-ref:

	rsem-generate-data-matrix line6u-single-rsem-cuffref.genes.results \
		line6u-paired-rsem-cuffref.genes.results line7i-single-rsem-cuffref.genes.results \
		line7i-paired-rsem-cuffref.genes.results > line6u_vs_i.gene.cuffref.counts.matrix
	rsem-run-ebseq line6u_vs_i.gene.cuffref.counts.matrix 2,2 line6u_vs_i.cuffref.degenes
	rsem-control-fdr line6u_vs_i.cuffref.degenes 0.05 line6u_vs_i.cuffref.degenes.fdr.05

ebseq-line7-models-ref:

	rsem-generate-data-matrix line7u-single-rsem-cuffref.genes.results \
		line7u-paired-rsem-cuffref.genes.results line7i-single-rsem-cuffref.genes.results \
		line7i-paired-rsem-cuffref.genes.results > line7u_vs_i.gene.cuffref.counts.matrix
	rsem-run-ebseq line7u_vs_i.gene.cuffref.counts.matrix 2,2 line7u_vs_i.cuffref.degenes
	rsem-control-fdr line7u_vs_i.cuffref.degenes 0.05 line7u_vs_i.cuffref.degenes.fdr.05

filter-low-isopct:

	python ~/mdv-protocol/filter-low-isopct.py 1.0 asm_cuff_models.bed *isoforms.results > asm_cuff_models.flt.bed

rsem-output-cuffref-to-fasta:

	python ~/mdv-protocol/rsem-output-to-fasta.py line6u_vs_i.cuffref.degenes.fdr.05 asm_cuff_ref_models.bed.fa > line6u_vs_i.cuffref.degenes.fdr.05.fa
	python ~/mdv-protocol/rsem-output-to-fasta.py line7u_vs_i.cuffref.degenes.fdr.05 asm_cuff_ref_models.bed.fa > line7u_vs_i.cuffref.degenes.fdr.05.fa

translate-degenes:

	estscan -t line6u_vs_i.degenes.fdr.05.fa.prot -M ~/mdv-protocol/gallus.hm line6u_vs_i.degenes.fdr.05.fa > line6u_vs_i.degenes.fdr.05.fa.nucl
	estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M ~/mdv-protocol/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

get-longest-sequences-cuffref:

	python ~/mdv-protocol/gene-rep.py line6u_vs_i.cuffref.degenes.fdr.05.fa > line6u_vs_i.cuffref.degenes.fdr.05.fa.longest
	#python ~/mdv-protocol/gene-rep.py line7u_vs_i.cuffref.degenes.fdr.05.fa > line7u_vs_i.cuffref.degenes.fdr.05.fa.longest
	#python ~/mdv-protocol/gene-rep.py line6u_vs_i.cuffref.degenes.fdr.05.fa.prot > line6u_vs_i.cuffref.degenes.fdr.05.fa.prot.longest
	#python ~/mdv-protocol/gene-rep.py line7u_vs_i.cuffref.degenes.fdr.05.fa.prot > line7u_vs_i.cuffref.degenes.fdr.05.fa.prot.longest

run-blast-models-ref:

	#qsub -v db="Gallus_prot",input="line6u_vs_i.cuffref.degenes.fdr.05.fa.prot.longest",program="blastp",output="line6u_vs_i.cuffref.degenes.fdr.05.fa.prot.longest.gallus.xml" ~/mdv-protocol/blast.sh
	#qsub -v db="Gallus_prot",input="line7u_vs_i.cuffref.degenes.fdr.05.fa.prot.longest",program="blastp",output="line7u_vs_i.cuffref.degenes.fdr.05.fa.prot.longest.gallus.xml" ~/mdv-protocol/blast.sh

	qsub -v db="Gallus_prot",input="line6u_vs_i.cuffref.degenes.fdr.05.fa.longest",program="blastx",output="line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.gallus.v2.xml" ~/mdv-protocol/blast.sh
	#qsub -v db="Gallus_prot",input="line7u_vs_i.cuffref.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.gallus.xml" ~/mdv-protocol/blast.sh

get-tophits:

	python protocol/get_top_hits.py line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.gallus.xml \
		> line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.xml

	#python protocol/get_top_hits.py line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.gallus.xml \
	#	> line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.xml

##TODO
