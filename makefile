local-assembly: run-tophat-pe run-tophat-se extract-reads
run-tophat-pe:
	cd tophat; qsub -v outdir="line6u_pe",index="gal4selected",left="../reads/line6u.1_trim1.fastq",right="../reads/line6u.1_trim2.fastq",unpaired="../reads/line6u.1_trim_unpaired.fastq" ../protocol/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line6i_pe",index="gal4selected",left="../reads/line6i.1_trim1.fastq",right="../reads/line6i.1_trim2.fastq",unpaired="../reads/line6i.1_trim_unpaired.fastq" ../protocol/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7u_pe",index="gal4selected",left="../reads/line7u.1_trim1.fastq",right="../reads/line7u.1_trim2.fastq",unpaired="../reads/line7u.1_trim_unpaired.fastq" ../protocol/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7i_pe",index="gal4selected",left="../reads/line7i.1_trim1.fastq",right="../reads/line7i.1_trim2.fastq",unpaired="../reads/line7i.1_trim_unpaired.fastq" ../protocol/tophat_pe_job.sh 

run-tophat-se:
	cd tophat; qsub -v outdir="line6u_se",index="gal4selected",input="../reads/line6u.fq_trim.fastq" ../protocol/tophat_se_job.sh
	cd tophat; qsub -v outdir="line6i_se",index="gal4selected",input="../reads/line6i.fq_trim.fastq" ../protocol/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7u_se",index="gal4selected",input="../reads/line7u.fq_trim.fastq" ../protocol/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7i_se",index="gal4selected",input="../reads/line7i.fq_trim.fastq" ../protocol/tophat_se_job.sh

extract-reads:
	cd tophat; for dir in line??_?e; \
		do ../protocol/extract_reads.sh $$dir/accepted_hits.bam ../protocol/chromosomes.txt; \
	done

merge-bams:
	cd tophat; \
	for chr in $(cat ../protocol/chromosomes.txt); do printf "merging %s..\n" "$chr";  \
		samtools merge -n merged/"$chr".bam \
		line6u_pe/"$chr".bam line6u_se/"$chr".bam \
		line6i_pe/"$chr".bam line6i_se/"$chr".bam \
		line7u_pe/"$chr".bam line7u_se/"$chr".bam \
		line7i_pe/"$chr".bam line7i_se/"$chr".bam; \
	done

run-velveth-local:
	cd tophat/merged; \
	for f in *.bam; \
		do qsub -v outdir=$$(basename "$$f" .bam),input="$$f" ../../protocol/velveth_local_job.sh; \
	done

run-velvetg-local:
	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" ../../protocol/velvetg_local_job.sh; \
	done

run-oases-local:
	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" ../../protocol/oases_local_job.sh; \
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

	#cat tophat/merged/local_merged.fa.clean assembly/global_merged.fa.clean >> all.fa.clean
	#qsub -v input="all.fa.clean",output="all.fa.clean.nr",c="1.0" protocol/cdhit_job.sh

	cd assembly; qsub -v input="global_merged.fa.clean",output="global_merged.fa.clean.nr",c="1.0" ../protocol/cdhit_job.sh

align-transcripts:

	python protocol/split-fa.py all.fa.clean.nr
	for f in subsets*.fa; do \
		qsub -v input="$$f" protocol/blat_job.sh; \
	done
	cat subsets*.fa.psl > all.fa.clean.nr.psl
	sort -k 10 all.fa.clean.nr.psl > all.fa.clean.nr.psl.sorted
	pslReps -nohead -singleHit all.fa.clean.nr.psl.sorted all.fa.clean.nr.psl.best info
	rm subsets*.fa.psl
	rm subsets*.fa

run-cufflinks:

	cd tophat; for d in line??_?e; do qsub -v outdir="$$d",input="$$d/accepted_hits.bam" \
		../protocol/cufflinks_job.sh; echo $$d; done

run-cuffmerge:

	cd tophat; cuffmerge -o merged_cuff_denovo -s gal4selected.fa -p 4 ../protocol/merge_list.txt

build-gene-models:

	qsub -v input="all.fa.clean.nr.psl.best",ref="tophat/gal4selected.fa" protocol/run_gimme.sh

build-gene-models-with-cufflinks:

	python ~/gimme/src/utils/gff2bed.py tophat/merged_cuff_denovo/transcripts.gtf > tophat/merged_cuff_denovo/transcripts.bed
	qsub -v input1="all.fa.clean.nr.psl.best",input2="tophat/merged_cuff_denovo/transcripts.bed",ref="tophat/gal4selected.fa" protocol/run_gimme2.sh

clean-gene-models:

	sed 1d all.fa.clean.nr.psl.best.bed > asm_models.bed
	sed 1d all.fa.clean.nr.psl.best.merged.bed > asm_cuff_models.bed
	python ~/gimme/src/utils/get_transcript_seq.py asm_models.bed tophat/gal4selected.fa | sed 1d > asm_models.bed.fa
	python ~/gimme/src/utils/get_transcript_seq.py asm_cuff_models.bed tophat/gal4selected.fa | sed 1d > asm_cuff_models.bed.fa
	qsub -v input="asm_models.bed.fa",output="asm_models.bed.fa.nr99",c="0.99" protocol/cdhit_job.sh
	qsub -v input="asm_cuff_models.bed.fa",output="asm_cuff_models.bed.fa.nr99",c="0.99" protocol/cdhit_job.sh
	python ~/gimme/src/utils/cdhit_transcript.py asm_models.bed asm_models.bed.fa.nr99 > asm_models.nr99.bed
	python ~/gimme/src/utils/cdhit_transcript.py asm_cuff_models.bed asm_cuff_models.bed.fa.nr99 > asm_cuff_models.nr99.bed

rsem-gimme-models:
	#cat asm_cuff_models.bed.fa | python protocol/fasta-to-gene-list.py > asm_cuff_models.txt
	#qsub -v list="asm_cuff_models.txt",input="asm_cuff_models.bed.fa",sample="asm_cuff_models_rsem" protocol/rsem_prepare_reference.sh

	#qsub -v input_read="reads/line6u.se.fq",sample_name="line6u-single-rsem-full",index="asm_cuff_models_rsem" \
	#	protocol/rsem_calculate_expr_single.sh
	#qsub -v input_read="reads/line6i.se.fq",sample_name="line6i-single-rsem-full",index="asm_cuff_models_rsem" \
	#	protocol/rsem_calculate_expr_single.sh
	#qsub -v input_read="reads/line7u.se.fq",sample_name="line7u-single-rsem-full",index="asm_cuff_models_rsem" \
	#	protocol/rsem_calculate_expr_single.sh
	#qsub -v input_read="reads/line7i.se.fq",sample_name="line7i-single-rsem-full",index="asm_cuff_models_rsem" \
	#	protocol/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line6u.pe.1",input_read2="reads/line6u.pe.2",sample_name="line6u-paired-rsem-full",index="asm_cuff_models_rsem" protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line6i.pe.1",input_read2="reads/line6i.pe.2",sample_name="line6i-paired-rsem-full",index="asm_cuff_models_rsem" protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7u.pe.1",input_read2="reads/line7u.pe.2",sample_name="line7u-paired-rsem-full",index="asm_cuff_models_rsem" protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7i.pe.1",input_read2="reads/line7i.pe.2",sample_name="line7i-paired-rsem-full",index="asm_cuff_models_rsem" protocol/rsem_calculate_expr_paired.sh

filter-low-isopct:

	python protocol/filter-low-isopct.py 1.0 asm_cuff_models.bed *isoforms.results > asm_cuff_models.flt.bed

ebseq-line6:

	rsem-generate-data-matrix line6u-single-rsem-full.genes.results \
		line6u-paired-rsem-full.genes.results line6i-single-rsem-full.genes.results \
		line6i-paired-rsem-full.genes.results > line6u_vs_i.gene.counts.matrix
	rsem-run-ebseq line6u_vs_i.gene.counts.matrix 2,2 line6u_vs_i.degenes
	rsem-control-fdr line6u_vs_i.degenes 0.05 line6u_vs_i.degenes.fdr.05

ebseq-line7:

	rsem-generate-data-matrix line7u-single-rsem-full.genes.results \
		line7u-paired-rsem-full.genes.results line7i-single-rsem-full.genes.results \
		line7i-paired-rsem-full.genes.results > line7u_vs_i.gene.counts.matrix
	rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

translate-degenes:

	#estscan -t line6u_vs_i.degenes.fdr.05.fa.prot -M protocol/gallus.hm line6u_vs_i.degenes.fdr.05.fa > line6u_vs_i.degenes.fdr.05.fa.nucl
	estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M protocol/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

get_longest_sequences:
	python protocol/gene-rep.py line6u_vs_i.degenes.fdr.05.fa > line6u_vs_i.degenes.fdr.05.fa.longest
	python protocol/gene-rep.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest
	python protocol/gene-rep.py line6u_vs_i.degenes.fdr.05.fa.prot > line6u_vs_i.degenes.fdr.05.fa.prot.longest
	python protocol/gene-rep.py line7u_vs_i.degenes.fdr.05.fa.prot > line7u_vs_i.degenes.fdr.05.fa.prot.longest

split_sequences:

	python protocol/split-fa.py line6u_vs_i.degenes.fdr.05.fa.prot.longest 100 line6u_vs_i.degenes.fdr.05.fa.prot.longest
	python protocol/split-fa.py line7u_vs_i.degenes.fdr.05.fa.prot.longest 100 line7u_vs_i.degenes.fdr.05.fa.prot.longest

run-blastp:

	for f in *prot.longest*.fa; do \
		qsub -v input="$$f",program="blastp" protocol/blast.sh; \
	done

run-blastx:
