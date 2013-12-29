local-assembly: run-tophat-pe run-tophat-se extract-reads
run-tophat-pe:
	cd tophat; qsub -v outdir="line6u_pe",index="gal4selected",left="../reads/line6u.1_trim1.fastq",right="../reads/line6u.1_trim2.fastq",unpaired="../reads/line6u.1_trim_unpaired.fastq" ../protocols/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line6i_pe",index="gal4selected",left="../reads/line6i.1_trim1.fastq",right="../reads/line6i.1_trim2.fastq",unpaired="../reads/line6i.1_trim_unpaired.fastq" ../protocols/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7u_pe",index="gal4selected",left="../reads/line7u.1_trim1.fastq",right="../reads/line7u.1_trim2.fastq",unpaired="../reads/line7u.1_trim_unpaired.fastq" ../protocols/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7i_pe",index="gal4selected",left="../reads/line7i.1_trim1.fastq",right="../reads/line7i.1_trim2.fastq",unpaired="../reads/line7i.1_trim_unpaired.fastq" ../protocols/tophat_pe_job.sh 

run-tophat-se:
	cd tophat; qsub -v outdir="line6u_se",index="gal4selected",input="../reads/line6u.fq_trim.fastq" ../protocols/tophat_se_job.sh
	cd tophat; qsub -v outdir="line6i_se",index="gal4selected",input="../reads/line6i.fq_trim.fastq" ../protocols/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7u_se",index="gal4selected",input="../reads/line7u.fq_trim.fastq" ../protocols/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7i_se",index="gal4selected",input="../reads/line7i.fq_trim.fastq" ../protocols/tophat_se_job.sh

extract-reads:
	cd tophat; for dir in line??_?e; \
		do ../protocols/extract_reads.sh $$dir/accepted_hits.bam chromosomes.txt; \
	done

merge-bams:
	for chr in $(cat chromosomes.txt); do printf "merging %s..\n" "$chr";  \
		samtools merge -n merged/"$chr".bam \
		line6u_pe/"$chr".bam line6u_se/"$chr".bam \
		line6i_pe/"$chr".bam line6i_se/"$chr".bam \
		line7u_pe/"$chr".bam line7u_se/"$chr".bam \
		line7i_pe/"$chr".bam line7i_se/"$chr".bam; \
	done

run-velveth-local:
	cd tophat/merged; \
	for f in *.bam; \
		do qsub -v outdir=$$(basename "$$f" .bam),input="$$f" ../../protocols/velveth_local_job.sh; \
	done

run-velvetg-local:
	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" ../../protocols/velvetg_local_job.sh; \
	done

run-oases-local:
	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" ../../protocols/oases_local_job.sh; \
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
	qsub protocols/cdhit.sh

align-transcripts:
	python protocols/split-fa.py all.fa.clean.nr
	for f in subsets*.fa; do \
		qsub -v input="$$f" protocols/blat_job.sh; \
	done
