run-tophat-se:

	# create a directory
	if [ ! -d miso/bam-data ]; then mkdir -p miso/bam-data; fi

	# Tophat 2.0.9 and Bowtie 2.1.0 are required
	cd miso/bam-data; qsub -v "outdir=line6u_se,index=../../gal4selected,\
		input=../../reads/line6u.se.fq" $(protocol)/tophat_se_job.sh
	cd miso/bam-data; qsub -v "outdir=line6i_se,index=../../gal4selected,\
		input=../../reads/line6i.se.fq" $(protocol)/tophat_se_job.sh
	cd miso/bam-data; qsub -v "outdir=line7u_se,index=../../gal4selected,\
		input=../../reads/line7u.se.fq" $(protocol)/tophat_se_job.sh
	cd miso/bam-data; qsub -v "outdir=line7i_se,index=../../gal4selected,\
		input=../../reads/line7i.se.fq" $(protocol)/tophat_se_job.sh

	cd miso/bam-data; qsub -v "outdir=line6u_pe,index=../../gal4selected,\
		input=../../reads/line6u.pe.1" $(protocol)/tophat_se_job.sh
	cd miso/bam-data; qsub -v "outdir=line6i_pe,index=../../gal4selected,\
		input=../../reads/line6i.pe.1" $(protocol)/tophat_se_job.sh
	cd miso/bam-data; qsub -v "outdir=line7u_pe,index=../../gal4selected,\
		input=../../reads/line7u.pe.1" $(protocol)/tophat_se_job.sh
	cd miso/bam-data; qsub -v "outdir=line7i_pe,index=../../gal4selected,\
		input=../../reads/line7i.pe.1" $(protocol)/tophat_se_job.sh

merge-bam-files:

	cd miso/bam-data; \
		samtools view -H line6u_se/accepted_hits.bam > inh.sam; \
		samtools index line6u.bam

	cd miso/bam-data; \
		samtools merge -h inh.sam line6u.bam \
		line6u_se/accepted_hits.bam line6u_pe/accepted_hits.bam; \
		samtools index line6u.bam

	cd miso/bam-data; \
		samtools merge -h inh.sam line6i.bam \
		line6i_se/accepted_hits.bam line6i_pe/accepted_hits.bam; \
		samtools index line6i.bam

	cd miso/bam-data; \
		samtools merge -h inh.sam line7u.bam \
		line7u_se/accepted_hits.bam line7u_pe/accepted_hits.bam; \
		samtools index line7u.bam

	cd miso/bam-data; \
		samtools merge -h inh.sam line7i.bam \
		line7i_se/accepted_hits.bam line7i_pe/accepted_hits.bam; \
		samtools index line7i.bam

filter-low-isopct:

	cd miso; \
	python $(protocol)/filter-low-isopct.py 1.0 \
		../gimme/gimme.bed ../gimme/*isoforms.results > gimme.flt.bed

build-se-models:

	cd miso; \
		python $(gimmedir)/src/utils/find_SE.py gimme.flt.bed > gimme.flt.SE.gff
	cd miso; \
		index_gff.py --index gimme.flt.SE.gff indexes/SE

build-a3ss-models:

	cd miso; \
		python $(gimmedir)/src/utils/find_A3SS.py gimme.flt.bed > gimme.flt.A3SS.gff
	cd miso; \
		index_gff.py --index gimme.flt.A3SS.gff indexes/A3SS

build-a5ss-models:

	cd miso; \
		python $(gimmedir)/src/utils/find_A5SS.py gimme.flt.bed > gimme.flt.A5SS.gff
	cd miso; \
		index_gff.py --index gimme.flt.A5SS.gff indexes/A5SS

run-miso-se:

	cd miso; \
		qsub -v "input_bam=bam-data/line6u.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/SE,output_dir=results/SE/line6u,\
		event=SE" $(protocol)/miso.sh

	cd miso; \
		qsub -v "input_bam=bam-data/line6i.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/SE,output_dir=results/SE/line6i,\
		event=SE" $(protocol)/miso.sh

	cd miso; \
		qsub -v "input_bam=bam-data/line7u.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/SE,output_dir=results/SE/line7u,\
		event=SE" $(protocol)/miso.sh

	cd miso; \
		qsub -v "input_bam=bam-data/line7i.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/SE,output_dir=results/SE/line7i,\
		event=SE" $(protocol)/miso.sh

run-miso-a3ss:

	cd miso; \
		qsub -v "input_bam=bam-data/line6u.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/A3SS,output_dir=results/A3SS/line6u,\
		event=A3SS" $(protocol)/miso.sh

	cd miso; \
		qsub -v "input_bam=bam-data/line6i.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/A3SS,output_dir=results/A3SS/line6i,\
		event=A3SS" $(protocol)/miso.sh

	cd miso; \
		qsub -v "input_bam=bam-data/line7u.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/A3SS,output_dir=results/A3SS/line7u,\
		event=A3SS" $(protocol)/miso.sh

	cd miso; \
		qsub -v "input_bam=bam-data/line7i.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/A3SS,output_dir=results/A3SS/line7i,\
		event=A3SS" $(protocol)/miso.sh

run-miso-a5ss:

	cd miso; \
		qsub -v "input_bam=bam-data/line6u.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/A5SS,output_dir=results/A5SS/line6u,\
		event=A5SS" $(protocol)/miso.sh

	cd miso; \
		qsub -v "input_bam=bam-data/line6i.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/A5SS,output_dir=results/A5SS/line6i,\
		event=A5SS" $(protocol)/miso.sh

	cd miso; \
		qsub -v "input_bam=bam-data/line7u.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/A5SS,output_dir=results/A5SS/line7u,\
		event=A5SS" $(protocol)/miso.sh

	cd miso; \
		qsub -v "input_bam=bam-data/line7i.bam,\
		setting_file=$(protocol)/miso_settings.txt,\
		index_dir=indexes/A5SS,output_dir=results/A5SS/line7i,\
		event=A5SS" $(protocol)/miso.sh

summarize-se:

	cd miso; run_miso.py --summarize-samples results/SE/line6u results/SE/line6u
	cd miso; run_miso.py --summarize-samples results/SE/line6i results/SE/line6i
	cd miso; run_miso.py --summarize-samples results/SE/line7u results/SE/line7u
	cd miso; run_miso.py --summarize-samples results/SE/line7i results/SE/line7i

summarize-a3ss:

	cd miso; run_miso.py --summarize-samples results/A3SS/line6u results/A3SS/line6u
	cd miso; run_miso.py --summarize-samples results/A3SS/line6i results/A3SS/line6i
	cd miso; run_miso.py --summarize-samples results/A3SS/line7u results/A3SS/line7u
	cd miso; run_miso.py --summarize-samples results/A3SS/line7i results/A3SS/line7i

summarize-a5ss:

	cd miso; run_miso.py --summarize-samples results/A5SS/line6u results/A5SS/line6u
	cd miso; run_miso.py --summarize-samples results/A5SS/line6i results/A5SS/line6i
	cd miso; run_miso.py --summarize-samples results/A5SS/line7u results/A5SS/line7u
	cd miso; run_miso.py --summarize-samples results/A5SS/line7i results/A5SS/line7i

compare-miso-se:

	cd miso; run_miso.py --compare-samples results/SE/line6u results/SE/line7u results/SE/comparisons
	cd miso; run_miso.py --compare-samples results/SE/line6i results/SE/line7i results/SE/comparisons

	cd miso; run_miso.py --compare-samples results/SE/line6u results/SE/line6i results/SE/comparisons
	cd miso; run_miso.py --compare-samples results/SE/line7u results/SE/line7i results/SE/comparisons

compare-miso-a3ss:

	cd miso; run_miso.py --compare-samples results/A3SS/line6u results/A3SS/line7u results/A3SS/comparisons
	cd miso; run_miso.py --compare-samples results/A3SS/line6i results/A3SS/line7i results/A3SS/comparisons
	cd miso; run_miso.py --compare-samples results/A3SS/line6u results/A3SS/line6i results/A3SS/comparisons
	cd miso; run_miso.py --compare-samples results/A3SS/line7u results/A3SS/line7i results/A3SS/comparisons

compare-miso-a5ss:

	cd miso; run_miso.py --compare-samples results/A5SS/line6u results/A5SS/line7u results/A5SS/comparisons
	cd miso; run_miso.py --compare-samples results/A5SS/line6i results/A5SS/line7i results/A5SS/comparisons
	cd miso; run_miso.py --compare-samples results/A5SS/line6u results/A5SS/line6i results/A5SS/comparisons
	cd miso; run_miso.py --compare-samples results/A5SS/line7u results/A5SS/line7i results/A5SS/comparisons

filter-miso-se:

	cd miso; \
		python $(protocol)/miso-filter.py \
		results/SE/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf \
		SE 0.20 0.20 10 2 > \
		results/SE/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt

	cd miso; \
		python $(protocol)/miso-filter.py \
		results/SE/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf \
		SE 0.20 0.20 10 2 > \
		results/SE/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt

	cd miso; \
		python $(protocol)/miso-filter.py \
		results/SE/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf \
		SE 0.20 0.20 10 2 > \
		results/SE/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt

	cd miso; \
		python $(protocol)/miso-filter.py \
		results/SE/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf \
		SE 0.20 0.20 10 2 > \
		results/SE/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt

filter-miso-a3ss:

	cd miso; python $(protocol)/miso-filter.py \
		results/A3SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf \
		ASS 0.20 0.20 10 > \
	results/A3SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt

	cd miso; python $(protocol)/miso-filter.py \
		results/A3SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf \
		ASS 0.20 0.20 10 > \
		results/A3SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt

	cd miso; python $(protocol)/miso-filter.py \
		results/A3SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf \
		ASS 0.20 0.20 10 > \
		results/A3SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt

	cd miso; python $(protocol)/miso-filter.py \
		results/A3SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf \
		ASS 0.20 0.20 10 > \
		results/A3SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt

filter-miso-a5ss:

	cd miso; python $(protocol)/miso-filter.py \
		results/A5SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf \
		ASS 0.20 0.20 10 > \
		results/A5SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt

	cd miso; python $(protocol)/miso-filter.py \
		results/A5SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf \
		ASS 0.20 0.20 10 > \
		results/A5SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt

	cd miso; python $(protocol)/miso-filter.py \
		results/A5SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf \
		ASS 0.20 0.20 10 > \
		results/A5SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt

	cd miso; python $(protocol)/miso-filter.py \
		results/A5SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf \
		ASS 0.20 0.20 10 > \
		results/A5SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt

miso-to-fa-se:

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/miso-to-fa.py line6i_vs_line7i.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line6i_vs_line7i.miso_bf.flt.fa

	cd miso/results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		python $(protocol)/miso-to-fa.py line6u_vs_line7u.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line6u_vs_line7u.miso_bf.flt.fa

	cd miso/results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		python $(protocol)/miso-to-fa.py line7u_vs_line7i.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line7u_vs_line7i.miso_bf.flt.fa

	cd miso/results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		python $(protocol)/miso-to-fa.py line6u_vs_line6i.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line6u_vs_line6i.miso_bf.flt.fa

miso-to-fa-a3ss:

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/miso-to-fa.py line6i_vs_line7i.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line6i_vs_line7i.miso_bf.flt.fa

	cd miso/results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python $(protocol)/miso-to-fa.py line6u_vs_line7u.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line6u_vs_line7u.miso_bf.flt.fa

	cd miso/results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python $(protocol)/miso-to-fa.py line7u_vs_line7i.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line7u_vs_line7i.miso_bf.flt.fa

	cd miso/results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python $(protocol)/miso-to-fa.py line6u_vs_line6i.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line6u_vs_line6i.miso_bf.flt.fa

miso-to-fa-a5ss:

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/miso-to-fa.py line6i_vs_line7i.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line6i_vs_line7i.miso_bf.flt.fa

	cd miso/results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python $(protocol)/miso-to-fa.py line6u_vs_line7u.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line6u_vs_line7u.miso_bf.flt.fa

	cd miso/results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python $(protocol)/miso-to-fa.py line7u_vs_line7i.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line7u_vs_line7i.miso_bf.flt.fa

	cd miso/results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python $(protocol)/miso-to-fa.py line6u_vs_line6i.miso_bf.flt \
		../../../../../../gimme/gimme.bed.fa > line6u_vs_line6i.miso_bf.flt.fa

blast-miso-se:

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
	qsub -v "db=Gallus_prot,input=line6i_vs_line7i.miso_bf.flt.fa,\
	program=blastx,output=line6i_vs_line7i.miso_bf.flt.fa.gallus.xml" \
	$(protocol)/blast.sh

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
	qsub -v "db=Human_prot,input=line6i_vs_line7i.miso_bf.flt.fa,\
	program=blastx,output=line6i_vs_line7i.miso_bf.flt.fa.human.xml" \
	$(protocol)/blast.sh

blast-miso-a3ss:

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
	qsub -v "db=Gallus_prot,input=line6i_vs_line7i.miso_bf.flt.fa,\
	program=blastx,output=line6i_vs_line7i.miso_bf.flt.fa.gallus.xml" \
	$(protocol)/blast.sh

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
	qsub -v "db=Human_prot,input=line6i_vs_line7i.miso_bf.flt.fa,\
	program=blastx,output=line6i_vs_line7i.miso_bf.flt.fa.human.xml" \
	$(protocol)/blast.sh

blast-miso-a5ss:

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
	qsub -v "db=Gallus_prot,input=line6i_vs_line7i.miso_bf.flt.fa,\
	program=blastx,output=line6i_vs_line7i.miso_bf.flt.fa.gallus.xml" \
	$(protocol)/blast.sh

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
	qsub -v "db=Human_prot,input=line6i_vs_line7i.miso_bf.flt.fa,\
	program=blastx,output=line6i_vs_line7i.miso_bf.flt.fa.human.xml" \
	$(protocol)/blast.sh

translate-isoforms-miso-se:

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		estscan -t line6i_vs_line7i.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line6i_vs_line7i.miso_bf.flt.fa > line6i_vs_line7i.miso_bf.flt.fna

	cd miso/results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		estscan -t line6u_vs_line7u.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line6u_vs_line7u.miso_bf.flt.fa > line6u_vs_line7u.miso_bf.flt.fna

	cd miso/results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		estscan -t line6u_vs_line6i.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line6u_vs_line6i.miso_bf.flt.fa > line6u_vs_line6i.miso_bf.flt.fna

	cd miso/results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		estscan -t line7u_vs_line7i.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line7u_vs_line7i.miso_bf.flt.fa > line7u_vs_line7i.miso_bf.flt.fna

translate-isoforms-miso-a3ss:

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		estscan -t line6i_vs_line7i.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line6i_vs_line7i.miso_bf.flt.fa > line6i_vs_line7i.miso_bf.flt.fna

	cd miso/results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		estscan -t line6u_vs_line7u.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line6u_vs_line7u.miso_bf.flt.fa > line6u_vs_line7u.miso_bf.flt.fna

	cd miso/results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		estscan -t line6u_vs_line6i.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line6u_vs_line6i.miso_bf.flt.fa > line6u_vs_line6i.miso_bf.flt.fna

	cd miso/results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		estscan -t line7u_vs_line7i.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line7u_vs_line7i.miso_bf.flt.fa > line7u_vs_line7i.miso_bf.flt.fna

translate-isoforms-miso-a5ss:

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		estscan -t line6i_vs_line7i.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line6i_vs_line7i.miso_bf.flt.fa > line6i_vs_line7i.miso_bf.flt.fna

	cd miso/results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		estscan -t line6u_vs_line7u.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line6u_vs_line7u.miso_bf.flt.fa > line6u_vs_line7u.miso_bf.flt.fna

	cd miso/results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		estscan -t line6u_vs_line6i.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line6u_vs_line6i.miso_bf.flt.fa > line6u_vs_line6i.miso_bf.flt.fna

	cd miso/results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		estscan -t line7u_vs_line7i.miso_bf.flt.faa -M $(protocol)/gallus.hm \
		line7u_vs_line7i.miso_bf.flt.fa > line7u_vs_line7i.miso_bf.flt.fna

interpro-isoforms-miso-se:

	cd miso/results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		qsub -v input="line6u_vs_7u.miso_bf.flt.faa" $(protocol)/iprscan.sh

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		qsub -v input="line6i_vs_7i.miso_bf.flt.faa" $(protocol)/iprscan.sh

	cd miso/results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		qsub -v input="line6u_vs_7u.miso_bf.flt.faa" $(protocol)/iprscan.sh

	cd miso/results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		qsub -v input="line7u_vs_7i.miso_bf.flt.faa" $(protocol)/iprscan.sh

interpro-isoforms-miso-a3ss:

	cd miso/results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		qsub -v input="line6u_vs_7u.miso_bf.flt.faa" $(protocol)/iprscan.sh

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		qsub -v input="line6i_vs_7i.miso_bf.flt.faa" $(protocol)/iprscan.sh

	cd miso/results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		qsub -v input="line6u_vs_7u.miso_bf.flt.faa" $(protocol)/iprscan.sh

	cd miso/results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		qsub -v input="line7u_vs_7i.miso_bf.flt.faa" $(protocol)/iprscan.sh

interpro-isoforms-miso-a5ss:

	cd miso/results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		qsub -v input="line6u_vs_7u.miso_bf.flt.faa" $(protocol)/iprscan.sh

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		qsub -v input="line6i_vs_7i.miso_bf.flt.faa" $(protocol)/iprscan.sh

	cd miso/results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		qsub -v input="line6u_vs_7u.miso_bf.flt.faa" $(protocol)/iprscan.sh

	cd miso/results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		qsub -v input="line7u_vs_7i.miso_bf.flt.faa" $(protocol)/iprscan.sh

blat-domains-se:

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/interpro_to_bed.py line6i_vs_line7i.miso_bf.flt.faa.tsv > \
		line6i_vs_line7i.miso_bf.flt.faa.bed

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		seqtk subseq line6i_vs_line7i.miso_bf.flt.faa \
		line6i_vs_line7i.miso_bf.flt.faa.bed > line6i_vs_line7i.miso_bf.flt.faa.domains

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/split-fa.py line6i_vs_line7i.miso_bf.flt.faa.domains \
		1000 line6i_vs_line7i.miso_bf.flt.faa.domains; \

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in *domains_*.fa; do \
			qsub -v "input=$$f,index=../../../../../../gal4selected.2bit" \
			$(protocol)/blat_domains.sh; \
		done

blat-domains-a3ss:

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/interpro_to_bed.py line6i_vs_line7i.miso_bf.flt.faa.tsv > \
		line6i_vs_line7i.miso_bf.flt.faa.bed

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		seqtk subseq line6i_vs_line7i.miso_bf.flt.faa \
		line6i_vs_line7i.miso_bf.flt.faa.bed > line6i_vs_line7i.miso_bf.flt.faa.domains

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/split-fa.py line6i_vs_line7i.miso_bf.flt.faa.domains \
		1000 line6i_vs_line7i.miso_bf.flt.faa.domains; \

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in *domains_*.fa; do \
			qsub -v "input=$$f,index=../../../../../../gal4selected.2bit" \
			$(protocol)/blat_domains.sh; \
		done

blat-domains-a5ss:

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/interpro_to_bed.py line6i_vs_line7i.miso_bf.flt.faa.tsv > \
		line6i_vs_line7i.miso_bf.flt.faa.bed

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		seqtk subseq line6i_vs_line7i.miso_bf.flt.faa \
		line6i_vs_line7i.miso_bf.flt.faa.bed > line6i_vs_line7i.miso_bf.flt.faa.domains

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/split-fa.py line6i_vs_line7i.miso_bf.flt.faa.domains \
		1000 line6i_vs_line7i.miso_bf.flt.faa.domains; \

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in *domains_*.fa; do \
			qsub -v "input=$$f,index=../../../../../../gal4selected.2bit" \
			$(protocol)/blat_domains.sh; \
		done

annotate-domains-se:

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *domains*fa.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k 10 line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl > \
		line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl \
		line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl info

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/merge_bed_gff3_interpro.py line6i_vs_line7i.miso_bf.flt.faa.bed \
		line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl > \
		line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed

annotate-domains-a3ss:

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *domains*fa.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k 10 line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl > \
		line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl \
		line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl info

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/merge_bed_gff3_interpro.py line6i_vs_line7i.miso_bf.flt.faa.bed \
		line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl > \
		line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed

annotate-domains-a5ss:

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *domains*fa.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k 10 line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl > \
		line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl \
		line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl info

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/merge_bed_gff3_interpro.py line6i_vs_line7i.miso_bf.flt.faa.bed \
		line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl > \
		line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed

domains-to-bigBed-se:

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k1,1 -k2,2n line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed \
		| uniq > line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		bedToBigBed line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed \
		$(protocol)/gal4.chrom.sizes \
		line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bigBed.bb

domains-to-bigBed-a3ss:

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k1,1 -k2,2n line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed \
		| uniq > line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		bedToBigBed line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed \
		$(protocol)/gal4.chrom.sizes \
		line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bigBed.bb

domains-to-bigBed-a5ss:

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k1,1 -k2,2n line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed \
		| uniq > line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		bedToBigBed line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed \
		$(protocol)/gal4.chrom.sizes \
		line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bigBed.bb

miso-to-kegg-se:

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.gallus.xml > \
		line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt; \
		python $(protocol)/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.human.xml > \
		line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt

miso-to-kegg-a3ss:

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.gallus.xml > \
		line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt; \
		python $(protocol)/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.human.xml > \
		line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt

miso-to-kegg-a5ss:

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.gallus.xml > \
		line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt; \
		python $(protocol)/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.human.xml > \
		line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt

interpro-isoforms-ensbl-se:

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cut -f 2 line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt \
		| sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt \
		| sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/ensbl_genes_to_isoforms.py \
		../../../../../..//Gallus_gallus.Galgal4.73.pep.all.fa \
		gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python $(protocol)/ensbl_genes_to_isoforms.py \
		../../../../../../Homo_sapiens.GRCh37.74.pep.all.fa \
		human.ensbl.genes.list.txt > human.ensbl.genes.faa

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		qsub -v input="gallus.ensbl.genes.faa" $(protocol)/iprscan.sh
	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		qsub -v input="human.ensbl.genes.faa" $(protocol)/iprscan.sh

miso-snps-se:

	cd miso; grep exon gimme.flt.SE.gff > gimme.flt.SE.exons.gff
	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/miso-call-snps.py line6i_vs_line7i.miso_bf.flt \
		../../../../../../gal4selected.fa line6i_vs_line7i.miso_bf.flt.line6.snps \
		../../../../../bam-data/line6u.bam ../../../../../bam-data/line6i.bam

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/miso-call-snps.py line6i_vs_line7i.miso_bf.flt \
		../../../../../../gal4selected.fa line7i_vs_line7i.miso_bf.flt.line7.snps \
		../../../../../bam-data/line7u.bam ../../../../../bam-data/line7i.bam

miso-snps-a3ss:

	cd miso; grep exon gimme.flt.A3SS.gff > gimme.flt.A3SS.exons.gff
	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/miso-call-snps.py line6i_vs_line7i.miso_bf.flt \
		../../../../../../gal4selected.fa line6i_vs_line7i.miso_bf.flt.line6.snps \
		../../../../../bam-data/line6u.bam ../../../../../bam-data/line6i.bam

	cd miso/results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/miso-call-snps.py line6i_vs_line7i.miso_bf.flt \
		../../../../../../gal4selected.fa line7i_vs_line7i.miso_bf.flt.line7.snps \
		../../../../../bam-data/line7u.bam ../../../../../bam-data/line7i.bam

miso-snps-a5ss:

	cd miso; grep exon gimme.flt.A5SS.gff > gimme.flt.A3SS.exons.gff
	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/miso-call-snps.py line6i_vs_line7i.miso_bf.flt \
		../../../../../../gal4selected.fa line6i_vs_line7i.miso_bf.flt.line6.snps \
		../../../../../bam-data/line6u.bam ../../../../../bam-data/line6i.bam

	cd miso/results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/miso-call-snps.py line6i_vs_line7i.miso_bf.flt \
		../../../../../../gal4selected.fa line7i_vs_line7i.miso_bf.flt.line7.snps \
		../../../../../bam-data/line7u.bam ../../../../../bam-data/line7i.bam

find-deu-snps-se:

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/select_events.py ../../../../../gimme.flt.SE.gff \
		line6i_vs_line7i.miso_bf.flt SE | sort -k1,1 -k4,5n > line6i_vs_line7i.miso_bf.flt.SE.gff

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *line6.snps_* | sort -k1,1 -k2,2n | awk '$$6>19' >> line6i_vs_line7i.miso_bf.flt.line6.snps.vcf; \
		cat *line7.snps_* | sort -k1,1 -k2,2n | awk '$$6>19' >> line6i_vs_line7i.miso_bf.flt.line7.snps.vcf

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat $(protocol)/vcf.header.txt > line6i_vs_line7i.miso_bf.flt.line6.SE.snps.vcf; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.SE.gff \
		-a line6i_vs_line7i.miso_bf.flt.line6.snps.vcf -wb \
		| sort -k1,1 -k2,2n >> line6i_vs_line7i.miso_bf.flt.line6.SE.snps.vcf

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat $(protocol)/vcf.header.txt > line6i_vs_line7i.miso_bf.flt.line7.SE.snps.vcf; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.SE.gff \
		-a line6i_vs_line7i.miso_bf.flt.line7.snps.vcf -wb \
		| sort -k1,1 -k2,2n >> line6i_vs_line7i.miso_bf.flt.line7.SE.snps.vcf

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		GATK -T LeftAlignAndTrimVariants -R $(projpath)/gal4selected.sorted.fa \
		-trim -V:testsnp,VCF line6i_vs_line7i.miso_bf.flt.line7.SE.snps.vcf \
		-o line6i_vs_line7i.miso_bf.flt.line7.SE.trimmed.snps.vcf; \
		GATK -T LeftAlignAndTrimVariants -R $(projpath)/gal4selected.sorted.fa \
		-trim -V:testsnp,VCF line6i_vs_line7i.miso_bf.flt.line6.SE.snps.vcf \
		-o line6i_vs_line7i.miso_bf.flt.line6.SE.trimmed.snps.vcf

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/get_variant_seq.py line6i_vs_line7i.miso_bf.flt.SE.gff \
		line6i_vs_line7i.miso_bf.flt.line6.SE.trimmed.snps.vcf line6.variant \
		$(projpath)/gal4selected.sorted.fa

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/get_variant_seq.py line6i_vs_line7i.miso_bf.flt.SE.gff \
		line6i_vs_line7i.miso_bf.flt.line7.SE.trimmed.snps.vcf line7.variant \
		$(projpath)/gal4selected.sorted.fa

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/check_strand.py line6i_vs_line7i.miso_bf.flt.SE.gff line7.variant; \
		python $(protocol)/check_strand.py line6i_vs_line7i.miso_bf.flt.SE.gff line6.variant

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/find_diff_snps.py \
		line6i_vs_line7i.miso_bf.flt.line6.SE.trimmed.snps.vcf \
		line6i_vs_line7i.miso_bf.flt.line7.SE.snps.vcf 20 | \
		sort -k1,1 -k2,2n  > line6i_vs_line7i.miso_bf.flt.line6.SE.diff.snps.vcf

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python $(protocol)/find_diff_snps.py \
		line6i_vs_line7i.miso_bf.flt.line7.SE.trimmed.snps.vcf \
		line6i_vs_line7i.miso_bf.flt.line6.SE.snps.vcf 20 | \
		sort -k1,1 -k2,2n  > line6i_vs_line7i.miso_bf.flt.line7.SE.diff.snps.vcf

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.SE.gff \
		-a line6i_vs_line7i.miso_bf.flt.line6.SE.diff.snps.vcf \
		-wb > line6i_vs_line7i.miso_bf.flt.line6.SE.diff.snps.exons.vcf

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.SE.gff \
		-a line6i_vs_line7i.miso_bf.flt.line7.SE.diff.snps.vcf \
		-wb > line6i_vs_line7i.miso_bf.flt.line7.SE.diff.snps.exons.vcf

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in line6.variant*.fa; do \
			blat $$f $$(echo $$f | sed s/line6/line7/) -out=blast $$f.blast; \
		done

	cd miso/results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in line6.variant*.fa.blast; do \
			echo "working on " $$f; \
			python $(protocol)/blast_to_snps.py $$f > $$f.snps; \
		done

convert-snps-to-csv:

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat line6i_vs_line7i.miso_bf.flt.line6.A5SS.diff.snps.exons.vcf | python ~/mdv-protocol/to_csv.py > line6i_vs_line7i.miso_bf.flt.line6.A5SS.diff.snps.exons.csv

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat line6i_vs_line7i.miso_bf.flt.line7.A5SS.diff.snps.exons.vcf | python ~/mdv-protocol/to_csv.py > line6i_vs_line7i.miso_bf.flt.line7.A5SS.diff.snps.exons.csv

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat line6i_vs_line7i.miso_bf.flt.line6.A3SS.diff.snps.exons.vcf | python ~/mdv-protocol/to_csv.py > line6i_vs_line7i.miso_bf.flt.line6.A3SS.diff.snps.exons.csv

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat line6i_vs_line7i.miso_bf.flt.line7.A3SS.diff.snps.exons.vcf | python ~/mdv-protocol/to_csv.py > line6i_vs_line7i.miso_bf.flt.line7.A3SS.diff.snps.exons.csv

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat line6i_vs_line7i.miso_bf.flt.line6.SE.diff.snps.exons.vcf | python ~/mdv-protocol/to_csv.py > line6i_vs_line7i.miso_bf.flt.line6.SE.diff.snps.exons.csv

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat line6i_vs_line7i.miso_bf.flt.line7.SE.diff.snps.exons.vcf | python ~/mdv-protocol/to_csv.py > line6i_vs_line7i.miso_bf.flt.line7.SE.diff.snps.exons.csv
