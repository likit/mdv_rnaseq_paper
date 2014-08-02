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

	cd miso; python ~/miso-protocol/miso-filter.py results/A3SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf ASS 0.20 0.20 10 > results/A3SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt
	cd miso; python ~/miso-protocol/miso-filter.py results/A3SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf ASS 0.20 0.20 10 > results/A3SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt

	cd miso; python ~/miso-protocol/miso-filter.py results/A3SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf ASS 0.20 0.20 10 > results/A3SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt
	cd miso; python ~/miso-protocol/miso-filter.py results/A3SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf ASS 0.20 0.20 10 > results/A3SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt

filter-miso-a5ss:

	cd miso; python ~/miso-protocol/miso-filter.py results/A5SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf ASS 0.20 0.20 10 > results/A5SS/comparisons/line6u_vs_line7u/bayes-factors/line6u_vs_line7u.miso_bf.flt
	cd miso; python ~/miso-protocol/miso-filter.py results/A5SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf ASS 0.20 0.20 10 > results/A5SS/comparisons/line6i_vs_line7i/bayes-factors/line6i_vs_line7i.miso_bf.flt

	cd miso; python ~/miso-protocol/miso-filter.py results/A5SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf ASS 0.20 0.20 10 > results/A5SS/comparisons/line6u_vs_line6i/bayes-factors/line6u_vs_line6i.miso_bf.flt
	cd miso; python ~/miso-protocol/miso-filter.py results/A5SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf ASS 0.20 0.20 10 > results/A5SS/comparisons/line7u_vs_line7i/bayes-factors/line7u_vs_line7i.miso_bf.flt

blast-miso-se:

	#cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
	#qsub -v db="nr",input="line6i_vs_line7i.miso_bf.flt.fa",program="blastx",output="line6i_vs_line7i.miso_bf.flt.fa.nr.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
	#qsub -v db="nr",input="line6u_vs_line7u.miso_bf.flt.fa",program="blastx",output="line6u_vs_line7u.miso_bf.flt.fa.nr.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
	#qsub -v db="Gallus_prot",input="line6i_vs_line7i.miso_bf.flt.fa",program="blastx",output="line6i_vs_line7i.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	#cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
	#qsub -v db="Gallus_prot",input="line6u_vs_line7u.miso_bf.flt.fa",program="blastx",output="line6u_vs_line7u.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	#cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
	#qsub -v db="Human_prot",input="line6i_vs_line7i.miso_bf.flt.fa",program="blastx",output="line6i_vs_line7i.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
	#qsub -v db="Human_prot",input="line6u_vs_line7u.miso_bf.flt.fa",program="blastx",output="line6u_vs_line7u.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
	qsub -v db="Gallus_prot",input="line7u_vs_line7i.miso_bf.flt.fa",program="blastx",output="line7u_vs_line7i.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
	qsub -v db="Gallus_prot",input="line6u_vs_line6i.miso_bf.flt.fa",program="blastx",output="line6u_vs_line6i.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
	qsub -v db="Human_prot",input="line7u_vs_line7i.miso_bf.flt.fa",program="blastx",output="line7u_vs_line7i.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
	qsub -v db="Human_prot",input="line6u_vs_line6i.miso_bf.flt.fa",program="blastx",output="line6u_vs_line6i.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

blast-miso-a3ss:

	#cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
	#qsub -v db="nr",input="line6i_vs_line7i.miso_bf.flt.fa",program="blastx",output="line6i_vs_line7i.miso_bf.flt.fa.nr.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
	#qsub -v db="nr",input="line6u_vs_line7u.miso_bf.flt.fa",program="blastx",output="line6u_vs_line7u.miso_bf.flt.fa.nr.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
	#qsub -v db="Gallus_prot",input="line6i_vs_line7i.miso_bf.flt.fa",program="blastx",output="line6i_vs_line7i.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	#cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
	#qsub -v db="Gallus_prot",input="line6u_vs_line7u.miso_bf.flt.fa",program="blastx",output="line6u_vs_line7u.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	#cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
	#qsub -v db="Human_prot",input="line6i_vs_line7i.miso_bf.flt.fa",program="blastx",output="line6i_vs_line7i.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
	#qsub -v db="Human_prot",input="line6u_vs_line7u.miso_bf.flt.fa",program="blastx",output="line6u_vs_line7u.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
	#qsub -v db="Gallus_prot",input="line7u_vs_line7i.miso_bf.flt.fa",program="blastx",output="line7u_vs_line7i.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
	qsub -v db="Gallus_prot",input="line6u_vs_line6i.miso_bf.flt.fa",program="blastx",output="line6u_vs_line6i.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	#cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
	#qsub -v db="Human_prot",input="line7u_vs_line7i.miso_bf.flt.fa",program="blastx",output="line7u_vs_line7i.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
	#qsub -v db="Human_prot",input="line6u_vs_line6i.miso_bf.flt.fa",program="blastx",output="line6u_vs_line6i.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

blast-miso-a5ss:

	#cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
	#qsub -v db="nr",input="line6i_vs_line7i.miso_bf.flt.fa",program="blastx",output="line6i_vs_line7i.miso_bf.flt.fa.nr.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
	#qsub -v db="nr",input="line6u_vs_line7u.miso_bf.flt.fa",program="blastx",output="line6u_vs_line7u.miso_bf.flt.fa.nr.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
	#qsub -v db="Gallus_prot",input="line6i_vs_line7i.miso_bf.flt.fa",program="blastx",output="line6i_vs_line7i.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	#cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
	#qsub -v db="Gallus_prot",input="line6u_vs_line7u.miso_bf.flt.fa",program="blastx",output="line6u_vs_line7u.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	#cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
	#qsub -v db="Human_prot",input="line6i_vs_line7i.miso_bf.flt.fa",program="blastx",output="line6i_vs_line7i.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
	#qsub -v db="Human_prot",input="line6u_vs_line7u.miso_bf.flt.fa",program="blastx",output="line6u_vs_line7u.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
	#qsub -v db="Gallus_prot",input="line7u_vs_line7i.miso_bf.flt.fa",program="blastx",output="line7u_vs_line7i.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
	qsub -v db="Gallus_prot",input="line6u_vs_line6i.miso_bf.flt.fa",program="blastx",output="line6u_vs_line6i.miso_bf.flt.fa.gallus.xml" ~/miso-protocol/blast.sh
	#cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
	#qsub -v db="Human_prot",input="line7u_vs_line7i.miso_bf.flt.fa",program="blastx",output="line7u_vs_line7i.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh

	#cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
	#qsub -v db="Human_prot",input="line6u_vs_line6i.miso_bf.flt.fa",program="blastx",output="line6u_vs_line6i.miso_bf.flt.fa.human.xml" ~/miso-protocol/blast.sh
