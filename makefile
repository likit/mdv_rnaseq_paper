###############################################
#####     Prepare reference genome       ######
###############################################

init:

	java -jar /opt/software/picardTools/1.113/CreateSequenceDictionary.jar R=gal4selected.sorted.fa O=gal4selected.sorted.dict
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

rsem-calc-gimme-models-ref-rspd:

	qsub -v input_read="reads/line6u.se.fq",sample_name="line6u-single-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line6i.se.fq",sample_name="line6i-single-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line6u.pe.1",input_read2="reads/line6u.pe.2",sample_name="line6u-paired-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line6i.pe.1",input_read2="reads/line6i.pe.2",sample_name="line6i-paired-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh

	qsub -v input_read="reads/line7u.se.fq",sample_name="line7u-single-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7i.se.fq",sample_name="line7i-single-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" \
		~/mdv-protocol/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line7u.pe.1",input_read2="reads/line7u.pe.2",sample_name="line7u-paired-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7i.pe.1",input_read2="reads/line7i.pe.2",sample_name="line7i-paired-rsem-cuffref-rspd",index="asm_cuff_ref_models_rsem" ~/mdv-protocol/rsem_calculate_expr_paired.sh

ebseq-line6-models-ref:

	rsem-generate-data-matrix line6u-single-rsem-cuffref.genes.results \
		line6u-paired-rsem-cuffref.genes.results line6i-single-rsem-cuffref.genes.results \
		line6i-paired-rsem-cuffref.genes.results > line6u_vs_i.gene.cuffref.counts.matrix
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

get-longest-sequences-cuffref:

	python ~/mdv-protocol/gene-rep.py line6u_vs_i.cuffref.degenes.fdr.05.fa > line6u_vs_i.cuffref.degenes.fdr.05.fa.longest
	python ~/mdv-protocol/gene-rep.py line7u_vs_i.cuffref.degenes.fdr.05.fa > line7u_vs_i.cuffref.degenes.fdr.05.fa.longest
	python ~/mdv-protocol/gene-rep.py line6u_vs_i.cuffref.degenes.fdr.05.fa.prot > line6u_vs_i.cuffref.degenes.fdr.05.fa.prot.longest
	python ~/mdv-protocol/gene-rep.py line7u_vs_i.cuffref.degenes.fdr.05.fa.prot > line7u_vs_i.cuffref.degenes.fdr.05.fa.prot.longest

run-blast-cuffref-gallus:

	qsub -v db="Gallus_prot",input="line6u_vs_i.cuffref.degenes.fdr.05.fa.longest",program="blastx",output="line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.gallus.xml" ~/mdv-protocol/blast.sh
	qsub -v db="Gallus_prot",input="line7u_vs_i.cuffref.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.gallus.xml" ~/mdv-protocol/blast.sh

run-blast-cuffref-human:

	qsub -v db="Human_prot",input="line6u_vs_i.cuffref.degenes.fdr.05.fa.longest",program="blastx",output="line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.human.xml" ~/mdv-protocol/blast.sh
	qsub -v db="Human_prot",input="line7u_vs_i.cuffref.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.human.xml" ~/mdv-protocol/blast.sh

get-tophits:

	python protocol/get_top_hits.py line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.gallus.xml \
		> line6u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.txt

	python protocol/get_top_hits.py line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.longest.gallus.xml \
		> line7u_vs_i.cuffref.degenes.fdr.05.fa.nucl.tophits.txt

miso-to-fa-se:

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line6i_vs_line7i.miso_bf.flt \
		$(path)/asm_cuff_ref_models.bed.fa > line6i_vs_line7i.miso_bf.flt.fa

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line6u_vs_line7u.miso_bf.flt \
		$(path)/asm_cuff_ref_models.bed.fa > line6u_vs_line7u.miso_bf.flt.fa

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line7u_vs_line7i.miso_bf.flt \
		$(path)/asm_cuff_ref_models.bed.fa > line7u_vs_line7i.miso_bf.flt.fa

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line6u_vs_line6i.miso_bf.flt \
		$(path)/asm_cuff_ref_models.bed.fa > line6u_vs_line6i.miso_bf.flt.fa


miso-to-fa-a3ss:

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line6i_vs_line7i.miso_bf.flt $(path)/asm_cuff_ref_models.bed.fa > line6i_vs_line7i.miso_bf.flt.fa

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line6u_vs_line7u.miso_bf.flt $(path)/asm_cuff_ref_models.bed.fa > line6u_vs_line7u.miso_bf.flt.fa

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line7u_vs_line7i.miso_bf.flt $(path)/asm_cuff_ref_models.bed.fa > line7u_vs_line7i.miso_bf.flt.fa

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line6u_vs_line6i.miso_bf.flt $(path)/asm_cuff_ref_models.bed.fa > line6u_vs_line6i.miso_bf.flt.fa

translate-isoforms-miso-se:

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		estscan -t line6i_vs_line7i.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line6i_vs_line7i.miso_bf.flt.fa > line6i_vs_line7i.miso_bf.flt.fna

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		estscan -t line6u_vs_line7u.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line6u_vs_line7u.miso_bf.flt.fa > line6u_vs_line7u.miso_bf.flt.fna

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		estscan -t line7u_vs_line7i.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line7u_vs_line7i.miso_bf.flt.fa > line7u_vs_line7i.miso_bf.flt.fna

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		estscan -t line6u_vs_line6i.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line6u_vs_line6i.miso_bf.flt.fa > line6u_vs_line6i.miso_bf.flt.fna

interpro-isoforms-miso-se:

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; qsub -v input="line6u_vs_7u.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh
	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; qsub -v input="line6i_vs_7i.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh
	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; qsub -v input="line6u_vs_7u.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh
	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; qsub -v input="line7u_vs_7i.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh

annotate-domains-se:

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line6u_vs_line7u.miso_bf.flt.faa.tsv > line6u_vs_line7u.miso_bf.flt.faa.bed

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		seqtk subseq line6u_vs_line7u.miso_bf.flt.faa line6u_vs_line7u.miso_bf.flt.faa.bed > line6u_vs_line7u.miso_bf.flt.faa.domains

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		for f in *domains_*.fa; do \
			qsub -v input=$$f ~/mdv-protocol/blat_domains.sh; \
		done

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		cat *domains*fa.psl > line6u_vs_line7u.miso_bf.flt.faa.domains.all.psl

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		sort -k 10 line6u_vs_line7u.miso_bf.flt.faa.domains.all.psl > line6u_vs_line7u.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6u_vs_line7u.miso_bf.flt.faa.domains.all.sorted.psl line6u_vs_line7u.miso_bf.flt.faa.domains.all.best.psl info;  \

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line6u_vs_line7u.miso_bf.flt.faa.bed line6u_vs_line7u.miso_bf.flt.faa.domains.all.best.psl > line6u_vs_line7u.miso_bf.flt.faa.domains.annots.bed

	#------------------------------
	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line6i_vs_line7i.miso_bf.flt.faa.tsv > line6i_vs_line7i.miso_bf.flt.faa.bed

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		seqtk subseq line6i_vs_line7i.miso_bf.flt.faa line6i_vs_line7i.miso_bf.flt.faa.bed > line6i_vs_line7i.miso_bf.flt.faa.domains

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/split-fa.py line6i_vs_line7i.miso_bf.flt.faa.domains 1000 line6i_vs_line7i.miso_bf.flt.faa.domains; \
		for f in *domains_*.fa; do \
			qsub -v input=$$f ~/mdv-protocol/blat_domains.sh; \
		done

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *domains*fa.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k 10 line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl info;  \

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line6i_vs_line7i.miso_bf.flt.faa.bed line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed

	#-------------------------------

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line7u_vs_line7i.miso_bf.flt.faa.tsv > line7u_vs_line7i.miso_bf.flt.faa.bed

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		seqtk subseq line7u_vs_line7i.miso_bf.flt.faa line7u_vs_line7i.miso_bf.flt.faa.bed > line7u_vs_line7i.miso_bf.flt.faa.domains

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/split-fa.py line7u_vs_line7i.miso_bf.flt.faa.domains 1000 line7u_vs_line7i.miso_bf.flt.faa.domains; \
		for f in *domains_*.fa; do \
			qsub -v input=$$f ~/mdv-protocol/blat_domains.sh; \
		done

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		cat *domains*fa.psl > line7u_vs_line7i.miso_bf.flt.faa.domains.all.psl

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		sort -k 10 line7u_vs_line7i.miso_bf.flt.faa.domains.all.psl > line7u_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line7u_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl line7u_vs_line7i.miso_bf.flt.faa.domains.all.best.psl info;  \

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line7u_vs_line7i.miso_bf.flt.faa.bed line7u_vs_line7i.miso_bf.flt.faa.domains.all.best.psl > line7u_vs_line7i.miso_bf.flt.faa.domains.annots.bed

	#-------------------------------

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line6u_vs_line6i.miso_bf.flt.faa.tsv > line6u_vs_line6i.miso_bf.flt.faa.bed

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		seqtk subseq line6u_vs_line6i.miso_bf.flt.faa line6u_vs_line6i.miso_bf.flt.faa.bed > line6u_vs_line6i.miso_bf.flt.faa.domains

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/split-fa.py line6u_vs_line6i.miso_bf.flt.faa.domains 1000 line6u_vs_line6i.miso_bf.flt.faa.domains; \
		for f in *domains_*.fa; do \
			qsub -v input=$$f ~/mdv-protocol/blat_domains.sh; \
		done

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		cat *domains*fa.psl > line6u_vs_line6i.miso_bf.flt.faa.domains.all.psl

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		sort -k 10 line6u_vs_line6i.miso_bf.flt.faa.domains.all.psl > line6u_vs_line6i.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6u_vs_line6i.miso_bf.flt.faa.domains.all.sorted.psl line6u_vs_line6i.miso_bf.flt.faa.domains.all.best.psl info;  \

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line6u_vs_line6i.miso_bf.flt.faa.bed line6u_vs_line6i.miso_bf.flt.faa.domains.all.best.psl > line6u_vs_line6i.miso_bf.flt.faa.domains.annots.bed

domains-to-bigBed-se:

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		sort -k1,1 -k2,2n line6u_vs_line7u.miso_bf.flt.faa.domains.annots.bed | uniq > line6u_vs_line7u.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line6u_vs_line7u.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line6u_vs_line7u.miso_bf.flt.faa.domains.annots.bigBed.bb

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		sort -k1,1 -k2,2n line6u_vs_line6i.miso_bf.flt.faa.domains.annots.bed | uniq > line6u_vs_line6i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line6u_vs_line6i.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line6u_vs_line6i.miso_bf.flt.faa.domains.annots.bigBed.bb

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		sort -k1,1 -k2,2n line7u_vs_line7i.miso_bf.flt.faa.domains.annots.bed | uniq > line7u_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line7u_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line7u_vs_line7i.miso_bf.flt.faa.domains.annots.bigBed.bb

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k1,1 -k2,2n line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed | uniq > line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bigBed.bb

domains-to-bigBed-a3ss:

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		sort -k1,1 -k2,2n line6u_vs_line7u.miso_bf.flt.faa.domains.annots.bed | uniq > line6u_vs_line7u.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line6u_vs_line7u.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line6u_vs_line7u.miso_bf.flt.faa.domains.annots.bigBed.bb

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		sort -k1,1 -k2,2n line6u_vs_line6i.miso_bf.flt.faa.domains.annots.bed | uniq > line6u_vs_line6i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line6u_vs_line6i.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line6u_vs_line6i.miso_bf.flt.faa.domains.annots.bigBed.bb

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		sort -k1,1 -k2,2n line7u_vs_line7i.miso_bf.flt.faa.domains.annots.bed | uniq > line7u_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line7u_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line7u_vs_line7i.miso_bf.flt.faa.domains.annots.bigBed.bb

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k1,1 -k2,2n line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed | uniq > line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bigBed.bb

domains-to-bigBed-a5ss:

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		sort -k1,1 -k2,2n line6u_vs_line7u.miso_bf.flt.faa.domains.annots.bed | uniq > line6u_vs_line7u.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line6u_vs_line7u.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line6u_vs_line7u.miso_bf.flt.faa.domains.annots.bigBed.bb

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		sort -k1,1 -k2,2n line6u_vs_line6i.miso_bf.flt.faa.domains.annots.bed | uniq > line6u_vs_line6i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line6u_vs_line6i.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line6u_vs_line6i.miso_bf.flt.faa.domains.annots.bigBed.bb

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		sort -k1,1 -k2,2n line7u_vs_line7i.miso_bf.flt.faa.domains.annots.bed | uniq > line7u_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line7u_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line7u_vs_line7i.miso_bf.flt.faa.domains.annots.bigBed.bb

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k1,1 -k2,2n line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed | uniq > line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed; \
		~/bedToBigBed line6i_vs_line7i.miso_bf.flt.faa.domains.annots.sorted.bed ~/mdv-protocol/gal4.chrom.sizes line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bigBed.bb

miso-to-kegg-se:

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line7u.miso_bf.flt.fa.gallus.xml > line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line7u.miso_bf.flt.fa.human.xml > line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.gallus.xml > line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.human.xml > line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line6i.miso_bf.flt.fa.gallus.xml > line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line6i.miso_bf.flt.fa.human.xml > line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line7u_vs_line7i.miso_bf.flt.fa.gallus.xml > line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line7u_vs_line7i.miso_bf.flt.fa.human.xml > line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt

miso-to-kegg-a3ss:

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line7u.miso_bf.flt.fa.gallus.xml > line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line7u.miso_bf.flt.fa.human.xml > line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line6i.miso_bf.flt.fa.gallus.xml > line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line6i.miso_bf.flt.fa.human.xml > line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line7u_vs_line7i.miso_bf.flt.fa.gallus.xml > line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line7u_vs_line7i.miso_bf.flt.fa.human.xml > line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.gallus.xml > line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.human.xml > line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt

miso-to-kegg-a5ss:

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line7u.miso_bf.flt.fa.gallus.xml > line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line7u.miso_bf.flt.fa.human.xml > line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line6i.miso_bf.flt.fa.gallus.xml > line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line6u_vs_line6i.miso_bf.flt.fa.human.xml > line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line7u_vs_line7i.miso_bf.flt.fa.gallus.xml > line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line7u_vs_line7i.miso_bf.flt.fa.human.xml > line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.gallus.xml > line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt; \
		python ~/mdv-protocol/get_top_hits.py line6i_vs_line7i.miso_bf.flt.fa.human.xml > line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt

translate-isoforms-miso-a3ss:

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		estscan -t line6i_vs_line7i.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line6i_vs_line7i.miso_bf.flt.fa > line6i_vs_line7i.miso_bf.flt.fna

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		estscan -t line6u_vs_line7u.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line6u_vs_line7u.miso_bf.flt.fa > line6u_vs_line7u.miso_bf.flt.fna

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		estscan -t line6u_vs_line6i.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line6u_vs_line6i.miso_bf.flt.fa > line6u_vs_line6i.miso_bf.flt.fna

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		estscan -t line7u_vs_line7i.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line7u_vs_line7i.miso_bf.flt.fa > line7u_vs_line7i.miso_bf.flt.fna

interpro-isoforms-miso-a3ss:

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; qsub -v input="line6u_vs_7u.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh
	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; qsub -v input="line6i_vs_7i.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh

	d miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; qsub -v input="line6u_vs_6i.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh
	d miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; qsub -v input="line7u_vs_7i.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh

anotate-domains-a3ss:

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line6u_vs_line7u.miso_bf.flt.faa.tsv > line6u_vs_line7u.miso_bf.flt.faa.bed

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		seqtk subseq line6u_vs_line7u.miso_bf.flt.faa line6u_vs_line7u.miso_bf.flt.faa.bed > line6u_vs_line7u.miso_bf.flt.faa.domains

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/split-fa.py line6u_vs_line7u.miso_bf.flt.faa.domains 400 line6u_vs_line7u.miso_bf.flt.faa.domains; \
		for f in *domains_*.fa; do \
			qsub -v input=$$f ~/mdv-protocol/blat_domains.sh; \
		done

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		cat *domains*fa.psl > line6u_vs_line7u.miso_bf.flt.faa.domains.all.psl

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		sort -k 10 line6u_vs_line7u.miso_bf.flt.faa.domains.all.psl > line6u_vs_line7u.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6u_vs_line7u.miso_bf.flt.faa.domains.all.sorted.psl line6u_vs_line7u.miso_bf.flt.faa.domains.all.best.psl info;  \

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line6u_vs_line7u.miso_bf.flt.faa.bed line6u_vs_line7u.miso_bf.flt.faa.domains.all.best.psl > line6u_vs_line7u.miso_bf.flt.faa.domains.annots.bed

	#----------------------------------

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line6i_vs_line7i.miso_bf.flt.faa.tsv > line6i_vs_line7i.miso_bf.flt.faa.bed

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		seqtk subseq line6i_vs_line7i.miso_bf.flt.faa line6i_vs_line7i.miso_bf.flt.faa.bed > line6i_vs_line7i.miso_bf.flt.faa.domains

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/split-fa.py line6i_vs_line7i.miso_bf.flt.faa.domains 1000 line6i_vs_line7i.miso_bf.flt.faa.domains; \
		for f in *domains_*.fa; do \
			qsub -v input=$$f ~/mdv-protocol/blat_domains.sh; \
		done

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *domains*fa.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k 10 line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl info;  \

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line6i_vs_line7i.miso_bf.flt.faa.bed line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed

	#----------------------------------

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line7u_vs_line7i.miso_bf.flt.faa.tsv > line7u_vs_line7i.miso_bf.flt.faa.bed

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		seqtk subseq line7u_vs_line7i.miso_bf.flt.faa line7u_vs_line7i.miso_bf.flt.faa.bed > line7u_vs_line7i.miso_bf.flt.faa.domains

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		qsub -v input="line7u_vs_line7i.miso_bf.flt.faa.domains" ~/mdv-protocol/blat_domains.sh

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		sort -k 10 line7u_vs_line7i.miso_bf.flt.faa.domains.psl > line7u_vs_line7i.miso_bf.flt.faa.domains.sorted.psl; \
		pslReps -nohead -singleHit line7u_vs_line7i.miso_bf.flt.faa.domains.sorted.psl line7u_vs_line7i.miso_bf.flt.faa.domains.best.psl info;  \

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line7u_vs_line7i.miso_bf.flt.faa.bed line7u_vs_line7i.miso_bf.flt.faa.domains.best.psl > line7u_vs_line7i.miso_bf.flt.faa.domains.annots.bed

	#----------------------------------

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line6u_vs_line6i.miso_bf.flt.faa.tsv > line6u_vs_line6i.miso_bf.flt.faa.bed

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		seqtk subseq line6u_vs_line6i.miso_bf.flt.faa line6u_vs_line6i.miso_bf.flt.faa.bed > line6u_vs_line6i.miso_bf.flt.faa.domains

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		qsub -v input="line6u_vs_line6i.miso_bf.flt.faa.domains" ~/mdv-protocol/blat_domains.sh

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		sort -k 10 line6u_vs_line6i.miso_bf.flt.faa.domains.psl > line6u_vs_line6i.miso_bf.flt.faa.domains.sorted.psl; \
		pslReps -nohead -singleHit line6u_vs_line6i.miso_bf.flt.faa.domains.sorted.psl line6u_vs_line6i.miso_bf.flt.faa.domains.best.psl info;  \

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line6u_vs_line6i.miso_bf.flt.faa.bed line6u_vs_line6i.miso_bf.flt.faa.domains.best.psl > line6u_vs_line6i.miso_bf.flt.faa.domains.annots.bed

miso-to-fa-a5ss:

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line6i_vs_line7i.miso_bf.flt $(path)/asm_cuff_ref_models.bed.fa > line6i_vs_line7i.miso_bf.flt.fa

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line6u_vs_line7u.miso_bf.flt $(path)/asm_cuff_ref_models.bed.fa > line6u_vs_line7u.miso_bf.flt.fa

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line7u_vs_line7i.miso_bf.flt $(path)/asm_cuff_ref_models.bed.fa > line7u_vs_line7i.miso_bf.flt.fa

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/miso-protocol/miso-to-fa.py line6u_vs_line6i.miso_bf.flt $(path)/asm_cuff_ref_models.bed.fa > line6u_vs_line6i.miso_bf.flt.fa

translate-isoforms-miso-a5ss:

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		estscan -t line6i_vs_line7i.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line6i_vs_line7i.miso_bf.flt.fa > line6i_vs_line7i.miso_bf.flt.fna

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		estscan -t line6u_vs_line7u.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line6u_vs_line7u.miso_bf.flt.fa > line6u_vs_line7u.miso_bf.flt.fna

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		estscan -t line6u_vs_line6i.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line6u_vs_line6i.miso_bf.flt.fa > line6u_vs_line6i.miso_bf.flt.fna

	cd /mnt/ls12/preeyanon/mdv-pipeline/miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		estscan -t line7u_vs_line7i.miso_bf.flt.faa -M ~/mdv-protocol/gallus.hm line7u_vs_line7i.miso_bf.flt.fa > line7u_vs_line7i.miso_bf.flt.fna

interpro-isoforms-miso-a5ss:

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; qsub -v input="line6u_vs_7u.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh
	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; qsub -v input="line6i_vs_7i.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; qsub -v input="line6u_vs_6i.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh
	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; qsub -v input="line7u_vs_7i.miso_bf.flt.faa" ~/mdv-protocol/iprscan.sh

annotate-domains-a5ss:

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line6u_vs_line7u.miso_bf.flt.faa.tsv > line6u_vs_line7u.miso_bf.flt.faa.bed

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		seqtk subseq line6u_vs_line7u.miso_bf.flt.faa line6u_vs_line7u.miso_bf.flt.faa.bed > line6u_vs_line7u.miso_bf.flt.faa.domains

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/split-fa.py line6u_vs_line7u.miso_bf.flt.faa.domains 400 line6u_vs_line7u.miso_bf.flt.faa.domains; \
		for f in *domains_*.fa; do \
			qsub -v input=$$f ~/mdv-protocol/blat_domains.sh; \
		done

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		cat *domains*fa.psl > line6u_vs_line7u.miso_bf.flt.faa.domains.all.psl

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		sort -k 10 line6u_vs_line7u.miso_bf.flt.faa.domains.all.psl > line6u_vs_line7u.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6u_vs_line7u.miso_bf.flt.faa.domains.all.sorted.psl line6u_vs_line7u.miso_bf.flt.faa.domains.all.best.psl info;  \

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line6u_vs_line7u.miso_bf.flt.faa.bed line6u_vs_line7u.miso_bf.flt.faa.domains.all.best.psl > line6u_vs_line7u.miso_bf.flt.faa.domains.annots.bed

	#----------------------------------

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line6i_vs_line7i.miso_bf.flt.faa.tsv > line6i_vs_line7i.miso_bf.flt.faa.bed

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		seqtk subseq line6i_vs_line7i.miso_bf.flt.faa line6i_vs_line7i.miso_bf.flt.faa.bed > line6i_vs_line7i.miso_bf.flt.faa.domains

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/split-fa.py line6i_vs_line7i.miso_bf.flt.faa.domains 400 line6i_vs_line7i.miso_bf.flt.faa.domains; \
		for f in *domains_*.fa; do \
			qsub -v input=$$f ~/mdv-protocol/blat_domains.sh; \
		done

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *domains*fa.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		sort -k 10 line6i_vs_line7i.miso_bf.flt.faa.domains.all.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl; \
		pslReps -nohead -singleHit line6i_vs_line7i.miso_bf.flt.faa.domains.all.sorted.psl line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl info;  \

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line6i_vs_line7i.miso_bf.flt.faa.bed line6i_vs_line7i.miso_bf.flt.faa.domains.all.best.psl > line6i_vs_line7i.miso_bf.flt.faa.domains.annots.bed

	#----------------------------------

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line7u_vs_line7i.miso_bf.flt.faa.tsv > line7u_vs_line7i.miso_bf.flt.faa.bed

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		seqtk subseq line7u_vs_line7i.miso_bf.flt.faa line7u_vs_line7i.miso_bf.flt.faa.bed > line7u_vs_line7i.miso_bf.flt.faa.domains

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		qsub -v input="line7u_vs_line7i.miso_bf.flt.faa.domains" ~/mdv-protocol/blat_domains.sh

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		sort -k 10 line7u_vs_line7i.miso_bf.flt.faa.domains.psl > line7u_vs_line7i.miso_bf.flt.faa.domains.sorted.psl; \
		pslReps -nohead -singleHit line7u_vs_line7i.miso_bf.flt.faa.domains.sorted.psl line7u_vs_line7i.miso_bf.flt.faa.domains.best.psl info;  \

	d miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
	 python ~/mdv-protocol/merge_bed_gff3_interpro.py line7u_vs_line7i.miso_bf.flt.faa.bed line7u_vs_line7i.miso_bf.flt.faa.domains.best.psl > line7u_vs_line7i.miso_bf.flt.faa.domains.annots.bed

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/interpro_to_bed.py line6u_vs_line6i.miso_bf.flt.faa.tsv > line6u_vs_line6i.miso_bf.flt.faa.bed

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		seqtk subseq line6u_vs_line6i.miso_bf.flt.faa line6u_vs_line6i.miso_bf.flt.faa.bed > line6u_vs_line6i.miso_bf.flt.faa.domains

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		qsub -v input="line6u_vs_line6i.miso_bf.flt.faa.domains" ~/mdv-protocol/blat_domains.sh

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		sort -k 10 line6u_vs_line6i.miso_bf.flt.faa.domains.psl > line6u_vs_line6i.miso_bf.flt.faa.domains.sorted.psl; \
		pslReps -nohead -singleHit line6u_vs_line6i.miso_bf.flt.faa.domains.sorted.psl line6u_vs_line6i.miso_bf.flt.faa.domains.best.psl info;  \

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/mdv-protocol/merge_bed_gff3_interpro.py line6u_vs_line6i.miso_bf.flt.faa.bed line6u_vs_line6i.miso_bf.flt.faa.domains.best.psl > line6u_vs_line6i.miso_bf.flt.faa.domains.annots.bed

interpro-isoforms-ensbl-se:

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		cut -f 2 line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py /mnt/ls12/preeyanon/mdv-pipeline/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py /mnt/ls12/preeyanon/mdv-pipeline/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cut -f 2 line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py /mnt/ls12/preeyanon/mdv-pipeline/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py /mnt/ls12/preeyanon/mdv-pipeline/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; qsub -v input="gallus.ensbl.genes.faa" ~/mdv-protocol/iprscan.sh
	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; qsub -v input="human.ensbl.genes.faa" ~/mdv-protocol/iprscan.sh

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; qsub -v input="gallus.ensbl.genes.faa" ~/mdv-protocol/iprscan.sh
	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; qsub -v input="human.ensbl.genes.faa" ~/mdv-protocol/iprscan.sh

interpro-isoforms-ensbl-a3ss:

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		cut -f 2 line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cut -f 2 line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		cut -f 2 line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		cut -f 2 line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

interpro-isoforms-ensbl-a5ss:

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		cut -f 2 line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cut -f 2 line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		cut -f 2 line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		cut -f 2 line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt | sort | uniq | grep -v NA | grep -v score > gallus.ensbl.genes.list.txt; \
		cut -f 2 line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt | sort | uniq | grep -v NA | grep -v score > human.ensbl.genes.list.txt; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Gallus_gallus.Galgal4.73.pep.all.fa gallus.ensbl.genes.list.txt > gallus.ensbl.genes.faa; \
		python ~/mdv-protocol/ensbl_genes_to_isoforms.py $(path)/Homo_sapiens.GRCh37.74.pep.all.fa human.ensbl.genes.list.txt > human.ensbl.genes.faa

miso-snps-se:

	cd miso; grep exon asm_cuff_ref_models.flt.SE.gff > asm_cuff_ref_models.flt.SE.exons.gff
	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line7u.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line7u.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/SE/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line7u.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line7u.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6i_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line6i_vs_line7i.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6i_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line6i_vs_line7i.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

	d miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
	 python ~/miso-protocol/miso-call-snps.py line6u_vs_line6i.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line6i.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	d miso/cuffref-results/SE/comparisons/line6u_vs_line6i/bayes-factors; \
	 python ~/miso-protocol/miso-call-snps.py line6u_vs_line6i.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line6i.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line7u_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line7u_vs_line7i.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/SE/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line7u_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line7u_vs_line7i.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

find-deu-snps-se:

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/select_events.py /mnt/ls12/preeyanon/mdv-pipeline/miso/asm_cuff_ref_models.flt.SE.gff \
		line6i_vs_line7i.miso_bf.flt SE | sort -k1,1 -k4,5n > line6i_vs_line7i.miso_bf.flt.SE.gff

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *line6.snps_* | sort -k1,1 -k2,2n | awk '$$6>19' >> line6i_vs_line7i.miso_bf.flt.line6.snps.vcf; \
		cat *line7.snps_* | sort -k1,1 -k2,2n | awk '$$6>19' >> line6i_vs_line7i.miso_bf.flt.line7.snps.vcf

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat ~/mdv-protocol/vcf.header.txt > line6i_vs_line7i.miso_bf.flt.line6.SE.snps.vcf; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.SE.gff -a line6i_vs_line7i.miso_bf.flt.line6.snps.vcf -wb \
		| sort -k1,1 -k2,2n >> line6i_vs_line7i.miso_bf.flt.line6.SE.snps.vcf

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		cat ~/mdv-protocol/vcf.header.txt > line6i_vs_line7i.miso_bf.flt.line7.SE.snps.vcf; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.SE.gff -a line6i_vs_line7i.miso_bf.flt.line7.snps.vcf -wb \
		| sort -k1,1 -k2,2n >> line6i_vs_line7i.miso_bf.flt.line7.SE.snps.vcf

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		GATK -T LeftAlignAndTrimVariants -R /mnt/ls12/preeyanon/mdv-pipeline/gal4selected.sorted.fa -trim -V:testsnp,VCF line6i_vs_line7i.miso_bf.flt.line7.SE.snps.vcf -o line6i_vs_line7i.miso_bf.flt.line7.SE.trimmed.snps.vcf; \
		GATK -T LeftAlignAndTrimVariants -R /mnt/ls12/preeyanon/mdv-pipeline/gal4selected.sorted.fa -trim -V:testsnp,VCF line6i_vs_line7i.miso_bf.flt.line6.SE.snps.vcf -o line6i_vs_line7i.miso_bf.flt.line6.SE.trimmed.snps.vcf

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/get_variant_seq.py line6i_vs_line7i.miso_bf.flt.SE.gff line6i_vs_line7i.miso_bf.flt.line6.SE.trimmed.snps.vcf line6.variant; \
		python ~/mdv-protocol/get_variant_seq.py line6i_vs_line7i.miso_bf.flt.SE.gff line6i_vs_line7i.miso_bf.flt.line7.SE.trimmed.snps.vcf line7.variant

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/check_strand.py line6i_vs_line7i.miso_bf.flt.SE.gff line7.variant; \
		python ~/mdv-protocol/check_strand.py line6i_vs_line7i.miso_bf.flt.SE.gff line6.variant

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/find_diff_snps.py line6i_vs_line7i.miso_bf.flt.line6.SE.trimmed.snps.vcf line6i_vs_line7i.miso_bf.flt.line7.SE.snps.vcf 20 | sort -k1,1 -k2,2n  > line6i_vs_line7i.miso_bf.flt.line6.SE.diff.snps.vcf

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/find_diff_snps.py line6i_vs_line7i.miso_bf.flt.line7.SE.trimmed.snps.vcf line6i_vs_line7i.miso_bf.flt.line6.SE.snps.vcf 20 | sort -k1,1 -k2,2n  > line6i_vs_line7i.miso_bf.flt.line7.SE.diff.snps.vcf

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.SE.gff -a line6i_vs_line7i.miso_bf.flt.line6.SE.diff.snps.vcf -wb > line6i_vs_line7i.miso_bf.flt.line6.SE.diff.snps.exons.vcf

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.SE.gff -a line6i_vs_line7i.miso_bf.flt.line7.SE.diff.snps.vcf -wb > line6i_vs_line7i.miso_bf.flt.line7.SE.diff.snps.exons.vcf

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in line6.variant*.fa; do \
			blat $$f $$(echo $$f | sed s/line6/line7/) -out=blast $$f.blast; \
		done

	cd miso/cuffref-results/SE/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in line6.variant*.fa.blast; do \
			echo "working on " $$f; \
			python ~/mdv-protocol/blast_to_snps.py $$f > $$f.snps; \
		done

miso-snps-a3ss:

	cd miso; grep exon asm_cuff_ref_models.flt.A3SS.gff > asm_cuff_ref_models.flt.A3SS.exons.gff
	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line7u.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line7u.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line7u.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line7u.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6i_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line6i_vs_line7i.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6i_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line6i_vs_line7i.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line6i.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line6i.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/A3SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line6i.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line6i.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line7u_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line7u_vs_line7i.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/A3SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line7u_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line7u_vs_line7i.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

find-deu-snps-a3ss:

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/select_events.py /mnt/ls12/preeyanon/mdv-pipeline/miso/asm_cuff_ref_models.flt.A3SS.gff \
		line6i_vs_line7i.miso_bf.flt A3SS | sort -k1,1 -k4,5n > line6i_vs_line7i.miso_bf.flt.A3SS.gff

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *line6.snps_* | sort -k1,1 -k2,2n | awk '$$6>19' >> line6i_vs_line7i.miso_bf.flt.line6.snps.vcf; \
		cat *line7.snps_* | sort -k1,1 -k2,2n | awk '$$6>19' >> line6i_vs_line7i.miso_bf.flt.line7.snps.vcf

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat ~/mdv-protocol/vcf.header.txt > line6i_vs_line7i.miso_bf.flt.line6.A3SS.snps.vcf; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.A3SS.gff -a line6i_vs_line7i.miso_bf.flt.line6.snps.vcf -wb \
		| sort -k1,1 -k2,2n >> line6i_vs_line7i.miso_bf.flt.line6.A3SS.snps.vcf

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat ~/mdv-protocol/vcf.header.txt > line6i_vs_line7i.miso_bf.flt.line7.A3SS.snps.vcf; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.A3SS.gff -a line6i_vs_line7i.miso_bf.flt.line7.snps.vcf -wb \
		| sort -k1,1 -k2,2n >> line6i_vs_line7i.miso_bf.flt.line7.A3SS.snps.vcf

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		GATK -T LeftAlignAndTrimVariants -R /mnt/ls12/preeyanon/mdv-pipeline/gal4selected.sorted.fa -trim -V:testsnp,VCF line6i_vs_line7i.miso_bf.flt.line7.A3SS.snps.vcf -o line6i_vs_line7i.miso_bf.flt.line7.A3SS.trimmed.snps.vcf; \
		GATK -T LeftAlignAndTrimVariants -R /mnt/ls12/preeyanon/mdv-pipeline/gal4selected.sorted.fa -trim -V:testsnp,VCF line6i_vs_line7i.miso_bf.flt.line6.A3SS.snps.vcf -o line6i_vs_line7i.miso_bf.flt.line6.A3SS.trimmed.snps.vcf

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/get_variant_seq.py line6i_vs_line7i.miso_bf.flt.A3SS.gff line6i_vs_line7i.miso_bf.flt.line6.A3SS.trimmed.snps.vcf line6.variant; \
		python ~/mdv-protocol/get_variant_seq.py line6i_vs_line7i.miso_bf.flt.A3SS.gff line6i_vs_line7i.miso_bf.flt.line7.A3SS.trimmed.snps.vcf line7.variant

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/check_strand.py line6i_vs_line7i.miso_bf.flt.A3SS.gff line7.variant; \
		python ~/mdv-protocol/check_strand.py line6i_vs_line7i.miso_bf.flt.A3SS.gff line6.variant

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/find_diff_snps.py line6i_vs_line7i.miso_bf.flt.line6.A3SS.trimmed.snps.vcf line6i_vs_line7i.miso_bf.flt.line7.A3SS.snps.vcf 20 | sort -k1,1 -k2,2n  > line6i_vs_line7i.miso_bf.flt.line6.A3SS.diff.snps.vcf

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/find_diff_snps.py line6i_vs_line7i.miso_bf.flt.line7.A3SS.trimmed.snps.vcf line6i_vs_line7i.miso_bf.flt.line6.A3SS.snps.vcf 20 | sort -k1,1 -k2,2n  > line6i_vs_line7i.miso_bf.flt.line7.A3SS.diff.snps.vcf

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.A3SS.gff -a line6i_vs_line7i.miso_bf.flt.line6.A3SS.diff.snps.vcf -wb > line6i_vs_line7i.miso_bf.flt.line6.A3SS.diff.snps.exons.vcf

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.A3SS.gff -a line6i_vs_line7i.miso_bf.flt.line7.A3SS.diff.snps.vcf -wb > line6i_vs_line7i.miso_bf.flt.line7.A3SS.diff.snps.exons.vcf

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in line6.variant*.fa; do \
			blat $$f $$(echo $$f | sed s/line6/line7/) -out=blast $$f.blast; \
		done

	cd miso/cuffref-results/A3SS/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in line6.variant*.fa.blast; do \
			echo "working on " $$f; \
			python ~/mdv-protocol/blast_to_snps.py $$f > $$f.snps; \
		done


miso-snps-a5ss:

	cd miso; grep exon asm_cuff_ref_models.flt.A5SS.gff > asm_cuff_ref_models.flt.A5SS.exons.gff
	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line7u.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line7u.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line7u/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line7u.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line7u.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6i_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line6i_vs_line7i.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6i_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line6i_vs_line7i.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line6i.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line6i.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/A5SS/comparisons/line6u_vs_line6i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line6u_vs_line6i.miso_bf.flt $(path)/gal4selected.fa line6u_vs_line6i.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line7u_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line7u_vs_line7i.miso_bf.flt.line6.snps $(path)/miso/bam-data/line6u.bam $(path)/miso/bam-data/line6i.bam

	cd miso/cuffref-results/A5SS/comparisons/line7u_vs_line7i/bayes-factors; \
		python ~/miso-protocol/miso-call-snps.py line7u_vs_line7i.miso_bf.flt $(path)/gal4selected.fa line7u_vs_line7i.miso_bf.flt.line7.snps $(path)/miso/bam-data/line7u.bam $(path)/miso/bam-data/line7i.bam

find-deu-snps-a5ss:

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/select_events.py /mnt/ls12/preeyanon/mdv-pipeline/miso/asm_cuff_ref_models.flt.A5SS.gff \
		line6i_vs_line7i.miso_bf.flt A5SS | sort -k1,1 -k4,5n > line6i_vs_line7i.miso_bf.flt.A5SS.gff

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat *line6.snps_* | sort -k1,1 -k2,2n | awk '$$6>19' >> line6i_vs_line7i.miso_bf.flt.line6.snps.vcf; \
		cat *line7.snps_* | sort -k1,1 -k2,2n | awk '$$6>19' >> line6i_vs_line7i.miso_bf.flt.line7.snps.vcf

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat ~/mdv-protocol/vcf.header.txt > line6i_vs_line7i.miso_bf.flt.line6.A5SS.snps.vcf; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.A5SS.gff -a line6i_vs_line7i.miso_bf.flt.line6.snps.vcf -wb \
		| sort -k1,1 -k2,2n >> line6i_vs_line7i.miso_bf.flt.line6.A5SS.snps.vcf

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		cat ~/mdv-protocol/vcf.header.txt > line6i_vs_line7i.miso_bf.flt.line7.A5SS.snps.vcf; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.A5SS.gff -a line6i_vs_line7i.miso_bf.flt.line7.snps.vcf -wb \
		| sort -k1,1 -k2,2n >> line6i_vs_line7i.miso_bf.flt.line7.A5SS.snps.vcf

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		GATK -T LeftAlignAndTrimVariants -R /mnt/ls12/preeyanon/mdv-pipeline/gal4selected.sorted.fa -trim -V:testsnp,VCF line6i_vs_line7i.miso_bf.flt.line7.A5SS.snps.vcf -o line6i_vs_line7i.miso_bf.flt.line7.A5SS.trimmed.snps.vcf; \
		GATK -T LeftAlignAndTrimVariants -R /mnt/ls12/preeyanon/mdv-pipeline/gal4selected.sorted.fa -trim -V:testsnp,VCF line6i_vs_line7i.miso_bf.flt.line6.A5SS.snps.vcf -o line6i_vs_line7i.miso_bf.flt.line6.A5SS.trimmed.snps.vcf

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/get_variant_seq.py line6i_vs_line7i.miso_bf.flt.A5SS.gff line6i_vs_line7i.miso_bf.flt.line6.A5SS.trimmed.snps.vcf line6.variant; \
		python ~/mdv-protocol/get_variant_seq.py line6i_vs_line7i.miso_bf.flt.A5SS.gff line6i_vs_line7i.miso_bf.flt.line7.A5SS.trimmed.snps.vcf line7.variant

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/mdv-protocol/check_strand.py line6i_vs_line7i.miso_bf.flt.A5SS.gff line7.variant; \
		python ~/mdv-protocol/check_strand.py line6i_vs_line7i.miso_bf.flt.A5SS.gff line6.variant

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/find_diff_snps.py line6i_vs_line7i.miso_bf.flt.line6.A5SS.trimmed.snps.vcf line6i_vs_line7i.miso_bf.flt.line7.A5SS.snps.vcf 20 | sort -k1,1 -k2,2n  > line6i_vs_line7i.miso_bf.flt.line6.A5SS.diff.snps.vcf

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		python ~/miso-protocol/find_diff_snps.py line6i_vs_line7i.miso_bf.flt.line7.A5SS.trimmed.snps.vcf line6i_vs_line7i.miso_bf.flt.line6.A5SS.snps.vcf 20 | sort -k1,1 -k2,2n  > line6i_vs_line7i.miso_bf.flt.line7.A5SS.diff.snps.vcf

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.A5SS.gff -a line6i_vs_line7i.miso_bf.flt.line6.A5SS.diff.snps.vcf -wb > line6i_vs_line7i.miso_bf.flt.line6.A5SS.diff.snps.exons.vcf

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		intersectBed -b line6i_vs_line7i.miso_bf.flt.A5SS.gff -a line6i_vs_line7i.miso_bf.flt.line7.A5SS.diff.snps.vcf -wb > line6i_vs_line7i.miso_bf.flt.line7.A5SS.diff.snps.exons.vcf

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in line6.variant*.fa; do \
			blat $$f $$(echo $$f | sed s/line6/line7/) -out=blast $$f.blast; \
		done

	cd miso/cuffref-results/A5SS/comparisons/line6i_vs_line7i/bayes-factors; \
		for f in line6.variant*.fa.blast; do \
			echo "working on " $$f; \
			python ~/mdv-protocol/blast_to_snps.py $$f > $$f.snps; \
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

degenes-annotation:

	Rscript protocol/degenes_annotation.R

degenes-annotation-human:

	Rscript protocol/degenes_annotation_human.R

goseq-kegg-degenes:

	Rscript protocol/goseq_line6.R
	Rscript protocol/goseq_line7.R
