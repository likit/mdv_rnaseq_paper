projdir=/mnt/ls12/preeyanon/mdv/
line6ui=/comparisons/line6u_vs_line6i/bayes-factors/
line7ui=/comparisons/line7u_vs_line7i/bayes-factors/
line67u=/comparisons/line6u_vs_line7u/bayes-factors/
line67i=/comparisons/line6i_vs_line7i/bayes-factors/

all: degenes miso misc

misc:

	scp -C preeyano@hpc.msu.edu:$(projdir)/gimme/gimme.bed results/

degenes:

	scp -C preeyano@hpc.msu.edu:$(projdir)/gimme/line6u_vs_i.degenes.fdr.05 results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/gimme/line7u_vs_i.degenes.fdr.05 results/

	scp -C preeyano@hpc.msu.edu:$(projdir)/gimme/line7u_vs_i.degenes.fdr.05.fa.tophits.txt results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/gimme/line6u_vs_i.degenes.fdr.05.fa.tophits.txt results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/gimme/line7u_vs_i.degenes.fdr.05.fa.longest.gallus.xml results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/gimme/line6u_vs_i.degenes.fdr.05.fa.longest.gallus.xml results/
	
	scp -C preeyano@hpc.msu.edu:$(projdir)/gimme/line6u_vs_i.degenes.fdr.05.fa.tophits.txt results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/gimme/line7u_vs_i.degenes.fdr.05.fa.tophits.txt results/

miso:

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/gimme.flt.SE.gff results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/gimme.flt.A3SS.gff results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/gimme.flt.A5SS.gff results/
	
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line6ui)/line6u_vs_line6i.miso_bf \
		results/line6u_vs_line6i.miso_bf.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line7ui)/line7u_vs_line7i.miso_bf \
		results/line7u_vs_line7i.miso_bf.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67u)/line6u_vs_line7u.miso_bf \
		results/line6u_vs_line7u.miso_bf.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67i)/line6i_vs_line7i.miso_bf \
		results/line6i_vs_line7i.miso_bf.se

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line6ui)/line6u_vs_line6i.miso_bf \
		results/line6u_vs_line6i.miso_bf.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line7ui)/line7u_vs_line7i.miso_bf \
		results/line7u_vs_line7i.miso_bf.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67u)/line6u_vs_line7u.miso_bf \
		results/line6u_vs_line7u.miso_bf.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67i)/line6i_vs_line7i.miso_bf \
		results/line6i_vs_line7i.miso_bf.a3ss

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line6ui)/line6u_vs_line6i.miso_bf \
		results/line6u_vs_line6i.miso_bf.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line7ui)/line7u_vs_line7i.miso_bf \
		results/line7u_vs_line7i.miso_bf.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67u)/line6u_vs_line7u.miso_bf \
		results/line6u_vs_line7u.miso_bf.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67i)/line6i_vs_line7i.miso_bf \
		results/line6i_vs_line7i.miso_bf.a5ss

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line6ui)/line6u_vs_line6i.miso_bf.flt \
		results/line6u_vs_line6i.miso_bf.flt.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line7ui)/line7u_vs_line7i.miso_bf.flt \
		results/line7u_vs_line7i.miso_bf.flt.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67u)/line6u_vs_line7u.miso_bf.flt \
		results/line6u_vs_line7u.miso_bf.flt.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67i)/line6i_vs_line7i.miso_bf.flt \
		results/line6i_vs_line7i.miso_bf.flt.se

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line6ui)/line6u_vs_line6i.miso_bf.flt \
		results/line6u_vs_line6i.miso_bf.flt.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line7ui)/line7u_vs_line7i.miso_bf.flt \
		results/line7u_vs_line7i.miso_bf.flt.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67u)/line6u_vs_line7u.miso_bf.flt \
		results/line6u_vs_line7u.miso_bf.flt.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67i)/line6i_vs_line7i.miso_bf.flt \
		results/line6i_vs_line7i.miso_bf.flt.a3ss

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line6ui)/line6u_vs_line6i.miso_bf.flt \
		results/line6u_vs_line6i.miso_bf.flt.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line7ui)/line7u_vs_line7i.miso_bf.flt \
		results/line7u_vs_line7i.miso_bf.flt.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67u)/line6u_vs_line7u.miso_bf.flt \
		results/line6u_vs_line7u.miso_bf.flt.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67i)/line6i_vs_line7i.miso_bf.flt \
		results/line6i_vs_line7i.miso_bf.flt.a5ss

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line6ui)/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt \
		results/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line7ui)/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt \
		results/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67u)/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt \
		results/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67i)/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt \
		results/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt.se

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line6ui)/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt \
		results/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line7ui)/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt \
		results/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67u)/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt \
		results/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67i)/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt \
		results/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt.a3ss

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line6ui)/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt \
		results/line6u_vs_line6i.miso_bf.flt.fa.gallus.tophits.txt.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line7ui)/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt \
		results/line7u_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67u)/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt \
		results/line6u_vs_line7u.miso_bf.flt.fa.gallus.tophits.txt.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67i)/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt \
		results/line6i_vs_line7i.miso_bf.flt.fa.gallus.tophits.txt.a5ss
	
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line6ui)/line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt \
		results/line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line7ui)/line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt \
		results/line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67u)/line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt \
		results/line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt.se
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67i)/line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt \
		results/line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt.se

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line6ui)/line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt \
		results/line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line7ui)/line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt \
		results/line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67u)/line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt \
		results/line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt.a3ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67i)/line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt \
		results/line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt.a3ss

	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line6ui)/line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt \
		results/line6u_vs_line6i.miso_bf.flt.fa.human.tophits.txt.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line7ui)/line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt \
		results/line7u_vs_line7i.miso_bf.flt.fa.human.tophits.txt.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67u)/line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt \
		results/line6u_vs_line7u.miso_bf.flt.fa.human.tophits.txt.a5ss
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67i)/line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt \
		results/line6i_vs_line7i.miso_bf.flt.fa.human.tophits.txt.a5ss
	
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67i)/line6i_vs_line7i.miso_bf.flt.line6.SE.diff.snps.exons.csv \
		results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67i)/line6i_vs_line7i.miso_bf.flt.line6.A3SS.diff.snps.exons.csv \
		results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67i)/line6i_vs_line7i.miso_bf.flt.line6.A5SS.diff.snps.exons.csv \
		results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/SE/$(line67i)/line6i_vs_line7i.miso_bf.flt.line7.SE.diff.snps.exons.csv \
		results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A3SS/$(line67i)/line6i_vs_line7i.miso_bf.flt.line7.A3SS.diff.snps.exons.csv \
		results/
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/A5SS/$(line67i)/line6i_vs_line7i.miso_bf.flt.line7.A5SS.diff.snps.exons.csv \
		results/
	
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/plots/chr9-25729.ev1.pdf \
		results/pfn2_miso.pdf
	scp -C preeyano@hpc.msu.edu:$(projdir)/miso/results/plots/chr7-24308.ev1.pdf \
		results/itgb2_miso.pdf
