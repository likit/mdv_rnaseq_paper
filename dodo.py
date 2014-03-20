import os
import doit

def task_ebseq_cuff_ref():
    '''run ebseq'''

    os.chdir(doit.initial_workdir)
    os.chdir('tophat/merged_cuff_ref')
    for sample in ['line6', 'line7']:
        cmd = 'rsem-generate-data-matrix %su-single-rsem-full.genes.results' \
                ' %su-paired-rsem-full.genes.results %si-single-rsem-full.genes.results' \
                ' %si-paired-rsem-full.genes.results > %su_vs_i.gene.counts.matrix'
        cmd1 = cmd % ((sample,) * 5)
        cmd2 = 'rsem-run-ebseq %su_vs_i.gene.counts.matrix 2,2 %su_vs_i.cuffref.degenes' % (sample, sample)
        cmd3 = 'rsem-control-fdr %su_vs_i.degenes 0.05 %su_vs_i.cuffref.degenes.fdr.05' % (sample, sample)
        yield {'name':sample,
                'actions':[cmd1, cmd2, cmd3]
            }

def task_prepare_for_blast2go_cuff_ref():
    '''prepare sequences and alignments for blast2go'''

    os.chdir(doit.initial_workdir)
    os.chdir('tophat/merged_cuff_ref')
    for sample in ['line6', 'line7']:
        cmd1 = "python ~/mdv-protocol/rsem-output-cufflinks-to-fasta.py " \
                "%su_vs_i.cuffref.degenes.fdr.05 merged.bed.fa " \
                "> %su_vs_i.cuffref.degenes.fdr.05.fa" % (sample, sample)

        cmd2 = "estscan -t %su_vs_i.cuffref.degenes.fdr.05.fa.prot " \
                "-M ~/mdv-protocol/gallus.hm %su_vs_i.cuffref.degenes.fdr.05.fa " \
                "> %su_vs_i.cuffref.degenes.fdr.05.fa.nucl" % ((sample,) * 3)

        cmd3 = "python ~/mdv-protocol/gene-rep.py %su_vs_i.cuffref.degenes.fdr.05.fa.prot " \
                "> %su_vs_i.cuffref.degenes.fdr.05.fa.prot.longest" % (sample, sample)

        cmd4 = "python ~/mdv-protocol/gene-rep.py %su_vs_i.cuffref.degenes.fdr.05.fa " \
                "> %su_vs_i.cuffref.degenes.fdr.05.fa.longest" % (sample, sample)

        yield {'name':sample,
                'actions':[cmd1, cmd2,
                            cmd3, cmd4]
            }


def task_blast_cuff_ref_human():
    '''run BLAST'''

    os.chdir(doit.initial_workdir)
    os.chdir('tophat/merged_cuff_ref')
    for sample in ['line6', 'line7']:
        cmd1 = 'qsub -v db="Human_prot",input="%su_vs_i.cuffref.degenes.fdr.05.fa.prot.longest"' \
                ',program="blastp",output="%su_vs_i.cuffref.degenes.fdr.05.fa.prot.longest.xml"' \
                ' ~/mdv-protocol/blast.sh' % (sample, sample)

        cmd2 = 'qsub -v db="Human_prot",input="%su_vs_i.cuffref.degenes.fdr.05.fa.longest"' \
                ',program="blastx",output="%su_vs_i.cuffref.degenes.fdr.05.fa.longest.gallus.xml"' \
                ' ~/mdv-protocol/blast.sh' % (sample, sample)

        yield {'name':sample,
                'actions':[cmd1, cmd2,]
            }

def task_blast_cuff_ref_gallus():
    '''run BLAST'''

    os.chdir(doit.initial_workdir)
    os.chdir('tophat/merged_cuff_ref')
    for sample in ['line6', 'line7']:
        cmd1 = 'qsub -v db="Gallus_prot",input="%su_vs_i.cuffref.degenes.fdr.05.fa.prot.longest"' \
                ',program="blastp",output="%su_vs_i.cuffref.degenes.fdr.05.fa.prot.longest.xml"' \
                ' ~/mdv-protocol/blast.sh' % (sample, sample)

        cmd2 = 'qsub -v db="Gallus_prot",input="%su_vs_i.cuffref.degenes.fdr.05.fa.longest"' \
                ',program="blastx",output="%su_vs_i.cuffref.degenes.fdr.05.fa.longest.gallus.xml"' \
                ' ~/mdv-protocol/blast.sh' % (sample, sample)

        yield {'name':sample,
                'actions':[cmd2,]
            }

def task_blast2go_cuff_ref():
    '''run blast2go'''

    os.chdir(doit.initial_workdir)
    os.chdir('tophat/merged_cuff_ref')
    for sample in ['line6', 'line7']:
        cmd1 = 'qsub -v input="%su_vs_i.cuffref.degenes.fdr.05.fa.prot.longest.xml"' \
                ',outdir="%su_cuffref_prot_blast2go_outdir"' \
                ' ~/mdv-protocol/b2g_job.sh' % (sample, sample)
        cmd2 = 'qsub -v input="%su_vs_i.cuffref.degenes.fdr.05.fa.longest.xml"' \
                ',outdir="%su_cuffref_nucl_blast2go_outdir"' \
                ' ~/mdv-protocol/b2g_job.sh' % (sample, sample)

        yield {'name':sample,
                'actions':[cmd1, cmd2,]
            }

def task_filter_low_isopct():
    os.chdir(doit.initial_workdir)
    cmd = 'python ~/mdv-protocol/filter-low-isopct.py 1.0 asm_cuff_ref_models.bed *cuffref*isoforms.results > asm_cuff_ref_models.flt.bed'
    return {
            'actions': [cmd],
            }
