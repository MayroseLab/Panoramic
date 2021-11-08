"""
Use EVM to annotate a genome.
The following steps are included:
1. Mask TEs in the input genome using EDTA (optional)
2. Liftover reference genes using Liftoff (optional)
3. Run three ab-initio predictors (all optional):
   - Augustus
   - GlimmerHMM
   - SNAP
4. Run PASA assembly to process transcript evidence (optional)
5. Run genomeThreader to process protein evidence (optional)
6. Use outputs from steps 2-5 as inputs to EVM
7. Detect and correct chimeric genes using chimeraBuster (optional)
8. Calculate AED scores for gene models and use them to filter low qulity predictions
9. Produce final outputs (gff and fasta files)
"""

import os
import re
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir) + '/util'
import sys
sys.path.append(utils_dir)
from snakemakeUtils import *

PIPELINE = 'EVM-annotation'

LOGS_DIR = config['out_dir'] + "/logs"
CONDA_ENV_DIR = os.path.dirname(pipeline_dir) + "/conda_env"

onstart:
    write_config_file(config)
    if not os.path.isdir(LOGS_DIR):
        os.mkdir(LOGS_DIR)

onsuccess:
    print("%s pipeline finished, no error" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

onerror:
    print("%s pipeline failed" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

#------------------------------------
#                RULES              |
#------------------------------------

localrules: all, skip_liftover, skip_augustus, skip_glimmerHmm, skip_snap, skip_genomeThreader, skip_PASA_assembly

rule all:
    input:
        os.path.join(config['out_dir'], 'EVM.filter.rename.gff3'),
        prot=os.path.join(config['out_dir'], 'EVM.filter.rename.prot.fasta'),
        trans=os.path.join(config['out_dir'], 'EVM.filter.rename.trans.fasta')

# determine chromosomes list and partitions
chr_list = {}
with open(config['input_genome']) as f:
    for line in f:
        if line.startswith('>'):
            chr_name = line.strip().split()[0][1:]
            chr_list[chr_name] = 0
        else:
            chr_list[chr_name] += len(line.strip())
    for chr_name in chr_list:
        if chr_list[chr_name] <= config['overlap_size']:
            chr_list[chr_name] = ['%s-%s' %(1,chr_list[chr_name])]
        else:
            chr_list[chr_name] = ["%s-%s" %(start, min(start+config['segment_size']-1, chr_list[chr_name])) for start in range(1, chr_list[chr_name], config['segment_size'] - config['overlap_size']) if (min(start+config['segment_size']-1, chr_list[chr_name]) - start) > config['overlap_size']]


if config['mask_genome'] == 1:
    genome_fasta = os.path.join(config['out_dir'], os.path.basename(config['input_genome'])+'.mod.MAKER.masked')
    if 'reference_cds' in config and config['reference_cds']:
        ref_cds = config['reference_cds']
    else:
        ref_cds = None
else:
    genome_fasta = config['input_genome']
    ref_cds = None

if ref_cds:
    rule mask_repeats:
        input:
            genome=config['input_genome'],
            cds=config['reference_cds']
        output:
            os.path.join(config['out_dir'], os.path.basename(config['input_genome'])+'.mod.MAKER.masked')
        params:
            sample=config['sample_name'],
            out_dir=config['out_dir'],
            queue=config['queue'],
            priority=config['priority'],
            ppn=config['ppn'],
            logs_dir=LOGS_DIR
        conda:
            CONDA_ENV_DIR + '/edta.yml'
        shell:
            """
            cd {params.out_dir}
            EDTA.pl --genome {input.genome} --cds {input.cds} --anno 1 --sensitive --threads {params.ppn} --force 1 || true
            """
else:
    rule mask_repeats:
        input:
            genome=config['input_genome'],
        output:
            os.path.join(config['out_dir'], os.path.basename(config['input_genome'])+'.mod.MAKER.masked')
        params:
            sample=config['sample_name'],
            out_dir=config['out_dir'],
            queue=config['queue'],
            priority=config['priority'],
            ppn=config['ppn'],
            logs_dir=LOGS_DIR
        conda:
            CONDA_ENV_DIR + '/edta.yml'
        shell:
            """
            cd {params.out_dir}
            EDTA.pl --genome {input.genome} --anno 1 --sensitive --threads {params.ppn} --force 1 || true
            """

if config['reference_liftover'] == 1:
    rule reference_liftover:
        """
        Perform liftover of reference
        genes on input genome
        """
        input:
            src_gff=config['reference_gff'],
            src_genome=config['reference_fasta'],
            target_genome=genome_fasta
        output:
            os.path.join(config['out_dir'], 'liftover.gff3')
        params:
            sample=config['sample_name'],
            queue=config['queue'],
            priority=config['priority'],
            ppn=config['ppn'],
            logs_dir=LOGS_DIR
        conda:
            CONDA_ENV_DIR + '/liftoff.yml'
        shell:
            """
            liftoff {input.target_genome} {input.src_genome} -g {input.src_gff} -o {output} -p {params.ppn}
            """
else:
    rule skip_liftover:
        """
        Create empty liftover result
        """
        output:
            os.path.join(config['out_dir'], 'liftover.gff3')
        params:
            sample=config['sample_name'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR
        shell:
            """
            touch {output}
            """
    config['liftover_weight'] = 0

rule split_genome_to_chr:
    """
    Split the input genome fasta to
    multiple fastas with one chromosome
    in each.
    """
    input:
        genome_fasta
    output:
        expand(os.path.join(config['out_dir'],'{CHR}','{CHR}.fa'), CHR=chr_list.keys())
    params:
        sample=config['sample_name'],
        split_script = os.path.join(utils_dir, 'split_multifasta.py'),
        out_dir=config['out_dir'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/biopython.yml'
    shell:
        """
        python {params.split_script} {input} {params.out_dir} True
        """

if config['augustus_species']:

    rule run_augustus:
        input:
            os.path.join(config['out_dir'],'{CHR}','{CHR}.fa')
        output:
            os.path.join(config['out_dir'],'{CHR}','augustus.out')
        params:
            sample=config['sample_name'],
            species=config['augustus_species'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        conda:
            CONDA_ENV_DIR + '/augustus.yml'
        shell:
            """
            augustus --gff3=on --species={params.species} {input} > {output}
            """
else:
    rule skip_augustus:
        output:
            os.path.join(config['out_dir'],'{CHR}','augustus.out')
        params:
            sample=config['sample_name'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        shell:
            """
            touch {output}
            """

if config['glimmerhmm_species']:
    rule run_glimmerHmm:
        input:
            os.path.join(config['out_dir'],'{CHR}','{CHR}.fa')
        output:
            os.path.join(config['out_dir'],'{CHR}','glimmerHmm.out')
        params:
            sample=config['sample_name'],
            species=config['glimmerhmm_species'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        conda:
            CONDA_ENV_DIR + '/glimmerHmm.yml'
        shell:
            """
            speciesDir="$CONDA_PREFIX/share/glimmerhmm/trained_dir/{params.species}/"
            glimmerhmm {input} $speciesDir -o {output} -g
            """
else:
    rule skip_glimmerHmm:
        output:
            os.path.join(config['out_dir'],'{CHR}','glimmerHmm.out')
        params:
            sample=config['sample_name'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        shell:
            """
            touch {output}
            """

if config['snap_species']:

    rule run_snap:
        input:
            os.path.join(config['out_dir'],'{CHR}','{CHR}.fa')
        output:
            os.path.join(config['out_dir'],'{CHR}','snap.out')
        params:
            sample=config['sample_name'],
            species=config['snap_species'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        conda:
            CONDA_ENV_DIR + '/snap.yml'
        shell:
            """
            if [ -f {params.species} ]; then
                speciesHmm={params.species}
            else
                speciesHmm="$CONDA_PREFIX/share/snap/HMM/{params.species}.hmm"
            fi
            snap $speciesHmm {input} > {output}
            """
else:
    rule skip_snap:
        output:
            os.path.join(config['out_dir'],'{CHR}','snap.out')
        params:
            sample=config['sample_name'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        shell:
            """
            touch {output}
            """

if config['proteins_fasta']:
    rule copy_proteins:
        """
        Copy input proteins to avoid
        clashing with other genomeThreader
        runs.
        """
        input:
            config['proteins_fasta']
        output:
            os.path.join(config['out_dir'],'{CHR}','proteins.fasta')
        params:
            sample=config['sample_name'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        shell:
            """
            cp {input} {output}
            """

    rule run_genomeThreader:
        input:
            chrom=os.path.join(config['out_dir'],'{CHR}','{CHR}.fa'),
            proteins=os.path.join(config['out_dir'],'{CHR}','proteins.fasta')
        output:
            os.path.join(config['out_dir'],'{CHR}','genomeThreader.out')
        params:
            sample=config['sample_name'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        conda:
            CONDA_ENV_DIR + '/genomethreader.yml'
        shell:
            """
            gth -genomic {input.chrom} -protein {input.proteins} -o {output} -gff3out -intermediate -fastdp
            """
else:
    rule skip_genomeThreader:
        """
        Create empty genomeThreader output
        """
        output:
            os.path.join(config['out_dir'],'{CHR}','genomeThreader.out')
        params:
            sample=config['sample_name'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        shell:
            """
            touch {output}
            """
    config['proteins_weight'] = 0


if config['transcripts_fasta']:
    rule prep_PASA_conf:
        """
        Prepare the configuration file for PASA assembly
        """
        input:
            fasta=os.path.join(config['out_dir'],'{CHR}','{CHR}.fa'),
            template=os.path.join(pipeline_dir, 'alignAssembly.config.template')
        output:
            os.path.join(config['out_dir'], '{CHR}','alignAssembly.config')
        params:
            sample=config['sample_name'],
            db_path=lambda wildcards: os.path.join(config['out_dir'], wildcards.CHR, 'pasa_db.sqlite'),
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR
        shell:
            """
            sed 's|<DBPATH>|{params.db_path}|' {input} > {output}
            """
    rule PASA_assembly:
        input:
            fasta=os.path.join(config['out_dir'],'{CHR}','{CHR}.fa'),
            transcripts=config['transcripts_fasta'],
            pasa_conf=os.path.join(config['out_dir'], '{CHR}','alignAssembly.config')
        output:
            os.path.join(config['out_dir'], '{CHR}', 'pasa_db.sqlite.pasa_assemblies.gff3')
        params:
            sample=config['sample_name'],
            exec_dir=lambda wildcards: os.path.join(config['out_dir'], wildcards.CHR),
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
            ppn=config['ppn']
        conda:
            CONDA_ENV_DIR + '/pasa.yml'
        shell:
            """
            cd {params.exec_dir}
            $CONDA_PREFIX/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c {input.pasa_conf} -C -R -g {input.fasta} -t {input.transcripts} --ALIGNERS blat --TRANSDECODER --CPU {params.ppn} || touch {output}
            """
else:
    rule skip_PASA_assembly:
        """
        Create empty PASA output
        """
        output:
            os.path.join(config['out_dir'], '{CHR}', 'pasa_db.sqlite.pasa_assemblies.gff3')
        params:
            sample=config['sample_name'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        shell:
            """
            touch {output}
            """
    config['transcripts_weight'] = 0

rule combine_chromosomes:
    input:
        augustus=expand(os.path.join(config['out_dir'],'{CHR}','augustus.out'), CHR=chr_list.keys()),
        glimmerhmm=expand(os.path.join(config['out_dir'],'{CHR}','glimmerHmm.out'), CHR=chr_list.keys()),
        snap=expand(os.path.join(config['out_dir'],'{CHR}','snap.out'), CHR=chr_list),
        genomethreader=expand(os.path.join(config['out_dir'],'{CHR}','genomeThreader.out'), CHR=chr_list.keys()),
        pasa=expand(os.path.join(config['out_dir'],'{CHR}','pasa_db.sqlite.pasa_assemblies.gff3'), CHR=chr_list.keys())
    output:
        all_augustus=os.path.join(config['out_dir'],'augustus.out'),
        all_glimmerhmm=os.path.join(config['out_dir'],'glimmerHmm.out'),
        all_snap=os.path.join(config['out_dir'],'snap.out'),
        all_genomethreader=os.path.join(config['out_dir'],'genomeThreader.out'),
        all_pasa=os.path.join(config['out_dir'],'pasa_db.sqlite.pasa_assemblies_all.gff3')
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cat {input.augustus} > {output.all_augustus}
        cat {input.glimmerhmm} > {output.all_glimmerhmm}
        cat {input.snap} > {output.all_snap}
        cat {input.genomethreader} > {output.all_genomethreader}
        cat {input.pasa} > {output.all_pasa}
        """

rule convert_augustus_to_EVM_gff3:
    input:
        os.path.join(config['out_dir'],'augustus.out')
    output:
        os.path.join(config['out_dir'],'augustus.EVM.gff3')
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/EVM.yml'
    shell:
        """
        $CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl {input} > {output}
        """

rule convert_snap_zff_to_gff3:
    input:
        os.path.join(config['out_dir'],'snap.out')
    output:
        os.path.join(config['out_dir'],'snap.gff3')
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snap.yml'
    shell:
        """
        zff2gff3.pl {input} > {output}
        """

rule convert_snap_gff3_to_EVM_gff3:
    input:
        os.path.join(config['out_dir'],'snap.gff3')
    output:
        os.path.join(config['out_dir'],'snap.EVM.gff3')
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/EVM.yml'
    shell:
        """
        $CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/misc/SNAP_to_GFF3.pl {input} > {output}
        """

rule convert_glimmerHmm_to_EVM_gff3:
    input:
        os.path.join(config['out_dir'],'glimmerHmm.out')
    output:
        os.path.join(config['out_dir'],'glimmerHmm.EVM.gff3')
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/EVM.yml'
    shell:
        """
        $CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/misc/glimmerHMM_to_GFF3.pl {input} > {output}
        """

rule convert_genomeThreader_to_EVM_gff3:
    input:
        os.path.join(config['out_dir'],'genomeThreader.out')
    output:
        os.path.join(config['out_dir'],'genomeThreader.EVM.gff3')
    params:
        sample=config['sample_name'],
        convert_script=os.path.join(pipeline_dir, 'genomeThreader_to_evm_gff3.pl'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/EVM.yml'
    shell:
        """
        {params.convert_script} {input} > {output}
        """

rule collect_predictions:
    """
    Concatenate Augustus, SNAP, GlimmerHmm
    and genomeThreader into one gff3 file
    """
    input:
        liftover=os.path.join(config['out_dir'], 'liftover.gff3'),
        augustus=os.path.join(config['out_dir'],'augustus.EVM.gff3'),
        snap=os.path.join(config['out_dir'],'snap.EVM.gff3'),
        glimmerHmm=os.path.join(config['out_dir'],'glimmerHmm.EVM.gff3'),
        genomeThreader=os.path.join(config['out_dir'],'genomeThreader.EVM.gff3')
    output:
        os.path.join(config['out_dir'],'all_predictions.EVM.gff3')
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cat {input} > {output}
        """

rule create_weights_file:
    output:
        os.path.join(config['out_dir'],'weights.tsv')
    params:
        sample=config['sample_name'],
        liftover_weight=config['liftover_weight'],
        ab_initio_weight=config['ab-initio_weight'],
        transcripts_weight=config['transcripts_weight'],
        proteins_weight=config['proteins_weight'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        echo -e "OTHER_PREDICTION\tLiftoff\t{params.liftover_weight}" >> {output}
        echo -e "ABINITIO_PREDICTION\tAugustus\t{params.ab_initio_weight}" >> {output}
        echo -e "ABINITIO_PREDICTION\tGlimmerHMM\t{params.ab_initio_weight}" >> {output}
        echo -e "ABINITIO_PREDICTION\tSNAP\t{params.ab_initio_weight}" >> {output}
        echo -e "TRANSCRIPT\tassembler-pasa_db.sqlite\t{params.transcripts_weight}" >> {output}
        echo -e "OTHER_PREDICTION\tgenomeThreader\t{params.proteins_weight}" >> {output}
        """

def partitions_files(chr_dict, base_dir, file_name):
    """
    Expects a dict like:
    {'1': ['1-500000','480001-980000',...], '2': ['1-500000','480001-980000',...]}
    """
    res = []
    for chrom, partitions in chr_dict.items():
        for p in partitions:
            p_file = os.path.join(base_dir, chrom, '%s_%s' %(chrom, p), file_name)
            res.append(p_file)
    return res

rule create_EVM_partitions_list:
    """
    Create a TSV file in the format
    expected by EVM, to be used as
    partitions list
    """
    input:
        genome_fasta
    output:
        os.path.join(config['out_dir'],'EVM_patitions.list')
    params:
        sample=config['sample_name'],
        partition_script=os.path.join(pipeline_dir, 'create_EVM_partitions_list.py'),
        segment_size=config['segment_size'],
        overlap_size=config['overlap_size'],
        out_dir=config['out_dir'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/biopython.yml'
    shell:
        """
        python {params.partition_script} {input} {params.segment_size} {params.overlap_size} {params.out_dir} {output}
        """

rule partition_genome:
    input:
        genome=genome_fasta,
        lst=os.path.join(config['out_dir'],'EVM_patitions.list')
    output:
        partitions_files(chr_list, config['out_dir'], os.path.basename(genome_fasta))
    params:
        sample=config['sample_name'],
        partition_script=os.path.join(pipeline_dir, 'partition_fasta.py'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/biopython.yml'
    shell:
        """
        python {params.partition_script} {input.genome} {input.lst}
        """

rule partition_predictions:
    input:
        gff=os.path.join(config['out_dir'],'all_predictions.EVM.gff3'),
        lst=os.path.join(config['out_dir'],'EVM_patitions.list')
    output:
        partitions_files(chr_list, config['out_dir'], 'all_predictions.EVM.gff3')
    params:
        sample=config['sample_name'],
        partition_script=os.path.join(pipeline_dir, 'partition_gff.py'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.partition_script} {input.gff} {input.lst}
        """

rule partition_PASA:
    input:
        gff=os.path.join(config['out_dir'],'pasa_db.sqlite.pasa_assemblies_all.gff3'),
        lst=os.path.join(config['out_dir'],'EVM_patitions.list')
    output:
        partitions_files(chr_list, config['out_dir'], 'pasa_db.sqlite.pasa_assemblies_all.gff3')
    params:
        sample=config['sample_name'],
        partition_script=os.path.join(pipeline_dir, 'partition_gff.py'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        if [ -s {input.gff} ]
        then
        python {params.partition_script} {input.gff} {input.lst}
        else
        touch {output}
        fi
        """

rule run_EVM:
    """
    Run EVM on each partition
    """
    input:
        weights=os.path.join(config['out_dir'],'weights.tsv'),
        pred=os.path.join(config['out_dir'],'{CHR}','{CHR}_{partition}','all_predictions.EVM.gff3'),
        trans=os.path.join(config['out_dir'],'{CHR}','{CHR}_{partition}','pasa_db.sqlite.pasa_assemblies_all.gff3'),
        genome=os.path.join(config['out_dir'],'{CHR}','{CHR}_{partition}',os.path.basename(genome_fasta))
    output:
        evm_out=os.path.join(config['out_dir'],'{CHR}','{CHR}_{partition}','evm.out'),
        evm_log=os.path.join(config['out_dir'],'{CHR}','{CHR}_{partition}','evm.out.log')
    params:
        sample=config['sample_name'],
        genome=os.path.basename(genome_fasta),
        pred='all_predictions.EVM.gff3',
        trans='pasa_db.sqlite.pasa_assemblies_all.gff3',
        exec_dir=lambda wildcards: os.path.join(config['out_dir'],wildcards.CHR,'%s_%s' %(wildcards.CHR, wildcards.partition)),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/EVM.yml'
    shell:
        """
        evidence_modeler.pl -G {params.genome} -g {params.pred} -w {input.weights} -e {params.trans} --exec_dir {params.exec_dir} > {output.evm_out} 2> {output.evm_log}
        """

rule recombine_EVM_partitions:
    """
    Collect results from EVM partitions
    using a dedicated script. This rule
    creates a evm.out file per chromosome
    """
    input:
        partitions_list=os.path.join(config['out_dir'],'EVM_patitions.list'),
        evm_res=partitions_files(chr_list, config['out_dir'], 'evm.out')
    output:
        expand(os.path.join(config['out_dir'], '{CHR}', 'evm.out'), CHR=chr_list.keys())
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/EVM.yml'
    shell:
        """
        $CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions {input.partitions_list} --output_file_name evm.out
        """

rule convert_EVM_to_gff3:
    input:
        evm_done=expand(os.path.join(config['out_dir'], '{CHR}', 'evm.out'), CHR=chr_list.keys()),
        partitions_list=os.path.join(config['out_dir'],'EVM_patitions.list'),
        genome=genome_fasta
    output:
        expand(os.path.join(config['out_dir'], '{CHR}', 'evm.out.gff3'), CHR=chr_list.keys())
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/EVM.yml'
    shell:
        """
        $CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions {input.partitions_list} --output evm.out --genome {input.genome}
        """

rule combine_EVM_gffs:
    input:
        expand(os.path.join(config['out_dir'], '{CHR}', 'evm.out.gff3'), CHR=chr_list.keys())
    output:
        os.path.join(config['out_dir'], 'evm.out.gff3')
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cat {input} > {output}
        """

if config['split_chimeras'] == 1 and config['transcripts_fasta']:
    evm_gff = os.path.join(config['out_dir'], 'evm.out.gff3.chimeraBuster.corrected.gff')
else:
    evm_gff = os.path.join(config['out_dir'], 'evm.out.gff3')

if not config['transcripts_fasta']:
    config['transcripts_fasta'] = ''

rule detect_chimeras:
    input:
        gff=os.path.join(config['out_dir'], 'evm.out.gff3'),
        genome=genome_fasta,
        transcripts=config['transcripts_fasta']
    output:
        os.path.join(config['out_dir'], 'chimeric_genes.list')
    params:
        sample=config['sample_name'],
        chimeraBuster_dir=config['chimeraBuster_dir'],
        out_dir=config['out_dir'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        os.path.join(config['chimeraBuster_dir'],'env.yml')
    shell:
        """
        python {params.chimeraBuster_dir}/detect_chimeric_genes.py -g {input.gff} -f {input.genome} -t {input.transcripts} -o {params.out_dir} -c {params.ppn}
        """

rule correct_chimeras:
    input:
        chimeras_list=os.path.join(config['out_dir'], 'chimeric_genes.list'),
        gff=os.path.join(config['out_dir'], 'evm.out.gff3'),
        genome=genome_fasta,
        transcripts=config['transcripts_fasta']
    output:
        evm_gff
    params:
        sample=config['sample_name'],
        chimeraBuster_dir=config['chimeraBuster_dir'],
        out_dir=config['out_dir'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        os.path.join(config['chimeraBuster_dir'],'env.yml')
    shell:
        """
        python {params.chimeraBuster_dir}/correct_chimeric_genes.py -g {input.gff} -f {input.genome} -t {input.transcripts} -l {input.chimeras_list} -o {params.out_dir} -c {params.ppn}
        """

rule calculate_AED:
    """
    Calculate and add AED scores for
    gene models predicted by EVM.
    """
    input:
        gff=evm_gff,
        pred=os.path.join(config['out_dir'],'all_predictions.EVM.gff3'),
        transcripts=os.path.join(config['out_dir'],'pasa_db.sqlite.pasa_assemblies_all.gff3'),
        weights=os.path.join(config['out_dir'],'weights.tsv'),
        feature_types=os.path.join(pipeline_dir, 'feature_types.tsv')
    output:
        os.path.join(config['out_dir'], 'EVM.AED.gff3')
    params:
        sample=config['sample_name'],
        add_aed_script=os.path.join(pipeline_dir, 'add_AED_to_gff3.py'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/add_AED.yml'
    shell:
        """
        python {params.add_aed_script} -g {input.gff} -t {input.transcripts} -a {input.pred} -w {input.weights} -l {input.feature_types} -i -o {output}
        """

rule filter_annotation:
    """
    Discard low-quality gene models
    based on AED and coding protein length
    """
    input:
        os.path.join(config['out_dir'], 'EVM.AED.gff3')
    output:
        os.path.join(config['out_dir'], 'EVM.filter.gff3')
    params:
        sample=config['sample_name'],
        filter_script=os.path.join(pipeline_dir, 'filter_gff.py'),
        max_AED=config['max_AED'],
        min_prot=config['min_protein'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.filter_script} {input} {params.max_AED} {params.min_prot} > {output}
        """

rule rename_features:
    """
    Add the sample name to IDs
    of all features.
    """
    input:
        os.path.join(config['out_dir'], 'EVM.filter.gff3')
    output:
        os.path.join(config['out_dir'], 'EVM.filter.rename.gff3')
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        sed -e 's/ID=/ID={params.sample}_/' -e 's/Parent=/Parent={params.sample}_/' {input} > {output}
        """

rule extract_fasta_sequences:
    """
    Extract protein and transcript
    sequences as FASTA
    """
    input:
        gff=os.path.join(config['out_dir'], 'EVM.filter.rename.gff3'),
        fasta=config['input_genome']
    output:
        prot=os.path.join(config['out_dir'], 'EVM.filter.rename.prot.fasta'),
        trans=os.path.join(config['out_dir'], 'EVM.filter.rename.trans.fasta')
    params:
        sample=config['sample_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/EVM.yml'
    shell:
        """
        $CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl {input.gff} {input.fasta} prot > {output.prot}
        $CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl {input.gff} {input.fasta} cDNA > {output.trans}
        """
