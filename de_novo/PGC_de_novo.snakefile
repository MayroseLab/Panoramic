"""
This pipeline constructs a pan genome
in the "de novo" approach. It consists
of the following general steps:
1. Download fastq files from ena
2. Assemble reads from each sample into
   contigs
3. Annotate each assembly
4. Perform QA and filtration on annotations
5. Cluster genes from all samples to detect
   orthologs
6. Summarize and create PAV matrix
"""
import os
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir) + '/util'
import sys
sys.path.append(utils_dir)
from snakemakeUtils import *

PIPELINE = 'De-novo assembly'

# print Welcome message
ver = get_git_commit()
print("#########################")
print("#### Panoramic %s ####" % ver)
if '-' in ver:
    print("WARNING: you are using an unstable release!\nYou may want to git-checkout to the latest tagged version")
print('Pan-genome construction using the %s pipeline' % PIPELINE)
print("#########################")

# get configfile path
if '--configfile' in sys.argv:
  i = sys.argv.index('--configfile')
elif '--configfiles' in sys.argv:
  i = sys.argv.index('--configfiles')
config_path = os.path.realpath(sys.argv[i+1])

# assert required params are in config
required = ['samples_info_file','hq_genomes_info_file','out_dir','reference_name',
            'reference_genome','reference_proteins','reference_annotation','trimming_modules',
            'merge_min_overlap','merge_max_mismatch_ratio','assembler','min_length','busco_set',
            'transcripts','proteins','min_protein','ppn']
for r in required:
  assert (r in config and config[r]), "Required argument %s is missing or empty in configuration file %s" %(r,config_path)
assemblers = ['spades','minia','megahit']
assert config['assembler'] in assemblers, "'assembler' must be one of: %s" % ', '.join(assemblers)

def init():
    #load_info_file
    config['samples_info'] = SampleInfoReader.sample_table_reader(filename=config['samples_info_file'],
                delimiter='\t', key_name='sample', col_names=['ena_ref'])
    # load HQ genomes info file
    config['hq_info'] = SampleInfoReader.sample_table_reader(filename=config['hq_genomes_info_file'],
                delimiter='\t', key_name='sample', col_names=['annotation_gff', 'genome_fasta', 'proteins_fasta'])
    # convert '_' to '-' in sample names
    config['samples_info'] = {s.replace('_','-'): config['samples_info'][s] for s in config['samples_info']}
    config['hq_info'] = {s.replace('_','-'): config['hq_info'][s] for s in config['hq_info']}
    # ensure not duplicate sample names exist
    all_names = list(config['samples_info'].keys()) + list(config['hq_info'].keys())
    assert len(all_names) == len(set(all_names)), "Can't use duplicate sample names!"
    #  split annotated and unannotated HQ
    config['hq_unannotated'] = {s: config['hq_info'][s] for s in config['hq_info'].keys() if config['hq_info'][s]['annotation_gff'] == '--' and config['hq_info'][s]['proteins_fasta'] == '--'}
    config['hq_info'] = {s: config['hq_info'][s] for s in config['hq_info'].keys() if config['hq_info'][s]['annotation_gff'] != '--' and config['hq_info'][s]['proteins_fasta'] != '--'}


init()
LOGS_DIR = config['out_dir'] + "/logs"
CONDA_ENV_DIR = os.path.dirname(pipeline_dir) + "/conda_env"
annotation_pipeline_dir = os.path.dirname(pipeline_dir) + '/EVM_annotation'
pan_genome_report_dir = os.path.dirname(pipeline_dir) + '/pan_genome_report'
annotation_templates_dir = annotation_pipeline_dir + "/annotation_templates/annotation"
if 'cluster_wrapper' in config and config['cluster_wrapper']:
    cluster_param = '--cluster \"%s\"' % config['cluster_wrapper']
else:
    cluster_param = ''

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

localrules: all_de_novo, prep_annotation_yaml 

rule all_de_novo:
    input:
        pav=config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv",
        cnv=config["out_dir"] + "/all_samples/pan_genome/pan_CNV.tsv",
        prot=config["out_dir"] + "/all_samples/pan_genome/pan_proteome.fasta",
        trans=config["out_dir"] + "/all_samples/pan_genome/pan_transcripts.fasta",
        report=config["out_dir"] + "/all_samples/stats/report.html"

def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['ena_ref']

def get_hq_sample_gff(wildcards):
    return config['hq_info'][wildcards.sample]['annotation_gff']
    

def get_hq_sample_genome(wildcards):
    return config['hq_info'][wildcards.sample]['genome_fasta']

def get_hq_sample_proteins(wildcards):
    return config['hq_info'][wildcards.sample]['proteins_fasta']
    

def get_unannotated_HQ_genome(wildcards):
    return config['hq_unannotated'][wildcards.sample]['genome_fasta']

import pandas as pd

def unannotated_HQ_list(hq_file):
    data = pd.read_csv(hq_file, sep="\t")
    samples = data[(data.annotation_gff == "--") & (data.proteins_fasta == "--")]['sample'].to_list()
    samples = [s.replace('_','-') for s in samples]
    return samples

def annotated_HQ_list(hq_file):
    data = pd.read_csv(hq_file, sep="\t")
    samples = data[~((data.annotation_gff == "--") & (data.proteins_fasta == "--"))]['sample'].to_list()
    samples = [s.replace('_','-') for s in samples]
    return samples

wildcard_constraints:
    sample="[^_]+"

"""
Download reads, preprocess, assemble,
and scaffold all input genomes using
a dedicated pipeline.
"""
if len(config['samples_info']) > 0:
    include: "../genome_assembly/genome_assembly.snakefile"
else:
    rule create_empty_assembly_stats:
        """
        In case only HQ samples are used,
        create an empty stats table
        """
        output:
            assembly_stats_tsv=config["out_dir"] + "/all_samples/stats/assembly_stats.tsv"
        params:
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
        shell:
            """
            echo -e "Assembly\t# contigs (>= 0 bp)\t# contigs (>= 1000 bp)\t# contigs (>= 5000 bp)\t# contigs (>= 10000 bp) # contigs (>= 25000 bp)\t# contigs (>= 50000 bp)\tTotal length (>= 0 bp)\tTotal length (>= 1000 bp)\tTotal length (>= 5000 bp)\tTotal length (>= 10000 bp)\tTotal length (>= 25000 bp)\tTotal length (>= 50000 bp)\t# contigs\tLargest contig\tTotal length\tGC (%%)\tN50\tN75\tL50\tL75\t# total reads\t# left\t# right Mapped (%%)\tProperly paired (%%)\tAvg. coverage depth\tCoverage >= 1x (%%)\t# N's per 100 kbp\t%% Complete BUSCOs\t%% unmapped (Chr0)\tQUAST report\tRead length (bp)" > {output}
            """
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
rule aggregate_annotation_evidence:
    """
    Collect all transcript and protein sequences
    into single files towards annotation
    """
    output:
        trans=config["out_dir"] + "/all_samples/all_transcripts.fasta",
        prot=config["out_dir"] + "/all_samples/all_proteins.fasta"
    params:
        trans_files=config['transcripts'].split(','),
        prot_files=config['proteins'].split(','),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cat {params.trans_files} > {output.trans}
        cat {params.prot_files} > {output.prot}
        """

rule prep_HQ_annotation_yaml:
    input:
        genome=get_unannotated_HQ_genome,
        ref_genome=config['reference_genome'],
        ref_gff=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans.gff',
        ref_cds=config['reference_cds'],
        transcripts=config["out_dir"] + "/all_samples/all_transcripts.fasta",
        proteins=config["out_dir"] + "/all_samples/all_proteins.fasta",
        yml_template=config['annotation_yml_template']
    output:
        config["out_dir"] + "/HQ_samples/{sample}/annotation/annotation.yml"
    params:
        annotation_dir=config["out_dir"] + "/HQ_samples/{sample}/annotation/",
        busco_set=config['busco_set'],
        queue=config['queue'],
        priority=config['priority'],
        ppn_=config['ppn'],
        max_ram_=config['max_ram'],
        max_jobs_=config['max_jobs'] / len(config['samples_info']),
        logs_dir=LOGS_DIR
    shell:
        """
        if [ -s {input.transcripts} ]
        then trans={input.transcripts}
        else trans=''
        fi
        if [ -s {input.proteins} ]
        then prot={input.proteins}
        else prot=''
        fi
        sed -e 's@<INPUT_GENOME>@{input.genome}@' -e 's@<SAMPLE_NAME>@{wildcards.sample}@' -e 's@<OUT_DIR>@{params.annotation_dir}@' -e 's@<REFERENCE_CDS>@{input.ref_cds}@' -e 's@<REFERENCE_LIFTOVER>@1@' -e 's@<REFERENCE_GFF>@{input.ref_gff}@' -e 's@<REFERENCE_FASTA>@{input.ref_genome}@' -e "s@<TRANSCRIPTS_FASTA>@$trans@" -e "s@<PROTEINS_FASTA>@$prot@" -e 's@<QUEUE>@{params.queue}@' -e 's@<PRIORITY>@{params.priority}@' -e 's@<PPN>@{params.ppn_}@' -e 's@<MAX_RAM>@{params.max_ram_}@' -e 's@<MAX_JOBS>@{params.max_jobs_}@' -e 's@<BUSCO_SET>@{params.busco_set}@'  {input.yml_template} > {output}
        """


rule LQ_prep_annotation_yaml:
    """
    Prepare yml config for annotation
    pipeline of each sample
    """
    input:
        genome=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffold.fasta",
        ref_genome=config['reference_genome'],
        ref_gff=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans.gff',
        ref_cds=config['reference_cds'],
        transcripts=config["out_dir"] + "/all_samples/all_transcripts.fasta",
        proteins=config["out_dir"] + "/all_samples/all_proteins.fasta",
        yml_template=config['annotation_yml_template']
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/annotation.yml"
    params:
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/",
        busco_set=config['busco_set'],
        queue=config['queue'],
        priority=config['priority'],
        ppn_=config['ppn'],
        max_ram_=config['max_ram'],
        max_jobs_=config['max_jobs']/len(config['samples_info']),
        logs_dir=LOGS_DIR
    shell:
        """
        if [ -s {input.transcripts} ]
        then trans={input.transcripts}
        else trans=''
        fi
        if [ -s {input.proteins} ]
        then prot={input.proteins}
        else prot=''
        fi
        sed -e 's@<INPUT_GENOME>@{input.genome}@' -e 's@<SAMPLE_NAME>@{wildcards.sample}@' -e 's@<OUT_DIR>@{params.annotation_dir}@' -e 's@<REFERENCE_CDS>@{input.ref_cds}@' -e 's@<REFERENCE_LIFTOVER>@1@' -e 's@<REFERENCE_GFF>@{input.ref_gff}@' -e 's@<REFERENCE_FASTA>@{input.ref_genome}@' -e "s@<TRANSCRIPTS_FASTA>@$trans@" -e "s@<PROTEINS_FASTA>@$prot@" -e 's@<QUEUE>@{params.queue}@' -e 's@<PRIORITY>@{params.priority}@' -e 's@<PPN>@{params.ppn_}@' -e 's@<MAX_RAM>@{params.max_ram_}@' -e 's@<MAX_JOBS>@{params.max_jobs_}@' -e 's@<BUSCO_SET>@{params.busco_set}@' {input.yml_template} > {output}
        """

random_sample = list(config['samples_info'].keys())[0]
random_ref = config['samples_info'][random_sample]['ena_ref']
if 'augustus_dir' not in config:
    config['augustus_dir'] = ''

rule install_EVM_dependencies:
    """
    Use conda/mamba to install
    all dependencies for the
    annotation pipeline. These
    will be used by all invocations
    of the pipeline.
    """
    input:
        config["out_dir"] + "/per_sample/%s/annotation_%s/annotation.yml" %(random_sample, random_ref)
    output:
        config["out_dir"] + '/EVM_dependencies.done'
    params:
        EVM_annotation_snakefile=annotation_pipeline_dir + '/EVM_annotation.snakefile',
        out_dir=config['out_dir'],
        augustus_dir=config['augustus_dir'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        cd {params.out_dir}
        snakemake -s {params.EVM_annotation_snakefile} --configfile {input} -j 1 --use-conda --conda-create-envs-only
        if [ ! -z "{params.augustus_dir}" ]; then
            augustus_conda=$(grep "name: augustus" {params.out_dir}/.snakemake/conda/*.yaml | sed 's/\.yaml.*//')
            cp -r {params.augustus_dir} $augustus_conda/config/species/
        fi
        touch {output}
        """

rule HQ_EVM_annotation:
    input:
        conf=config["out_dir"] + "/HQ_samples/{sample}/annotation/annotation.yml",
        dep=config["out_dir"] + '/EVM_dependencies.done'
    output:
        config["out_dir"] + "/HQ_samples/{sample}/annotation/EVM.filter.rename.gff3",
        config["out_dir"] + "/HQ_samples/{sample}/annotation/EVM.filter.rename.prot.fasta",
        config["out_dir"] + "/HQ_samples/{sample}/annotation/EVM.filter.rename.trans.fasta"
    params:
        EVM_annotation_snakefile=annotation_pipeline_dir + '/EVM_annotation.snakefile',
        queue=config['queue'],
        jobs=int(config['max_jobs']/len(config['samples_info'])),
        annotation_dir=config["out_dir"] + "/HQ_samples/{sample}/annotation/",
        cluster_param=cluster_param,
        priority=config['priority'],
        jobscript=utils_dir + '/jobscript.sh',
        dependencies_dir=config['out_dir'] + '/.snakemake/conda',
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        cd {params.annotation_dir}
        snakemake -s {params.EVM_annotation_snakefile} --configfile {input.conf} {params.cluster_param} -j {params.jobs} --latency-wait 60 --restart-times 3 --jobscript {params.jobscript} --keep-going -p --use-conda --ignore-incomplete --conda-prefix {params.dependencies_dir}
        """

rule LQ_EVM_annotation:
    """
    Run the EVM annotation pipeline,
    including liftover, ab-initio and
    evidence-based annotation
    """
    input:
        conf=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/annotation.yml",
        dep=config["out_dir"] + '/EVM_dependencies.done'
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/EVM.filter.rename.gff3",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/EVM.filter.rename.prot.fasta",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/EVM.filter.rename.trans.fasta"
    params:
        EVM_annotation_snakefile=annotation_pipeline_dir + '/EVM_annotation.snakefile',
        queue=config['queue'],
        jobs=int(config['max_jobs']/len(config['samples_info'])),
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}",
        cluster_param=cluster_param,
        priority=config['priority'],
        jobscript=utils_dir + '/jobscript.sh',
        dependencies_dir=config['out_dir'] + '/.snakemake/conda',
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        cd {params.annotation_dir}
        snakemake -s {params.EVM_annotation_snakefile} --configfile {input.conf} {params.cluster_param} -j {params.jobs} --latency-wait 60 --restart-times 3 --jobscript {params.jobscript} --keep-going -p --use-conda --ignore-incomplete --conda-prefix {params.dependencies_dir}
        """

rule remove_ref_alt_splicing:
    """
    In case the reference annotation
    contains genes with multiple mRNAs,
    only keep the longest transcript.
    """
    input:
        ref_gff=config['reference_annotation']
    output:
        config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans.gff'
    params:
        longest_trans_script=utils_dir + '/remove_alt_splicing_from_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.longest_trans_script} {input.ref_gff} {output}
        """

rule get_ref_proteins:
    """
    Filter reference proteins according
    to filtered gff and put the new file
    in orthofinder dir.
    """
    input:
        trans=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans.fasta',
        fasta=config['reference_proteins'],
        gff=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans.gff'
    output:
        config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta'
    params:
        filter_fasta_script=utils_dir + '/filter_fasta_by_gff.py',
        min_protein=config['min_protein'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.filter_fasta_script} {input.gff} {input.fasta} {output} mRNA ID --min_len {params.min_protein}
        """

rule make_ref_blast_db:
    """
    Create blast DB of ref proteins
    """
    input:
        config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta'
    output:
        config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta.pin',
        config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta.phr',
        config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta.psq'
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/GAWN.yml'
    shell:
        """
        makeblastdb -in {input} -input_type fasta -dbtype prot
        """

rule remove_hq_sample_alt_splicing:
    """
    For HQ genome samples, keep only
    the longest mRNA per gene
    """
    input:
        get_hq_sample_gff
    output:
        config["out_dir"] + "/HQ_samples/{sample}/{sample}_longest_trans.gff"
    params:
        longest_trans_script=utils_dir + '/remove_alt_splicing_from_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.longest_trans_script} {input} {output}
        """

rule get_hq_proteins:
    """
    Filter HQ samples proteins according
    to filtered gff and put the new file
    in orthofinder dir.
    """
    input:
        fasta=get_hq_sample_proteins,
        gff=config["out_dir"] + "/HQ_samples/{sample}/{sample}_longest_trans.gff"
    output:
        config["out_dir"] + "/HQ_samples/{sample}/{sample}_HQ.fasta"
    params:
        filter_fasta_script=utils_dir + '/filter_fasta_by_gff.py',
        min_protein=config['min_protein'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.filter_fasta_script} {input.gff} {input.fasta} {output} mRNA ID --min_len {params.min_protein}
        """

rule get_hq_sample_transcripts:
    """
    Get HQ samples transcripts according
    to filtered gff.
    """
    input:
        fasta=get_hq_sample_genome,
        gff=config["out_dir"] + "/HQ_samples/{sample}/{sample}_longest_trans.gff"
    output:
        config["out_dir"] + "/HQ_samples/{sample}/{sample}_longest_trans.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffread.yml'
    shell:
        """
        gffread {input.gff} -g {input.fasta} -w {output}
        """

rule get_ref_transcripts:
    """
    Filter reference transcripts according
    to filtered gff. These transcripts will
    be used for liftover.
    """
    input:
        fasta=config['reference_transcripts'],
        gff=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans.gff'
    output:
        config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans.fasta'
    params:
        filter_fasta_script=utils_dir + '/filter_fasta_by_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.filter_fasta_script} {input.gff} {input.fasta} {output} mRNA ID
        """

rule copy_LQ_to_orthofinder:
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/EVM.filter.rename.prot.fasta"
    output:
        config["out_dir"] + "/all_samples/orthofinder/{sample}_{ena_ref}.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        ln -s {input} {output}
        """

rule copy_annotated_HQ_to_orthofinder:
    input:
        config["out_dir"] + "/HQ_samples/{sample}/{sample}_HQ.fasta"
    output:
        config["out_dir"] + "/all_samples/orthofinder/{sample}.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        ln -s {input} {output}
        """

rule copy_unannotated_HQ_to_orthofinder:
    input:
        config["out_dir"] + "/HQ_samples/{sample}/annotation/EVM.filter.rename.prot.fasta" 
    output:
        config["out_dir"] + "/all_samples/orthofinder/{sample}.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        ln -s {input} {output}
        """

rule orthofinder:
    """
    Run OrthoFinder2 on all proteins
    from all annotated genomes to get
    initial orthogroups
    """
    input:
        lq = expand(config["out_dir"] + "/all_samples/orthofinder/{sample}_{ena_ref}.fasta",zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        annotated_hq = expand(config["out_dir"] + "/all_samples/orthofinder/{sample}.fasta",sample=annotated_HQ_list(config['hq_genomes_info_file'])),
        unannotated_hq = expand(config["out_dir"] + "/all_samples/orthofinder/{sample}.fasta",sample=unannotated_HQ_list(config['hq_genomes_info_file'])),
        ref=config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta'
    output:
        config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups/Orthogroups.tsv"
    params:
        orthofinder_dir=config["out_dir"] + "/all_samples/orthofinder",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/orthofinder2.yml'
    shell:
        """
        if [ -d {params.orthofinder_dir}/OrthoFinder/Results_orthofinder/  ]
        then
            rm -rf {params.orthofinder_dir}/OrthoFinder/Results_orthofinder/
        fi
        orthofinder -t {params.ppn} -a {params.ppn} -S mmseqs -n orthofinder -f {params.orthofinder_dir}
        """

rule break_orthogroups_MWOP:
    """
    Use orthofinder's ouput and further
    break orthogroups into smaller
    clusters representing single genes.
    This is done using the Maximum Weight
    Orthogonal Partitions (MWOP) algorithm.
    """
    input:
        config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups/Orthogroups.tsv"
    output:
        config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups_break_MWOP.tsv"
    params:
        orthofinder_dir=config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder",
        mwop_script=os.path.join(pipeline_dir,"break_OrthoFinder_clusters.py"),
        ref_genome_name=config['reference_name'] + "_REF",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/break_orthogroups.yml'
    shell:
        """
        python {params.mwop_script} {params.orthofinder_dir} bitscore --allow_gene_copies exclude_ref --ref_genome_name {params.ref_genome_name} {output} --cpus {params.ppn}
        """

rule create_PAV_matrix:
    """
    Create the final PAV and CNV matrices
    """
    input:
        config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups_break_MWOP.tsv"
    output:
        pav=config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv",
        cnv=config["out_dir"] + "/all_samples/pan_genome/pan_CNV.tsv",
        mapping=config["out_dir"] + "/all_samples/pan_genome/OG_to_gene_names.tsv"
    params:
        create_pav_mat_script=os.path.join(pipeline_dir,"create_PAV_matrix.py"),
        ref_name=config['reference_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/break_orthogroups.yml'
    shell:
        """
        python {params.create_pav_mat_script} {input} {params.ref_name} {output.pav} {output.cnv} {output.mapping}
        """

rule create_all_proteins_fasta:
    """
    Aggregate all proteins from all samples
    into one fasta file
    """
    input:
        lq=expand(config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/EVM.filter.rename.prot.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        hq=expand(config["out_dir"] + "/HQ_samples/{sample}/{sample}_HQ.fasta", sample=config['hq_info'].keys()),
        ref=config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta'
    output:
        config["out_dir"] + "/all_samples/pan_genome/all_proteins.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cat {input.ref} {input.hq} {input.lq} | tr ':' '_' > {output}
        """

rule create_all_transcripts_fasta:
    """
    Aggregate all transcripts from all samples
    into one fasta file
    """
    input:
        lq=expand(config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/EVM.filter.rename.trans.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        hq=expand(config["out_dir"] + "/HQ_samples/{sample}/{sample}_longest_trans.fasta", sample=config['hq_info'].keys()),
        ref=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans.fasta'
    output:
        config["out_dir"] + "/all_samples/pan_genome/all_transcripts.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cat {input.ref} {input.hq} {input.lq} | tr ':' '_' > {output}
        """

rule choose_representatives:
    """
    Select one representative gene
    from each orthogroup
    """
    input:
        mapping=config["out_dir"] + "/all_samples/pan_genome/OG_to_gene_names.tsv",
        mwop=config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups_break_MWOP.tsv"
    output:
        config["out_dir"] + "/all_samples/pan_genome/OG_to_gene_names_with_representatives.tsv"
    params:
        choose_representatives_script=os.path.join(pipeline_dir,"choose_representatives.py"),
        og_seq_dir=config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroup_Sequences",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/biopython.yml'
    shell:
        """
        python {params.choose_representatives_script} {params.og_seq_dir} {input.mapping} {input.mwop} {output}
        """

rule create_pan_proteome_fasta:
    """
    Create a fasta file with one
    representative protein sequence
    per pan gene.
    """
    input:
        represent=config["out_dir"] + "/all_samples/pan_genome/OG_to_gene_names_with_representatives.tsv",
        prot=config["out_dir"] + "/all_samples/pan_genome/all_proteins.fasta"
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_proteome.fasta"
    params:
        filter_script=utils_dir + '/filter_fasta_by_id_list.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/biopython.yml'
    shell:
        """
        awk '{{print $4"\t"$2}}' {input.represent} | python {params.filter_script} {input.prot} > {output}
        """  

rule create_pan_transcripts_fasta:
    """
    Create a fasta file with one
    representative transcript sequence
    per pan gene.
    """
    input:
        represent=config["out_dir"] + "/all_samples/pan_genome/OG_to_gene_names_with_representatives.tsv",
        trans=config["out_dir"] + "/all_samples/pan_genome/all_transcripts.fasta"
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_transcripts.fasta"
    params:
        filter_script=utils_dir + '/filter_fasta_by_id_list.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/biopython.yml'
    shell:
        """
        awk '{{print $4"\t"$2}}' {input.represent} | python {params.filter_script} {input.trans} > {output}
        """

rule calculate_stepwise_stats:
    """
    Perform genome stepwise addition
    analysis and generate stats
    """
    input:
        config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv"
    output:
        config["out_dir"] + "/all_samples/stats/stepwise_stats.tsv"
    params:
        stepwise_script=pan_genome_report_dir + '/calc_stepwise_stats.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/jupyter.yml'
    shell:
        """
        python {params.stepwise_script} {input} 100 {output}
        """

rule create_report_notebook:
    """
    Create pan genome report Jupyter notebook
    """
    input:
        pav_tsv=config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv",
        stepwise_tsv=config["out_dir"] + "/all_samples/stats/stepwise_stats.tsv",
        assembly_stats_tsv=config["out_dir"] + "/all_samples/stats/assembly_stats.tsv",
        proteins_fasta=config["out_dir"] + "/all_samples/pan_genome/pan_proteome.fasta"
    output:
        config["out_dir"] + "/all_samples/stats/report.ipynb"
    params:
        ref_name=config['reference_name'],
        conf=config_path,
        nb_template=pan_genome_report_dir + '/report_template.ipynb',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        sed -e 's|<PIPELINE>|Panoramic de-novo|' -e 's|<CONF>|{params.conf}|' -e 's|<PAV_TSV>|{input.pav_tsv}|' -e 's|<SYEPWISE_TSV>|{input.stepwise_tsv}|' -e 's|<REF_NAME>|{params.ref_name}|' -e 's|<PROT_FASTA>|{input.proteins_fasta}|' -e 's|<STATS_TSV>|{input.assembly_stats_tsv}|' {params.nb_template} > {output}
        """


rule create_report_html:
    """
    Run notebook and create pan
    genome report HTML
    """
    input:
        config["out_dir"] + "/all_samples/stats/report.ipynb"
    output:
        config["out_dir"] + "/all_samples/stats/report.html"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/jupyter.yml'
    shell:
        """
        jupyter nbconvert {input} --output {output} --to html --no-prompt --no-input --execute --NotebookClient.timeout=-1 --ExecutePreprocessor.timeout=-1
        """

