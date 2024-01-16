"""
This pipeline constructs a pan genome
in the "map-to-pan" approach. It consists
of the following general steps:
1. Download fastq files from ena
2. Assemble reads from each sample into
   contigs
3. Iteratively map contigs to reference
   and adding novel sequences to create the
   non-reference section of the pan genome
4. Annotate the non-reference contigs
5. Align reads from each sample to reference
   + non-reference genes to determine gene
   presence/absence
6. Summarize and create PAV matrix
"""

import os
mtp_pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(mtp_pipeline_dir) + '/util'
import sys
sys.path.append(utils_dir)
from snakemakeUtils import *
from collections import OrderedDict

PIPELINE = 'Map-to-pan'

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
required = ["samples_info_file","hq_genomes_info_file","out_dir","reference_name","reference_genome",
            "reference_annotation","reference_proteins","id_simplify_function","trimming_modules",
            "merge_min_overlap","merge_max_mismatch_ratio","assembler","min_length","busco_set",
            "transcripts","proteins","min_protein",
            "similarity_threshold_proteins","HQ_min_cov","LQ_min_cov","min_read_depth","ppn"]
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
                delimiter='\t', key_name='sample', col_names=['annotation_gff','genome_fasta', 'proteins_fasta'])
    # convert '_' to '-' in sample names
    config['samples_info'] = {s.replace('_','-'): config['samples_info'][s] for s in config['samples_info']}
    config['hq_info'] = {s.replace('_','-'): config['hq_info'][s] for s in config['hq_info']}
    # ensure not duplicate sample names exist
    all_names = list(config['samples_info'].keys()) + list(config['hq_info'].keys())
    assert len(all_names) == len(set(all_names)), "Can't use duplicate sample names!"
    
init()
config['samples_info'] = OrderedDict(config['samples_info'])
LOGS_DIR = config['out_dir'] + "/logs"
CONDA_ENV_DIR = os.path.dirname(mtp_pipeline_dir) + '/conda_env'
annotation_mtp_pipeline_dir = os.path.dirname(mtp_pipeline_dir) + '/EVM_annotation'
pan_genome_report_dir = os.path.dirname(mtp_pipeline_dir) + '/pan_genome_report'
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

localrules: all, prep_annotation_chunks_tsv, prep_annotation_yaml, calculate_n_chunks

n_samples = len(config['samples_info'])
last_sample = list(config['samples_info'].keys())[-1]
last_sample_ena = config['samples_info'][last_sample]['ena_ref']
rule all_map_to_pan:
    input:
        config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv",
        config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta",
        config["out_dir"] + "/all_samples/pan_genome/pan_proteome.fasta",
        config["out_dir"] + "/all_samples/pan_genome/pan_genes.gff",
        config["out_dir"] + "/all_samples/pan_genome/pan_transcripts.fasta",
        config["out_dir"] + "/all_samples/stats/report.html"

def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['ena_ref']

def get_hq_sample_gff(wildcards):
    return config['hq_info'][wildcards.sample]['annotation_gff']

def get_hq_sample_genome(wildcards):
    return config['hq_info'][wildcards.sample]['genome_fasta']

def get_hq_sample_proteins(wildcards):
    return config['hq_info'][wildcards.sample]['proteins_fasta']

wildcard_constraints:
    sample="[^_]+"

"""
Download reads, preprocess, assemble,
and scaffold all input genomes using
a dedicated pipeline.
"""
include: "../genome_assembly/genome_assembly.snakefile"

rule simplify_ref_gff_ID:
    """
    Edit ID (and Parent) attributes of ref gff
    To make them simpler (if needed)
    """
    input:
        config['reference_annotation']
    output:
        config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_simp.gff'
    params:
        simplify_gff_script=utils_dir + '/simplify_gff_id.py',
        simp_function=config['id_simplify_function'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.simplify_gff_script} {input} {output} "{params.simp_function}"
        """

rule remove_ref_alt_splicing:
    """
    In case the reference annotation
    contains genes with multiple mRNAs,
    only keep the longest transcript.
    Also, remove short proteins.
    """
    input:
        gff=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_simp.gff',
        prot_fasta=config['reference_proteins']
    output:
        gff=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff',
        table=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff.gene_to_mRNA'
    params:
        longest_trans_script=utils_dir + '/remove_alt_splicing_from_gff.py',
        min_protein=config['min_protein'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.longest_trans_script} {input.gff} {output.gff} {input.prot_fasta} {params.min_protein} ID
        """

rule iterative_map_to_pan_HQ:
    """
    Iteratively create HQ pan genome
    by mapping HQ samples to ref and
    adding novel sequences and genes.
    """
    input:
       samples=config['hq_genomes_info_file'],
       ref_genome=config['reference_genome'],
       ref_gff=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff',
       ref_proteins=config['reference_proteins']
    output:
       config["out_dir"] + "/HQ_samples/HQ_pan/pan_genome.fasta",
       config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes.gff",
       config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins.fasta"
    params:
        map_to_pan_script=mtp_pipeline_dir + '/iterative_map_to_pan.py',
        min_len=config['min_length'],
        min_protein=config['min_protein'],
        out_dir=config["out_dir"] + "/HQ_samples/HQ_pan",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/iterative_map_to_pan.yml'
    shell:
        """
        python {params.map_to_pan_script} {input.ref_genome} {input.ref_gff} {input.ref_proteins} {input.samples} {params.out_dir} --cpus {params.ppn} --min_len {params.min_len} --min_protein {params.min_protein}
        """

rule remove_alt_splicing_from_HQ_pan:
    """
    Remove splice variants from HQ pan gff.
    Only keep longest variant. Then match
    the proteins fasta
    """
    input:
        gff=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes.gff",
        prot=config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins.fasta"
    output:
        gff=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes_no_alt_splicing.gff3",
        prot=config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins_no_alt_splicing.fasta",
        mrna_to_gene=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes_no_alt_splicing.gff3" + '.gene_to_mRNA'
    params:
        longest_trans_script=utils_dir + '/remove_alt_splicing_from_gff.py',
        filter_fasta_script=utils_dir + '/filter_fasta_by_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.longest_trans_script} {input.gff} {output.gff}
        python {params.filter_fasta_script} {output.gff} {input.prot} {output.prot} mRNA ID
        """

rule prep_tsv_for_LQ_samples:
    """
    Prepare the TSV file required as
    input for creating the pan genome
    from assembled LQ samples
    """
    input:
        expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffold.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])
    output:
        config["out_dir"] + "/all_samples/pan_genome/samples.tsv"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        echo "sample\tgenome_fasta" > {output}
        echo "{input}" | tr ' ' '\n' | awk '{{split($0,a,"/"); print a[length(a)-3]"\t"$0}}' >> {output}
        """

rule iterative_map_to_pan_LQ:
    """
    Iteratively create LQ pan genome
    by mapping LQ samples to ref+HQ and
    adding novel sequences.
    """
    input:
       samples=config["out_dir"] + "/all_samples/pan_genome/samples.tsv",
       ref_genome=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genome.fasta",
       ref_gff=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes_no_alt_splicing.gff3",
       ref_proteins=config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins_no_alt_splicing.fasta"
    output:
       config["out_dir"] + "/all_samples/pan_genome/all_novel.fasta",
       config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta"
    params:
        map_to_pan_script=mtp_pipeline_dir + '/iterative_map_to_pan.py',
        min_len=config['min_length'],
        out_dir=config["out_dir"] + "/all_samples/pan_genome/",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/iterative_map_to_pan.yml'
    shell:
        """
         python {params.map_to_pan_script} {input.ref_genome} {input.ref_gff} {input.ref_proteins} {input.samples} {params.out_dir} --cpus {params.ppn} --min_len {params.min_len} --min_gap_break 100
        """

rule novel_to_chunks:
    """
    Divide novel sequences to
    chunks of concatenated
    contigs
    """
    input:
        config["out_dir"] + "/all_samples/pan_genome/all_novel.fasta"
    output:
        fasta=config["out_dir"] + "/all_samples/pan_genome/all_novel_chunks.fasta",
        bed=config["out_dir"] + "/all_samples/pan_genome/all_novel_chunks.fasta.chunks.bed"
    params:
        chunks_script=utils_dir + '/fasta_to_chunks.py',
        chunk_size=config['chunk_size'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/biopython.yml'
    shell:
        """
        python {params.chunks_script} {input} {output.fasta} {params.chunk_size}
        """


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

rule prep_annotation_yaml:
    """
    Prepare yml config for annotation
    of the novel sequences
    """
    input:
        novel=config["out_dir"] + "/all_samples/pan_genome/all_novel_chunks.fasta",
        transcripts=config["out_dir"] + "/all_samples/all_transcripts.fasta",
        proteins=config["out_dir"] + "/all_samples/all_proteins.fasta",
        yml_template=config['annotation_yml_template']
    output:
        config["out_dir"] + "/all_samples/annotation/annotation.yml"
    params:
        annotation_dir=config["out_dir"] + "/all_samples/annotation/",
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
        sed -e 's@<INPUT_GENOME>@{input.novel}@' -e 's@<SAMPLE_NAME>@novel_sequences@' -e 's@<OUT_DIR>@{params.annotation_dir}@' -e 's@<REFERENCE_CDS>@@' -e 's@<REFERENCE_LIFTOVER>@0@' -e 's@<REFERENCE_GFF>@@' -e 's@<REFERENCE_FASTA>@@' -e "s@<TRANSCRIPTS_FASTA>@$trans@" -e "s@<PROTEINS_FASTA>@$prot@" -e 's@<QUEUE>@{params.queue}@' -e 's@<PRIORITY>@{params.priority}@' -e 's@<PPN>@{params.ppn_}@' -e 's@<MAX_RAM>@{params.max_ram_}@' -e 's@<MAX_JOBS>@{params.max_jobs_}@' {input.yml_template} > {output}
        """
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
        config["out_dir"] + "/all_samples/annotation/annotation.yml"
    output:
        config["out_dir"] + '/EVM_dependencies.done'
    params:
        EVM_annotation_snakefile=annotation_mtp_pipeline_dir + '/EVM_annotation.snakefile',
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

rule EVM_annotation:
    """
    Run the EVM annotation pipeline
    on novel sequences, including
    ab-initio and evidence-based annotation.
    """
    input:
        conf=config["out_dir"] + "/all_samples/annotation/annotation.yml",
        dep=config["out_dir"] + '/EVM_dependencies.done'
    output:
        config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.gff3",
        config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.prot.fasta",
        config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.trans.fasta"
    params:
        EVM_annotation_snakefile=annotation_mtp_pipeline_dir + '/EVM_annotation.snakefile',
        queue=config['queue'],
        jobs=int(config['max_jobs']/len(config['samples_info'])),
        annotation_dir=config["out_dir"] + "/all_samples/annotation",
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

rule transform_coordinates:
    """
    Transform annotation GFF from chunk
    coordinates to original novel contig
    coordinates.
    This also discards gene models spanning
    multiple concatenated contigs.
    """
    input:
        gff=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.gff3",
        bed=config["out_dir"] + "/all_samples/pan_genome/all_novel_chunks.fasta.chunks.bed"
    output:
        config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.transform.gff3"
    params:
        transform_script=utils_dir + '/transform_gff_chunk_coordinates.py',
        jobscript=utils_dir + '/jobscript.sh',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/iterative_map_to_pan.yml'
    shell:
        """
        python {params.transform_script} {input.gff} {input.bed} {output}
        """

rule match_proteins:
    """
    Remove proteins discarded
    during coordinate transformation
    """
    input:
        prot=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.prot.fasta",
        gff=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.transform.gff3"
    output:
        config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.prot.no_conc.fasta"
    params:
        filter_fasta_script=utils_dir + '/filter_fasta_by_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.filter_fasta_script} {input.gff} {input.prot} {output} mRNA ID 
        """


rule remove_redundant_proteins:
    """
    Remove redundant non-ref proteins
    by clustering with CD-HIT and taking
    longest sequence from each cluster.
    """
    input:
        config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.prot.no_conc.fasta"
    output:
        config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.prot.no_redun.fasta"
    params:
        similarity_threshold=config['similarity_threshold_proteins'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn'],
    conda:
        CONDA_ENV_DIR + '/cd-hit.yml'
    shell:
        """
        cd-hit -i {input} -o {output} -c {params.similarity_threshold} -n 5 -M 0 -d 0 -T {params.ppn}
        """

rule match_gff:
    """
    Filter genes gff according to
    the genes remaining in proteins
    fasta after all filtrations.
    """
    input:
        prot=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.prot.no_redun.fasta",
        gff=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.transform.gff3"
    output:
        gff=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.no_redun.gff3",
        filter_list=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.no_redun.list"
    params:
        filter_gff_script=utils_dir + '/filter_gff_by_id_list.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        if [ -s {input.prot} ]
        then
            grep '>' {input.prot} | tr -d '>' | awk '{{print $1}}' > {output.filter_list}
            grep '>' {input.prot} | tr -d '>' | awk '{{print $1}}' | sed 's/novel_sequences_evm.model/novel_sequences_evm.TU/' >> {output.filter_list}
            python {params.filter_gff_script} {input.gff} {output.filter_list} {output.gff}
        else
            touch {output.filter_list}
            echo '##gff-version 3' > {output.gff}
       fi
        """

rule rename_pan_genes:
    """
    Assign PanGene names to
    gene models in novel sequences
    """
    input:
        prot=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.prot.no_redun.fasta",
        gff=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.no_redun.gff3"
    output:
        prot=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.prot.no_redun.PanGene.fasta",
        gff=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.no_redun.PanGene.gff3",
        id_map=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.no_redun.PanGene.gff3.id_map.tsv"
    params:
        rename_gff_script=utils_dir + '/rename_gff_features.py',
        rename_fasta_script=utils_dir + '/rename_fasta_records.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.rename_gff_script} {input.gff} {output.gff} PanGene
        python {params.rename_fasta_script} {input.prot} {output.id_map} {output.prot}
        """

rule create_pan_proteome:
    """
    Create pan genome proteins fasta
    """
    input:
        non_ref_proteins=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.prot.no_redun.PanGene.fasta",
        ref_plus_hq_proteins=config["out_dir"] + "/HQ_samples/HQ_pan/pan_proteins_no_alt_splicing.fasta"
    output:
        pan_proteome=config["out_dir"] + "/all_samples/pan_genome/pan_proteome.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        cat {input.ref_plus_hq_proteins} {input.non_ref_proteins} > {output.pan_proteome}
        """
        
rule create_pan_annotation:
    """
    Create pan genome annotation gff
    """
    input:
        non_ref_gff=config["out_dir"] + "/all_samples/annotation/EVM.filter.rename.no_redun.PanGene.gff3",
        ref_plus_hq_gff=config["out_dir"] + "/HQ_samples/HQ_pan/pan_genes_no_alt_splicing.gff3"
    output:
        pan_genes=config["out_dir"] + "/all_samples/pan_genome/pan_genes.gff"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        cat {input.ref_plus_hq_gff} {input.non_ref_gff} > {output.pan_genes}
        """

rule index_pan_genome:
    """
    Index pan genome for BWA runs
    """
    input:
        config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta"
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta.bwt"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/bwa.yml'
    shell:
        """
        bwa index {input}
        """

rule map_reads_to_pan:
    """
    Map reads from each sample to
    the pan genome.
    """
    input:
        r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz",
        r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz",
        pan_genome=config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta",
        pan_genome_index=config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta.bwt"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sam"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/bwa.yml'
    shell:
        """
        bwa mem -t {params.ppn} {input.pan_genome} {input.r1_paired} {input.r2_paired} > {output}
        """

rule sam_to_sorted_bam:
    """
    Convert SAM output to sorted BAM
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sam"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools view -@ {params.ppn} -bh {input} | samtools sort -@ {params.ppn} - -o {output}
        """

rule index_bam:
    input:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam.bai"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools index {input}
        """

rule map_HQ_genome_to_pan:
    """
    Map HQ genome sequence to pan
    genome in order to detect gene PAV
    """
    input:
        pan_genome=config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta",
        hq_genome=get_hq_sample_genome
    output:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sam"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/minimap2.yml'
    shell:
        """
        minimap2 -ax asm5 -t {params.ppn} {input.pan_genome} {input.hq_genome} -L > {output}
        """

rule HQ_sam_to_sorted_bam:
    """
    Convert SAM output for HQ samples to sorted BAM
    """
    input:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sam" 
    output:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sort.bam"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools view -@ {params.ppn} -bh {input} | samtools sort -@ {params.ppn} - -o {output}
        """

rule HQ_index_bam:
    input:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sort.bam"
    output:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sort.bam.bai"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools index {input}
        """

rule pan_genes_gff_to_bed:
    """
    Convert pan mRNA features from GFF3
    to bed for easy handling
    """
    input:
        config["out_dir"] + "/all_samples/pan_genome/pan_genes.gff"
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_mRNA.bed"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/bedops.yml'
    shell:
        """
        awk '$3 == "mRNA"' {input} | gff2bed | cut -f1,2,3,4 |  sort -k1,1 -k2,2n > {output}
        """

rule calculate_HQ_mRNA_cov:
    """
    Calculate coverage of HQ genome
    sequences on pan genes.
    """
    input:
        bed=config["out_dir"] + "/all_samples/pan_genome/pan_mRNA.bed",
        bam=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.sort.bam"
    output:
        order=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}.order",
        sorted_bed=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/pan_mRNA.sort.bed",
        filtered_bam=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.filter.sort.bam",
        cov=config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.bedCovHist"
    params:
        sort_script=utils_dir + "/custom_sort_bed.py",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/bedtools.yml'
    shell:
        """
        samtools view {input.bam} | cut -f3 | uniq > {output.order}
        python {params.sort_script} {input.bed} {output.order} > {output.sorted_bed}
        bedtools intersect -a {input.bam} -b {input.bed} -wa > {output.filtered_bam}
        bedtools coverage -a {output.sorted_bed} -b {output.filtered_bam} -hist -sorted > {output.cov}
        """

rule calculate_LQ_mRNA_cov:
    """
    Calculate coverage of LQ reads
    on pan genes.
    """
    input:
        bed=config["out_dir"] + "/all_samples/pan_genome/pan_mRNA.bed",
        bam=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.sort.bam"
    output:
        order=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}.order",
        sorted_bed=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/pan_mRNA.sort.bed",
        filtered_bam=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.filter.sort.bam",
        cov=config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.bedCovHist"
    params:
        sort_script=utils_dir + "/custom_sort_bed.py",
        queue=config['queue'],
        priority=config['priority'],
        ppn=2,
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/bedtools.yml'
    shell:
        """
        samtools view {input.bam} | cut -f3 | uniq > {output.order}
        python {params.sort_script} {input.bed} {output.order} > {output.sorted_bed}
        bedtools intersect -a {input.bam} -b {input.bed} -wa > {output.filtered_bam}
        bedtools coverage -a {output.sorted_bed} -b {output.filtered_bam} -hist -sorted > {output.cov}
        """

rule detect_HQ_PAV:
    """
    Determine per gene per sample PA
    based on sequence coverage
    """
    input:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_vs_pan.bedCovHist"
    output:
        config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_PAV.tsv"
    params:
        detect_pav_script=mtp_pipeline_dir + '/detect_gene_PAV_from_bed_cov.py',
        min_frac_covered=config['HQ_min_cov'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.detect_pav_script} {input} 1 {params.min_frac_covered} {wildcards.sample} > {output}
        """

rule detect_LQ_PAV:
    """
    Determine per gene per sample PA
    based on sequence coverage
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_map_to_pan.bedCovHist"
    output:
        config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_PAV.tsv"
    params:
        detect_pav_script=mtp_pipeline_dir + '/detect_gene_PAV_from_bed_cov.py',
        min_depth=config['min_read_depth'],
        min_frac_covered=config['LQ_min_cov'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.detect_pav_script} {input} {params.min_depth} {params.min_frac_covered} {wildcards.sample}_{wildcards.ena_ref} > {output}
        """

rule create_ref_genes_list:
    """
    Create a file containing the
    reference mRNAs list
    """
    input:
        config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff.gene_to_mRNA'
    output:
        config["out_dir"] + "/all_samples/ref/ref_genes.list"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cut -f2 {input} > {output}
        """  

rule create_pan_PAV:
    """
    Combine HQ and LQ PAV tables
    and add ref column to create
    the final PAV matrix
    """
    input:
        lq=expand(config["out_dir"] + "/per_sample/{sample}/map_to_pan_{ena_ref}/{ena_ref}_PAV.tsv", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        hq=expand(config["out_dir"] + "/HQ_samples/{sample}/map_to_pan/{sample}_PAV.tsv", sample=config['hq_info'].keys()),
        ref_genes=config["out_dir"] + "/all_samples/ref/ref_genes.list",
        names_sub=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans_simp.gff.gene_to_mRNA'
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv"
    params:
        create_PAV_matrix_script=mtp_pipeline_dir + '/create_PAV_matrix.py',
        ref_name=config['reference_name'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/pandas.yml'
    shell:
        """
        python {params.create_PAV_matrix_script} {input.hq} {input.lq} {params.ref_name} {input.ref_genes} {input.names_sub} {output}
        """

rule create_pan_transcripts_fasta:
    """
    Create a fasta file with
    pan-genome transcripts
    """
    input:
        fasta=config["out_dir"] + "/all_samples/pan_genome/pan_genome.fasta",
        gff=config["out_dir"] + "/all_samples/pan_genome/pan_genes.gff"
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_transcripts.fasta"
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
        CONDA_ENV_DIR + '/pandas.yml'
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
        proteins_fasta=config["out_dir"] + "/all_samples/pan_genome/pan_proteome.fasta",
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
        sed -e 's|<PIPELINE>|Panoramic map-to-pan|' -e 's|<CONF>|{params.conf}|' -e 's|<PAV_TSV>|{input.pav_tsv}|' -e 's|<SYEPWISE_TSV>|{input.stepwise_tsv}|' -e 's|<REF_NAME>|{params.ref_name}|' -e 's|<PROT_FASTA>|{input.proteins_fasta}|' -e 's|<STATS_TSV>|{input.assembly_stats_tsv}|' {params.nb_template} > {output}
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
