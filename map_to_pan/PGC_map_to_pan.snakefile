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
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir) + '/util'
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
i = sys.argv.index('--configfile')
config_path = os.path.realpath(sys.argv[i+1])

# assert required params are in config
required = ["samples_info_file","hq_genomes_info_file","out_dir","reference_name","reference_genome",
            "reference_annotation","reference_proteins","id_simplify_function","trimming_modules",
            "merge_min_overlap","merge_max_mismatch_ratio","assembler","min_length","busco_set",
            "repeats_library","transcripts","proteins","augustus_species","min_protein","max_aed",
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
CONDA_ENV_DIR = os.path.dirname(pipeline_dir) + '/conda_env'
annotation_pipeline_dir = os.path.dirname(pipeline_dir) + '/genome_annotation'
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
        map_to_pan_script=pipeline_dir + '/iterative_map_to_pan.py',
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
        expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffolds.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])
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
        map_to_pan_script=pipeline_dir + '/iterative_map_to_pan.py',
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
         python {params.map_to_pan_script} {input.ref_genome} {input.ref_gff} {input.ref_proteins} {input.samples} {params.out_dir} --cpus {params.ppn} --min_len {params.min_len}
        """

rule calculate_n_chunks:
    """
    Calculate the optimal number of
    chunks to split the novel sequences
    to, towards annotation.
    """
    input:
        config["out_dir"] + "/all_samples/pan_genome/all_novel.fasta"
    output:
        config["out_dir"] + "/all_samples/non_ref/chunks/calc_n_chunks.done"
    params:
        max_jobs=config['max_jobs'] - 1,
        logs_dir=LOGS_DIR
    run:
        # get total size
        s = 0
        with open(input[0]) as f:
            for line in f:
                if line.startswith('>'):
                    continue
                s += len(line.strip())
        n_chunks =  min(s//100000, params.max_jobs)
        with open(output[0], 'w') as fo:
            print(n_chunks, file=fo)

rule prep_annotation_chunks:
    """
    Divide non-ref contigs into chunks for efficient parallel analysis
    """
    input:
        d=config["out_dir"] + "/all_samples/non_ref/chunks/calc_n_chunks.done",
        fasta=config["out_dir"] + "/all_samples/pan_genome/all_novel.fasta"
    output:
        config["out_dir"] + "/all_samples/non_ref/chunks/chunks.done"
    params:
        out_dir=config["out_dir"] + "/all_samples/non_ref/chunks",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/faSplit.yml'
    shell:
        """
        n_chunks=`cat {input.d}`
        faSplit sequence {input.fasta} $n_chunks {params.out_dir}/chunk
        touch {output}
        """

rule prep_annotation_chunks_tsv:
    """
    Prepare TSV config for liftover run
    """
    input:
        config["out_dir"] + "/all_samples/non_ref/chunks/chunks.done"
    output:
        config["out_dir"] + "/all_samples/non_ref/chunks/chunks.tsv"
    params:
        chunks_dir=config["out_dir"] + "/all_samples/non_ref/chunks",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        echo "chunk\tpath" > {output}
        realpath {params.chunks_dir}/*.fa | awk '{{n=split($0,a,"/"); sub(".yml","",a[n]); print a[n]"\t"$0}}' >> {output}
        """

rule prep_annotation_yaml:
    """
    Prepare yml config for annotation run
    """
    input:
        chunks_tsv=config["out_dir"] + "/all_samples/non_ref/chunks/chunks.tsv"
    output:
        config["out_dir"] + "/all_samples/annotation/annotation.yml"
    params:
        annotation_dir=config["out_dir"] + "/all_samples/annotation",
        templates_dir=annotation_templates_dir,
        transcripts=config['transcripts'],
        proteins=config['proteins'],
        repeats_library=config['repeats_library'],
        augustus_species=config['augustus_species'],
        min_protein=config['min_protein'],
        maker_load=config['maker_load'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        echo "name: MAKER_wrapper" >> {output}
        echo "chunks_info_file: {input.chunks_tsv}" >> {output}
        echo "out_dir: {params.annotation_dir}" >> {output}
        echo "config_templates: {params.templates_dir}" >> {output}
        echo "queue: {params.queue}" >> {output}
        echo "priority: {params.priority}" >> {output}
        echo "sample: non_ref_contigs" >> {output}
        echo "logs_dir: {params.logs_dir}" >> {output}
        echo "maker_load: {params.maker_load}" >> {output}
        echo config_kv_pairs: est={params.transcripts} protein={params.proteins} rmlib={params.repeats_library} augustus_species={params.augustus_species} min_protein={params.min_protein} >> {output}
        """

rule maker_annotation:
    """
    Run MAKER on non-ref contigs
    """
    input:
        config["out_dir"] + "/all_samples/annotation/annotation.yml"
    output:
        config["out_dir"] + "/all_samples/annotation/maker.proteins.fasta",
        config["out_dir"] + "/all_samples/annotation/maker.genes.gff"
    params:
        run_maker_in_chunks_snakefile=annotation_pipeline_dir + '/run_MAKER_in_chunks.snakefile',
        queue=config['queue'],
        jobs=config['max_jobs'],
        annotation_dir=config["out_dir"] + "/all_samples/annotation",
        cluster_param=cluster_param,
        priority=config['priority'],
        jobscript=utils_dir + '/jobscript.sh',
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        cd {params.annotation_dir}
        snakemake -s {params.run_maker_in_chunks_snakefile} --configfile {input} {params.cluster_param} -j {params.jobs} --latency-wait 60 --restart-times 3 --jobscript {params.jobscript} --keep-going -p
        """

rule rename_genes:
    """
    Assign genes short, unique names (gff and fasta).
    Names consist of the genome name and a unique ID.
    """
    input:
        fasta=config["out_dir"] + "/all_samples/annotation/maker.proteins.fasta",
        gff=config["out_dir"] + "/all_samples/annotation/maker.genes.gff"
    output:
        fasta=config["out_dir"] + "/all_samples/annotation/maker.proteins.rename.fasta",
        gff=config["out_dir"] + "/all_samples/annotation/maker.genes.rename.gff",
        gff_map=config["out_dir"] + "/all_samples/annotation/gff.map"
    params:
        maker_load=config['maker_load'],
        out_dir=config["out_dir"] + "/all_samples/annotation/",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        {params.maker_load}
        maker_map_ids --prefix PanGene_ --justify 1 --iterate 1 {input.gff} > {output.gff_map}
        cp {input.gff} {output.gff}
        map_gff_ids {output.gff_map} {output.gff}
        cp {input.fasta} {output.fasta}
        map_fasta_ids {output.gff_map} {output.fasta}
        """

rule filter_annotation:
    """
    Remove unreliable genes from gff,
    having AED > X (set by user)
    """
    input:
        gff=config["out_dir"] + "/all_samples/annotation/maker.genes.rename.gff"
    output:
        lst=config["out_dir"] + "/all_samples/annotation/maker.genes.rename.filter.list",
        gff=config["out_dir"] + "/all_samples/annotation/maker.genes.rename.filter.gff"
    params:
        max_aed=config['max_aed'],
        filter_gff_script=utils_dir + '/filter_gff_by_id_list.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        awk '{{split($9,a,";"); split(a[1],b,"="); split(a[2],c,"="); split(a[5],d,"=")}} $3 == "mRNA" && d[2] <= {params.max_aed} {{print(b[2]"\\n"c[2])}}' {input.gff} > {output.lst}
        python {params.filter_gff_script} {input.gff} {output.lst} {output.gff}
        """

rule filter_proteins:
    """
    Filter proteins fasta according to filtered gff
    """
    input:
        gff=config["out_dir"] + "/all_samples/annotation/maker.genes.rename.filter.gff",
        fasta=config["out_dir"] + "/all_samples/annotation/maker.proteins.rename.fasta"
    output:
        config["out_dir"] + "/all_samples/annotation/maker.proteins.rename.filter.fasta"
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

rule prevent_duplicate_names:
    """
    If for any reason the filtered
    proteins contain duplicate names,
    prevent this by renaming. This
    helps prevent problems in next steps
    """
    input:
        config["out_dir"] + "/all_samples/annotation/maker.proteins.rename.filter.fasta"
    output:
        config["out_dir"] + "/all_samples/annotation/maker.proteins.rename.filter.no_dupl.fasta"
    params:
        dupl_script=utils_dir + '/prevent_duplicate_names.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.dupl_script} {input} {output}
        """

rule remove_redundant_proteins:
    """
    Remove redundant non-ref proteins
    by clustering with CD-HIT and taking
    longest sequence from each cluster.
    """
    input:
        config["out_dir"] + "/all_samples/annotation/maker.proteins.rename.filter.no_dupl.fasta"
    output:
        config["out_dir"] + "/all_samples/annotation/maker.proteins.rename.filter.no_dupl.no_redun.fasta"
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
        prot=config["out_dir"] + "/all_samples/annotation/maker.proteins.rename.filter.no_dupl.no_redun.fasta",
        gff=config["out_dir"] + "/all_samples/annotation/maker.genes.rename.filter.gff"
    output:
        gff=config["out_dir"] + "/all_samples/annotation/maker.genes.rename.filter.no_dupl.no_redun.gff",
        filter_list=config["out_dir"] + "/all_samples/annotation/maker.genes.rename.filter.no_dupl.no_redun.list"
    params:
        filter_gff_script=utils_dir + '/filter_gff_by_id_list.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        grep '>' {input.prot} | sed 's/>\(.*\)-R.*/\\1/' > {output.filter_list}
        grep '>' {input.prot} | sed 's/>\(.*-R[1-9]*\).*/\\1/' >> {output.filter_list}
        python {params.filter_gff_script} {input.gff} {output.filter_list} {output.gff}
        """

rule create_pan_proteome:
    """
    Create pan genome proteins fasta
    """
    input:
        non_ref_proteins=config["out_dir"] + "/all_samples/annotation/maker.proteins.rename.filter.no_dupl.no_redun.fasta",
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
        non_ref_gff=config["out_dir"] + "/all_samples/annotation/maker.genes.rename.filter.no_dupl.no_redun.gff",
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
        detect_pav_script=pipeline_dir + '/detect_gene_PAV_from_bed_cov.py',
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
        detect_pav_script=pipeline_dir + '/detect_gene_PAV_from_bed_cov.py',
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
        create_PAV_matrix_script=pipeline_dir + '/create_PAV_matrix.py',
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
