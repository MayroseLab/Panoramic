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
i = sys.argv.index('--configfile')
config_path = os.path.realpath(sys.argv[i+1])

# assert required params are in config
required = ['samples_info_file','hq_genomes_info_file','out_dir','reference_name',
            'reference_genome','reference_proteins','reference_annotation','trimming_modules',
            'merge_min_overlap','merge_max_mismatch_ratio','assembler','min_length','busco_set',
            'liftover_transcripts','min_identity','max_ratio_diff','repeats_library',
            'additional_transcripts','proteins','augustus_species','min_protein','max_aed','ppn']
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

init()
LOGS_DIR = config['out_dir'] + "/logs"
CONDA_ENV_DIR = os.path.dirname(pipeline_dir) + "/conda_env"
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

localrules: all_de_novo, prep_annotation_chunks_tsv, prep_annotation_yaml 

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

rule prep_liftover:
    """
    Prepare directory for transcripts
    liftover using GAWN
    """
    input:
        genome=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffolds.fasta",
        liftover_transcripts=config["out_dir"] + "/all_samples/ref/" + config['reference_name'] + '_longest_trans.fasta'
    output:
        genome=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/03_data/genome.fasta",
        liftover_transcripts=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/03_data/transcriptome.fasta"
    params:
        gawn_dir=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        rm -rf {params.gawn_dir}
        git clone https://github.com/enormandeau/gawn.git {params.gawn_dir}
        cp {input.genome} {output.genome}
        cp {input.liftover_transcripts} {output.liftover_transcripts}
        """

rule GAWN_liftover:
    """
    Use GAWN to perform liftover
    """
    input:
        genome=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/03_data/genome.fasta",
        liftover_transcripts=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/03_data/transcriptome.fasta"
    output:
        gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/genome.gff3"
    params:
        gawn_dir=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/GAWN.yml'
    shell:
        """
        cd {params.gawn_dir}
        sed -i 's/NCPUS=10/NCPUS={params.ppn}/' 02_infos/gawn_config.sh
        ./gawn 02_infos/gawn_config.sh
        """

rule create_liftover_transcripts_fasta:
    """
    Use liftover gff and genome fasta
    to create transcripts fasta file.
    (cDNA and protein fasta files are also
    created but are not used)
    """
    input:
        gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/genome.gff3",
        genome=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/03_data/genome.fasta"
    output:
        trans=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta",
        cdna=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_cDNA.fasta",
        prot=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_proteins.fasta"
    params:
        extract_script=utils_dir + '/extract_transcripts_and_proteins_from_gff.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.extract_script} {input.gff} {input.genome} {output.trans} {output.cdna} {output.prot}
        """

rule predict_liftover_proteins:
    """
    Use TransDecoder to predict proteins
    from liftover transcripts
    """
    input:
        fasta=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta",
        ref=config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta',
        ref_db=config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta.psq'
    output:
        gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta.transdecoder.gff3",
        pep=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta.transdecoder.pep"
    params:
        wdir=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/",
        min_protein=max(50,config['min_protein']),
        queue=config['queue'],
        priority=config['priority'],
        ppn=config['ppn'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/TransDecoder.yml'
    shell:
        """
        cd {params.wdir}
        TransDecoder.LongOrfs -t {input.fasta} -m {params.min_protein}
        blastp -query liftover_transcripts.fasta.transdecoder_dir/longest_orfs.pep -db {input.ref} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads {params.ppn} > blastp.outfmt6
        TransDecoder.Predict -t {input.fasta} --single_best_only --no_refine_starts --retain_blastp_hits blastp.outfmt6
        """

rule simplify_transdecoder_protein_names:
    """
    Remove unwanted parts from fasta headers
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta.transdecoder.pep"
    output:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_TransDecoder_proteins.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        sed 's/>\(.*\)\.p[0-9]* GENE.*/>\\1/' {input} > {output}
        """

rule blast_liftover_proteins:
    """
    Run Blastp of liftover proteins
    against the reference proteins
    """
    input:
        liftover=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_TransDecoder_proteins.fasta",
        ref=config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta',
        ref_db=config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta.psq'
    output:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta.transdecoder.pep_vs_ref.blast6"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/GAWN.yml'
    shell:
        """
        blastp -query {input.liftover} -db {input.ref} -out {output} -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -num_threads {params.ppn}
        """

# find min and max gene length in ref gff
def min_max_gene_len(gff):
    min_len = 999999
    max_len = 0
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
            if fields[2] == "gene":
                gene_len = int(fields[4]) - int(fields[3])
                if gene_len < min_len:
                    min_len = gene_len
                elif gene_len > max_len:
                    max_len = gene_len
    return min_len, max_len

ref_min_gene_len, ref_max_gene_len = min_max_gene_len(config['reference_annotation'])
        
rule improve_liftover_result:
    """
    Remove low quality predictions
    from liftover results and ensure
    one gene with one mRNA per
    ref transcript. Also rewrite GFF
    CDS features
    """
    input:
        gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/genome.gff3",
        transdecoder_gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta.transdecoder.gff3",
        blastp_res=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta.transdecoder.pep_vs_ref.blast6"
    output:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/genome_improve.gff3"
    params:
        improve_script=os.path.join(pipeline_dir,"improve_GAWN_liftover.py"),
        min_len=int(ref_min_gene_len*0.8),
        max_len=int(ref_max_gene_len*1.2),
        min_identity=config['min_identity'],
        max_ratio_diff=config['max_ratio_diff'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.improve_script} {input.gff} {params.min_len} {params.max_len} {input.blastp_res} {input.transdecoder_gff} {params.min_identity} {params.max_ratio_diff} {output}
        """

rule get_liftover_proteins:
    """
    Fetch the final set of filtered liftover proteins
    """
    input:
        gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/genome_improve.gff3",
        fasta=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_TransDecoder_proteins.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_TransDecoder_proteins_improve.fasta"
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

rule get_liftover_transcripts:
    """
    Fetch the final set of filtered liftover transcripts
    """
    input:
        gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/genome_improve.gff3",
        fasta=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts_improve.fasta"
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

def calc_n_chunks(max_jobs, n_samples):
    if n_samples < 1:
        return 0
    return max_jobs//n_samples - 1

rule prep_chunks:
    """
    Divide assembly into chunks for efficient parallel analysis
    """
    input:
        fasta=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffolds.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunks.lft"
    params:
        n_chunks=calc_n_chunks(config['max_jobs'],len(config['samples_info'])),
        out_pref=config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunk",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/faSplit.yml'
    shell:
        """
        chunkSize=`expr $(grep -v '>' {input} | tr -d '\n' | wc | awk '{{print $3}}') / {params.n_chunks}`
        faSplit gap {input} $chunkSize {params.out_pref} -minGapSize=10 -lift={output}
        """ 

rule prep_annotation_chunks_tsv:
    """
    Prepare TSV config for annotation run
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunks.lft" 
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/chunks.tsv"
    params:
        chunks_dir=config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        echo "chunk\tpath" > {output}
        realpath {params.chunks_dir}/*.fa | awk '{{n=split($0,a,"/"); print a[n]"\t"$0}}' >> {output}
        """

rule prep_annotation_yaml:
    """
    Prepare yml config for annotation run
    """
    input:
        chunks_tsv=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/chunks.tsv",
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/annotation.yml"
    params:
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}",
        templates_dir=annotation_templates_dir,
        liftover_transcripts=config['liftover_transcripts'],
        additional_transcripts=config['additional_transcripts'],
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
        echo "sample: {wildcards.sample}" >> {output}
        echo "logs_dir: {params.logs_dir}" >> {output}
        echo "maker_load: {params.maker_load}" >> {output}
        echo config_kv_pairs: est={params.liftover_transcripts},{params.additional_transcripts} protein={params.proteins} rmlib={params.repeats_library} augustus_species={params.augustus_species} min_protein={params.min_protein} >> {output}
        """

rule maker_annotation:
    """
    Run MAKER using evidence and ab-initio
    to obtain gene models on top of the liftover
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/annotation.yml"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.fasta",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts.fasta",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.all.gff"
    params:
        run_maker_in_chunks_snakefile=annotation_pipeline_dir + '/run_MAKER_in_chunks.snakefile',
        queue=config['queue'],
        jobs=calc_n_chunks(config['max_jobs'],len(config['samples_info'])),
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}",
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

rule make_chunks_bed:
    """
    Create a bed file containing chunks borders.
    Useful as part of annotation evidence collection.
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunks.lft"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/chunks.bed"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        awk '{{print $4"\t"$1"\t"$1+$3"\t"$2}}' {input} > {output}
        """

rule convert_chunks_to_chromosomes:
    """
    Convert gff coordinates and sequence names
    to transform from chunks to chromosomes
    """
    input:
        genes=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.gff",
        all_=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.all.gff",
        bed=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/chunks.bed"
    output:
        genes_out=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.gff",
        all_out=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.all.chr.gff"
    params:
        convert_script=os.path.join(utils_dir,"transform_gff_coordinates.py"),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.convert_script} {input.genes} {input.bed} {output.genes_out}
        python {params.convert_script} {input.all_} {input.bed} {output.all_out}
        """

rule get_novel_genes:
    """
    Fetch only genes not overlapping with
    liftover genes
    """
    input:
        liftover_gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/genome_improve.gff3",
        annotation_gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.gff"
    output:
        not_liftover_list=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/not_liftover.list",
        not_liftover_gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.not_liftover.gff"
    params:
        filter_script=os.path.join(utils_dir,"filter_gff_by_id_list.py"),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/index_gff.yml'
    shell:
        """
        bedtools intersect -a {input.annotation_gff} -b {input.liftover_gff} -v -f 0.1 | awk '$3 == "gene" || $3 == "mRNA" {{split($9,a,";"); split(a[1],b,"="); print b[2]}}' > {output.not_liftover_list}
        python {params.filter_script} {input.annotation_gff} {output.not_liftover_list} {output.not_liftover_gff}
        """

rule get_novel_proteins:
    """
    Fetch proteins selected in the novel gff
    """
    input:
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.fasta",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.not_liftover.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.not_liftover.fasta"
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

rule get_novel_transcripts:
    """
    Fetch transcripts selected in the novel gff
    """
    input:
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts.fasta",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.not_liftover.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts.not_liftover.fasta"
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

rule combine_liftover_with_novel_gff:
    """
    Combine liftover genes with novel genes
    """
    input:
        liftover_gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/genome_improve.gff3",
        novel_gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.not_liftover.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.gff"
    params:
        maker_load=config['maker_load'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        {params.maker_load}
        gff3_merge {input.liftover_gff} {input.novel_gff} -s > {output}
        """

rule combine_liftover_with_novel_proteins:
    """
    Combine liftover proteins with novel proteins
    """
    input:
        liftover_fasta=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_TransDecoder_proteins_improve.fasta",
        novel_fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.not_liftover.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.combine.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cat {input.liftover_fasta} {input.novel_fasta} > {output}
        """

rule combine_liftover_with_novel_transcripts:
    """
    Combine liftover proteins with novel transcripts
    """
    input:
        liftover_fasta=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts_improve.fasta",
        novel_fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts.not_liftover.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts.combine.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        cat {input.liftover_fasta} {input.novel_fasta} > {output}
        """

rule rename_genes:
    """
    Assign genes short, unique names (gff and fasta).
    Names consist of the genome name and a unique ID.
    """
    input:
        prot=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.combine.fasta",
        trans=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts.combine.fasta",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.gff"
    output:
        prot=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.combine.rename.fasta",
        trans=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts.combine.rename.fasta",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.rename.gff",
        gff_map=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/gff.map",
    params:
        maker_load=config['maker_load'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        {params.maker_load}
        maker_map_ids --prefix {wildcards.sample}_ --justify 1 --iterate 1 {input.gff} > {output.gff_map}
        cp {input.gff} {output.gff}
        map_gff_ids {output.gff_map} {output.gff}
        cp {input.prot} {output.prot}
        map_fasta_ids {output.gff_map} {output.prot}
        cp {input.trans} {output.trans}
        map_fasta_ids {output.gff_map} {output.trans}
        """

rule make_evidence_gffs:
    """
    Extract evidence features from MAKER gff.
    Useful as part of annotation evidence collection.
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.all.chr.gff"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.augustus.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.blastn.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.blastx.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.est2genome.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.pred_gff:maker.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.protein2genome.gff"
    params:
        split_gff_script=os.path.join(utils_dir,"split_gff_by_source.py"),
        out_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.split_gff_script} {input} augustus,blastn,blastx,est2genome,pred_gff:maker,protein2genome {params.out_dir} 
        """

rule make_contigs_bed:
    """
    Create a bed file with original contig borders.
    Useful as part of annotation evidence collection.
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffolds.agp"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/contigs.bed"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        awk '$5 == "W" {{print $1"\t"$2"\t"$3"\t"$6}}' {input} > {output}
        """

rule index_evidence_gff:
    """
    Sort, compress and index (tabix)
    evidence gffs for easy display in IGV
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/contigs.bed",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.augustus.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.blastn.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.blastx.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.est2genome.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.pred_gff:maker.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/maker.all.chr.protein2genome.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/chunks.bed",
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/done"
    params:
        evidence_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/index_gff.yml'
    shell:
        """
        cd {params.evidence_dir}
        for x in `ls -1 *.gff`; do srt=`echo $x | sed 's/\.gff/\.sort\.gff/'`; bedtools sort -i $x > $srt; bgzip $srt; tabix $srt.gz; done
        touch {output}
        """

rule filter_annotation:
    """
    Remove unreliable genes from gff.
    These are genes that do NOT come from
    lift-over AND have AED > X (set by user)
    """
    input:
        gff_map=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/gff.map",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.rename.gff"
    output:
        lst=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.rename.filter.list",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.rename.filter.gff"
    params:
        max_aed=config['max_aed'],
        filter_gff_script=utils_dir + '/filter_gff_by_id_list.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        awk '{{split($9,a,";"); split(a[1],b,"="); split(a[2],c,"="); split(a[5],d,"=")}} $3 == "mRNA" && ($2 == "indexed_genome" || d[2] <= {params.max_aed}) {{print(b[2]"\\n"c[2])}}' {input.gff} > {output.lst}
        python {params.filter_gff_script} {input.gff} {output.lst} {output.gff}
        """

rule filter_proteins:
    """
    Filter proteins fasta according to filtered gff  
    """
    input:
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.rename.filter.gff",
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.combine.rename.fasta",
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter.fasta"
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

rule filter_transcripts:
    """
    Filter transcripts fasta according to filtered gff
    """
    input:
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.rename.filter.gff",
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts.combine.rename.fasta",
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts_filter.fasta"
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
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter_nodupl.fasta"
    params:
        dupl_script=utils_dir + '/prevent_duplicate_names.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.dupl_script} {input} {output}
        """

rule annotation_busco:
    """
    Run BUSCO on filtered annotation proteins
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter_nodupl.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/run_BUSCO/short_summary_BUSCO.txt"
    params:
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}",
        busco_set=config['busco_set'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/busco.yml'
    shell:
        """
        cd {params.annotation_dir}
        busco -i {input} -o BUSCO -m proteins -l {params.busco_set} -c {params.ppn} -f
        """

rule prep_for_orthofinder:
    """
    Prepare orthofinder input - simplify
    fasta record names and put all fasta
    files into one dir with file names
    matching genome names.
    """
    input:
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter_nodupl.fasta",
        evidence=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/done"
    output:
        config["out_dir"] + "/all_samples/orthofinder/{sample}_{ena_ref}_LQ.fasta"
    params:
        of_dir=config["out_dir"] + "/all_samples/orthofinder",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        sed 's/ protein .*//' {input.fasta} > {output}
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

rule get_hq_sample_proteins:
    """
    Filter HQ samples proteins according
    to filtered gff and put the new file
    in orthofinder dir.
    """
    input:
        fasta=get_hq_sample_proteins,
        gff=config["out_dir"] + "/HQ_samples/{sample}/{sample}_longest_trans.gff"
    output:
        config["out_dir"] + "/all_samples/orthofinder/{sample}_HQ.fasta"
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
        fasta=config['liftover_transcripts'],
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

rule orthofinder:
    """
    Run OrthoFinder2 on all proteins
    from all annotated genomes to get
    initial orthogroups
    """
    input:
        expand(config["out_dir"] + "/all_samples/orthofinder/{sample}_{ena_ref}_LQ.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        config["out_dir"] + "/all_samples/orthofinder/" + config['reference_name'] + '_REF.fasta',
        expand(config["out_dir"] + "/all_samples/orthofinder/{sample}_HQ.fasta", sample=config['hq_info'].keys())
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
        lq=expand(config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter_nodupl.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        hq=expand(config["out_dir"] + "/all_samples/orthofinder/{sample}_HQ.fasta", sample=config['hq_info'].keys()),
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
        lq=expand(config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.transcripts_filter.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
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
