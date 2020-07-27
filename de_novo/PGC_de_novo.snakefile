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
print(sys.version)
sys.path.append(utils_dir)
from snakemakeUtils import *

def init():
    #load_info_file
    config['samples_info'] = SampleInfoReader.sample_table_reader(filename=config['samples_info_file'],
                delimiter='\t', key_name='sample', col_names=['ena_ref'])
    # load HQ genomes info file
    config['hq_info'] = SampleInfoReader.sample_table_reader(filename=config['hq_genomes_info_file'],
                delimiter='\t', key_name='sample', col_names=['annotation_gff','proteins_fasta'])
    # ensure not duplicate sample names exist
    all_names = list(config['samples_info'].keys()) + list(config['hq_info'].keys())
    assert len(all_names) == len(set(all_names)), "Can't use duplicate sample names!"

init()
LOGS_DIR = config['out_dir'] + "/logs"
CONDA_ENV_DIR = os.path.dirname(pipeline_dir) + "/conda_env"
annotation_pipeline_dir = os.path.dirname(pipeline_dir) + '/genome_annotation'
pan_genome_report_dir = os.path.dirname(pipeline_dir) + '/pan_genome_report'

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

localrules: all, prep_liftover_chunks_tsv, prep_annotation_chunks_tsv, prep_liftover_yaml, prep_annotation_yaml, require_evidence

rule all:
    input:
        pav=config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv",
        cnv=config["out_dir"] + "/all_samples/pan_genome/pan_CNV.tsv",
        prot=config["out_dir"] + "/all_samples/pan_genome/pan_proteins.fasta",
        report=config["out_dir"] + "/all_samples/stats/report.html"

def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['ena_ref']

def get_hq_sample_gff(wildcards):
    return config['hq_info'][wildcards.sample]['annotation_gff']

def get_hq_sample_proteins(wildcards):
    return config['hq_info'][wildcards.sample]['proteins_fasta']

rule download_fastq:
    """
    Download reads data from ENA
    """
    output:
        config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_1.fastq.gz",
        config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_2.fastq.gz"
    params:
        sample_out_dir=config["out_dir"] + "/per_sample/{sample}/data",
        ena_ref=get_sample,
        download_script=utils_dir + '/ena-fast-download.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    #conda:
    #    CONDA_ENV_DIR + '/ena_download.yml'
    shell:
        """
        module load curl
        python {params.download_script} {params.ena_ref} --output_directory {params.sample_out_dir}
        """

rule quality_trimming:
    """
    Trim/remove low quality reads
    """
    input:
        r1=config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_1.fastq.gz",
        r2=config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_2.fastq.gz"
    output:
        r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz",
        r1_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_unpaired.fastq.gz",
        r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz",
        r2_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_unpaired.fastq.gz"
    params:
        trimming_modules=config['trimming_modules'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/trimmomatic.yml'
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} {params.trimming_modules} -threads {params.ppn}
        """

rule merge_reads:
    """
    Merge read pairs to create long fragments
    """
    input:
        r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz",
        r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.extendedFrags.fastq.gz",
        config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_1.fastq.gz",
        config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_2.fastq.gz"
    params:
        merge_out_dir=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}",
        merge_min_overlap=config['merge_min_overlap'],
        merge_max_mismatch_ratio=config['merge_max_mismatch_ratio'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/flash.yml'
    shell:
        """
        flash {input.r1_paired} {input.r2_paired} -d {params.merge_out_dir} -m {params.merge_min_overlap} -x {params.merge_max_mismatch_ratio} -z -t {params.ppn} -o {wildcards.ena_ref}
        """

rule combine_unpaired:
    """
    Combine unpaired R1 and R2 into one file
    """
    input:
        r1_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_unpaired.fastq.gz",
        r2_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_unpaired.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_clean_unpaired.fastq.gz"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        cat {input.r1_unpaired} {input.r2_unpaired} > {output}
        """

rule genome_assembly:
    """
    De novo assembly of reads into contigs
    """
    input:
        r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_1.fastq.gz",
        r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_2.fastq.gz",
        unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_clean_unpaired.fastq.gz",
        merged=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.extendedFrags.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs.fasta"
    params:
        out_dir=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}",
        ppn=config['ppn'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/spades.yml'
    shell:
        """
        spades.py -o {params.out_dir} --pe1-1 {input.r1_paired} --pe1-2 {input.r2_paired} --pe1-m {input.merged} --pe1-s {input.unpaired} --threads {params.ppn}
        """

rule filter_contigs:
    """
    Discard contigs shorter than L or with coverage lower than C (L, C given by user)
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta"
    params:
        filter_script=utils_dir + '/filter_contigs.py',
        min_length=config['min_length'],
        min_coverage=config['min_coverage'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.filter_script} {input} {params.min_length} {params.min_coverage} {output}
        """

rule index_contigs_fasta:
    """
    Index contigs fasta. For later use by RaGOO.
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta.fai"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/samtools.yml'
    shell:
        """
        samtools faidx {input}
        """

rule assembly_quast:
    """
    Run QUAST on assembly to get assembly stats and QA
    """
    input:
        contigs=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta",
        r1=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz",
        r2=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/QUAST/report.html"
    params:
        out_dir=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/QUAST",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/quast.yml'
    shell:
        """
        quast {input.contigs} -o {params.out_dir} -t {params.ppn} -1 {input.r1} -2 {input.r2}
        """

rule ref_guided_assembly:
    """
    Assemble contigs into pseudomolecules
    by mapping to the reference genome and
    using reference-guided assembly
    """
    input:
        contigs=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta",
        ref_genome=config['ref_genome'],
    output:
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/ragoo.fasta",
        directory(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/orderings")
    params:
        ragoo_script=config['ragoo_script'],
        gap_size=config['gap_size'],
        out_dir=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}",
        queue=config['queue'],
        priority=config['priority'],
        ppn=config['ppn'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/RaGOO.yml'
    shell:
        """
        cd {params.out_dir}
        ln -f {input.contigs} contigs.fasta
        ln -f {input.ref_genome} ref.fasta
        python {params.ragoo_script} contigs.fasta ref.fasta -t {params.ppn} -g {params.gap_size}
        """

rule assembly_busco:
    """
    Run BUSCO on assembly
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/ragoo.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/run_BUSCO/short_summary_BUSCO.txt"
    params:
        assembly_dir=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output",
        busco_set=config['busco_set'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/busco.yml'
    shell:
        """
        cd {params.assembly_dir}
        run_busco -i {input} -o BUSCO -m genome -l {params.busco_set} -c {params.ppn} -f
        """

rule prep_liftover:
    """
    Prepare directory for transcripts
    liftover using GAWN
    """
    input:
        genome=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/ragoo.fasta",
        liftover_transcripts=config["out_dir"] + "/all_samples/ref/" + config['ref_name'] + '_longest_trans.fasta'
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
        config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta"
    output:
        gff=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta.transdecoder.gff3",
        pep=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/liftover_transcripts.fasta.transdecoder.pep"
    params:
        wdir=config["out_dir"] + "/per_sample/{sample}/liftover_{ena_ref}/gawn/05_results/",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/TransDecoder.yml'
    shell:
        """
        cd {params.wdir}
        TransDecoder.LongOrfs -t {input} -m 1
        TransDecoder.Predict -t {input} --single_best_only
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
        ref=config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '_REF.fasta',
        ref_db=config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '_REF.fasta.psq'
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
        blastp -query {input.liftover} -db {input.ref} -out {output} -max_target_seqs 1 -outfmt 6 -num_threads {params.ppn}
        """

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
        min_identity=config['min_identity'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.improve_script} {input.gff} {input.blastp_res} {input.transdecoder_gff} {params.min_identity} {output}
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

rule prep_chunks:
    """
    Divide assembly into chunks for efficient parallel analysis
    """
    input:
        fasta=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/ragoo.fasta",
    output:
        config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunks.lft"
    params:
        n_chunks=config['max_jobs']//len(config['samples_info']) - 1,
        out_pref=config["out_dir"] + "/per_sample/{sample}/chunks_{ena_ref}/chunk",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/faSplit.yml'
    shell:
        """
        chunkSize=`expr $(grep -v '>' {input} | tr -d '\n' | wc | awk '{{print $3}}') / {params.n_chunks}`
        faSplit gap {input} $chunkSize {params.out_pref} -noGapDrops -minGapSize=10 -lift={output}
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
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
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
        templates_dir=config["annotation_config_templates"],
        liftover_transcripts=config['liftover_transcripts'],
        additional_transcripts=config['additional_transcripts'],
        proteins=config['proteins'],
        repeats_library=config['repeats_library'],
        augustus_species=config['augustus_species'],
        min_protein=config['min_protein'],
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
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.gff",
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.all.gff"
    params:
        run_maker_in_chunks_snakefile=annotation_pipeline_dir + '/run_MAKER_in_chunks.snakefile',
        queue=config['queue'],
        jobs=config['max_jobs']//len(config['samples_info']),
        annotation_dir=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}",
        qsub_wrapper_script=utils_dir + '/pbs_qsub_snakemake_wrapper.py',
        priority=config['priority'],
        jobscript=utils_dir + '/jobscript.sh',
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        cd {params.annotation_dir}
        snakemake -s {params.run_maker_in_chunks_snakefile} --configfile {input} --cluster "python {params.qsub_wrapper_script}" -j {params.jobs} --latency-wait 60 --restart-times 3 --jobscript {params.jobscript}
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
        bedtools intersect -a {input.annotation_gff} -b {input.liftover_gff} -v | awk '$3 == "gene" || $3 == "mRNA" {{split($9,a,";"); split(a[1],b,"="); print b[2]}}' > {output.not_liftover_list}
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
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
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

rule rename_genes:
    """
    Assign genes short, unique names (gff and fasta).
    Names consist of the genome name and a unique ID.
    """
    input:
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.combine.fasta",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.gff"
    output:
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins.combine.rename.fasta",
        gff=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.genes.chr.combine.rename.gff",
        gff_map=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/gff.map",
        #fasta_map=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/fasta.map"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        maker_map_ids --prefix {wildcards.sample}_ --justify 1 --iterate 1 {input.gff} > {output.gff_map}
        cp {input.gff} {output.gff}
        map_gff_ids {output.gff_map} {output.gff}
        cp {input.fasta} {output.fasta}
        map_fasta_ids {output.gff_map} {output.fasta}
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
        ord_dir=directory(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/orderings"),
        faidx=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta.fai"
    output:
        config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/contigs.bed"
    params:
        ragoo_contigs_script=os.path.join(pipeline_dir,"ragoo_ordering_to_bed.py"),
        gap_size=config['gap_size'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        python {params.ragoo_contigs_script} {input.ord_dir} {input.faidx} {params.gap_size} {output}
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
        config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/QUAST/report.html",
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/run_BUSCO/short_summary_BUSCO.txt"
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
        #fasta_map=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/fasta.map"
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
        run_busco -i {input} -o BUSCO -m proteins -l {params.busco_set} -c {params.ppn} -f
        """

rule prep_for_orthofinder:
    """
    Prepare orthofinder input - simplify
    fasta record names and put all fasta
    files into one dir with file names
    matching genome names.
    """
    input:
        #ev=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/evidence/done",
        fasta=config["out_dir"] + "/per_sample/{sample}/annotation_{ena_ref}/maker.proteins_filter_nodupl.fasta"
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
        ref_gff=config['ref_annotation']
    output:
        config["out_dir"] + "/all_samples/ref/" + config['ref_name'] + '_longest_trans.gff'
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
        trans=config["out_dir"] + "/all_samples/ref/" + config['ref_name'] + '_longest_trans.fasta',
        fasta=config['ref_proteins'],
        gff=config["out_dir"] + "/all_samples/ref/" + config['ref_name'] + '_longest_trans.gff'
    output:
        config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '_REF.fasta'
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

rule make_ref_blast_db:
    """
    Create blast DB of ref proteins
    """
    input:
        config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '_REF.fasta'
    output:
        config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '_REF.fasta.pin',
        config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '_REF.fasta.phr',
        config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '_REF.fasta.psq'
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
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/gffutils.yml'
    shell:
        """
        python {params.filter_fasta_script} {input.gff} {input.fasta} {output} mRNA ID
        """

rule get_ref_transcripts:
    """
    Filter reference transcripts according
    to filtered gff. These transcripts will
    be used for liftover.
    """
    input:
        fasta=config['liftover_transcripts'],
        gff=config["out_dir"] + "/all_samples/ref/" + config['ref_name'] + '_longest_trans.gff'
    output:
        config["out_dir"] + "/all_samples/ref/" + config['ref_name'] + '_longest_trans.fasta'
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
        config["out_dir"] + "/all_samples/orthofinder/" + config['ref_name'] + '_REF.fasta',
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
        orthofinder -t {params.ppn} -a {params.ppn} -S diamond -n orthofinder -f {params.orthofinder_dir}
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
        ref_genome_name=config['ref_name'] + "_REF",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/break_orthogroups.yml'
    shell:
        """
        python {params.mwop_script} {params.orthofinder_dir} bitscore --allow_gene_copies exclude_ref --ref_genome_name {params.ref_genome_name} {output}
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
        ref_name=config['ref_name'] + "_REF",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/break_orthogroups.yml'
    shell:
        """
        python {params.create_pav_mat_script} {input} {params.ref_name} {output.pav} {output.cnv} {output.mapping}
        """

rule create_pan_proteins_fasta:
    """
    Create a fasta file with one
    representative protein sequence
    per pan gene.
    """
    input:
        mapping=config["out_dir"] + "/all_samples/pan_genome/OG_to_gene_names.tsv",
        mwop=config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroups_break_MWOP.tsv"
    output:
        config["out_dir"] + "/all_samples/pan_genome/pan_proteins.fasta"
    params:
        create_pan_prot_fasta_script=os.path.join(pipeline_dir,"create_pan_proteins_fasta.py"),
        og_seq_dir=config["out_dir"] + "/all_samples/orthofinder/OrthoFinder/Results_orthofinder/Orthogroup_Sequences",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/snakemake.yml'
    shell:
        """
        python {params.create_pan_prot_fasta_script} {params.og_seq_dir} {input.mapping} {input.mwop} {output}
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

rule prep_for_collect_stats:
    """
    Prepare the TSV required for
    collecting assembly stats
    """
    input:
        quast=expand(config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/QUAST/report.tsv", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        busco=expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/run_BUSCO/short_summary_BUSCO.txt", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        ragoo=expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragoo_output/ragoo.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])
    output:
        config["out_dir"] + "/all_samples/stats/assembly_stats_files.tsv"
    params:
        samples=' '.join(config['samples_info'].keys()),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        paste <(echo {params.samples} | tr ' ' '\n') <(echo {input.quast} | tr ' ' '\n') <(echo {input.busco} | tr ' ' '\n') <(echo {input.ragoo} | tr ' ' '\n') > {output}
        """

rule collect_assembly_stats:
    """
    Collect QUAST, BUSCO and RaGOO
    stats for LQ samples
    """
    input:
        config["out_dir"] + "/all_samples/stats/assembly_stats_files.tsv"
    output:
        config["out_dir"] + "/all_samples/stats/assembly_stats.tsv"
    params:
        collect_script=pan_genome_report_dir + '/collect_stats.py',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/jupyter.yml'
    shell:
        """
        python {params.collect_script} {input} {output}
        """

rule create_report_notebook:
    """
    Create pan genome report Jupyter notebook
    """
    input:
        pav_tsv=config["out_dir"] + "/all_samples/pan_genome/pan_PAV.tsv",
        stepwise_tsv=config["out_dir"] + "/all_samples/stats/stepwise_stats.tsv",
        assembly_stats_tsv=config["out_dir"] + "/all_samples/stats/assembly_stats.tsv",
        proteins_fasta=config["out_dir"] + "/all_samples/pan_genome/pan_proteins.fasta"
    output:
        config["out_dir"] + "/all_samples/stats/report.ipynb"
    params:
        ref_name=config['ref_name'] + "_REF",
        nb_template=pan_genome_report_dir + 'report_template.ipynb',
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        sed -e 's|<PAV_TSV>|{input.pav_tsv}|' -e 's|<SYEPWISE_TSV>|{input.stepwise_tsv}|' -e 's|<REF_NAME>|{params.ref_name}|' -e 's|<PROT_FASTA>|{input.proteins_fasta}|' -e 's|<STATS_TSV>|{input.assembly_stats_tsv}|' {params.nb_template} > {output}
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
        jupyter nbconvert {input} --output {output} --no-prompt --no-input --execute --NotebookClient.timeout=-1
        """
