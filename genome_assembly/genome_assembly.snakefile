"""
This pipeline assembles multiple
genomes from short reads.
The main steps are:
1. Download raw reads from ENA
2. Preprocess reads (quality trimming
   + PE merging)
3. Assemble contigs
4. filter out contamination
5. Reference-guided correction and
   scaffolding into pseudomolecules
6. Assembly quality assessment
"""

import os
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir) + '/util'
genome_assembly_dir = os.path.dirname(pipeline_dir) + '/genome_assembly'
import sys
sys.path.append(utils_dir)
from snakemakeUtils import *

PIPELINE = 'Genome-assembly'

def init():
    #load_info_file
    config['samples_info'] = SampleInfoReader.sample_table_reader(filename=config['samples_info_file'],
                delimiter='\t', key_name='sample', col_names=['ena_ref'])
    # convert '_' to '-' in sample names
    config['samples_info'] = {s.replace('_','-'): config['samples_info'][s] for s in config['samples_info']}
    # ensure not duplicate sample names exist
    all_names = list(config['samples_info'].keys())
    assert len(all_names) == len(set(all_names)), "Can't use duplicate sample names!"

init()
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

localrules: all

rule all:
    input:
        assembly_stats_tsv=config["out_dir"] + "/all_samples/stats/assembly_stats.tsv",
        assemblies=expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffold.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        r1=expand(config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        r2=expand(config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])


def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['ena_ref']

kingfisher_git_url = "https://github.com/wwood/kingfisher-download"
kingfisher_git_stable_commit = "cd7b2ed0c2488f10b91a1cf26ad3728ca26eba09"

wildcard_constraints:
    sample="[^_]+"

rule fetch_kingfisher:
    """
    Get kingfisher code
    """
    output:
        config["out_dir"] + "/kingfisher-download/bin/kingfisher"
    params:
        kingfisher_git_url=kingfisher_git_url,
        out_dir=config["out_dir"] + "/kingfisher-download/",
        kingfisher_git_stable_commit=kingfisher_git_stable_commit,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        rm -rf {params.out_dir}
        git clone {params.kingfisher_git_url} {params.out_dir}
        cd {params.out_dir}
        git checkout {params.kingfisher_git_stable_commit}
        """

rule download_fastq:
    """
    Download reads data from ENA
    """
    input:
        exe=config["out_dir"] + "/kingfisher-download/bin/kingfisher"
    output:
        config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_1.fastq.gz",
        config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_2.fastq.gz"
    params:
        sample_out_dir=config["out_dir"] + "/per_sample/{sample}/data",
        ena_ref=get_sample,
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/kingfisher.yml'
    shell:
        """
        cd {params.sample_out_dir}
        {input.exe} get -m ena-ascp -r {params.ena_ref}
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

if config['assembler'] == 'spades':
    rule genome_assembly:
        """
        De novo assembly of reads into contigs using SPAdes
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

elif config['assembler'] == 'minia':
    rule fetch_minia:
        """
        Download the gatb-minia-pipeline,
        then checkout a specific version
        and replace the main script with
        Panoramic's modified version.
        """
        output:
            config["out_dir"] + "/gatb-minia-pipeline/gatb"
        params:
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
            out_dir=config["out_dir"],
            gatb_modified=os.path.join(genome_assembly_dir,'gatb')
        shell:
            """
            cd {params.out_dir}
            git clone --recursive https://github.com/GATB/gatb-minia-pipeline
            cd {params.out_dir}/gatb-minia-pipeline
            git checkout 831ba4e
            cp {params.gatb_modified} ./
            """

    rule create_single_reads_list:
        """
        Create a file with merged and unpaired files for Minia
        """
        input:
            unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_clean_unpaired.fastq.gz",
            merged=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.extendedFrags.fastq.gz"
        output:
            config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/single_reads.list"
        params:
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR
        shell:
            """
            echo "{input.unpaired}" > {output}
            echo "{input.merged}" >> {output}
            """

    rule genome_assembly:
        """
        De novo assembly of reads into contigs using Minia
        """
        input:
            minia=config["out_dir"] + "/gatb-minia-pipeline/gatb",
            r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_1.fastq.gz",
            r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_2.fastq.gz",
            single_reads_list=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/single_reads.list"
        output:
            config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs.fasta"
        params:
            out_dir=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}",
            ppn=config['ppn'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR
        shell:
            """
            {input.minia} -1 {input.r1_paired} -2 {input.r2_paired} -s {input.single_reads_list} --nb-cores {params.ppn} --no-scaffolding -o {params.out_dir}/assembly --cleanup
            ln {params.out_dir}/assembly_final.contigs.fa {output}
            """

elif config['assembler'] == 'megahit':
    rule genome_assembly:
        """
        De novo assembly of reads into contigs using MEGAHIT
        """
        input:
            r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_1.fastq.gz",
            r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_2.fastq.gz",
            unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_clean_unpaired.fastq.gz",
            merged=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.extendedFrags.fastq.gz"
        output:
            config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs.fasta"
        conda:
            CONDA_ENV_DIR + '/megahit.yml'
        params:
            out_dir=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/megahit_out",
            ppn=config['ppn'],
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR
        shell:
            """
            if [ -d "{params.out_dir}" ]
            then
                rm -rf {params.out_dir}
            fi
            megahit -1 {input.r1_paired} -2 {input.r2_paired} -r {input.merged},{input.unpaired} -t {params.ppn} -o {params.out_dir} --min-contig-len 1
            ln {params.out_dir}/final.contigs.fa {output}
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
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/biopython.yml'
    shell:
        """
        python {params.filter_script} {input} {params.min_length} {output}
        """

rule download_db:
    output:
        db = directory(config["out_dir"] + "/all_samples/kraken-db")
    params:
       db_location=config["out_dir"] + '/all_samples/kraken-db',
       db=config['db'],
       queue=config['queue'],
       priority=config['priority'],
       logs_dir=LOGS_DIR

    shell:
        '''
        mkdir -p {params.db_location}
        wget -O db.tar.gz {params.db}
        tar -xvzf db.tar.gz -C {params.db_location}
        rm db.tar.gz
        '''

rule run_kraken:
    input:
        db=config["out_dir"] + '/all_samples/kraken-db',
        assemblies=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta"
    output:
        classification=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/classification",
        report=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/report",
        summary=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/summary"
    params:
        confidence=config['confidence'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ram=config['max_ram']
    conda:
        CONDA_ENV_DIR + '/kraken2.yml'
    shell:
        "kraken2 --use-names --threads 12 --confidence {params.confidence} --db {input.db} --input {input.assemblies} --output {output.classification} --report {output.report} &>{output.summary}"


rule filter_contamination:
    input:
        classification=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/classification",
        report=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/report",
        assemblies=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta"
    output:
        filtered_assemblies=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_no_contamination.fasta",
        contaminations_precent=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contaminations_precent"
    params:
        script=utils_dir + "/filter_contaminations.py",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/kraken2.yml'
    shell:
        "python {params.script} {input.report} {input.classification} {input.assemblies} {output.filtered_assemblies} {output.contaminations_precent}"


rule ref_guided_assembly:
    """
    Assemble contigs into pseudomolecules
    by mapping to the reference genome,
    breaking chimeric contigs and then scaffolding
    """
    input:
#        contigs=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_no_contamination.fasta",
        contigs=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta",
        ref_genome=config['reference_genome']
    output:
        corrected=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.correct.fasta",
        pm=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffold.fasta",
        pm_agp=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffold.agp"
    params:
        out_dir=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}",
        queue=config['queue'],
        priority=config['priority'],
        ppn=config['ppn'],
        logs_dir=LOGS_DIR
    conda:
        CONDA_ENV_DIR + '/RagTag.yml'
    shell:
        """
        cd {params.out_dir}
        ragtag.py correct {input.ref_genome} {input.contigs} -b 100
        ragtag.py scaffold {input.ref_genome} {output.corrected} -C -r -g 10 -t {params.ppn}
        """

rule assembly_busco:
    """
    Run BUSCO on assembly
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffold.fasta"
    output:
       config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/BUSCO/short_summary.BUSCO.txt"
    params:
        assembly_dir=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output",
        busco_set=config['busco_set'],
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn'],
        cpus=config['ppn']-1
    conda:
        CONDA_ENV_DIR + '/busco.yml'
    shell:
        """
        cd {params.assembly_dir}
        busco -i {input} -o BUSCO -m genome -l {params.busco_set} -c {params.cpus} -f --limit 2
        cp {params.assembly_dir}/BUSCO/short_summary.specific.{params.busco_set}.BUSCO.txt {output}
        """

rule assembly_quast:
    """
    Run QUAST on filtered assembly to get assembly stats and QA
    """
    input:
        contigs=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.correct.fasta"
    output:
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/QUAST/report.html",
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/QUAST/report.tsv"
    params:
        out_dir=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/QUAST",
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=config['ppn']
    conda:
        CONDA_ENV_DIR + '/quast.yml'
    shell:
        """
        quast {input.contigs} -o {params.out_dir} -t {params.ppn}
        """

rule get_data_stats:
    """
    Calculate raw data and RPP stats
    """
    input:
        raw_r1=config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_1.fastq.gz",
        raw_r2=config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_2.fastq.gz",
        clean_r1_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_paired.fastq.gz",
        clean_r2_paired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz",
        clean_r1_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_1_clean_unpaired.fastq.gz",
        clean_r2_unpaired=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}_2_clean_paired.fastq.gz",
        merged=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.extendedFrags.fastq.gz",
        unmerged_r1=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_1.fastq.gz",
        unmerged_r2=config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.notCombined_2.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.read_stats.tsv"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
        ppn=4
    shell:
        """
        set +o pipefail;
        echo 1
        zcat {input.raw_r1} | head -2 | tail -1 | awk '{{print "Read length\t"length($0)}}' > {output}
        echo 2
        zcat {input.raw_r1} {input.raw_r2} | sed -n '2~4p' | tr -d '\n' | wc | awk '{{print "Input bases\t"$3}}' >> {output}
        echo 3
        zcat {input.clean_r1_paired} {input.clean_r2_paired} {input.clean_r1_unpaired} {input.clean_r2_unpaired} | sed -n '2~4p' | tr -d '\n' | wc | awk '{{print "Clean bases\t"$3}}' >> {output}
        echo 4
        zcat {input.merged} | sed -n '2~4p' | tr -d '\n' | wc | awk '{{print "Bases in merged reads\t"$3}}' >> {output}
        """

rule prep_for_collect_stats:
    """
    Prepare the TSV required for
    collecting assembly stats
    """
    input:
        quast=expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/QUAST/report.tsv", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        busco=expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/BUSCO/short_summary.BUSCO.txt", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        ragtag=expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffold.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        data_stats=expand(config["out_dir"] + "/per_sample/{sample}/RPP_{ena_ref}/{ena_ref}.read_stats.tsv", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        contaminations_precent=expand(config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contaminations_precent", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])
    output:
        config["out_dir"] + "/all_samples/stats/assembly_stats_files.tsv"
    params:
        samples=' '.join(config['samples_info'].keys()),
        n_samples=len(config['samples_info']),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        paste <(echo {params.samples} | tr ' ' "\\n") <(echo {input.quast} | tr ' ' "\\n") <(echo {input.busco} | tr ' ' "\\n") <(echo {input.ragtag} | tr ' ' "\\n") <(echo {input.data_stats} | tr ' ' "\\n") <(echo {input.contaminations_precent} | tr ' ' "\\n") > {output}
        """

rule collect_assembly_stats:
    """
    Collect QUAST, BUSCO and RagTag
    stats for LQ samples
    """
    input:
        config["out_dir"] + "/all_samples/stats/assembly_stats_files.tsv"
    output:
        config["out_dir"] + "/all_samples/stats/assembly_stats.tsv"
    params:
        collect_script=os.path.join(genome_assembly_dir, 'collect_stats.py'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/jupyter.yml'
    shell:
        """
        python {params.collect_script} {input} {output}
        """
