"""
This pipeline assembles multiple
genomes from short reads.
The main steps are:
1. Download raw reads from ENA
2. Preprocess reads (quality trimming
   + PE merging)
3. Assemble contigs
4. Reference-guided correction and
   scaffolding into pseudomolecules
5. Assembly quality assessment
"""

import os
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir) + '/util'
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
        config["out_dir"] + "/all_samples/stats/assembly_stats.tsv"

def get_sample(wildcards):
    return config['samples_info'][wildcards.sample]['ena_ref']

ena_fast_download_url = "https://raw.githubusercontent.com/wwood/ena-fast-download/master/ena-fast-download.py"

wildcard_constraints:
    sample="[^_]+"

rule fetch_ena_fast_download_script:
    """
    Download latest version on ena-fast-download.py
    """
    output:
        config["out_dir"] + "/ena-fast-download.py"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    run:
        shell("wget %s -O {output}" % ena_fast_download_url)

rule download_fastq:
    """
    Download reads data from ENA
    """
    input:
        config["out_dir"] + "/ena-fast-download.py"
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
        CONDA_ENV_DIR + '/ena_download.yml'
    shell:
        """
        # find ssh key in conda env
        env=`grep -l aspera ./.snakemake/conda/*.yaml | xargs basename | sed 's/\.yaml//'`
        ssh="./.snakemake/conda/$env/etc/asperaweb_id_dsa.openssh"
        # download (retry 3 times)
        n=0
        until [ "$n" -ge 3 ]
        do
            python {input} {params.ena_ref} --output_directory {params.sample_out_dir} --ssh-key $ssh && break
            n=$((n+1)) 
            echo "Download failed! Retry in 30 sec ($n/3)..."
            sleep 30
        done
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
        Download latest version of gat-minia-pipeline
        """
        output:
            config["out_dir"] + "/gatb-minia-pipeline/gatb"
        params:
            queue=config['queue'],
            priority=config['priority'],
            logs_dir=LOGS_DIR,
            out_dir=config["out_dir"]
        shell:
            """
            cd {params.out_dir}
            git clone --recursive https://github.com/GATB/gatb-minia-pipeline
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
            {input.minia} -1 {input.r1_paired} -2 {input.r2_paired} -s {input.single_reads_list} --nb-cores {params.ppn} --no-scaffolding -o {params.out_dir}/assembly
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

rule ref_guided_assembly:
    """
    Assemble contigs into pseudomolecules
    by mapping to the reference genome,
    breaking chimeric contigs and then scaffolding
    """
    input:
        contigs=config["out_dir"] + "/per_sample/{sample}/assembly_{ena_ref}/contigs_filter.fasta",
        ref_genome=config['reference_genome'],
    output:
        corrected=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/contigs_filter.corrected.fasta",
        pm=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffolds.fasta",
        pm_agp=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffolds.agp"
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
        config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffolds.fasta"
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
        busco -i {input} -o BUSCO -m genome -l {params.busco_set} -c {params.cpus} -f
        cp {params.assembly_dir}/BUSCO/short_summary.specific.{params.busco_set}.BUSCO.txt {output}
        """

rule assembly_quast:
    """
    Run QUAST on filtered assembly to get assembly stats and QA
    """
    input:
        contigs=config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/contigs_filter.corrected.fasta",
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

rule get_read_length:
    """
    Find raw data read length for stats report
    """
    input:
        config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}_1.fastq.gz"
    output:
        config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}.read_length"
    params:
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR
    shell:
        """
        set +o pipefail;
        zcat {input} | head -2 | tail -1 | wc | awk '{{print $3}}' > {output}
        """

rule prep_for_collect_stats:
    """
    Prepare the TSV required for
    collecting assembly stats
    """
    input:
        quast=expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/QUAST/report.tsv", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        busco=expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/BUSCO/short_summary.BUSCO.txt", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        ragtag=expand(config["out_dir"] + "/per_sample/{sample}/RG_assembly_{ena_ref}/ragtag_output/ragtag.scaffolds.fasta", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()]),
        read_length=expand(config["out_dir"] + "/per_sample/{sample}/data/{ena_ref}.read_length", zip, sample=config['samples_info'].keys(),ena_ref=[x['ena_ref'] for x in config['samples_info'].values()])
    output:
        config["out_dir"] + "/all_samples/stats/assembly_stats_files.tsv"
    params:
        samples=' '.join(config['samples_info'].keys()),
        n_samples=len(config['samples_info']),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    shell:
        """
        set -e
        paste <(echo {params.samples} | tr ' ' "\\n") <(echo {input.quast} | tr ' ' "\\n") <(echo {input.busco} | tr ' ' "\\n") <(echo {input.ragtag} | tr ' ' "\\n") <(echo {input.read_length} | tr ' ' "\\n") > {output}
        exit 0
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
        collect_script=os.path.join(pipeline_dir, 'genome_assembly',  'collect_stats.py'),
        queue=config['queue'],
        priority=config['priority'],
        logs_dir=LOGS_DIR,
    conda:
        CONDA_ENV_DIR + '/jupyter.yml'
    shell:
        """
        python {params.collect_script} {input} {output}
        """
