import os
pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))
utils_dir = os.path.dirname(pipeline_dir)
import sys
sys.path.append(utils_dir)
from snakemakeUtils import *
from time import time

def init(): 
    #load_info_file
    config['chunks_info'] = SampleInfoReader.sample_table_reader(filename=config['chunks_info_file'], 
                delimiter='\t', key_name='chunk', col_names=['path'], opt_col_names=['pred_gff'])

init()

pipeline_dir = os.path.dirname(os.path.realpath(workflow.snakefile))

onstart:
    write_config_file(config)

onsuccess:
    print("%s pipeline finished, no error" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

onerror:
    print("%s pipeline failed" % config['name'])
    shell("cat {log} >> %s/run_log.txt" % config["out_dir"])

#------------------------------------
#                RULES              |
#------------------------------------

localrules: all, prep_maker_configs
all_chunks = config['chunks_info'].keys()

rule all:
    input:
        config["out_dir"] + "/maker.all.gff",
        config["out_dir"] + "/maker.genes.gff",
        config["out_dir"] + "/maker.transcripts.fasta",
        config["out_dir"] + "/maker.proteins.fasta",
        expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.rename.gff", chunk=all_chunks),
        expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.genes.rename.gff", chunk=all_chunks)
	
def get_chunk(wildcards):
    return config['chunks_info'][wildcards.chunk]

rule prep_maker_configs:
    input:
        unpack(get_chunk)
    output:
        bopts=config["out_dir"] + "/chunks/{chunk}/maker_bopts.ctl",
        opts=config["out_dir"] + "/chunks/{chunk}/maker_opts.ctl",
        exe=config["out_dir"] + "/chunks/{chunk}/maker_exe.ctl"
    params:
        templates=config["config_templates"],
        config_edit_script=pipeline_dir + '/edit_maker_conf.py',
        config_kv_pairs=config["config_kv_pairs"]
    run:
        shell("cp {params.templates}/maker_bopts.ctl {output.bopts}")
        shell("cp {params.templates}/maker_exe.ctl {output.exe}")
        if hasattr(input,'pred_gff'):
            shell("python {params.config_edit_script} {params.templates}/maker_opts.ctl {output.opts} --edits genome={input.path} pred_gff={input.pred_gff} {params.config_kv_pairs}")
        else:
            shell("python {params.config_edit_script} {params.templates}/maker_opts.ctl {output.opts} --edits genome={input.path} {params.config_kv_pairs}")
        

rule run_maker:
    input:
        bopts=config["out_dir"] + "/chunks/{chunk}/maker_bopts.ctl",
        opts=config["out_dir"] + "/chunks/{chunk}/maker_opts.ctl",
        exe=config["out_dir"] + "/chunks/{chunk}/maker_exe.ctl"
    output:
        config["out_dir"] + "/chunks/{chunk}/maker.done"
    log:
        index=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk_master_datastore_index.log"
    params:
        run_dir=config["out_dir"] + "/chunks/{chunk}",
        queue=config['queue'],
        priority=config['priority'],
        sample=config['sample'],
        logs_dir=config['logs_dir']
    shell:
        """
        cd {params.run_dir}
        module load miniconda/miniconda2-4.5.4-MakerMPI
        maker -b chunk -fix_nucleotides
        if [ -f {log.index} ] && (( `grep STARTED {log.index} | wc -l` <= `grep FINISHED {log.index} | wc -l` )); then touch {output}; fi
        """

rule create_full_gff:
    input:
       done=config["out_dir"] + "/chunks/{chunk}/maker.done",
       index=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk_master_datastore_index.log" 
    output:
        config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.gff"
    params:
        queue=config['queue'],
        priority=config['priority'],
        sample=config['sample'],
        logs_dir=config['logs_dir']
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        gff3_merge -d {input.index} -n -s > {output}
        """

rule create_genes_gff:
    input:
       done=config["out_dir"] + "/chunks/{chunk}/maker.done",
       index=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk_master_datastore_index.log"
    output:
        config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.genes.gff"
    params:
        queue=config['queue'],
        priority=config['priority'],
        sample=config['sample'],
        logs_dir=config['logs_dir']
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        gff3_merge -d {input.index} -n -g -s > {output}
        """

rule rename_gff_features:
    input:
        config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.{type}.gff"
    output:
        config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.{type}.rename.gff"
    params:
        queue=config['queue'],
        priority=config['priority'],
        sample=config['sample'],
        logs_dir=config['logs_dir']
    shell:
        """
        sed 's/ID=\\([^;]\\+\\)\\(.*\\)Name=\\([^;]\\+\\)/ID=\\1\\2Name=\\3__\\1/' {input} > {output}
        """

rule merge_full_gff:
    input:
       expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.gff", chunk=all_chunks)
    output:
        config["out_dir"] + "/maker.all.gff"
    params:
        queue=config['queue'],
        priority=config['priority'],
        sample=config['sample'],
        logs_dir=config['logs_dir']
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        gff3_merge -n -s {input} | grep -v '###' > {output}
        """

rule merge_genes_gff:
    input:
       expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.genes.gff", chunk=all_chunks)
    output:
        config["out_dir"] + "/maker.genes.gff"
    params:
        queue=config['queue'],
        priority=config['priority'],
        sample=config['sample'],
        logs_dir=config['logs_dir']
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        gff3_merge -n -s {input} | grep -v '###' > {output}
        """

#rule convert_gff_coords:
#    input:
#        config["out_dir"] + "/maker.{type}.gff"
#    output:
#        config["out_dir"] + "/maker.{type}.convert.gff"
#    params:
#        coord_conversion_script = config["coord_conversion_script"],
#        queue=config['queue'],
#        priority=config['priority'],
#        sample=config['sample'],
#        logs_dir=config['logs_dir']
#    shell:
#        "python {params.coord_conversion_script} {input} {output}"

rule create_fasta:
    input:
       done=config["out_dir"] + "/chunks/{chunk}/maker.done",
       index=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk_master_datastore_index.log",
       genes_gff=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.genes.gff"
    output:
        trans=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.maker.transcripts.fasta",
        prot=config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.maker.proteins.fasta"
    params:
        out_dir = config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/",
        queue=config['queue'],
        priority=config['priority'],
        sample=config['sample'],
        logs_dir=config['logs_dir']
    shell:
        """
        module load miniconda/miniconda2-4.5.4-MakerMPI
        cd {params.out_dir}
        fasta_merge -d {input.index}
        if [ ! -f {output.prot} ] && [ `grep -v '#' {input.genes_gff} | wc -l` == 0 ]; then touch {output.prot}; fi
        if [ ! -f {output.trans} ] && [ `grep -v '#' {input.genes_gff} | wc -l` == 0 ]; then touch {output.trans}; fi
        """

rule merge_transcripts_fasta:
    input:
       expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.maker.transcripts.fasta", chunk=all_chunks)
    output:
        config["out_dir"] + "/maker.transcripts.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        sample=config['sample'],
        logs_dir=config['logs_dir']
    shell:
        """
        cat {input} > {output}
        """
        
rule merge_proteins_fasta:
    input:
       expand(config["out_dir"] + "/chunks/{chunk}/chunk.maker.output/chunk.all.maker.proteins.fasta", chunk=all_chunks)
    output:
        config["out_dir"] + "/maker.proteins.fasta"
    params:
        queue=config['queue'],
        priority=config['priority'],
        sample=config['sample'],
        logs_dir=config['logs_dir']
    shell:
        """
        cat {input} > {output}
        """
