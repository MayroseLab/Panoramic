#!/bin/bash
#PBS -S /bin/bash
#PBS -N test_envs
#PBS -r y
#PBS -q itaym
#PBS -V
#PBS -e /groups/itay_mayrose_nosnap/galtol/orig/Panoramic/conda_env/error_logs/ERR
#PBS -o /groups/itay_mayrose_nosnap/galtol/orig/Panoramic/conda_env/error_logs/OUT
#PBS -p 3
#PBS -l nodes=1

source ~/.bashrc
hostname
conda activate snakemake-panoramic
export PATH=$CONDA_PREFIX/bin:$PATH

cd /groups/itay_mayrose_nosnap/galtol/orig/Panoramic/conda_env

for file in *.yml; do
    env_name=$(grep '^name:' $file | cut -d ' ' -f 2)
    if conda info --envs | grep -q $env_name; then
        conda env remove --name $env_name
    fi

    mamba env create -f $file
done

echo "Done!"

