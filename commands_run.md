todo: Need to get sampling data!!!! JDS HERE

#1 - Set-up directory with eDNA barcoding bioinformatic analysis and results

Before makeing directory clone starting repo and then populate with eDNA base
```
git clone git@github.com:tamucc-gcl/wrkshp_edna_metabarcoding.git
bash gcl_bioinformatic_tools/initialize/initialize_metabarcoding.sh wrkshp_edna_metabarcoding
cd wrkshp_edna_metabarcoding
```

Create some helpful directories to seperate large sequencing data from other data that may be useful in downstream analyses/classwork
```
mkdir data/fastq_files processed_data
```

Copy Data & YAML Used in Analysis from Kevin's Directory
```
#Metadata
cp /scratch/group/p.bio240270.000/kll_marb6596_dis_edna_prj_simons_gut-content/data/simons_2018/rb_input_metadata.txt ./data

#Barcodes
cp /scratch/group/p.bio240270.000/kll_marb6596_dis_edna_prj_simons_gut-content/data/simons_2018/rb_input_barcode.txt ./data

#Sample Map
cp /scratch/group/p.bio240270.000/kll_marb6596_dis_edna_prj_simons_gut-content/data/simons_2018/rb_input_sample_map.txt ./data

#Reads
cp /scratch/group/p.bio240270.000/kll_marb6596_dis_edna_prj_simons_gut-content/data/simons_2018/*fastq.gz ./data/fastq_files

#YAML 
cp /scratch/group/p.bio240270.000/kll_marb6596_dis_edna_prj_simons_gut-content/output/simons_2018/metabarcode_rainbowbridge_paired.yml ./
mv metabarcode_rainbowbridge_paired.yml rb_param.yml
```

Set-up github & DVC 
```
module load Anaconda3; source activate dvc
dvc init
git add .dvc
git commit -m "Initialize DVC"

dvc add data/fastq_files
echo "intermediate_files/" >> .gitignore
echo "*tar.xz" >> .gitignore

git add data.dvc .gitignore
git commit -m "Track data/fastq_files with DVC"
git add .
git commit -m "Add project code and outputs"
git push
```

Set-up Remote Storage
```
dvc remote add -d crest ssh://jselwyn@crest-files.tamucc.edu:/work/birdlab/dvc_storage/wrkshp_edna_metabarcoding
dvc push
conda deactivate

git add .
git commit -m "update Remote storage"
git push
```

Manually modify `rb_param.yml` to include the needed paths/directories

Run Rainbow Bridge Analysis
```
sbatch scripts/run_rainbowBridge.sh rb_param.yml intermediate_files final
```