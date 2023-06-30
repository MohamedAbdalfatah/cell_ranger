# Cell ranger for no permession 

In this Documentation I will update some script related to cell ranger falvours for people they doesn't have permession to the fastqs.

## Cell Multiplexing with cellranger multi

Here we will check how to do the preprocession analysis of Cell Multiplexing data, which is you will have different samples in the same library, you can know more information from here [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi]

# Steps 

1- Get lims **1-lims.sh** script

2- Copy files to accecable directory **copy_fastqs.sh** 

3- Generate fastq path **2-write_fastq_paths.py** script 

4- Create metadata **Create_Metadata.R**

5- Copy fastqs and Create cell ranger job **3-cp_cell.py**

6- Create *config.csv* file **create_cmo_config.py**

First I assume we are in **/home/groups/singlecell/mabdalfttah/projects**, we will need to create the subproject directory and jobs directory to save everything in those directories 

```{r}
mkdir subproject
mkdir subproject/jobs
cd subproject
```

**Step 1** LIMS information 

In this step we are calling all of the information related to the subproject to use it in the next steps, the only differencess between this script and other scritis the name of generated file **"lims_info_"$1".txt"**, this will generate "lims_info_subproject.txt" instead of "lims_info.txt"

```{r}
nano projects/scripts/1-lims.sh
```

```{r}
#!/bin/bash

# Get information for each library (flow cell, lane, sample id, etc.)
# $1  needs to be the name of the project
/home/groups/singlecell/mabdalfttah/projects/scripts/limsq.py -sp $1 | sed 's/;/\t/g' > "lims_info_"$1".txt"

echo "Created LIMS information file: lims_info.txt"
```

```{r}
chmod +x 1-lims.sh
./../scripts/1-lims.sh subproject
# Remove any fail file in LanePassFail column 
awk -F'\t' -v column="LanePassFail" 'BEGIN {OFS=FS} NR==1 {for (i=1; i<=NF; i++) if ($i == column) col=i} $col != "fail"' lims_info_DOLSORI_06.txt > tmp_file && mv tmp_file lims_info_DOLSORI_06.txt
```

**Step 2** Copy files to accecable directory
This is very important to people don't have access to read the FASTQ files in proudaction team directory, this file is generating the original fastq path in proudaction and copy them to accecable directory you select it when you run the script. the script in **/scratch_isilon/groups/singlecell/shared/projects/copy_files/copy_fastqs.sh** and this is what it is contain:

```{}
#!/bin/bash

# Step 1: Execute the 1-lims.sh script and save the output to lims_info.txt
sh /scratch_isilon/groups/singlecell/shared/projects/copy_files/scripts/1-lims.sh "$1"

# Step 2: Activate the desired conda environment
source /scratch/groups/hheyn/software/anaconda3/bin/activate cellranger

# Step 3: Run the 2-write_fastq_paths.py script and generate the fastq_paths.tab file
python /scratch_isilon/groups/singlecell/shared/projects/copy_files/scripts/2-write_fastq_paths.py --info_file lims_${1}.txt --subproject "$1"

# Step 3: Specify the target directory where you want to copy the files
[ -w "$2" ] && echo "You have permission to target directory" || echo "You do not have permission to target directory"
target_directory="$2"

# Step 4: Check if fastq_paths.tab file exists and copy the files
if [ -f "./fastq_paths.tab" ]; then
    while IFS=$'\t' read -r -a fields; do
        # Extract the path from the second column (index 1)
        path=${fields[1]}

        # Copy the file to the target directory using rsync
        rsync -avL "$path" "$target_directory"
    done < "./fastq_paths.tab"
else
    echo "fastq_paths.tab file not found."
fi

# Step 5: Change permissions of the target directory
chmod g+rwx "$2"

# Step 6: Deactivate conda env
source /scratch/groups/hheyn/software/anaconda3/bin/deactivate
```{}

To run this script you just need to write this:

```{}
sh /scratch/devel/pnieto/scripts/mo_copy.sh DOLSORI_05 /scratch_isilon/groups/singlecell/shared/projects/copy_files/fastq_dir
```{}

**DOLSORI_05** here is the subproject and **/scratch_isilon/groups/singlecell/shared/projects/copy_files/fastq_dir** is the target directory to copy the files too 

**Step 3** Generate fastq path 

Perfect! since we have nowall FASTQs in **fastq_dir** now we will generate the FASTQ path pf the new directory for our subproject, we will use this path for next script to create symbolic link in jobs directory for each sample . Since we changed the structure of the path here there is some changes in the script, first we changed where is the fastqs and second we changed the structure of creationg the path, and that is it:
 
1- fastq_path = "/scratch_isilon/groups/singlecell/shared/projects/copy_files/fastq_dir" 
change the path from **production dir** which we don't have permession to **fastq_dir** which is we have permession

2 - fastq_path_r1 = "{}/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, index)
    fastq_path_r2 = "{}/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, index)

This is the script:

```{}
#!/usr/bin/env python

# Writes fastq path by arranging proper flowcell, lane, index and read for a set of libraries

# Load packages
import numpy as np
import pandas as pd
import os
import argparse


# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to transfer feature-barcodes matrices from cluster to lcal")
parser.add_argument("--subproject",
                    dest = "subproject",
                    action = "store",
                    default = None,
                    help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--info_file",
                    dest = "info_file",
                    action = "store",
                    default = None,
                    help = "Tab-delimited file with the information of Illumina sequence of libraries for that subproject")


options = parser.parse_args()
subproject = options.subproject
info_file = options.info_file

# Read file
lims = pd.read_csv(info_file, sep = "\t", header = 0)

# Assemble fastq paths combining flowcell, lane and index
fastq_path = "/scratch_isilon/groups/singlecell/shared/projects/copy_files/fastq_dir"
fastq_path_list_r1 = []
fastq_path_list_r2 = []
for idx in lims.index:
    fc = lims.loc[idx, "flowcell"]
    lane = lims.loc[idx, "lane"]
    index = lims.loc[idx, "index"]
    fastq_path_r1 = "{}/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, index)
    fastq_path_r2 = "{}/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, index)
    fastq_path_list_r1.append(fastq_path_r1)
    fastq_path_list_r2.append(fastq_path_r2)
library_id_l = list(lims["id"].append(lims["id"]))
p_l = "P" * len(fastq_path_list_r1)
indx_l = list(range(1, len(fastq_path_list_r1) + 1))
pair_id = [p_l[x] + str(indx_l[x]) for x in range(len(indx_l))]
fastq_path_list_r1.extend(fastq_path_list_r2)
pair_id.extend(pair_id)
fastq_path_l = fastq_path_list_r1
read_l = (["R1"] * lims.shape[0]) + (["R2"] * lims.shape[0])
fastq_dict = {"library_id":library_id_l, "fastq_path":fastq_path_l, "read":read_l, "pair_id":pair_id}
fastq_df = pd.DataFrame(fastq_dict)

fastq_df.to_csv("fastq_paths.tab".format(subproject), header = True, index = False, sep="\t")

```

To run this script you need to pass the lims info from **STEP1** and the subproject

```{}
python ../scripts/2-write_fastq_paths.py --subproject  DOLSORI_05 --info_file lims_info_DOLSORI_05.txt
```


```{}
**Step 4** Cretae metadata

```{}
Path = "../Downloads/"
library(tidyverse)
#===============================================================================
Files = list.files(paste0(Path), pattern = "lims_info_DOLSORI")
All_Files = list()
metadata = list()
for (i in seq_along(Files)) {
  All_Files[[i]] = read.table(paste0(Path, Files[i]), sep = "\t", header = T)
  metadata[[i]] = data.frame(subproject = All_Files[[i]]$subproject, gem_id = All_Files[[i]]$SampleName,
                             library_id = All_Files[[i]]$id, library = All_Files[[i]]$library,
                             type = "not_hashed",donor_id = All_Files[[i]]$SampleName, flowcell = All_Files[[i]]$flowcell,
                             lane = All_Files[[i]]$lane, index = All_Files[[i]]$index)
  metadata[[i]]$gem_id = str_replace_all(string = metadata[[i]]$gem_id, pattern = "\\.", replacement = "_")
  
}
write.csv(metadata[[1]],paste0("../Downloads/DOLSORI_05.csv"), row.names = F)
```

**STEP 5** symlinks and Cellranger creation 

This script is create the subdirectory for each gem_id, create  the symbolic links for fastqs for each gem_id and create the job file for multi cellranger , most important thing changed in this script is the job mode templet to run cell ranger to work with the new cluster which is this one **/scratch/groups/singlecell/software/cellranger/6.1.1/external/martian/jobmanagers/new_cluster/slurm.template**

```{r}
# This script initializes the filesystem of this project:
# It creates a "jobs" folder which contains as many subdirectories as samples it has
# For each sample directory, it creates the following files/folders:
# 1. fastq: dir with the symlinks pointing to the fastq files
# 2. log: dir which contains standard error and output of cellranger
# 3. (sample_id).cmd: job script to compute the features-barcode matrix using cellranger


# Import required packages
import numpy as np
import pandas as pd
import os
import argparse
import subprocess
import re
import sys
import config_vars as cfg
from utils import *


# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to initialize the filesystem and scripts of this project")
parser.add_argument("--subproject",
                    dest = "subproject",
                    action = "store",
                    default = None,
                    help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--gem_id",
                    dest = "gem_id",
                    action = "store",
                    default = None,
                    help = "Gel Beads in Emulsion id")
parser.add_argument("--verbose",
                    dest = "verbose",
                    action = "store_true",
                    default = False,
                    help = "Print log in standard error")
parser.add_argument("--metadata",
                    dest = "metadata",
                    action = "store",
                    default = None,
                    help = "Metadata csv file for the tonsil atlas project")
parser.add_argument("--fastq_paths",
                    dest = "fastq_paths",
                    action = "store",
                    default = None,
                    help = "File that contains the paths of the fastqs for the subproject libraries")


def create_fastq_symlink_nh(gem_id, fastq_path_df, symlink_path):
    """Creates a symbolic link pointing to a fastq file using cellranger notation

    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      fastq_path_df: pandas dataframe with the fastq paths for that gem_id
      symlink_path: string specifying where to create the symlinks

    Returns:
      None
    """
    pair_ids = np.unique(fastq_path_df["pair_id"])
    for i in range(len(pair_ids)):
        filt = (fastq_path_df["pair_id"] == pair_ids[i])
        pair_df = fastq_path_df.loc[filt, :]
        for j in pair_df.index:
            fastq_path = pair_df.loc[j, "fastq_path"]
            lane = str(i + 1)
            read = pair_df.loc[j, "read"]
            read = read.replace("R", "")
            subprocess.run(["ln", "-s", fastq_path, "{}/{}_S1_L00{}_R{}_001.fastq.gz".format(symlink_path, gem_id, lane, read)])


def make_cellranger_nh(gem_id, jobscript_path, fastq_path, expected_cells):
    """Creates a cellranger script for a non-hashed GEM well

    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      jobscript_path: path to save the jobscript
      fastq_path: path to the fastq files
      expected_cells: expected number of high-quality cells in this experiment

    Returns:
      None
    """
    job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
    job_script = """#!/bin/bash
#SBATCH --job-name={}

#SBATCH --mail-type=all        # send email when job begins, ends, or fails
#SBATCH --mail-user=mohamed.abdalfttah@cnag.crg.eu

#SBATCH --output=%x.slurm.%J.out        # define where our output and error from the job will be stored
#SBATCH --error=%x.slurm.%J.err

#SBATCH --time=23:00:00 # set a maximum time that the job will take HH:MM:SS (process will be terminated after this is reached)

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --partition=genD
#SBATCH --mem=32G

echo [`date "+%Y-%m-%d %T"`] started job on $HOSTNAME

export TENX_IGNORE_DEPRECATED_OS=1
export HDF5_USE_FILE_LOCKING=FALSE


{} multi --id {} --csv config.csv  --jobmode /scratch/groups/singlecell/software/cellranger/6.1.1/external/martian/jobmanagers/new_cluster/slurm.template;
""".format(gem_id, cfg.cellranger_path, gem_id)
    job_script_file.write(job_script)
    job_script_file.close()


options = parser.parse_args()
subproject = options.subproject
gem_id = options.gem_id
metadata_path = options.metadata
fastq_paths = options.fastq_paths


# Read data
project_dir = "/home/groups/singlecell/mabdalfttah/projects/{}".format(subproject)
fastq_path_df = pd.read_csv(fastq_paths, sep = "\t", header = 0)
metadata_df = pd.read_csv(metadata_path, sep = ",", header = 0)
if options.verbose:
    sys.stderr.write("Files read successfully!\n")


# For each sample, create directories and jobscript
if not os.path.exists("{}/jobs".format(project_dir)):
    os.mkdir("{}/jobs".format(project_dir))
filt = (metadata_df["gem_id"] == gem_id)
metadata_df = metadata_df.loc[filt]


# Create directories
subproject_dir = "{}/jobs/{}".format(project_dir, gem_id)
fastq_dir = "{}/fastq".format(subproject_dir)
log_dir = "{}/log".format(subproject_dir)
for direct in [subproject_dir, fastq_dir, log_dir]:
    if not os.path.exists(direct):
        os.mkdir(direct)


# Define variables and subset dataframes
library_id = metadata_df.loc[filt, "library_id"]
fastq_sub_df = fastq_path_df.loc[fastq_path_df["library_id"].isin(library_id), :]
type = metadata_df["type"]
type = type.values[0]

# Create symmlinks to fastq files
create_fastq_symlink_nh(gem_id, fastq_sub_df, fastq_dir)
# Create cellranger script
fastq_dir_abs = "{}{}".format(cfg.project_path, fastq_dir)
make_cellranger_nh(gem_id, subproject_dir, fastq_dir_abs, 5000)
```

```{r}
python ../scripts/cp_cell.py  --subproject DOLSORI_05 --fastq_paths fastq_paths.tab --metadata DOLSORI_05.csv --gem_id Plex6_1
```

**STEP 5** Create config file

This is very important fro *Cell Multiplexing* and contains all of the information cellranger need to work, to create ths file we ned to have file like ths correlate the gem_ids we have with the CMOs in the library, the file should be looks like this:
```{}
subproject_folder,gem_id,GEX_library_name,CellPlex_library_name,sample_id,CMO_id
DOLSORI_05_06,Plex5_1,Plex5_GEX_1,Plex5_CellPlex_1,1956028,CMO301
DOLSORI_05_06,Plex5_1,Plex5_GEX_1,Plex5_CellPlex_1,1864260-181121,CMO302
DOLSORI_05_06,Plex5_1,Plex5_GEX_1,Plex5_CellPlex_1,1853936,CMO303
DOLSORI_05_06,Plex5_1,Plex5_GEX_1,Plex5_CellPlex_1,1864260-270521,CMO304
DOLSORI_05_06,Plex5_1,Plex5_GEX_1,Plex5_CellPlex_1,VN00213,CMO305
DOLSORI_05_06,Plex5_1,Plex5_GEX_1,Plex5_CellPlex_1,1893906,CMO307
DOLSORI_05_06,Plex5_1,Plex5_GEX_1,Plex5_CellPlex_1,VN00264,CMO308
DOLSORI_05_06,Plex5_2,Plex5_GEX_2,Plex5_CellPlex_2,1956028,CMO301
DOLSORI_05_06,Plex5_2,Plex5_GEX_2,Plex5_CellPlex_2,1864260-181121,CMO302
```

The most important column is CMO_id, sample_id and gem_ids, for each experemnets wet lab people should share with you the information of CMO_ids and sample_ids and after that you can create this file to pass it to this script.

```{}
#!/usr/bin/env python

import sys
import pandas as pd
import csv

# Get the input file and gem_id from command-line arguments
input_file = sys.argv[1]
gem_id = sys.argv[2]

# Specify the output file path and name
output_file_path = "/home/groups/singlecell/mabdalfttah/projects/DOLSORI_05/jobs/{}/config.csv".format(gem_id)

# Read the CSV file into a DataFrame
df = pd.read_csv(input_file)

# Filter the DataFrame based on the gem_id value
filtered_df = df[df['gem_id'] == gem_id]

# Update the 'fastq_id' column with the gem_id value
filtered_df['fastq_id'] = gem_id

# Data for the [gene-expression] section
gene_expression_data = [
    ['[gene-expression]'],
    ['reference', '/scratch/groups/singlecell/data/reference/refdata-gex-GRCh38-2020-A'],
    ['cmo-set', '/home/groups/singlecell/mabdalfttah/projects/data/CMO_reference.csv'],
    ['expect-cells', '30000'],
    ['chemistry', 'SC3Pv3'],
    ['no-secondary', 'true'],
    ['no-bam', 'false']
]

# Data for the [libraries] section
libraries_data = [
    ['[libraries]'],
    ['fastq_id', 'fastqs', 'feature_types'],
    [gem_id, '/home/groups/singlecell/mabdalfttah/projects/DOLSORI_06/jobs/{}/fastq'.format(gem_id), 'Multiplexing Capture'],
    [gem_id, '/home/groups/singlecell/mabdalfttah/projects/DOLSORI_05/jobs/{}/fastq'.format(gem_id), 'Gene Expression']
]

# Data for the [samples] section
samples_data = [
    ['[samples]'],
    ['sample_id', 'cmo_ids'],
] + filtered_df[['sample_id', 'CMO_id']].values.tolist()

# Open the file in write mode
with open(output_file_path, mode='w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)

    # Write the [gene-expression] section to the CSV file
    writer.writerows(gene_expression_data)

    # Add an empty line between sections
    writer.writerow([])

    # Write the [libraries] section to the CSV file
    writer.writerows(libraries_data)

    # Add an empty line between sections
    writer.writerow([])

    # Write the [samples] section to the CSV file
    writer.writerows(samples_data)

print('CSV file created successfully at {}'.format(output_file_path))
```

This script generate file called config.csv for each gm_id in it's subdirctory, it takes the CMO information and the gem_id and save the csv file as config.csv in the job directory. 

**NOTE**: the config file should  contain CMO_refernce.csv file, this file is refernce of each CMO in any experements we have in our lab, since it is same for all experements we put it in the config file becouse cellranger will need it, you don't need to do anything just change the path of this file in the script to your path of the file. In my case it is here **/home/groups/singlecell/mabdalfttah/projects/data/CMO_reference.csv**
in your case change just the path, but the file should be the same, please don't confuse beteen this fle and DOLSORI_05_06_CMO.csv which is should be uiqe for each experement and CMO_reference.csv is the same for all experements

Run the script
It takes gem_id adn CMO File
```{}
python ../scripts/create_cmo_config.py ../data/DOLSORI_05_06_CMO.csv Plex6_1
```

Perfect
