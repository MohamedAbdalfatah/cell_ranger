# Cell ranger for no permession 

## We have done:
1- lims 

2- fastq path

3- metadata

4- copy fastqs and create directory of gem_id

5- create cell ranger job

## We need to do:

1- Modify cellranger job to be multi and add config file:
    a) cellranger multi comand 
    b) Sample (gem_id)
    c) config_file.csv

2- create confi file manually or by script 
Here all of what we need to put this:

[gene-expression]

reference,/scratch/groups/singlecell/data/reference/refdata-gex-GRCh38-2020-A (No changed need)

cmo-set,/scratch/devel/ljimenez/projects/DOLSORI/01_cellranger_mapping/data/CMO_reference.csv (File avilable in Laura chat but we need to change the path)

expect-cells,30000 (No changed need)
chemistry,SC3Pv3 (No changed need)
no-secondary,true (No changed need)
no-bam,false (No changed need)
[libraries] (Change the path of the fastq files, we have directory for CMO and cDNA)
fastq_id,fastqs,feature_types 
Plex5_2_CMO,/scratch/devel/ljimenez/projects/DOLSORI/01_cellranger_mapping/subprojects/DOLSORI_05_06/jobs/Plex5_2/fastq/CMO,Multiplexing Capture
Plex5_2,/scratch/devel/ljimenez/projects/DOLSORI/01_cellranger_mapping/subprojects/DOLSORI_05_06/jobs/Plex5_2/fastq/cDNA,Gene Expression

[samples] (We need to write a script to generate this part, we can generate it from wetlab refernce )
sample_id,cmo_ids
1956028,CMO301
1864260-181121,CMO302
1853936,CMO303
1864260-270521,CMO304
VN00213,CMO305
1893906,CMO307
VN00264,CMO308

3- understand the structure of the project 

In this file I'm going to update the cellranger files for people dosn't have permession 

```{r}
mkdir subproject
mkdir subproject/jobs
cd subproject
```
**Step 1** LIMS information 
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
# Copy the feature file
cp ../cellranger_mapping/feature_reference.csv .
```

**Step 2** Write fastq files 

here there is some changes: 
1- fastq_path = "/scratch_isilon/groups/singlecell/shared/projects/copy_files/fastq_dir" 
change the path from pridaction which we don't have permession to fastq_dir which is we have permession
2 - fastq_path_r1 = "{}/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, index)
    fastq_path_r2 = "{}/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, index)


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

```{}
python ../scripts/2-write_fastq_paths.py --subproject  DOLSORI_05 --info_file lims_info_DOLSORI_05.txt
```


```{}
**Step 3** Cretae metadata

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

**STEP 4** symlinks and Cellranger creation 

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

#SBATCH --cpus-per-task=1
#SBATCH --partition=genB,main
#SBATCH --ntasks=1
#SBATCH --qos=normal

export TENX_IGNORE_DEPRECATED_OS=1
export HDF5_USE_FILE_LOCKING=FALSE


{} multi --id {}  --jobmode /scratch/groups/hheyn/software/cellranger/6.1.1/external/martian/jobmanagers/slurm.template;
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
python ../scripts/cp_cell.py  --subproject DOLSORI_05 --fastq_paths fastq_paths.tab --metadata DOLSORI_05.csv --gem_id Plex5_1
```

**STEP 5** Create config file

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

Run the script
it takes the CMO information and the gem_id and save the csv file as config.csv in the job directory 

```{}
python ../scripts/create_cmo_config.py ../data/DOLSORI_05_06_CMO.csv Plex5_1
```
