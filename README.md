# Cell ranger for no permession 

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

**STEP 4** Cellranger creation 

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
parser.add_argument("--feat_ref",
                    dest = "feat_ref",
                    action = "store",
                    default = None,
                    help = "Feature hashtag reference correspondence path")
parser.add_argument("--fastq_paths",
                    dest = "fastq_paths",
                    action = "store",
                    default = None,
                    help = "File that contains the paths of the fastqs for the subproject libraries")
options = parser.parse_args()
subproject = options.subproject
gem_id = options.gem_id
metadata_path = options.metadata
feat_ref = options.feat_ref
fastq_paths = options.fastq_paths


# Read data
project_dir = "projects/{}".format(subproject)
fastq_path_df = pd.read_csv(fastq_paths, sep = "\t", header = 0)
metadata_df = pd.read_csv(metadata_path, sep = ",", header = 0)
feature_reference_df = pd.read_csv(feat_ref, sep = ",", header = 0)
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

if type == "not_hashed":
    # Create symmlinks to fastq files
    create_fastq_symlink_nh(gem_id, fastq_sub_df, fastq_dir)


    # Create cellranger script
    fastq_dir_abs = "{}{}".format(cfg.project_path, fastq_dir)
    make_cellranger_nh(gem_id, subproject_dir, fastq_dir_abs, 5000)

else:
    if not os.path.exists("{}/cDNA".format(fastq_dir)):
        os.mkdir("{}/cDNA".format(fastq_dir))
    if not os.path.exists("{}/HTO".format(fastq_dir)):
        os.mkdir("{}/HTO".format(fastq_dir))
    for lib in library_id:
        # Create symmlinks to fastq files
        lib_type = metadata_df.loc[metadata_df["library_id"] == lib, "type"]
        lib_type = lib_type.values[0]
        if lib_type == "hashed_cdna":
            symlink_path = "{}/cDNA".format(fastq_dir)
        elif lib_type == "hashed_hto":
            symlink_path = "{}/HTO".format(fastq_dir)
        create_fastq_symlink_h(gem_id, lib, lib_type, fastq_sub_df, symlink_path)


    # Create libraries.csv: file indicating fastq path and type of reads (HTO/cDNA)
    write_libraries_csv(gem_id, os.path.abspath(subproject_dir))


    # Create feature_reference.csv: file indicating correspondence between HTO and sample id
    filt_row = ((feature_reference_df["subproject"] == subproject) & (feature_reference_df["gem_id"] == gem_id))
    filt_col = np.invert(feature_reference_df.columns.isin(["subproject", "gem_id"]))
    feature_ref_sub = feature_reference_df.loc[filt_row, filt_col]
    feature_ref_sub.to_csv("{}/feature_reference.csv".format(subproject_dir), header = True, index = False)


    # Create cellranger script
    make_cellranger_h(gem_id, subproject_dir, 20000)

```

```{r}
python ../scripts/3-make_cellranger.py  --subproject DOLSORI_06 --fastq_paths fastq_paths.tab --metadata DOLSORI_06.csv --gem_id Plex5_CellPlex_1 --feat_ref feature_reference.csv
```
