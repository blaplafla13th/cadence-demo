#!/bin/bash
micromamba install bioconda::bcftools bioconda::tabix
bash gip.sh -data_path . -data_name chr6_327k.all.0.25.missing -num_al 3 -num_gpus 1