#!/bin/bash
set -e

dx-set-timeout 10h

dx download -f project-J04Bkq8J60bFVqf61yx508gk:file-J6yfGx8J60b5Z2Yvg91v3YqG # custom data
dx download -f project-J04Bkq8J60bFVqf61yx508gk:file-J7PkbQ0J60b6V7BvKGfvfZJg # fgr.R
dx download -f project-J04Bkq8J60bFVqf61yx508gk:file-J7PpQ00J60bJZPkVb83XX283 # fgr_main.r
dx download -f project-J04Bkq8J60bFVqf61yx508gk:file-J7PkbQ8J60b6V7BvKGfvfZJk # formatting.r

mkdir -p R_libs

echo "Running fgr_main.R..."

Rscript fgr_main.R > run_report.txt 2>&1

dx upload "AD_fgr.rds" --destination project-J04Bkq8J60bFVqf61yx508gk:/dharding/results/
dx upload "VD_fgr.rds" --destination project-J04Bkq8J60bFVqf61yx508gk:/dharding/results/
dx upload "VRD_fgr.rds" --destination project-J04Bkq8J60bFVqf61yx508gk:/dharding/results/
dx upload "run_report.txt" --destination project-J04Bkq8J60bFVqf61yx508gk:/dharding/results/

echo "Finished. Shutting down..."

dx-set-timeout 20s