library(cluster)
library(cowplot)
library(lattice)
library(oro.nifti)

target <- readNIfTI("M:/Documents/scans/segmentation_results/m2m_sub-003/NBM_L.nii", reorient=F)
scalp  <- readNIfTI("M:/Documents/scans/segmentation_results/m2m_sub-003/scalp.nii", reorient=F)
exclude <- readNIfTI("M:/Documents/scans/segmentation_results/m2m_sub-003/ears.nii", reorient=F)

setwd("M:/Documents/repos/PRESTUS_forked/toolboxes/TUS_entry_tool")
source("TUS_entry.R")

transducer_pos = TUS_entry(target=target, scalp=scalp, maximal_distance = 10, exclu=exclude, transducer_size=10)

cat(transducer_pos$report)
