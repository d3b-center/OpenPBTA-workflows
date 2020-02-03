#!/bin/bash
module load matlab

export mcr_installation_dir=/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/install/GISTIC_2_0_23/
export maindir=/mnt/isilon/maris_lab/target_nbl_ngs/PBTA/gistic/2020-01-30/
export basedir=/mnt/isilon/maris_lab/target_nbl_ngs/PBTA/gistic/2020-01-30/pbta-cnv-consensus-gistic/
mkdir $maindir
mkdir $basedir
echo --- running GISTIC ---

#specify segfile
export segfile=/mnt/isilon/maris_lab/target_nbl_ngs/PBTA/v14-pbta-cnv-consensus.seg

#optional array file specifies samples in seg file on which to perform analysis
#export arrayfile=/mnt/isilon/maris_lab/target_nbl_ngs/PPTC-PDX-genomics/CHOP-processed/snp/arrays/osteo-array.txt

#run command
#export refgenefile=/mnt/isilon/diskin_lab/apps_ZV/gistic/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat 
export refgenefile=/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/install/GISTIC_2_0_23/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/install/GISTIC_2_0_23/gp_gistic2_from_seg -v 30 -b $basedir -alf $arrayfile -seg $segfile -m $markerfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -twoside 1 -brlen 0.98 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme -js 2 -rx 0
