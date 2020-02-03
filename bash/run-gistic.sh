#!/bin/bash
module load matlab

export mcr_installation_dir=/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/install/GISTIC_2_0_23/
export maindir=/mnt/isilon/maris_lab/target_nbl_ngs/PBTA/gistic/2020-01-30/
export basedir=/mnt/isilon/maris_lab/target_nbl_ngs/PBTA/gistic/2020-01-30/pbta-cnv-cnvkit-gistic/
mkdir $maindir
mkdir $basedir
echo --- running GISTIC ---

#specify segfile
export segfile=/mnt/isilon/maris_lab/target_nbl_ngs/PBTA/v14-pbta-cnv-cnvkit.seg

#run command
export refgenefile=/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/install/GISTIC_2_0_23/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/install/GISTIC_2_0_23/gp_gistic2_from_seg -v 30 -b $basedir -alf $arrayfile -seg $segfile -m $markerfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -twoside 1 -brlen 0.98 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme -js 2 -rx 0
