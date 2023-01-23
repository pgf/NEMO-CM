#!/bin/sh -l
#BSUB -J rebuild_nemo          # Name of the job.
#BSUB -o logs/postrun_%J.out   # Appends std output to file %J.out.
#BSUB -e logs/postrun_%J.err   # Appends std error to  file %J.out.
#BSUB -q s_long             # queue
#BSUB -n 4                 # Number of CPUs
#BSUB -R "span[ptile=4]"   #
#BSUB -N 
#BSUB -W 24:00
#BSUB -L /bin/sh
#BSUB -P R000

#######################################################################
. /users_home/csp/sp1/SPS/CMCC-SPS3.5/src/scripts_oper/descr_SPS3.5.sh
. /users_home/csp/sp1/SPS/CMCC-SPS3.5/src/templates/libraries4nemo_reb.sh

set -evxu
echo "-----------STARTING postrun-------- "`date`

#######################################################################


#HERE=$PWD
HERE=/users_home/csp/as34319/CESM2/work/rebuild_nemo
REBUILD="$HERE/rebuild_nemo"
REBUILDEXE="$HERE/rebuild_nemo.exe"

#######################################################################

#######################################################################

NCPU=4

outdir=/work/csp/dp16116/CESM2/2000_cam6-nemo4_025deg_tc01b/run
ln -sf $REBUILD $outdir
ln -sf $REBUILDEXE $outdir
cd $outdir
root=2000_cam6-nemo4_025deg_tc01_00104832_restart
if [ -f ${root}_0000.nc ]; then
      npes=`ls ${root}_[0-9][0-9][0-9][0-9].nc | wc -w`
      ${REBUILD} ${root} ${npes} #&& \
#      rm -f ${root}_[0-9][0-9][0-9][0-9].nc || exit 1 &
fi

