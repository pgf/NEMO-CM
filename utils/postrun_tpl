#!/bin/bash
#BSUB -J JOBNAME       # Name of the job.
#BSUB -o JOBNAME.out   # Appends std output to file %J.out.
#BSUB -e JOBNAME.out   # Appends std error to  file %J.out.
#BSUB -q serial_6h             # queue
#BSUB -n REBPROC                 # Number of CPUs
#BSUB -R "span[ptile=16]"   #
#BSUB -N 
#BSUB -W 01:30
#BSUB -L /bin/bash

#######################################################################

set -xv

#######################################################################
# LOAD environemnt
#######################################################################

which module
stat=$?
set -exv
if [ ${stat} -ne 0 ]; then
  source /usr/share/Modules/init/sh
fi
compver="COMPILERVER"

module purge
module purge   # Needed! Not an error!

if [ "20${compver}" == "2013" ] ; then
   # IFORT 2013
   module load SZIP/szip-2.1
   module load CMAKE/cmake-3.3.0-rc1
   module load UDUNITS/udunits-2.1.24
   module load UDUNITS/udunits-1.12.11
   module load GRIB_API/grib_api-1.10.0
   module load NETCDF/netcdf-4.3
   module load HDF5/hdf5-1.8.11
   module load ESMF/esmf-6.3.0rp1-mpiuni-64-O
   module load MAGICS/Magics-2.18.12
   module load CDO/cdo-1.6.4
   module load NCO/nco-4.2.5
   module unload INTEL/intel_xe_2013
   module unload INTEL/intel_xe_2015.3.187
   module unload IMPI/intel_mpi_5.0.3.048
   module load INTEL/intel_xe_2013.5.192

elif [ "20${compver}" == "2015" ] ; then
   # IFORT 2015
   module load SZIP/szip-2.1_int15
   module load UDUNITS/udunits-1.12.11_int15
   module load UDUNITS/udunits-2.2.19
   module load GRIB_API/grib_api-1.13.1
   module load HDF5/hdf5-1.8.15-patch1
   module load NETCDF/netcdf-C_4.3.3.1-F_4.4.2_C++_4.2.1
   module load ESMF/esmf-6.3.0rp1-mpiuni-64-O_int15
   module load MAGICS/magics-2.24.7
   module load EMOS/emos_000392
   module load EMOS/libemos-4.0.6
   module load CDO/cdo-1.7.0rc2
   module load NCO/nco-4.4.9
   module unload INTEL/intel_xe_2013
   module unload INTEL/intel_xe_2013.5.192
   module load INTEL/intel_xe_2015.3.187
   module load IMPI/intel_mpi_5.0.3.048

else
   exit "BAD COMPILER VERSION: ${compver}"
fi

#######################################################################
# General settings
#######################################################################

OUTDIR="DOUT_S_ROOT"
AUXDIR=`echo ${OUTDIR} | sed -e "s;\/hist;\/grid;"`
EXPID="CASE"
CONF="ORCA1"   # redefined later
restart="RESTART"

CDO="cdo"

REBUILD="./rebuild_nemo"

#######################################################################
# show environment settings
#######################################################################

module list
printenv | sort

#######################################################################
#           REBUILD single PEs output files if present                #
#######################################################################

NCPU=`echo ${LSB_HOSTS} | wc -w`

y1="YYYYMMDD1"
y2="YYYYMMDD2"

nproc=0
flist=`ls ${EXPID}_*_${y1}_${y2}_*grid*0000.nc 2> /dev/null | wc -w`

if [ ${flist} -gt 0 ] ; then

   flist=`ls ${EXPID}_*_${y1}_${y2}_*0000.nc`

   for ff in ${flist} ; do
   
       thisfile=`echo ${ff} |  sed -e "s/_0000.nc//"`
        
       nproc=`ps -f -u ${USER} | grep "${REBUILD}" | grep ${y1} | wc -l`
       while [ ${nproc} -ge ${NCPU} ]; do
         sleep 10
         nproc=`ps -f -u ${USER} | grep "${REBUILD}" | grep ${y1} | wc -l`
       done

       if [ -f ${thisfile}.nc -a -f ${thisfile}_0000.nc ]; then
         rm -f ${thisfile}.nc
       fi
  
       # scalar files only renamed
       isscal=0
       if [[ "${thisfile}" == *"scalar"* ]] ; then
         isscal=1
         mv ${thisfile}_0000.nc ${thisfile}.nc
       fi

       # standard file
       if [ -f ${thisfile}_0000.nc -a ${isscal} -eq 0 ]; then
         npes=`ls ${thisfile}_[0-9][0-9][0-9][0-9].nc | wc -w`
           ${REBUILD} ${thisfile} ${npes} && \
         rm -f ${thisfile}_[0-9][0-9][0-9][0-9].nc || exit 1 &
       fi
   done

   if [ ! -d ${OUTDIR} ] ; then  mkdir -p ${OUTDIR} ; fi ;

   mv ${EXPID}_*${y1}_${y2}*[z-aA-Z].nc ${OUTDIR}/ || exit 1

else
   echo "No files found for rebuild process."
fi

wait

#######################################################################
#   Create metrics for post-processing (only once if not restart)     #
#######################################################################
if [ "${restart}" == "FALSE" ] ; then

#
# 1. Check coherence of model horizontal grid resolution
#
if [ -f ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_T.nc ]; then
  fname="${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_T.nc"
elif [ -f ${OUTDIR}/${EXPID}_1d_${y1}_${y2}_grid_T.nc ]; then
  fname="${OUTDIR}/${EXPID}_1d_${y1}_${y2}_grid_T.nc"
elif [ -f ${OUTDIR}/${EXPID}_1y_${y1}_${y2}_grid_T.nc ]; then
  fname="${OUTDIR}/${EXPID}_1y_${y1}_${y2}_grid_T.nc"
else
  echo "ERROR: unable to find an output grid_T file to verify this model resolution"
  exit 1
fi
nlon=`${CDO} griddes ${fname} | grep -i xsize | cut -d'=' -f2 | tr -d '[:blank:]'`

case ${nlon} in
182)
  CONF="ORCA2"
  ;;
362)
  CONF="ORCA1"
  ;;
1442)
  CONF="ORCA025"
  ;;
*)
  echo "ERROR: unknown resolution! nlon = ${nlon}"
  exit 1
  ;;
esac

#
# 2. Create CDO grid description files
#
# Here we used monthly files by default

if [ ! -f ${CONF}_grid_T_griddes.txt -a -f ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_T.nc ]; then
    ${CDO} -f nc griddes -selname,tos \
    ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_T.nc >${CONF}_grid_T_griddes.txt
fi
if [ ! -f ${CONF}_grid_U_griddes.txt -a -f ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_U.nc ]; then
    ${CDO} -f nc griddes -selname,uo \
    ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_U.nc >${CONF}_grid_U_griddes.txt
fi
if [ ! -f ${CONF}_grid_V_griddes.txt -a -f ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_V.nc ]; then
    ${CDO} -f nc griddes -selname,vo \
    ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_V.nc >${CONF}_grid_V_griddes.txt
fi
if [ ! -f ${CONF}_grid_T_zaxisdes.txt -a -f ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_T.nc ]; then
    ${CDO} -f nc zaxisdes -selname,thetao \
    ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_T.nc | \
    sed -e "s/generic/depth_below_sea/" | sed -n "/bounds/q;p" >${CONF}_grid_T_zaxisdes.txt
fi
if [ ! -f ${CONF}_grid_W_zaxisdes.txt -a -f ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_W.nc ]; then
    ${CDO} -f nc zaxisdes -selname,wo \
    ${OUTDIR}/${EXPID}_1m_${y1}_${y2}_grid_W.nc | \
    sed -e "s/generic/depth_below_sea/" | sed -n "/bounds/q;p" >${CONF}_grid_W_zaxisdes.txt
fi

#
# 3. Rebuild mesh_mask and create grid metrics
#

meshmask="${EXPID}_mesh_mask.nc"

if [ ! -f ${meshmask} ]; then
  if [ -f mesh_mask_0000.nc ]; then
    npes=`ls mesh_mask_[0-9][0-9][0-9][0-9].nc | wc -w`
    ${REBUILD} mesh_mask ${npes} && \
    mv mesh_mask.nc ${meshmask} && \
    rm -f mesh_mask_[0-9][0-9][0-9][0-9].nc || exit 1
  fi
fi

for v in tmaskutil umaskutil vmaskutil ; do
  if [ ! -f ${EXPID}_${v}.nc -a -f ${meshmask} ]; then
    grid=`echo ${v} | cut -c1 | tr '[:lower:]' '[:upper:]'`
    ${CDO} -b F32 setgrid,${CONF}_grid_${grid}_griddes.txt -selname,${v} ${meshmask} ${EXPID}_${v}.nc
  fi
done

for v in tmask umask vmask ; do
  if [ ! -f ${EXPID}_${v}.nc -a -f ${meshmask} ]; then
    grid=`echo ${v} | cut -c1 | tr '[:lower:]' '[:upper:]'`
    ${CDO} -f nc setzaxis,${CONF}_grid_T_zaxisdes.txt -setgrid,${CONF}_grid_${grid}_griddes.txt \
      -mul -selname,${v} ${meshmask} ${EXPID}_${v}util.nc ${EXPID}_${v}.nc
  fi
done

for v in t u v ; do
  g=`echo ${v} | tr '[:lower:]' '[:upper:]'`
  if [ ! -f ${EXPID}_${v}area.nc -a -f ${meshmask} ]; then
    ${CDO} selname,e1${v} ${meshmask} ${EXPID}_e1${v}.nc
    ${CDO} setmissval,-1.e30 -ifthen ${EXPID}_${v}maskutil.nc ${EXPID}_e1${v}.nc ${EXPID}_e1${v}_masked.nc
    ${CDO} selname,e2${v} ${meshmask} ${EXPID}_e2${v}.nc
    ${CDO} setmissval,-1.e30 -ifthen ${EXPID}_${v}maskutil.nc ${EXPID}_e2${v}.nc ${EXPID}_e2${v}_masked.nc
    ${CDO} -f nc -setgrid,${CONF}_grid_${g}_griddes.txt -chname,e1${v},cell_area -mul \
      ${EXPID}_e1${v}.nc ${EXPID}_e2${v}.nc ${EXPID}_${v}area.nc
    ${CDO} -f nc setmissval,-1.e30 -ifthen ${EXPID}_${v}maskutil.nc \
      ${EXPID}_${v}area.nc ${EXPID}_${v}area_masked.nc
    ${CDO} fldsum ${EXPID}_${v}area_masked.nc ${EXPID}_${v}area_int.nc
  fi

  if [ ! -f ${EXPID}_${v}volume.nc -a -f ${meshmask} ]; then
    ${CDO} -chname,e3${v}_0,e3${v} -selname,e3${v}_0 ${meshmask} ${EXPID}_e3${v}.nc
    ${CDO} -setmissval,-1.e30 -ifthen ${EXPID}_${v}maskutil.nc ${EXPID}_e3${v}.nc ${EXPID}_e3${v}_masked.nc
    ${CDO} -f nc -setzaxis,${CONF}_grid_T_zaxisdes.txt -setgrid,${CONF}_grid_${g}_griddes.txt \
      -chname,e3${v},cell_volume -mul ${EXPID}_e3${v}.nc ${EXPID}_${v}area.nc ${EXPID}_${v}volume.nc
    ${CDO} -f nc setmissval,-1.e30 -ifthen ${EXPID}_${v}mask.nc \
      ${EXPID}_${v}volume.nc ${EXPID}_${v}volume_masked.nc
    ${CDO} vertsum -fldsum ${EXPID}_${v}volume_masked.nc ${EXPID}_${v}volume_int.nc
  fi
done

# 
# Store metrics in archive ocean grid folder ${AUXDIR}
#

if [ ! -d ${AUXDIR} ]; then
  mkdir -p ${AUXDIR}
fi

for v in T_griddes U_griddes V_griddes T_zaxisdes W_zaxisdes ; do
  if [ ! -f ${AUXDIR}/${CONF}_grid_${v}.txt ]; then
    mv ${CONF}_grid_${v}.txt ${AUXDIR}/
  fi
done
#
if [ ! -f ${AUXDIR}/${meshmask} ]; then
  mv ${meshmask} ${AUXDIR}/
fi
#
for v in tmaskutil umaskutil vmaskutil tmask umask vmask ; do
  if [ ! -f ${AUXDIR}/${EXPID}_${v}.nc ]; then
    mv ${EXPID}_${v}.nc ${AUXDIR}/
  fi
done
#
for v in t u v ; do
  if [ ! -f ${AUXDIR}/${EXPID}_${v}area.nc -a -f ${EXPID}_${v}area.nc ]; then
    mv ${EXPID}_${v}area.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}area_masked.nc -a -f ${EXPID}_${v}area_masked.nc ]; then
    mv ${EXPID}_${v}area_masked.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}area_int.nc -a -f ${EXPID}_${v}area_int.nc ]; then
    mv ${EXPID}_${v}area_int.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}volume.nc -a -f ${EXPID}_${v}volume.nc ]; then
    mv ${EXPID}_${v}volume.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}volume_masked.nc -a -f ${EXPID}_${v}volume_masked.nc ]; then
    mv ${EXPID}_${v}volume_masked.nc ${AUXDIR}/
  fi
  if [ ! -f ${AUXDIR}/${EXPID}_${v}volume_int.nc -a -f ${EXPID}_${v}volume_int.nc ]; then
    mv ${EXPID}_${v}volume_int.nc ${AUXDIR}/
  fi
  for m in e1 e2 e3 ; do
    if [ ! -f ${AUXDIR}/${EXPID}_${m}${v}.nc -a -f ${EXPID}_${m}${v}.nc ]; then
      mv ${EXPID}_${m}${v}.nc ${AUXDIR}/
    fi
    if [ ! -f ${AUXDIR}/${EXPID}_${m}${v}_masked.nc -a -f ${EXPID}_${m}${v}_masked.nc ]; then
      mv ${EXPID}_${m}${v}_masked.nc ${AUXDIR}/
    fi
  done
done
#

fi # control on restart
