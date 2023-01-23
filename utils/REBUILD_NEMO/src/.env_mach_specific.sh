# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
source /usr/share/Modules/init/sh
module purge 
module load ANACONDA2/python2.7 INTEL/intel_xe_2015.3.187 SZIP/szip-2.1_int15 ESMF/esmf-6.3.0rp1-mpiuni-64-O_int15 HDF5/hdf5-1.8.15-patch1 NETCDF/netcdf-C_4.3.3.1-F_4.4.2_C++_4.2.1 CMAKE/cmake-3.3.0-rc1
module unload INTEL/intel_xe_2013.5.192 INTEL/intel_xe_2013 HDF5/hdf5-1.8.10-patch1
module load INTEL/intel_xe_2015.3.187
export OMP_STACKSIZE=256M
export I_MPI_EXTRA_FILESYSTEM_LIST=gpfs
export I_MPI_EXTRA_FILESYSTEM=on
export I_MPI_PLATFORM=snb
export I_MPI_HYDRA_BOOTSTRAP=lsf
export I_MPI_LSF_USE_COLLECTIVE_LAUNCH=1
export I_MPI_DAPL_UD=on
export I_MPI_DAPL_SCALABLE_PROGRESS=on
export COMPILER=intel
export MPILIB=mpi-serial
export DEBUG=FALSE
export OS=LINUX
