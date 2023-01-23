#!/bin/bash -exv
#
# Download NEMO model files in the local path

target_dir="nemocore"

if [ ! -d ${target_dir} ] ; then
   #mkdir -p ${target_dir}/src ${target_dir}/ext ${target_dir}/tools ${target_dir}/cfgs
   mkdir -p ${target_dir}
fi

# clean up
if [ -d tmp ] ; then  rm -rf tmp ; fi

# get code core file
branch="branch_4.2"
commit="05678c0d"
gitrepo="https://forge.nemo-ocean.eu/nemo/nemo.git"

echo
echo "Get code from branch ${branch} at commit ${commit}"
echo "Remote origin : ${gitrepo} "
echo

git -c http.sslVerify=false clone --branch ${branch} ${gitrepo} tmp

# check if cloning was successful
if [ ! -d tmp ] ; then echo "Cloning of remote repo failed." ; exit 1 ; fi

# switch to the requested commit
echo
echo "Checkout requested commit"
cd tmp
git reset --hard ${commit}
echo

# cherrypick folders
target_code="src/OCE src/TOP cfgs/SHARED ext/IOIPSL tools/REBUILD_NEMO"

for code in ${target_code} ; do
    echo "Getting data from ${code}"
    if [ -d ../../${target_dir}/${code} ] ; then rm -rf ../../${target_dir}/${code} ; fi
    mkdir -p ../../${target_dir}/${code}
    mv ${code}/* ../../${target_dir}/${code}
done

# copy license file
echo "Copy license file"
mv LICENSE ../../${target_dir}

# Clean up 
cd ..
rm -rf tmp
if [ -d ${target_dir} ] ; then rm -rf ${target_dir} ; fi

# End
echo
echo "Operation completed"
