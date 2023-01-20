# NEMO-CM
Repository for the NEMO model used in CMCC-CM

The original NEMO code is available at https://forge.nemo-ocean.eu/nemo/nemo.git

To update the NEMO code you need to:

1) add the gitlab repository as a remote 
git remote add upstream https://forge.nemo-ocean.eu/nemo/nemo.git

2) retrieve the gitlab 
git fetch upstream

3) merge the actual code with the updated one available on gitlab
git merge upstream/branch_4.2

4) correct the merge conflicts

5) delete the folders not needed in the CMCC-CM version of NEMO4.2
CMCC-CM needs only
	ext/IOIPSL
	src/OCE

6) push back to the repository the updated version
git push 

