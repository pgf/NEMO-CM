Repository for the NEMO model used in CMCC-CM

The original NEMO code is available at https://forge.nemo-ocean.eu/nemo/nemo.git

This repository is structured as follow:
1) utils     -> contains a set of scripts needed to manage NEMO code inside CMCC-CM 
2) drivers   -> contains the NUOPC cap code
3) cfgs      -> contains the configurations file for each allowed NEMO grid
4) interface -> contains the NEMO code modified to run in CMCC-CM. This overwrites the original NEMO code files.
5) nemocore  -> contains the original NEMO code required in CMCC-CM and set to a specific commit


Update NEMOCORE

The script utils/download_nemo.sh
allows to update the nemocore folder to a specific NEMO commit.
This operation leads to an updated version of the code stored in the nemocore folder.
Note that the .h90 files need to be the same in the nemocore and interface folders.
Any specific change required in the interface folder needs to be done manually.

