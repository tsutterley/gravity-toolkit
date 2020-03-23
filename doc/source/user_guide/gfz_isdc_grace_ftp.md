gfz_isdc_grace_ftp.py
=====================

 - Syncs GRACE/GRACE-FO and auxiliary data from the [GFZ Information System and Data Center (ISDC)](http://isdc.gfz-potsdam.de/grace-isdc/)
 - Syncs CSR/GFZ/JPL files for RL04/RL05/RL06 GAA/GAB/GAC/GAD/GSM (GAA and GAB are GFZ/JPL only)
 - Creates an index file for each data product

#### Calling Sequence
```bash
python gfz_isdc_grace_ftp.py --directory=<path_to_grace_directory> --release=RL06
```

#### Command Line Options
 - `-D X`, `--directory=X`: Working Data Directory
 - `-C X`, `--center=X`: GRACE Processing Center (CSR,GFZ,JPL)
 - `-R X`, `--release=X`: GRACE data releases to sync (RL05,RL06)
 - `--newsletters`: Sync GRACE Newsletters
 - `-L`, `--list`: Only print files that are to be transferred
 - `--clobber`: Overwrite existing data in transfer
 - `--checksum`: Compare hashes to check if overwriting existing data
 - `-M X`, `--mode=X`: Permission mode of directories and files synced
 - `-l`, `--log`: Output log file
