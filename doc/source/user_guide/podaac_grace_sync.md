podaac_grace_sync.py
====================

 - Syncs GRACE/GRACE-FO and auxiliary data from the [NASA JPL PO.DAAC Drive Server](https://podaac-tools.jpl.nasa.gov/drive)  
 - Syncs CSR/GFZ/JPL files for RL04/RL05/RL06 GAA/GAB/GAC/GAD/GSM (GAA and GAB are GFZ/JPL only)
 - Gets the latest technical note (TN) files
 - Gets the monthly GRACE/GRACE-FO newsletters
 - Creates an index file for each data product

#### Calling Sequence
```bash
python podaac_grace_sync.py --user=<username> --directory=<path_to_grace_directory> --release=RL06
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/podaac_grace_sync.py)

#### Command Line Options
 - `-U X`, `--user=X`: Username for NASA Earthdata Login
 - `-N X`, `--netrc=X`: Path to .netrc file for authentication
 - `-D X`, `--directory=X`: Working Data Directory
 - `-C X`, `--center=X`: GRACE Processing Center (CSR,GFZ,JPL)
 - `-R X`, `--release=X`: GRACE data releases to sync (RL05,RL06)
 - `--newsletters`: Sync GRACE Newsletters
 - `-L`, `--list`: Only print files that are to be transferred
 - `--clobber`: Overwrite existing data in transfer
 - `--checksum`: Compare hashes to check if overwriting existing data
 - `-M X`, `--mode=X`: Permission mode of directories and files synced
 - `-l`, `--log`: Output log file
