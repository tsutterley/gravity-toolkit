gfz_isdc_grace_ftp.py
=====================

- Syncs GRACE/GRACE-FO and auxiliary data from the [GFZ Information System and Data Center (ISDC)](http://isdc.gfz-potsdam.de/grace-isdc/)
- Syncs CSR/GFZ/JPL files for RL04/RL05/RL06 GAA/GAB/GAC/GAD/GSM (GAA and GAB are GFZ/JPL only)
- Creates an index file for each data product

#### Calling Sequence
```bash
python gfz_isdc_grace_ftp.py --directory <path_to_grace_directory> --release RL06
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/gfz_isdc_grace_ftp.py)

#### Command Line Options
- `-D X`, `--directory X`: Working Data Directory
- `-m X`, `--mission X`: Mission to sync between GRACE and GRACE-FO
   * `'grace'`
   * `'grace-fo'`
- `-c X`, `--center X`: GRACE/GRACE-FO Processing Center (CSR,GFZ,JPL)
- `-r X`, `--release X`: GRACE/GRACE-FO Data Releases to sync (RL05,RL06)
- `--newsletters`: Sync GRACE Newsletters
- `-L`, `--list`: Only print files that are to be transferred
- `-C`, `--clobber`: Overwrite existing data in transfer
- `--checksum`: Compare hashes to check if overwriting existing data
- `-M X`, `--mode X`: Permission mode of directories and files synced
- `-l`, `--log`: Output log file
