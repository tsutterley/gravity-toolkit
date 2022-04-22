podaac_grace_sync.py
====================

- Syncs GRACE/GRACE-FO and auxiliary data from the [NASA JPL PO.DAAC Drive Server](https://podaac-tools.jpl.nasa.gov/drive)
- Syncs CSR/GFZ/JPL files for RL04/RL05/RL06 GAA/GAB/GAC/GAD/GSM (GAA and GAB are GFZ/JPL only)
- Gets the latest technical note (TN) files
- Gets the monthly GRACE/GRACE-FO newsletters
- Creates an index file for each data product

#### Calling Sequence
```bash
python podaac_grace_sync.py --user <username> --directory <path_to_grace_directory> --release RL06
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/podaac_grace_sync.py)

#### Command Line Options
- `-U X`, `--user X`: Username for NASA Earthdata Login
- `-W X`, `--webdav X`: WebDAV Password for JPL PO.DAAC Drive Login
- `-N X`, `--netrc X`: Path to .netrc file for authentication
- `-D X`, `--directory X`: Working Data Directory
- `-m X`, `--mission X`: Mission to sync between GRACE and GRACE-FO
   * `'grace'`
   * `'grace-fo'`
- `-c X`, `--center X`: GRACE/GRACE-FO Processing Center (CSR,GFZ,JPL)
- `-r X`, `--release X`: GRACE/GRACE-FO Data Releases to sync (RL06)
- `-v X`, `--version X`: GRACE/GRACE-FO Level-2 Data Version to sync (0,1)
- `--newsletters`: Sync GRACE Newsletters
- `-t X`, `--timeout X`: Timeout in seconds for blocking operations
- `-L`, `--list`: Only print files that are to be transferred
- `-C`, `--clobber`: Overwrite existing data in transfer
- `--checksum`: Compare hashes to check if overwriting existing data
- `-M X`, `--mode X`: Permission mode of directories and files synced
- `-l`, `--log`: Output log file
