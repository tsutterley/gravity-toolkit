podaac_cumulus.py
=================

- Syncs GRACE/GRACE-FO data from [NASA JPL PO.DAAC Cumulus AWS S3 bucket](https://podaac.jpl.nasa.gov/cloud-datasets/about)
- S3 Cumulus syncs are only available in AWS instances in `us-west-2`
- Creates an index file for each data product

#### Calling Sequence
```bash
python podaac_cumulus.py --user <username> --directory <path_to_grace_directory> --release RL06
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/podaac_cumulus.py)

#### Command Line Options
- `-U X`, `--user X`: Username for NASA Earthdata Login
- `-W X`, `--password X`: Password for NASA Earthdata Login
- `-N X`, `--netrc X`: Path to .netrc file for authentication
- `-D X`, `--directory X`: Working Data Directory
- `-c X`, `--center X`: GRACE/GRACE-FO Processing Center (CSR,GFZ,JPL)
- `-r X`, `--release X`: GRACE/GRACE-FO Data Releases to sync (RL06)
- `-v X`, `--version X`: GRACE/GRACE-FO Level-2 Data Version to sync (0,1)
- `-t X`, `--timeout X`: Timeout in seconds for blocking operations
- `-C`, `--clobber`: Overwrite existing data in transfer
- `-M X`, `--mode X`: Permission mode of directories and files synced
- `-l`, `--log`: Output log file
