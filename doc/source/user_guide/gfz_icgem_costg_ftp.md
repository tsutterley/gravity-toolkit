gfz_icgem_costg_ftp.py
======================

- Syncs GRACE/GRACE-FO/Swarm COST-G data from the [GFZ International Centre for Global Earth Models (ICGEM)](http://icgem.gfz-potsdam.de/)
- Creates an index file for each data product

#### Calling Sequence
```bash
python gfz_icgem_costg_ftp.py --directory <path_to_grace_directory> \
    --mission Grace Grace-FO Swarm
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/gfz_icgem_costg_ftp.py)

#### Command Line Options
- `-D X`, `--directory X`: Working Data Directory
- `-m X`, `--mission X`: Mission to sync between GRACE, GRACE-FO and Swarm
    * `'Grace'`
    * `'Grace-FO'`
    * `'Swarm'`
- `-t X`, `--timeout X`: Timeout in seconds for blocking operations
- `-l`, `--log`: Output log file
- `-L`, `--list`: Only print files that are to be transferred
- `-C`, `--clobber`: Overwrite existing data in transfer
- `--checksum`: Compare hashes to check if overwriting existing data
- `-M X`, `--mode X`: Permission mode of directories and files synced
