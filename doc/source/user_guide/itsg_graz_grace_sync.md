itsg_graz_grace_sync.py
=======================

- Syncs GRACE/GRACE-FO and auxiliary data from the [ITSG GRAZ server](https://www.tugraz.at/institute/ifg/downloads/gravity-field-models)
- Creates an index file for each data product

#### Calling Sequence
```bash
python itsg_graz_grace_sync.py --directory <path_to_grace_directory> \
    --release Grace2018 Grace_operational --lmax 60
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/itsg_graz_grace_sync.py)

#### Command Line Options
- `-D X`, `--directory X`: Working Data Directory
- `-r X`, `--release X`: GRAZ Data Releases to sync
    * `'Grace2014'`
    * `'Grace2016'`
    * `'Grace2018'`
    * `'Grace_operational'`
- `--lmax X`: Maximum degree and order of GRAZ products
    * `60`
    * `96`
    * `120`
- `-t X`, `--timeout X`: Timeout in seconds for blocking operations
- `-l`, `--log`: Output log file
- `-L`, `--list`: Only print files that are to be transferred
- `-C`, `--clobber`: Overwrite existing data in transfer
- `--checksum`: Compare hashes to check if overwriting existing data
- `-M X`, `--mode X`: Permission mode of directories and files synced
