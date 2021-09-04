esa_costg_swarm_sync.py
=======================

- Syncs Swarm gravity field products from the [ESA Swarm Science Server](https://earth.esa.int/eogateway/missions/swarm/data)
- Creates an index file for each data product

#### Calling Sequence
```bash
python esa_costg_swarm_sync.py --directory <path_to_grace_directory>
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/esa_costg_swarm_sync.py)

#### Command Line Options
- `-D X`, `--directory X`: Working Data Directory
- `-t X`, `--timeout X`: Timeout in seconds for blocking operations
- `-l`, `--log`: Output log file
- `-L`, `--list`: Only print files that are to be transferred
- `-C`, `--clobber`: Overwrite existing data in transfer
- `--checksum`: Compare hashes to check if overwriting existing data
- `-M X`, `--mode X`: Permission mode of directories and files synced
