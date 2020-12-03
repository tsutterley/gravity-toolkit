cnes_grace_sync.py
==================

 - Syncs GRACE/GRACE-FO and auxiliary data from the [CNES Server](http://grgs.obs-mip.fr/grace)  
 - Syncs CNES/GRGS files for RL01/RL02/RL03/RL04/RL05
 - Creates an index file for each data product

#### Calling Sequence
```bash
python cnes_grace_sync.py --directory <path_to_grace_directory> --release RL05
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/cnes_grace_sync.py)

#### Command Line Options
 - `-D X`, `--directory X`: Working Data Directory
 - `-r X`, `--release X`: GRACE/GRACE-FO Data Releases to sync (RL04,RL05)
 - `-C`, `--clobber`: Overwrite existing data in transfer
 - `-M X`, `--mode X`: Permission mode of directories and files synced
 - `-l`, `--log`: Output log file
