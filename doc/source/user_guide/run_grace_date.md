run_grace_date.py
=================

 - Wrapper program for running GRACE date and months programs
 - Reads GRACE/GRACE-FO index files from `podaac_grace_sync.py` or `gfz_isdc_grace_ftp.py`  
 - Reads dates of each GRACE/GRACE-FO file and assigns the month number  
 - Creates an index of dates for GRACE/GRACE-FO files  
 - Creates an index of dates for all GRACE/GRACE-FO processing centers  

#### Calling Sequence
```bash
python run_grace_date.py --release RL06 --mode 0o775
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/run_grace_date.py)

#### Command Line Options  
 - `-D X`, `--directory X`: GRACE/GRACE-FO working data directory
 - `-c X`, `--center X`: GRACE/GRACE-FO Processing Center (CSR,GFZ,JPL) 
    * `'CSR'`: University of Texas Center for Space Research  
    * `'GFZ'`: German Research Centre for Geosciences (GeoForschungsZentrum)
    * `'JPL'`: Jet Propulsion Laboratory    
 - `-r X`, `--release X`: GRACE/GRACE-FO data release (RL04,RL05,RL06)
 - `-V`, `--verbose`: verbose output of program progress
 - `-M X`, `--mode X`: permissions mode of output file
