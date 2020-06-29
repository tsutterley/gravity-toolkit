podaac_webdav.py
====================

 - Retrieves and prints a user's PO.DAAC WebDAV credentials  
 - WebDAV credentials can be used in the [PO.DAAC sync program](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/podaac_grace_sync.py)  

#### Calling Sequence
```bash
python podaac_webdav.py --user=<username> 
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/scripts/podaac_webdav.py)

#### Command Line Options
 - `-U X`, `--user=X`: Username for NASA Earthdata Login
 - `-N X`, `--netrc=X`: Path to .netrc file for authentication
 - `-A`, `--append`: Append .netrc file instead of printing