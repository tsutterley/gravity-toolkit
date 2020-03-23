NASA Earthdata
==============

#### NASA Data Distribution Centers  
The NASA Earth Science Data Information Systems Project funds and operates [12 Distributed Active Archive Centers (DAACs)](https://earthdata.nasa.gov/about/daacs) throughout the United States.  These centers have recently transitioned from ftp to https servers.
The https updates are designed to increase performance and improve security during data retrieval. NASA Earthdata uses [OAuth2](https://wiki.earthdata.nasa.gov/pages/viewpage.action?pageId=71700485), an approach to authentication that protects your personal information.  
- https://urs.earthdata.nasa.gov/documentation  
- https://wiki.earthdata.nasa.gov/display/EL/Knowledge+Base  

#### PO.DAAC
The [Physical Oceanography Distributed Active Archive Center (PO.DAAC)](https://podaac.jpl.nasa.gov/) provides data and related information pertaining to the physical processes and conditions of the global oceans, including measurements of ocean winds, temperature, topography, salinity, circulation and currents, and sea ice.  If any problems contact JPL PO.DAAC support at [podaac@podaac.jpl.nasa.gov](mailto:podaac@podaac.jpl.nasa.gov) or the NASA EOSDIS support team [support@earthdata.nasa.gov](mailto:support@earthdata.nasa.gov).  

#### WebDAV  
PO.DAAC uses passwords generated using the Web Distributed Authoring and Versioning (WebDAV) API.  This password is created at the [PO.DAAC Drive](https://podaac-tools.jpl.nasa.gov/drive) website.  Use this password rather than your Earthdata password when retrieving data from PO.DAAC.  [More information](https://podaac-tools.jpl.nasa.gov/drive/help).

#### Steps to Sync from PO.DAAC
1. [Register with NASA Earthdata Login system](https://urs.earthdata.nasa.gov/users/new)  
2. [After registering, login to the system](https://urs.earthdata.nasa.gov/home)
3. Add `PO.DAAC Drive OPS` [application to Earthdata](https://wiki.earthdata.nasa.gov/display/EL/How+To+Pre-authorize+an+application)  
4. Retrieve [WebDAV password](https://podaac-tools.jpl.nasa.gov/drive/) to access PO.DAAC servers and sync data

#### Other Data Access Examples   
-  [Curl and Wget](https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget)   
-  [Python](https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python)    
