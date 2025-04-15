# PDKE

Pacific Drought Knowledge Exchange  
Repository created for sharing code related to PDKE products.

## Setup:

### Setting up a new instance on Jetstream2:

- Login to Jetstream2: <https://jetstream2.exosphere.app/exosphere/projects>
- Use "ACCESS CI (XSEDE)" for identity provider, not "University of Hawaii"
- Click on existing allocation ('RAPID: Tuning and Assessing Lahaina Wildfire Models with AI Enhanced Data' or 'CIS240457: AI Agents on Jetstream2 Training for UH')
- Click on red “Create” button on the right.
- Select “Instance”
- Select the latest Ubuntu
- Enter the instance name
- Select “Flavor” (m3.medium for PDKE)
- For “Enable Web Desktop?” Select: “Yes”
- Create

> sudo apt-get update  
> sudo apt-get install r-base  
> sudo apt install r-base-core  
> sudo apt install r-cran-devtools  
> sudo apt-get install libtiff-dev  

> sudo R  
> install.packages(c("devtools", "rmarkdown", "leaflet", "shiny", "Rcpp", "leaflet.extras", "jsonlite", "here", "gstat", "raster", "sp", "maptools", "RColorBrewer", "gridExtra", "ggplot2", "grid", "data.table", "DescTools", "lubridate", "latticeExtra", "rasterVis", "plotrix", "plyr", "dplyr", "xts", "timeSeries", "ggfortify", "changepoint", "scales", "reshape2", "hydroTSM", "lmomco", "parallel", "SPEI", "ggpubr", "terrainr", "ggmap", "ggthemes", "zoo", "classInt", "magrittr", "tidyverse", "rvg", "knitr", "xtable", "flextable", "officer", "mschart", "purrr", "pdftools", "httr", "stringr", "zip", "magick", "sf"), dependencies = TRUE)  
> install.packages("tiff", repos="https://packagemanager.posit.co/cran/2023-10-13", dep = TRUE)

Alternatively, install R packages separately:
> sudo R  
> install.packages("devtools", dep = TRUE)  
> install.packages("rmarkdown", dep = TRUE)  
> install.packages("leaflet", dep = TRUE)  
> install.packages("shiny", dep = TRUE)  
> install.packages("leaflet.extras", dep = TRUE)  
> install.packages("jsonlite", dep = TRUE)  
> install.packages("here", dep = TRUE)  
> install.packages("gstat", dep = TRUE)  
> install.packages("raster", dep = TRUE)  
> install.packages("sp", dep = TRUE)   
> install.packages("maptools", dep = TRUE)  
> install.packages("RColorBrewer", dep = TRUE)  
> install.packages("gridExtra", dep = TRUE)  
> install.packages("ggplot2", dep = TRUE)  
> install.packages("grid", dep = TRUE)  
> install.packages("data.table", dep = TRUE)  
> install.packages("DescTools", dep = TRUE)  
> install.packages("lubridate", dep = TRUE)  
> install.packages("latticeExtra", dep = TRUE)  
> install.packages("rasterVis", dep = TRUE)  
> install.packages("plotrix", dep = TRUE)  
> install.packages("plyr", dep = TRUE)  
> install.packages("dplyr", dep = TRUE)  
> install.packages("xts", dep = TRUE)  
> install.packages("timeSeries", dep = TRUE)  
> install.packages("ggfortify", dep = TRUE)  
> install.packages("changepoint", dep = TRUE)  
> install.packages("scales", dep = TRUE)  
> install.packages("reshape2", dep = TRUE)  
> install.packages("hydroTSM", dep = TRUE)  
> install.packages("lmomco", dep = TRUE)  
> install.packages("parallel", dep = TRUE)  
> install.packages("SPEI", dep = TRUE)  
> install.packages("ggpubr", dep = TRUE)  
> install.packages("terrainr", dep = TRUE)  
> install.packages("ggmap", dep = TRUE)  
> install.packages("ggthemes", dep = TRUE)  
> install.packages("zoo", dep = TRUE)  
> install.packages("classInt", dep = TRUE)  
> install.packages("magrittr", dep = TRUE)  
> install.packages("tidyverse", dep = TRUE)  
> install.packages("rvg", dep = TRUE)  
> install.packages("knitr", dep = TRUE)  
> install.packages("xtable", dep = TRUE)  
> install.packages("flextable", dep = TRUE)  
> install.packages("officer", dep = TRUE)  
> install.packages("mschart", dep = TRUE)  
> install.packages("purrr", dep = TRUE)  
> install.packages("pdftools", dep = TRUE)  
> install.packages("httr", dep = TRUE)  
> install.packages("stringr", dep = TRUE)  
> install.packages("zip", dep = TRUE) 
> install.packages("magick", dep = TRUE)  
> install.packages("tiff", repos="https://packagemanager.posit.co/cran/2023-10-13", dep = TRUE)  
> install.packages("sf", dep = TRUE)  

- install markdown: 
> sudo su - -c "R -e \"library(devtools); install_github('rstudio/rmarkdown')\""  
- install shiny r package: 
> sudo R   
> install.packages('shiny', repos='https://cran.rstudio.com/')  
- install shiny server: <https://posit.co/download/shiny-server/>
> sudo apt-get install gdebi-core  
> sudo wget https://download3.rstudio.org/ubuntu-20.04/x86_64/shiny-server-1.5.23.1030-amd64.deb  
> sudo gdebi shiny-server-1.5.23.1030-amd64.deb  
- Clone repo into /srv/shiny-server  
> cd /srv/shiny-server  
> sudo git clone https://github.com/landyachtz05/PDKE.git  
> sudo chown -R exouser:exouser PDKE  
- edit the r server config file, should only have to edit /etc/shiny-server/shiny-server.conf but I had to edit /opt/shiny-server/config/default.config before it would stick.
> listen 80;  
> location / {  
>   site_dir /srv/shiny-server/PDKE/PDKESite;  
>   log_dir /var/log/shiny-server;  
> }


### Set up the folder structure  
> mkdir /srv/shiny-server/PDKE/CCVD  
> cd /srv/shiny-server/PDKE/CCVD  
> mkdir /srv/shiny-server/PDKE/CCVD/CCVD_OUTPUTS  
> sudo chmod -R 777 CCVD_OUTPUTS  
> mkdir /srv/shiny-server/PDKE/CCVD/MINI_PPT  
> sudo chmod -R 777 MINI_PPT  
> mkdir /srv/shiny-server/PDKE/CCVD/MINI_Phase2  
> sudo chmod -R 777 MINI_Phase2  

### Download files from Google Drive (ask Admin for access)  

<https://drive.google.com/drive/u/1/folders/1QwQU3ulooAeSRmm6rNC8JMGo22aXdgzN>  
- Download and unzip all the folders from the drive, the following assumes they were unzipped into /home/exouser/Downloads  
- Go into each all the folders except for the second CCVD_INPUTS folder and move their contents one level up into the /home/exouser/Downloads directory.  Delete all the old, now-empty directories and zip files (except for the second CCVD_INPUTS file).  
- Once that is done, you'll do the following, the actual commands are just below.  
- Move the Shapefiles folder into the PDKESite directory, include the folder so the structure is PDKESite/Shapefiles.  
- Move the data folder into the CCVD folder so the structure is CCVD/data.  
- Move the CCVD_INPUTS folder, and it's contents, into the CCVD folder so the structure is CCVD/CCVD_INPUTS.  
- Move the IMAGE folder into the CCVD folder so the structure is CCVD/IMAGE.  
- Move the NEW_RF_MAPS folder into the CCVD folder so the structure is CCVD/NEW_RF_MAPS.  
- Make links to the Shapefiles and MINI_PPT folders so they can be accessed from the web interface.  
> mv /home/exouser/Downloads/Shapefiles /srv/shiny-server/PDKE/PDKESite/  
> mv /home/exouser/Downloads/data /srv/shiny-server/PDKE/CCVD/  
> mv /home/exouser/Downloads/CCVD_INPUTS /srv/shiny-server/PDKE/CCVD/  
> mv /home/exouser/Downloads/IMAGE /srv/shiny-server/PDKE/CCVD/  
> mv /home/exouser/Downloads/NEW_RF_MAPS /srv/shiny-server/PDKE/CCVD/  
> sudo ln -s /srv/shiny-server/PDKE/PDKESite/Shapefiles /srv/shiny-server/PDKE/PDKESite/www/shapefile   
> sudo ln -s /srv/shiny-server/PDKE/CCVD/MINI_PPT /srv/shiny-server/PDKE/PDKESite/www/results   
> sudo chown -R shiny:shiny /srv/shiny-server/PDKE/PDKESite/Shapefiles  
> sudo chown -R shiny:shiny /srv/shiny-server/PDKE/CCVD/data  
> sudo chown -R shiny:shiny /srv/shiny-server/PDKE/CCVD/CCVD_INPUTS  
> sudo chown -R shiny:shiny /srv/shiny-server/PDKE/CCVD/IMAGE  
> sudo chown -R shiny:shiny /srv/shiny-server/PDKE/CCVD/NEW_RF_MAPS  
- Make an empty file which will be used for logging  
> sudo vi /srv/shiny-server/PDKE/locations.csv   
- Save the empty file and exit vi.  
> sudo chown shiny:shiny locations.csv  
- Unzip the second CCVD_INPUTS folder and copy the inner folder of the new folder up one level into the Downloads directory.  Delete the now-empty CCVD_INPUTS_....002 folder and the zip file.  
> sudo cp -R ~/Downloads/CCVD_INPUTS/* .  
> rm -rf ~/Downloads/CCVD_INPUTS  

### Setup your credentials file  

Edit the credentials.json file in the project's root directory.  It should look something like this:  
> {  
>   "ENV_TYPE": "linux",  
>   "PROJ_LIB_VAL": "/software/r9/matlab/R2021a/toolbox/map/mapproj/projdata/",  
>   "RSCRIPT_PATH": "/usr/bin/Rscript",  
>   "bearer": ""  
> }  
Set ENV_TYPE to either "linux" or "windows" depending on the box type you are running this on.  
Edit the PROJ_LIB_VAL to be the path to wherever your proj.db file is (do a search, 'sudo find / -name "proj.db"', or get it via installing proj at: <https://proj.org/en/stable/index.html>).  
Edit RSCRIPT_PATH to point at your Rscript (try whereis Rscript). 
Ask the PDKE administrator for the bearer value.

### Setup the cronjobs

The cronjobs below will update rainfall (GetRainfallDataFromHCDP.py) and delete files in CCVD_OUTPUTS & MINI_PPT older than 7 days (delete_week_old_files.py)

Enter cronjobs editor with:  
> sudo crontab -e

Select 1 for editor options and paste in the following:
> 0 0 1 * * /usr/bin/env python3 /srv/shiny-server/PDKE/GetRainfallDataFromHCDP.py >> /var/log/rainfall_update.log 2>&1  
> 0 2 * * * /usr/bin/env python3 /srv/shiny-server/PDKE/delete_week_old_files.py >> /var/log/delete_week_old_files.log 2>&1

Change cronjob timing as need, the above uses:
- 0 0 1 * * : occurs once every month 
- 0 2 * * * : occurs everyday at 2 AM

Use the following to check logs:
> cat /var/log/rainfall_update.log  
> cat /var/log/delete_week_old_files.log

### if on an old jetstream instance that doesn't allow download, need to copy up the files:  
> scp -r /Users/jgeis/PDKE/PDKESite/Shapefiles exouser@149.165.154.114:/srv/shiny-server/PDKE/PDKESite/  
> scp -r /Users/jgeis/PDKE/CCVD/CCVD_INPUTS exouser@149.165.154.114:/srv/shiny-server/PDKE/CCVD/  
> scp -r /Users/jgeis/PDKE/CCVD/data exouser@149.165.154.114:/srv/shiny-server/PDKE/CCVD/  
> scp -r /Users/jgeis/PDKE/CCVD/IMAGE exouser@149.165.154.114:/srv/shiny-server/PDKE/CCVD/  
> scp -r /Users/jgeis/PDKE/CCVD/NEW_RF_MAPS exouser@149.165.154.114:/srv/shiny-server/PDKE/CCVD/  

<https://docs.posit.co/shiny-server/#stopping-and-starting>  
> sudo systemctl start shiny-server  
> sudo systemctl stop shiny-server  

You can restart the server with:  
> sudo systemctl restart shiny-server  
This command will shutdown all running Shiny processes, disconnect all open connections, and re-initialize the server.  

If you wish to reload the configuration but keep the server and all Shiny processes running without interruption, you can use the systemctl command to send a SIGHUP signal:  
> sudo systemctl kill -s HUP --kill-who=main shiny-server  
This will cause the server to re-initialize, but will not interrupt the current processes or any of the open connections to the server.  

You can check the status of the shiny-server service using:  
> sudo systemctl status shiny-server  

### If Apache is running, kill it  
sudo /etc/init.d/apache2 stop  

TODO: Possibly, hopefully, no longer valid  

#### Running/Updating on Jetstream:  

Site URL: <http://149.165.154.114:3838/PDKESite/>. Note: I expect this to change.    
GitHub: <https://github.com/landyachtz05/PDKE> 

To Access for editing:  
Login (2 methods): 
- ssh <username>@149.165.154.114  
- Jetstream site (hoping this will change if we can git clone in /srv/shiny-server:  
<https://jetstream2.exosphere.app/exosphere/projects>  
  - Use "ACCESS CI (XSEDE)" for identity provider, not "University of Hawaii"  
  - Click on existing allocation ('RAPID: Tuning and Assessing Lahaina Wildfire Models with AI Enhanced Data' or 'CIS240457: AI Agents on Jetstream2 Training for UH')  
  - Click on Instances  
  - Click on PDKE-Shiny  
  - Click on Web Desktop  
Once in a terminal in either method:  
- cd /home/exouser/workflow  
- git fetch  
- git status  
- git pull    

Shiny Server: <https://support.posit.co/hc/en-us/articles/218499158-Shiny-Server-Quick-Start-Guides>  
<https://posit.co/download/shiny-server/>  

config file values:  
- /etc/shiny-server/shiny-server.conf   
- listen 3838 (change to 8080)  
- location / 
- site_dir /srv/shiny-server/PDKE/PDKESite;  
- log_dir /var/log/shiny-server;  

/opt/shiny-server  

/var/log/shiny-server.log  
/var/log/shiny-server/PDKESite-shiny-20250331-222658-41981.log  

Code is found at:  
/srv/shiny-server/PDKESite  
/home/exouser/workflow  

Tried using the workflow dir as it was supposed to let me automate github updates, but the automated github updates didn’t work (this was attempted via binders <https://docs.jetstream-cloud.org/ui/exo/binder/>).  Then I tried to point the shiny-server config file to this directory, but there were permission issues.  Finally had to resort to doing git updates in the workflow dir, then doing a copy of the affected files (typically just CCVD_portfolio_content.R, CCVD_portfolio_ppt.R, and the 3 r files in the PDKESite dir).  

Had to do an ln for shapefile access while in the main scripts when we started from PDKESite.  
/srv/shiny-server/shapefile -> PDKESite/Shapefiles/  
ln -s /srv/shiny-server/PDKE/PDKESite/Shapefiles /srv/shiny-server/PDKE/shapefile 

Shiny Server will use the [default configuration](https://github.com/rstudio/shiny-server/blob/master/config/default.config) unless an alternate configuration is provided at `/etc/shiny-server/shiny-server.conf`. Using the default configuration, Shiny Server will look for Shiny apps in `/srv/shiny-server/` and host them on port 3838.  

sudo cp -R ~/MY-APP /srv/shiny-server/  

---

### Using R server (website)  

The R server is hosted on Jetstream and is currently called PDKE\_Shiny.  To access it, go to the following link:  
<https://jetstream2.exosphere.app/exosphere/projects/ae2152821a6a4d5d866e10698d616466/regions/IU/servers/c1f9d4fa-8ff0-4e43-b15a-326f30f19fa6>   

Note, as of 3/10/2025 we have a new Jetstream VM that we are going to use for setting up PDKE:  
CIS240457: AI Agents on Jetstream2 Training for UH  
<https://jetstream2.exosphere.app/exosphere/projects/fa8f488bff214b96baf622f56818cace/regions/IU/overview>  

---
### Running program via RStudio on your local machine  

- Install r and RStudio   
- Install r shiny: install.packages("shiny")  
- Go down to 'Setup your credentials file' and follow the instructions there, then come back here.  
- Open server.R and click the run button.  The GUI opens and you either select a pre-existing shapefile or create a new one by clicking on the map.   
- Make sure the 'Name' and 'Short Name' fields are filled in.   
- Enter your email address and click "Generate Data" button.  
- This starts a process that takes about 40 minutes.  When it's finished running, look in the MINI_PPT folder and there will be a new pdf named something like: \<Name\>\_\<YYYY\>\_\<MM\>\_\<DD\>\_CCVD_Portfolio.pdf.  

---


DEREK:   
- Which version of R? Try set up anaconda env to match Derek's R  
2.) Koheo 1-2_2025_03_13_02_39_25  PDKE:  define Mean_CLIM   
  Koheo 1-2_2025_03_13_02_39_25  PDKE:  INPUTS_FOLDER: /srv/shiny-server/PDKE//CCVD/CCVD_INPUTS/   
  Error in (function (file = if (onefile) "Rplots.pdf" else "Rplot%03d.pdf",  :   
    cannot open file 'Rplots.pdf'   
  Calls: plot ... .plotraster2 -> .rasterImagePlot -> <Anonymous> -> <Anonymous>  
  Execution halted  
  Solved by setting PDKE permissions to 777, can't let this continue, need to figure out where the permissions issue is.  It's in MINI_Phase2  
3.) email address for "contact us"  
4.) Forest Reserves and Sanctuaries isn't showing anything.  

----------------------------
Trying to figure out stack overflow issue on jetstream.

pj_obj_create
ulimit -s 16384

Error: C stack usage  7969476 is too close to the limit

7969812
7972820 (3008)
7976932 (4112)
7969476

exouser@pdke-shiny:/srv/shiny-server/PDKE$ ulimit -a
real-time non-blocking time  (microseconds, -R) unlimited
core file size              (blocks, -c) 0
data seg size               (kbytes, -d) unlimited
scheduling priority                 (-e) 0
file size                   (blocks, -f) unlimited
pending signals                     (-i) 120101
max locked memory           (kbytes, -l) 3848984
max memory size             (kbytes, -m) unlimited
open files                          (-n) 1024
pipe size                (512 bytes, -p) 8
POSIX message queues         (bytes, -q) 819200
real-time priority                  (-r) 0
stack size                  (kbytes, -s) unlimited
cpu time                   (seconds, -t) unlimited
max user processes                  (-u) 120101
virtual memory              (kbytes, -v) unlimited
file locks                          (-x) unlimited


* soft nofile 65536
* hard nofile 65536
* soft nproc 65536
* hard nproc 65536





           
------------

Haleaha_2025_03_28_02_20_21  PDKE:  Frazier et al 2016 
Warning: sf layer has inconsistent datum (+proj=longlat +datum=NAD83 +no_defs).
Need '+proj=longlat +datum=WGS84'





If cairo freetype2 are already installed, check that 'pkg-config' is in your
PATH and PKG_CONFIG_PATH contains a cairo freetype2.pc file. If pkg-config
is unavailable you can set INCLUDE_DIR and LIB_DIR manually via:
R CMD INSTALL --configure-vars='INCLUDE_DIR=... LIB_DIR=...'






The program is likely using a significant amount of disk space due to the following reasons:

1. Large Raster Files
The program processes numerous raster files (e.g., climate data, elevation, rainfall, etc.) and creates intermediate raster outputs (e.g., masked, cropped, transformed rasters). These intermediate files can consume a lot of disk space, especially if they are not deleted after use.
Solution:

Use in-memory processing where possible by setting rasterOptions(tmpdir = "RAM") or explicitly specifying options(rasterTmpDir = "/path/to/temp") to control where temporary files are stored.
Delete intermediate raster files after they are no longer needed using removeTmpFiles() or unlink().
Example:

2. Generated Plots
The program generates numerous plots (e.g., maps, graphs) and saves them as PNG files. These files can take up significant space, especially if they are high-resolution.
Solution:

Reduce the resolution of the plots by lowering the dpi value in png() calls.
Save only essential plots and skip unnecessary ones.
Example:

3. Duplicated Data
The program creates multiple versions of the same data (e.g., monthly, seasonal, and annual rasters). These duplicates can significantly increase disk usage.
Solution:

Avoid saving intermediate results unless absolutely necessary.
Use in-memory operations for temporary calculations.
Example:

4. Temporary Files
The raster and rgdal packages often create temporary files during processing. These files are stored in the system's temporary directory and can accumulate over time.
Solution:

Set a custom temporary directory and clean it up after the script finishes.
Use removeTmpFiles() to delete temporary files created by the raster package.
Example:

5. Large Output Files
The program writes CSV files, shapefiles, and other outputs for each analysis unit. These files can accumulate and take up significant space.
Solution:

Compress output files where possible (e.g., use .zip or .gz for CSV files).
Save only essential outputs and avoid redundant files.
Example:

6. High-Resolution Data
The program may be working with high-resolution raster data, which increases file sizes.
Solution:

Downsample raster data to a lower resolution if high precision is not required.
Example:

7. Unused Variables
The program creates many variables (e.g., monthly rasters, cropped rasters) that are not reused later.
Solution:

Remove unused variables from memory using rm() and call gc() to free up memory.
Example:

Summary of Recommendations:
Use in-memory processing and clean up temporary files (removeTmpFiles()).
Reduce plot resolution and save only essential plots.
Avoid saving intermediate results unless necessary.
Compress output files (e.g., CSV, shapefiles).
Downsample raster data if high resolution is not required.
Remove unused variables and call gc() to free memory.
By implementing these changes, you can significantly reduce the disk space usage of the program.


- 'alpha' must be between 0 and 1.
geom_polygon(data = Coast, aes(x = long, y = lat, group = group), col = "black", fill = "aliceblue", alpha = 0.5, size = 1)



PDKE: runs, but has 'error' status on instance.  Also maxes out disk space with 2 runs.  Each run takes 15G
PDKE2: Error: node stack overflow
PDKE3: (clone of PDKE on new instance) runs, but maxes out CPUs and uses 15G of disk space per run.


Unzip them all (except CCVD_INPUTS #2 
ls air_temp/data_map_newer/1990





1: packages ?grid?, ?parallel? are base packages, and should not be updated
2: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?tiff? had non-zero exit status
3: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?gifski? had non-zero exit status
4: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?webp? had non-zero exit status
5: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?xslt? had non-zero exit status
6: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?terra? had non-zero exit status
7: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?odbc? had non-zero exit status
8: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?pdftools? had non-zero exit status
9: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?magick? had non-zero exit status
10: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?sf? had non-zero exit status
11: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?gdtools? had non-zero exit status
12: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?equatags? had non-zero exit status
13: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?tesseract? had non-zero exit status
14: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?lwgeom? had non-zero exit status
15: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?raster? had non-zero exit status
16: In install.packages(c("devtools", "rmarkdown", "leaflet",  ... :
  installation of package ?flextable? had non-zero exit status

