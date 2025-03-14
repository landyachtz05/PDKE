# PDKE

Pacific Drought Knowledge Exchange  
Repository created for sharing code related to PDKE products.

## Setup:
### Download and install R and RStudio (free version)

git clone https://github.com/landyachtz05/PDKE.git
R: <https://cran.r-project.org/>  
RStudio: <https://posit.co/products/open-source/rstudio/>

### Setup your credentials file

Create a credentials.json file in the project's root directory.  It should look something like this:  
> {  
>   "ENV_TYPE": "linux",  
>   "PROJ_LIB_VAL": "/software/r9/matlab/R2021a/toolbox/map/mapproj/projdata/proj.db",  
>   "RSCRIPT_PATH": "/usr/bin/Rscript",  
>   "bearer": ""  
> }  
Set ENV_TYPE to either "linux" or "windows" depending on the box type you are running this on.  
Edit the PROJ_LIB_VAL to be the path to wherever your proj.db file is (do a search, 'sudo find / -name "proj.db"', or get it via installing proj at: <https://proj.org/en/stable/index.html>).  
Edit RSCRIPT_PATH to point at your Rscript (try whereis Rscript). 
Ask the PDKE administrator for the bearer value.  

### Set up the folder structure

- Make a folder called results in the project's root directory. 
- Make a folder called CCVD in the project's root directory. 
- In the CCVD directory, make a new folder called CCVD_OUTPUTS.  The structure should be CCVD/CCVD_OUTPUTS. 
- In the CCVD directory, make a new folder called MINI_PPT.  The structure should be CCVD/MINI_PPT.  
- In the CCVD directory, make a new folder called MINI_Phase2.  The structure should be CCVD/MINI_Phase2.


### Download files from Google Drive (ask Admin for access)

<https://drive.google.com/drive/u/1/folders/1QwQU3ulooAeSRmm6rNC8JMGo22aXdgzN>
- Download the Shapefiles folder into the PDKESite directory, include the folder so the structure is PDKESite/Shapefiles.
- In Shapefiles folder, make a folder called SelectedPolygon, so the structure is PDKESite/Shapefiles/SelectedPolygon.
- In Shapefiles folder, make a folder called UserDefinedPolygon, so the structure is PDKESite/Shapefiles/UserDefinedPolygon.
- Download the data folder, and it's contents, into the CCVD folder so the structure is CCVD/data.
- Download the CCVD_INPUTS folder, and it's contents, into the CCVD folder so the structure is CCVD/CCVD_INPUTS.
- Download the IMAGE folder, and it's contents, into the CCVD folder so the structure is CCVD/IMAGE.
- Download the NEW_RF_MAPS folder, and it's contents, into the CCVD folder so the structure is CCVD/NEW_RF_MAPS.
- ln -s /srv/shiny-server/PDKE/PDKESite/Shapefiles /srv/shiny-server/PDKE/shapefile 


### if on an old jetstream instance that doesn't allow download, need to copy up the files:
> scp -r /Users/jgeis/PDKE_testing/PDKESite/Shapefiles exouser@149.165.154.114:/srv/shiny-server/PDKE/PDKESite/
> scp -r /Users/jgeis/PDKE_testing/CCVD/CCVD_INPUTS exouser@149.165.154.114:/srv/shiny-server/PDKE/CCVD/
> scp -r /Users/jgeis/PDKE_testing/CCVD/data exouser@149.165.154.114:/srv/shiny-server/PDKE/CCVD/
> scp -r /Users/jgeis/PDKE_testing/CCVD/IMAGE exouser@149.165.154.114:/srv/shiny-server/PDKE/CCVD/
> scp -r /Users/jgeis/PDKE_testing/CCVD/NEW_RF_MAPS exouser@149.165.154.114:/srv/shiny-server/PDKE/CCVD/

## To run the program

### Using RStudio on your local machine

- Install r and RStudio 
- Install r shiny: install.packages("shiny")
- Open server.R and click the run button.  The GUI opens and you either select a pre-existing shapefile or create a new one by clicking on the map.  
- Make sure the 'Name' and 'Short Name' fields are filled in.  
- Enter your email address and click "Generate Data" button.  
- This starts a process that takes about 40 minutes.  When it's finished running, look in the MINI_PPT folder and there will be a new pdf named something like: \<Name\>\_\<YYYY\>\_\<MM\>\_\<DD\>\_CCVD_Portfolio.pdf.  

### Using R server (website)

The R server is hosted on Jetstream and is currently called PDKE\_Shiny.  To access it, go to the following link:  
<https://jetstream2.exosphere.app/exosphere/projects/ae2152821a6a4d5d866e10698d616466/regions/IU/servers/c1f9d4fa-8ff0-4e43-b15a-326f30f19fa6>   

Note, as of 3/10/2025 we have a new Jetstream VM that we are going to use for setting up PDKE:
CIS240457: AI Agents on Jetstream2 Training for UH
<https://jetstream2.exosphere.app/exosphere/projects/fa8f488bff214b96baf622f56818cace/regions/IU/overview>

#### Setting up a new instance on Jetstream2:

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
- RStudio is installed by default
- sudo apt-get update
- sudo apt-get install r-base 
- RStudio is installed by default 
- install devtools: sudo apt install r-cran-devtools

> sudo R, 
> install.packages("rmarkdown", dep = TRUE)
> install.packages("leaflet", dep = TRUE)
> install.packages("shiny", dep = TRUE)
> install.packages("leaflet.extras", dep = TRUE)
> install.packages("sf", dep = TRUE)
> install.packages("jsonlite", dep = TRUE)
> install.packages("here", dep = TRUE)
>
> install.packages("gstat", dep = TRUE)
> install.packages("raster", dep = TRUE)
> install.packages("sp", dep = TRUE) 
> install.packages("maptools", dep = TRUE)
> install.packages("rgdal", dep = TRUE)
> install.packages("RColorBrewer", dep = TRUE)
> install.packages("gridExtra", dep = TRUE)
> install.packages("ggplot2", dep = TRUE)
> install.packages("grid", dep = TRUE)
> install.packages("data.table", dep = TRUE)
> install.packages("devtools", dep = TRUE)
> install.packages("DescTools", dep = TRUE)
> install.packages("lubridate", dep = TRUE)
> install.packages("rgeos", dep = TRUE)
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
> install.packages("tiff", dep = TRUE)
> install.packages("lmomco", dep = TRUE)
> install.packages("parallel", dep = TRUE)
> install.packages("SPEI", dep = TRUE)
> install.packages("sf", dep = TRUE)
> install.packages("ggpubr", dep = TRUE)
> install.packages("terrainr", dep = TRUE)
> install.packages("ggmap", dep = TRUE)
> install.packages("ggthemes", dep = TRUE)
> install.packages("zoo", dep = TRUE)
> install.packages("classInt", dep = TRUE)
> install.packages("jsonlite", dep = TRUE)
>
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
> install.packages("magick", dep = TRUE)
> install.packages("httr", dep = TRUE)
> install.packages("stringr", dep = TRUE)
> install.packages("zip", dep = TRUE)
> install.packages("tiff", repos="https://packagemanager.posit.co/cran/2023-10-13", dep = TRUE)

- install markdown: sudo su - -c "R -e \"library(devtools); install_github('rstudio/rmarkdown')\""
- install shiny r package: sudo su - -c "R -e \"install.packages('shiny', repos='https://cran.rstudio.com/')\""
- install shiny server: <https://posit.co/download/shiny-server/>
  - sudo apt-get install gdebi-core
  - wget https://download3.rstudio.org/ubuntu-20.04/x86_64/shiny-server-1.5.23.1030-amd64.deb
  - sudo gdebi shiny-server-1.5.23.1030-amd64.deb
- Clone repo into /srv/shiny-server
- edit the r server config file, should use /etc/shiny-server/shiny-server.conf but had to edit /opt/shiny-server/config/default.config before it would stick.
> listen 80;  
> location / {  
>   site_dir /srv/shiny-server/PDKE/PDKESite;  
>   log_dir /var/log/shiny-server;  
> }

<https://docs.posit.co/shiny-server/#stopping-and-starting>  
sudo systemctl start shiny-server  
sudo systemctl stop shiny-server  

You can restart the server with:  
sudo systemctl restart shiny-server  
This command will shutdown all running Shiny processes, disconnect all open connections, and re-initialize the server.

If you wish to reload the configuration but keep the server and all Shiny processes running without interruption, you can use the systemctl command to send a SIGHUP signal:  
sudo systemctl kill -s HUP --kill-who=main shiny-server  
This will cause the server to re-initialize, but will not interrupt the current processes or any of the open connections to the server.

You can check the status of the shiny-server service using:  
sudo systemctl status shiny-server



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

To see apache error log:
sudo less /var/log/apache2/error.log
Apache: /var/lib/apache2, /etc/apache2

config file values:
- /etc/shiny-server/shiny-server.conf  
- listen 3838 (change to 8080)  
- location / 
- site_dir /srv/shiny-server/PDKE/PDKESite;
- log_dir /var/log/shiny-server;  

/opt/shiny-server  
/var/logs/shiny-server.log  

Code is found at:  
/srv/shiny-server/PDKESite  
/home/exouser/workflow  

Tried using the workflow dir as it was supposed to let me automate github updates, but the automated github updates didn’t work (this was attempted via binders <https://docs.jetstream-cloud.org/ui/exo/binder/>).  Then I tried to point the shiny-server config file to this directory, but there were permission issues.  Finally had to resort to doing git updates in the workflow dir, then doing a copy of the affected files (typically just CCVD_portfolio_content.R, CCVD_portfolio_ppt.R, and the 3 r files in the PDKESite dir).  

Had to do an ln for shapefile access while in the main scripts when we started from PDKESite.  
/srv/shiny-server/shapefile -> PDKESite/Shapefiles/  
ln -s /srv/shiny-server/PDKE/PDKESite/Shapefiles /srv/shiny-server/PDKE/shapefile 

Shiny Server will use the [default configuration](https://github.com/rstudio/shiny-server/blob/master/config/default.config) unless an alternate configuration is provided at `/etc/shiny-server/shiny-server.conf`. Using the default configuration, Shiny Server will look for Shiny apps in `/srv/shiny-server/` and host them on port 3838.  

sudo cp -R ~/MY-APP /srv/shiny-server/



DEREK:
- Which version of R? Try set up anaconda env to match Derek's R
1.) rgdal and rgeos - talk to Matt

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


exouser@pdke:/var/log$ sudo find / -name sp.h
/software/u22/r/4.4.1-old/lib/R/library/sp/include/sp.h
/usr/lib/R/site-library/sp/include/sp.h

# this worked to install rgdal
download: 
https://cran.r-project.org/src/contrib/Archive/rgdal/
https://cran.r-project.org/src/contrib/Archive/rgeos/
sudo R
set working directory to where it was downloaded using setwd('/path/to/source/code/downloaded')
setwd('/home/exouser/Downloads')
install.packages("sp")
install.packages("rgdal_1.6-7.tar.gz",
                 repos=NULL,
                 type = 'source',
                 configure.args="--with-proj-include=/software/u22/r/4.4.1-old/lib/R/library/sp/include/")
# this doesn't work for rgeos
install.packages("rgeos_0.6-4.tar.gz",
                 repos=NULL,
                 type = 'source',
                 configure.args="--with-proj-include=/usr/lib/R/site-library/sp/include/")
                 
                 https://cran.r-project.org/src/contrib/Archive/rgeos/rgeos_0.6-4.tar.gz