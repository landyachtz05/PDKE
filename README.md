# PDKE

Pacific Drought Knowledge Exchange  
Repository created for sharing code related to PDKE products.

## Setup:
### Download and install R and RStudio (free version)

R: <https://cran.r-project.org/>  
RStudio: <https://posit.co/products/open-source/rstudio/>

### Setup your credentials file

Create a credentials.json file in the project's root directory.  It should look something like this:  
> {  
>   "PROJ_LIB_VAL": "/software/r9/matlab/R2021a/toolbox/map/mapproj/projdata/proj.db",  
>   "RSCRIPT_PATH": "/usr/bin/Rscript",  
>   "bearer": ""  
> }  
Edit the PROJ_LIB_VAL to be the path to wherever your proj.db file is (do a search, 'sudo find / -name "proj.db"', or get it via installing proj at: <https://proj.org/en/stable/index.html>).  
Edit RSCRIPT_PATH to point at your Rscript (try whereis Rscript). 
Ask the PDKE administrator for the bearer value.  

### Set up the folder structure

- Make a folder called CCVD in the project's root directory. 
- In the CCVD directory, make a new folder called CCVD_OUTPUTS.  The structure should be CCVD/CCVD_OUTPUTS. 
- In the CCVD directory, make a new folder called MINI_PPT.  The structure should be CCVD/MINI_PPT.  
- In the CCVD directory, make a new folder called MINI_Phase2.  The structure should be CCVD/MINI_Phase2.  

### Download files from Google Drive (ask Admin for access)

<https://drive.google.com/drive/u/1/folders/1QwQU3ulooAeSRmm6rNC8JMGo22aXdgzN>
- Download the Shapefiles folder into the PDKESite directory, include the folder so the structure is PDKESite/Shapefiles.
- Download the data folder, and it's contents, into the CCVD folder so the structure is CCVD/data.
- Download the CCVD_INPUTS folder, and it's contents, into the CCVD folder so the structure is CCVD/CCVD_INPUTS.
- Download the IMAGE folder, and it's contents, into the CCVD folder so the structure is CCVD/IMAGE.
- Download the NEW_RF_MAPS folder, and it's contents, into the CCVD folder so the structure is CCVD/NEW_RF_MAPS.

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
- sudo R, 
  - install.packages("rmarkdown", dep = TRUE)
  - install.packages("leaflet", dep = TRUE)
  - install.packages("shiny", dep = TRUE)
  - install.packages("leaflet.extras", dep = TRUE)
  - install.packages("sf", dep = TRUE)
  - install.packages("jsonlite", dep = TRUE)
  - install.packages("here", dep = TRUE)
- install markdown: sudo su - -c "R -e \"library(devtools); install_github('rstudio/rmarkdown')\""
- install shiny r package: sudo su - -c "R -e \"install.packages('shiny', repos='https://cran.rstudio.com/')\""
- install shiny server: <https://posit.co/download/shiny-server/>
  - sudo apt-get install gdebi-core
  - wget https://download3.rstudio.org/ubuntu-20.04/x86_64/shiny-server-1.5.23.1030-amd64.deb
  - sudo gdebi shiny-server-1.5.23.1030-amd64.deb
- Clone repo into /srv/shiny-server
- edit the r server config file, should use </etc/shiny-server/shiny-server.conf> but had to edit cd /opt/shiny-server/config/default.config before it would stick.
  - listen 80;  
  - location / {  
    site_dir /srv/shiny-server/PDKE/PDKESite;  
    log_dir /var/log/shiny-server;  
  - }


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
- listen 3838 (to do, change to 8080)  
- location / (to do, try point at the PDKESite dir so can get rid of the subdir in the url, otherwise, move the shiny files into the project's root dir)  
- site_dir /srv/shiny-server;  
- log_dir /var/log/shiny-server;  

/opt/shiny-server  
/var/logs/shiny-server.log  

Code is found at:  
/srv/shiny-server/PDKESite  
/home/exouser/workflow  

Tried using the workflow dir as it was supposed to let me automate github updates, but the automated github updates didn’t work (this was attempted via binders <https://docs.jetstream-cloud.org/ui/exo/binder/>).  Then I tried to point the shiny-server config file to this directory, but there were permission issues.  Finally had to resort to doing git updates in the workflow dir, then doing a copy of the affected files (typically just CCVD_portfolio_content.R, CCVD_portfolio_ppt.R, and the 3 r files in the PDKESite dir).  

Had to do an ln for shapefile access while in the main scripts when we started from PDKESite.  
/srv/shiny-server/shapefile -> PDKESite/Shapefiles/  

Shiny Server will use the [default configuration](https://github.com/rstudio/shiny-server/blob/master/config/default.config) unless an alternate configuration is provided at `/etc/shiny-server/shiny-server.conf`. Using the default configuration, Shiny Server will look for Shiny apps in `/srv/shiny-server/` and host them on port 3838.  

sudo cp -R ~/MY-APP /srv/shiny-server/

