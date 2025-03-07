# PDKE

Pacific Drought Knowledge Exchange  
Repository created for sharing code related to PDKE products.

## Setup:
### Download and install R and RStudio (free version)

R: <https://cran.r-project.org/>  
RStudio: <https://posit.co/products/open-source/rstudio/>

### Setup your credentials file

Create a credentials.json file in the project's root directory.  It should look like this:  
> {  
>   "PROJ_LIB_VAL": "/opt/anaconda3/share/proj/",  
>   "RSCRIPT_PATH": "/Library/Frameworks/R.framework/Resources/bin/Rscript",  
>   "bearer": ""  
> }  
Edit the PROJ_LIB_VAL to be the path to wherever your proj.db file is (do a search).  
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
- Open server.R and click the run button.  The GUI opens and you either select a pre-existing shapefile or create a new one by clicking on the map.  
- Make sure the 'Name' and 'Short Name' fields are filled in.  
- Enter your email address and click "Generate Data" button.  
- This starts a process that takes about 40 minutes.  When it's finished running, look in the MINI_PPT folder and there will be a new pdf named something like: \<Name\>\_\<YYYY\>\_\<MM\>\_\<DD\>\_CCVD_Portfolio.pdf.  

### Using R server (website)

The R server is hosted on Jetstream and is currently called PDKE\_Shiny.  To access it, go to the following link:  
<https://jetstream2.exosphere.app/exosphere/projects/ae2152821a6a4d5d866e10698d616466/regions/IU/servers/c1f9d4fa-8ff0-4e43-b15a-326f30f19fa6>   

### Setting up a new instance on Jetstream2:

- Login to Jetstream2: <https://jetstream2.exosphere.app/exosphere/projects>
- Use "ACCESS CI (XSEDE)" for identity provider, not "University of Hawaii"
- Click on existing allocation (RAPID: Tuning and Assessing Lahaina Wildfire Models with AI Enhanced Data)
- Click on red “Create” button on the right.
- Select “Instance”
- Select the latest Ubuntu
- Enter the instance name
- Select “Flavor” (most of the small sites I’ve been making, I’ve used m3.small)
- Add to storage size if needed
- For “Enable Web Desktop?” Select: “Yes”
- Create
- Once created and ready, open web desktop, then open terminal and type the following
  - sudo apt install apache2
  - sudo adduser <username>
  - sudo passwd <username>
  - sudo groupadd <groupname>
  - sudo usermod -a -G <groupname> <username>
- If repo is private, generate a personal access token following instructions here: <https://stackoverflow.com/questions/2505096/clone-a-private-repository-github>  
- Figure out where you want to put the repo, then clone it.  Ideally it will be cloned into the shiny-server directory, but I had permission issues with that.  I ended up cloning it into the workflow directory, then copying the files over to the shiny-server directory.  Want to take another try at this.

### Running/Updating on Jetstream:

Site URL: <http://149.165.154.114:3838/PDKESite/>. Note: I expect this to change.    
GitHub: <https://github.com/landyachtz05/PDKE> 

To Access for editing: 
Login (2 methods): 
- ssh <username>@149.165.154.114
- Jetstream site (hoping this will change if we can git clone in /srv/shiny-server:
<https://jetstream2.exosphere.app/exosphere/projects>
  - Use "ACCESS CI (XSEDE)" for identity provider, not "University of Hawaii"
  - Click on existing allocation (RAPID: Tuning and Assessing Lahaina Wildfire Models with AI Enhanced Data)
  - Click on Instances
  - Click on PDKE-Shiny
  - Click on Web Desktop
Once in a terminal in either method:
- cd /home/exouser/workflow
- git fetch
- git status
- git pull  

Shiny Server: <https://support.posit.co/hc/en-us/articles/218499158-Shiny-Server-Quick-Start-Guides>

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







