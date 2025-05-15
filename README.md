# PDKE

Pacific Drought Knowledge Exchange  
Repository created for sharing code related to PDKE products.


### Using R server (website)  

The R server is hosted on Jetstream and is currently called PDKE\_Shiny.  To access it, go to the following link:  
<https://jetstream2.exosphere.app/exosphere/projects/ae2152821a6a4d5d866e10698d616466/regions/IU/servers/6d871a80-9fe0-418f-aea4-b7278daa85ac>   

## To setup from scratch:

### Jetstream2 Setup:

- Login to Jetstream2: 
  - <https://jetstream2.exosphere.app/exosphere/projects>
- Use "ACCESS CI (XSEDE)" for identity provider, not "University of Hawaii"
- Select an allocation (e.g. 'RAPID: Tuning and Assessing Lahaina Wildfire Models with AI Enhanced Data' or 'CIS240457: AI Agents on Jetstream2 Training for UH')
- Create a New Instance:
  - Click on the red “Create” button on the right.
  - Select “Instance”
  - Choose the latest Ubuntu image. 
  - Enter the instance name
  - Select “Flavor” (m3.large for PDKE)
  - Select custom disk size and set to 500 GB 
  - For “Enable Web Desktop?” Select: “Yes”
  - Click "Create." 
- Access the Instance using the Web Desktop or SSH into the instance (e.g., ssh exouser@<instance-ip>) 

#### Install required libraries 
Open a terminal window and run the following commands.  The r-base and gdebi-core are loaded separately so you can move on to the next step while the install is running. 
```
sudo apt-get update
wget -P /home/exouser/Downloads https://cran.r-project.org/src/contrib/Archive/maptools/maptools_1.1-8.tar.gz
sudo apt-get install -y r-base gdebi-core
sudo apt-get install -y r-cran-devtools libgdal-dev git libproj-dev libgeos-dev libudunits2-dev cargo libxslt1-dev libpoppler-cpp-dev tesseract-ocr-eng libtesseract-dev libmagick++-dev
```

#### Install R dependencies 
Open another terminal window and run the following commands: 
```
sudo R
install.packages(c('devtools', "rmarkdown", "leaflet", "shiny", "Rcpp", "leaflet.extras", "jsonlite", "here", "gstat", "raster", "sp", "RColorBrewer", "gridExtra", "ggplot2", "grid", "data.table", "DescTools", "lubridate", "latticeExtra", "rasterVis", "plotrix", "plyr", "dplyr", "xts", "timeSeries", "ggfortify", "changepoint", "scales", "reshape2", "hydroTSM", "lmomco", "parallel", "SPEI", "ggpubr", "terrainr", "ggmap", "ggthemes", "zoo", "classInt", "magrittr", "tidyverse", "rvg", "knitr", "xtable", "flextable", "officer", "mschart", "purrr", "pdftools", "httr", "stringr", "zip", "magick", "sf"), dependencies = TRUE)
```

## Proxy and SSL Setup
#### While waiting - the previous command takes a while - you can continue with setting up the proxy and ssl 
Instructions for this were obtained from the following links:

- https://www.digitalocean.com/community/tutorials/how-to-install-nginx-on-ubuntu-20-04
- https://www.digitalocean.com/community/tutorials/how-to-use-certbot-standalone-mode-to-retrieve-let-s-encrypt-ssl-certificates-on-ubuntu-20-04
- https://www.digitalocean.com/community/tutorials/how-to-set-up-shiny-server-on-ubuntu-20-04

### Setup nginx:
```
sudo apt update;
sudo apt install nginx;

sudo ufw allow ssh;
sudo ufw allow 49528/tcp;
sudo ufw allow from any to 172.17.0.1 port 5901 proto tcp;
sudo ufw allow 'Nginx Full';
sudo ufw allow 3838;
sudo ufw enable;
sudo service nginx start;
```
Create the directory for your_domain as follows, using the -p flag to create any necessary parent directories:
```
sudo mkdir -p /var/www/ccvd.manoa.hawaii.edu/html
```
Assign ownership of the directory with the $USER environment variable:
```
sudo chown -R $USER:$USER /var/www/ccvd.manoa.hawaii.edu/html
```
Set permissions
```
sudo chmod -R 755 /var/www/ccvd.manoa.hawaii.edu
```
Create a sample index.html page using nano or your favorite editor:
```
sudo vi /var/www/ccvd.manoa.hawaii.edu/html/index.html
```
Add the following lines to the file:
```
<html>
    <head>
        <title>Welcome to your_domain!</title>
    </head>
    <body>
        <h1>Success!  The your_domain server block is working!</h1>
    </body>
</html>
```
Create a server block
```
sudo vi /etc/nginx/sites-available/ccvd.manoa.hawaii.edu
```
Paste the following configuration block
```
server {
        listen 80;
        listen [::]:80;
        root /var/www/ccvd.manoa.hawaii.edu/html;
        index index.html index.htm index.nginx-debian.html;
        server_name ccvd.manoa.hawaii.edu www.ccvd.manoa.hawaii.edu;
        location / {
                try_files $uri $uri/ =404;
        }
}
```
Enable the file by creating a link from it to the sites-enabled directory, which Nginx reads from during startup:
```
sudo ln -s /etc/nginx/sites-available/ccvd.manoa.hawaii.edu /etc/nginx/sites-enabled/
```
Edit nginx.conf.  Find the server_names_hash_bucket_size directive and remove the # symbol to uncomment the line. If you are using nano, you can quickly search for words in the file by pressing CTRL and w.
```
sudo vi /etc/nginx/nginx.conf
```
Verify no errors
```
sudo nginx -t
```
If the command fails with 'could not build the server_names_hash, you should increase server_names_hash_bucket_size: 64'.  Solution: increase hash by a factor or 2, so 32 gets changed to 64, 64 goes to 128.  Then reload nginx with the following:
```
sudo systemctl restart nginx
```

###  Get SSL certificate
```
sudo certbot --nginx -d ccvd.manoa.hawaii.edu
```

###  Setup shiny:
```
sudo su - -c "R -e \"install.packages('shiny', repos='http://cran.rstudio.com/')\""
sudo wget https://download3.rstudio.org/ubuntu-20.04/x86_64/shiny-server-1.5.23.1030-amd64.deb
sudo gdebi shiny-server-1.5.23.1030-amd64.deb
```
Edit nginx.conf
```
sudo vi /etc/nginx/nginx.conf
```
Add the following to the http block:
```
    # Map proxy settings for RStudio
    map $http_upgrade $connection_upgrade {
        default upgrade;
        '' close;
    }
```
Open your domain’s server block:
```
sudo vi /etc/nginx/sites-available/ccvd.manoa.hawaii.edu
```
Find and edit the 'location' block to be only the following:
```
location / {
       proxy_pass http://your_server_ip:3838;
       proxy_redirect http://your_server_ip:3838/ https://$host/;
       proxy_http_version 1.1;
       proxy_set_header Upgrade $http_upgrade;
       proxy_set_header Connection $connection_upgrade;
       proxy_read_timeout 20d;
}
```
Test the file for errors
```
sudo nginx -t
sudo unlink /etc/nginx/sites-enabled/default
sudo systemctl restart nginx
```

### PDKE & Shiny Server Setup 
Open another terminal and run the following command to install the required packages:  
```
sudo R
install.packages("tiff", repos="https://packagemanager.posit.co/cran/2023-10-13", dep = TRUE)
install.packages("/home/exouser/Downloads/maptools_1.1-8.tar.gz", repos = NULL, type = "source")
```

###  Install and Setup Shiny Server   
Back in the original terminal use: 
```
sudo wget https://download3.rstudio.org/ubuntu-20.04/x86_64/shiny-server-1.5.23.1030-amd64.deb
sudo gdebi shiny-server-1.5.23.1030-amd64.deb

cd /srv/shiny-server
sudo git clone https://github.com/landyachtz05/PDKE.git
sudo chown -R exouser:exouser /srv/shiny-server

cd /srv/shiny-server/PDKE
mkdir -p CCVD/CCVD_OUTPUTS CCVD/MINI_PPT CCVD/MINI_Phase2 Shapefiles/SelectedPolygon Shapefiles/UserDefinedPolygon
sudo chmod -R 777 CCVD/CCVD_OUTPUTS CCVD/MINI_PPT CCVD/MINI_Phase2
```

## Download files 

#### (preferred method) Download files from Google Drive (ask Admin for access to the Google Drive)  

<https://drive.google.com/drive/u/1/folders/1QwQU3ulooAeSRmm6rNC8JMGo22aXdgzN>  

- Download and unzip all the folders from the drive, the following assumes they were unzipped into /home/exouser/Downloads  
- Go into each all the folders except for the second CCVD_INPUTS folder and move their contents one level up into the /home/exouser/Downloads directory.  Delete all the old, now-empty directories and zip files (except for the second CCVD_INPUTS file).  
- Once that is done, you'll do the following, the actual commands are just below:
  - Move the Shapefiles folder into the PDKESite directory, include the folder so the structure is PDKESite/Shapefiles.
  - Move the data folder into the CCVD folder so the structure is CCVD/data.
  - Move the CCVD_INPUTS folder, and it's contents, into the CCVD folder so the structure is CCVD/CCVD_INPUTS.
  - Move the IMAGE folder into the CCVD folder so the structure is CCVD/IMAGE.
  - Move the NEW_RF_MAPS folder into the CCVD folder so the structure is CCVD/NEW_RF_MAPS.
  - Make links to the Shapefiles and MINI_PPT folders so they can be accessed from the web interface.
```
mv /home/exouser/Downloads/data /home/exouser/Downloads/CCVD_INPUTS /home/exouser/Downloads/IMAGE /home/exouser/Downloads/NEW_RF_MAPS /srv/shiny-server/PDKE/CCVD/  
mv /home/exouser/Downloads/Shapefiles /srv/shiny-server/PDKE/PDKESite/  
mv /home/exouser/Downloads/credentials.json /srv/shiny-server/PDKE/  
```
Unzip the second CCVD_INPUTS folder and copy the inner folder of the new folder up one level into the Downloads directory.  Delete the now-empty CCVD_INPUTS_....002 folder and the zip file.  
```
sudo cp -R /home/exouser//Downloads/CCVD_INPUTS/* /srv/shiny-server/PDKE/CCVD/CCVD_INPUTS/ 
rm -rf /home/exouser//Downloads/CCVD_INPUTS  
```
Make links to the Shapefiles and MINI_PPT folders so they can be accessed from the web interface.
```
sudo ln -s /srv/shiny-server/PDKE/PDKESite/Shapefiles /srv/shiny-server/PDKE/PDKESite/www/shapefile   
sudo ln -s /srv/shiny-server/PDKE/CCVD/MINI_PPT /srv/shiny-server/PDKE/PDKESite/www/results   
```
Make an empty file which will be used for logging.  Open it, save the empty file, and exit vi.
```
sudo vi /srv/shiny-server/PDKE/locations.csv
```
Fix permissions
```
sudo chown -R shiny:shiny /srv/shiny-server/PDKE/PDKESite /srv/shiny-server/PDKE/CCVD/CCVD_INPUTS /srv/shiny-server/PDKE/CCVD/CCVD_OUTPUTS /srv/shiny-server/PDKE/CCVD/MINI_PPT /srv/shiny-server/PDKE/CCVD/MINI_Phase2 /srv/shiny-server/PDKE/PDKESite/Shapefiles /srv/shiny-server/PDKE/CCVD/data /srv/shiny-server/PDKE/CCVD/IMAGE /srv/shiny-server/PDKE/CCVD/NEW_RF_MAPS 
sudo chown shiny:shiny locations.csv  
```

### Setup your credentials file  
Edit the credentials.json file in the project's root directory.  
```
sudo vi /srv/shiny-server/PDKE/credentials.json
```
It should look something like this:  
```
{  
  "ENV_TYPE": "linux",  
  "PROJ_LIB_VAL": "/usr/share/proj/",  
  "RSCRIPT_PATH": "/usr/bin/Rscript",  
  "bearer": "ask admin for this"  
}
```
- Set ENV_TYPE to either "linux" or "windows" depending on the box type you are running this on.  
- Edit the PROJ_LIB_VAL to be the path to wherever your proj.db file is (do a search, 'sudo find / -name "proj.db"', or get it via installing proj at: <https://proj.org/en/stable/index.html>).  
- Edit RSCRIPT_PATH to point at your Rscript (try whereis Rscript). 
- Ask the PDKE administrator for the bearer value.

### Setup the cronjobs

The cronjobs below will update rainfall (GetRainfallDataFromHCDP.py) and delete files in CCVD_OUTPUTS & MINI_PPT older than 7 days (delete_week_old_files.py)

Enter cronjobs editor with:  

```
sudo crontab -e
```

Select 1 for editor options and paste in the following:

```
0 0 1 * * /usr/bin/env python3 /srv/shiny-server/PDKE/GetRainfallDataFromHCDP.py >> /var/log/rainfall_update.log 2>&1  
0 2 * * * /usr/bin/env python3 /srv/shiny-server/PDKE/delete_week_old_files.py >> /var/log/delete_week_old_files.log 2>&1
```

Change cronjob timing as need, the above uses:
- 0 0 1 * * : occurs once every month 
- 0 2 * * * : occurs everyday at 2 AM

Use the following to check logs:
```
cat /var/log/rainfall_update.log  
cat /var/log/delete_week_old_files.log
```

## Useful information:

### Shiny Server commands
<https://docs.posit.co/shiny-server/#stopping-and-starting>  
```
sudo systemctl start shiny-server  
sudo systemctl stop shiny-server  
sudo systemctl restart shiny-server
```
If you wish to reload the configuration but keep the server and all Shiny processes running without interruption, you can use the systemctl command to send a SIGHUP signal.  This will cause the server to re-initialize, but will not interrupt the current processes or any of the open connections to the server.  
:  
```
sudo systemctl kill -s HUP --kill-who=main shiny-server
```
You can check the status of the shiny-server service using:  
```
sudo systemctl status shiny-server
```
View the logs at:
```
sudo less /var/log/shiny-server.log 
sudo tail -f /var/log/shiny-server/PDKESite-shiny-<date-time>.log
```

### UFW info:
```
sudo ufw app list;
sudo ufw status;
```

### Nginx info:
File locations and uses:  

- /var/www/html: The actual web content, which by default only consists of the default Nginx page you saw earlier, is served out of the /var/www/html directory. This can be changed by altering Nginx configuration files.  
- /etc/nginx: The Nginx configuration directory. All of the Nginx configuration files reside here.  
- /etc/nginx/nginx.conf: The main Nginx configuration file. This can be modified to make changes to the Nginx global configuration.  
- /etc/nginx/sites-available/: The directory where per-site server blocks can be stored. Nginx will not use the configuration files found in this directory unless they are linked to the sites-enabled directory. Typically, all server block configuration is done in this directory, and then enabled by linking to the other directory.  
- /etc/nginx/sites-enabled/: The directory where enabled per-site server blocks are stored. Typically, these are created by linking to configuration files found in the sites-available directory.  
- /etc/nginx/snippets: This directory contains configuration fragments that can be included elsewhere in the Nginx configuration. Potentially repeatable configuration segments are good candidates for refactoring into snippets.  
- /var/log/nginx/access.log: Every request to your web server is recorded in this log file unless Nginx is configured to do otherwise.  
- /var/log/nginx/error.log: Any Nginx errors will be recorded in this log.  

Commands:
```
sudo apt install net-tools
sudo systemctl stop nginx
sudo systemctl restart nginx
sudo systemctl start nginx
sudo systemctl reload nginx
sudo systemctl disable nginx
sudo systemctl enable nginx
sudo systemctl status nginx;
sudo service nginx status
```

### PDKE Debugging:
Error:
```
Error in (function (file = if (onefile) "Rplots.pdf" else "Rplot%03d.pdf",  : 
  cannot open file 'Rplots.pdf'
Calls: plot ... plot.magick-image -> plot -> plot.raster -> plot.new -> <Anonymous>
Execution halted
```
Solution:   
- You have a write permission issue somewhere.  Set directories (likely CCVD_OUPUT, CCVD_
MINI, and/or CCVD_Phase2) to 777.  This is a temporary fix, but it works for now.

### Sean Umida's PDKE setup instructions:
https://docs.google.com/document/d/1OXuDyXvY6pl52HOKsZnMY2UvB8dJUcW6UoqyAZ75R1Q/edit?tab=t.0  
