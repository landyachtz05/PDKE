# script intended to be run as a monthly cron job to get statewide rainfall
# and temperature data

# Rainfall:
# 1.) Look in ~/pdke/ccvd/new_rf_maps/statewide_rf_mm/rf_mm and get the 
# latest year/month by parsing the newest file (last one in a sorted list)
# 2.) Get today's year and month.
# 3.) Make calls to the HCDP API to get all the files between the most recent and now
#     ex: https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/production/rainfall/new/month/statewide/data_map/2024/rainfall_new_month_statewide_data_map_2024_05.tif
#
# Temperature:
# 1.) Look in ~/pdke/ccvd/ccvd_inputs/air_temp/data_map_newer and open the folder
#     for the latest year.
# 2.) Get the most recent month by parsing it out of the name of the newest file (last one in a sorted list)
# 3.) Get today's year and month.
# 4.) Make calls to the HCDP API to get all the files between the most recent and now
#     ex: https:/ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/production/temperature/mean/month/statewide/data_map/2023/temperature_mean_month_statewide_data_map_2023_12.tif
#
# To run:
# python GetRainfallDataFromHCDP.py

import os
from datetime import datetime
from pathlib import Path
import requests

base_dir = str(Path(__file__).resolve().parent)

# Specify the directory to find existing rainfall files and save any newly downloaded rainfall files
rainfall_directory_path = base_dir + "/CCVD/NEW_RF_MAPS/statewide_rf_mm/rf_mm/"
#rainfall_directory_path = "/Users/jgeis/Work/PDKE/CCVD/NEW_RF_MAPS/statewide_rf_mm/rf_mm/"

# Specify the directory to find existing temperature files and save any newly downloaded temperature files
temperature_directory_path = base_dir + "/CCVD/CCVD_INPUTS/air_temp/data_map_newer/"
#temperature_directory_path = "/Users/jgeis/Work/PDKE/CCVD/CCVD_INPUTS/air_temp/data_map_newer/"

def list_files_sorted(directory):
  # Check if the directory exists
  if not os.path.exists(directory):
      return "Directory not found"
  # Get a list of all files in the directory
  files = os.listdir(directory)
  # Sort the list of files alphabetically
  sorted_files = sorted(files)
  # Return the last entry of the sorted list
  return sorted_files[-1]

# file name format: rainfall_new_month_statewide_data_map_yyyy_mm.tif
def get_year_from_filename(filename):
  file_parts = filename.split('_')
  file_year = int(file_parts[-2])  # Extract yyyy
  #print("File year:", file_year)  
  return file_year

# file name format: rainfall_new_month_statewide_data_map_yyyy_mm.tif
def get_month_from_filename(filename):
  file_parts = filename.split('_')
  file_month = int(file_parts[-1].split('.')[0])  # Extract mm 
  #print("File month:", file_month)  
  return file_month

# downlad either temperature or rainfall data to the given directory
# directory = the location of the existing rainfall or temperature data.
# filename = the name of the last file in the given directory
# is_rain = boolean, true if supposed to get rainfall data, false if getting temperature data.
def download_data(directory, filename, is_rain):
  file_year = get_year_from_filename(filename)
  file_month = get_month_from_filename(filename)
  #print("File year:", file_year)
  #print("File Month:", file_month)
  
  # Get current year and month
  current_date = datetime.now()
  current_year = int(current_date.strftime("%Y"))
  current_month = int(current_date.strftime("%m"))
  #print("Current year:", current_year)
  #print("Current month:", current_month)
  
  # Loop through each month between file's year and month and current year and month
  for year in range(file_year, current_year + 1):
    print("year: ", year)
    #start_month = file_month if year == file_year else 1
    start_month = file_month + 1 if year == file_year else 1  # Start from the month after file's month
    end_month = current_month if year == current_year else 12
    #print("start_month: ", start_month)
    #print("end_month: ", end_month)
    for month in range(start_month, end_month + 1):
        #print("Year:", year, "Month:", month)
  
        # Construct URL with year and month values
        url = f"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/production/rainfall/new/month/statewide/data_map/{year}/rainfall_new_month_statewide_data_map_{year}_{month:02d}.tif"
        if (not is_rain):
          url = f"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/production/temperature/mean/month/statewide/data_map/{year}/temperature_mean_month_statewide_data_map_{year}_{month:02d}.tif"

        
        # Print URL for verification
        print("URL:", url)
  
        # # Here you can add code to download the file using the URL
        response = requests.get(url)
        if response.status_code == 200:
            # Write the downloaded content to a file in the specified directory
            if (is_rain):
              with open(os.path.join(directory, f"rainfall_new_month_statewide_data_map_{year}_{month:02d}.tif"), "wb") as f:
                  f.write(response.content)
              print("File downloaded successfully")
            else: 
              with open(os.path.join(directory, f"temperature_mean_month_statewide_data_map_{year}_{month:02d}.tif"), "wb") as f:
                  f.write(response.content)
        else:
            print("Failed to download file")

def get_last_sorted_directory(directory_path):
    # Check if the directory exists
    if not os.path.exists(directory_path):
        return "Directory not found"

    # Get a list of all folders in the directory
    folders = [folder for folder in os.listdir(directory_path) if os.path.isdir(os.path.join(directory_path, folder))]

    # Sort the list of folders alphabetically
    sorted_folders = sorted(folders)

    # Return the last directory after sorting
    if sorted_folders:
        return sorted_folders[-1]
    else:
        return "No folders found"

  
#print("temperature_directory_path: ", temperature_directory_path)
#print("rainfall_directory_path: ", rainfall_directory_path)

# get the most recent rainfall file
last_rain_file = list_files_sorted(rainfall_directory_path)
#print("Last file in the directory:", last_rain_file)

# get the most recent temperature file
last_temp_dir = get_last_sorted_directory(temperature_directory_path)
#print("Last sorted directory:", last_temp_dir)
temp_dir = temperature_directory_path + last_temp_dir
last_temp_file = list_files_sorted(temp_dir)
#print("Last file in the directory:", last_temp_file)

# download any missing rainfall files
download_data(rainfall_directory_path, last_rain_file, True)
# download any missing temperature files
download_data(temp_dir, last_temp_file, False)


