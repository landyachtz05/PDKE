# script intended to be run as a monthly cron job to get statewide rainfall
# and temperature data

# 1.) Look in ~/pdke/ccvd/new_rf_maps/statewide_rf_mm/rf_mm and get the 
# latest year/month by parsing the newest file (last one in a sorted list)
# 2.) Get today's year and month.
# 3.) Make calls to the HCDP API to get all the files between the most recent and now
#     ex: https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/production/rainfall/new/month/statewide/data_map/2024/rainfall_new_month_statewide_data_map_2024_05.tif

import os
from datetime import datetime
import requests

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

# Example usage
# Specify the directory to find existing rainfall files and save any newly downloaded rainfall files
directory_path = "/Users/jgeis/Work/PDKE/CCVD/NEW_RF_MAPS/statewide_rf_mm/rf_mm/"
last_file = list_files_sorted(directory_path)
print("Last file in the directory:", last_file)

# get the year and month of the most recent download
# file name format: rainfall_new_month_statewide_data_map_yyyy_mm.tif
file_parts = last_file.split('_')
file_year = int(file_parts[-2])  # Extract yyyy
file_month = int(file_parts[-1].split('.')[0])  # Extract mm 
print("File year:", file_year)
print("File Month:", file_month)

# Get current year and month
current_date = datetime.now()
current_year = int(current_date.strftime("%Y"))
current_month = int(current_date.strftime("%m"))
print("Current year:", current_year)
print("Current month:", current_month)


# Loop through each month between file's year and month and current year and month
for year in range(file_year, current_year + 1):
  print("year: ", year)
  #start_month = file_month if year == file_year else 1
  start_month = file_month + 1 if year == file_year else 1  # Start from the month after file's month
  end_month = current_month if year == current_year else 12
  print("start_month: ", start_month)
  print("end_month: ", end_month)
  for month in range(start_month, end_month + 1):
      print("Year:", year, "Month:", month)

      # Construct URL with year and month values
      url = f"https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/production/rainfall/new/month/statewide/data_map/{year}/rainfall_new_month_statewide_data_map_{year}_{month:02d}.tif"
      
      # Print URL for verification
      print("URL:", url)

      # # Here you can add code to download the file using the URL
      response = requests.get(url)
      if response.status_code == 200:
          # Write the downloaded content to a file in the specified directory
          with open(os.path.join(directory_path, f"rainfall_new_month_statewide_data_map_{year}_{month:02d}.tif"), "wb") as f:
              f.write(response.content)
          print("File downloaded successfully")
      else:
          print("Failed to download file")


