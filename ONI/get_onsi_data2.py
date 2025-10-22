# ======================================================================
# ENSO ONI Data Fetcher
# ======================================================================
#
# This script fetches the latest Oceanic Niño Index (ONI) data from the 
# NOAA Climate Prediction Center website and saves it to a CSV file.
#
# The ONI is one of the primary indices used to monitor El Niño and 
# La Niña events, calculated as the 3-month running mean of SST anomalies 
# in the Niño 3.4 region (5°N-5°S, 120°-170°W).
#
# ======================================================================
# Usage:
# ======================================================================
#
# Run the script from the command line:
#   python3 get_onsi_data2.py
#
# Note: If running in VSCode, use a separate terminal window rather than
# the integrated terminal as it was flaky in my environment.  Maybe it will be better for you.
#
# ======================================================================
# Requirements:
# ======================================================================
#
# - Python 3.x
# - Required packages: pandas, requests
#   Install with: pip install pandas requests
#
# ======================================================================
# Output:
# ======================================================================
#
# The script will update the file: <project_path>/PDKE/ONI/enso_oni_raw.csv
# with the latest ONI values from 1950 to the present.
#
# Data format: Year, followed by 12 columns of 3-month running means
# (DJF, JFM, FMA, MAM, AMJ, MJJ, JJA, JAS, ASO, SON, OND, NDJ)
#
# ======================================================================
# Maintainer: Jennifer Geis and Derek Ford
# Last Updated: October 22, 2025
# ======================================================================
import pandas as pd
import requests
import os
import io

def fetch_oni_data():
    """Fetch the latest ONI data (monthly) from NOAA PSL CSV."""
    url = "https://psl.noaa.gov/data/correlation/oni.csv"
    
    print("Fetching URL:", url)
    response = requests.get(url)
    response.raise_for_status()

    # Read CSV directly from the response text
    df = pd.read_csv(io.StringIO(response.text))
    print(f"Loaded {len(df)} rows from the new ONI CSV")

    # Clean up column names (in case of extra spaces)
    df.columns = [col.strip() for col in df.columns]
    
    # Rename for consistency (optional)
    if "ONI from CPC" in df.columns:
        df.rename(columns={"ONI from CPC": "ONI"}, inplace=True)

    # Convert Date to datetime
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")

    # Replace missing values (-9999) with NaN
    df["ONI"].replace(-9999, pd.NA, inplace=True)

    return df

def main():
    # Get the directory containing the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)

    # Create paths relative to project root
    csv_path = os.path.join(project_root, "ONI", "enso_oni_raw.csv")

    try:
        print("Fetching latest ONI data...")
        new_data = fetch_oni_data()

        print("\nSample of data to be written:")
        print(new_data.head())
        print(f"Total rows: {len(new_data)}")

        # Save to CSV
        print(f"\nOverwriting file: {csv_path}")
        os.makedirs(os.path.dirname(csv_path), exist_ok=True)
        new_data.to_csv(csv_path, index=False)
        
        print("CSV file updated successfully!")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
