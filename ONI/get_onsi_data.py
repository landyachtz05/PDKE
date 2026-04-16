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
#   python3 get_onsi_data.py
#
# Note: If running in VSCode, use a separate terminal window rather than
# the integrated terminal as it was flaky in my environment.  Maybe it will be better for you.
#
# ======================================================================
# Requirements:
# ======================================================================
#
# - Python 3.x
# - Required packages: pandas, requests, beautifulsoup4
#   Install with: pip install pandas requests beautifulsoup4
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
# Last Updated: May 21, 2025
# ======================================================================
import pandas as pd
import requests
import os
import io

def fetch_oni_data():
    """Fetch the latest ONI data (monthly) from NOAA CPC ASCII text file."""
    url = "https://www.cpc.ncep.noaa.gov/data/indices/Rnino34.ascii.txt"
    
    print("Fetching URL:", url)
    response = requests.get(url)
    response.raise_for_status()

    # Read the text file using whitespace as the delimiter
    df = pd.read_csv(io.StringIO(response.text), sep=r'\s+')
    print(f"Loaded {len(df)} rows from the ASCII text file")

    # Clean up column names just in case of trailing spaces
    df.columns = [col.strip() for col in df.columns]
    print("Columns found:", df.columns.tolist())

    # Ensure required columns are present
    required_cols = ["YR", "MTH", "ANOM"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Could not find the expected '{col}' column in the text file.")

    # Create the 'Date' column by combining YR and MTH, adding '01' for the day
    df["Date"] = pd.to_datetime(
        df["YR"].astype(str) + "-" + df["MTH"].astype(str).str.zfill(2) + "-01", 
        errors="coerce"
    )

    # Rename ANOM to ONI
    df.rename(columns={"ANOM": "ONI"}, inplace=True)

    # Keep only the requested columns
    df = df[["Date", "ONI"]]

    # Drop any rows where Date failed to parse (e.g., text footers in NOAA files)
    df.dropna(subset=["Date"], inplace=True)

    # Replace missing values (-9999 or -99.99) with NaN
    df["ONI"] = df["ONI"].replace([-9999, -99.99], pd.NA)

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
