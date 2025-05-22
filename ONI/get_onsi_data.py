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
from bs4 import BeautifulSoup
import io
import csv
import os

def fetch_oni_data():
    """Fetch the latest ONI data from NOAA website."""
    url = "https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php"
    
    temp_html_file = "oni_page.html"
    
    try:
        # Get the webpage content
        print("Fetching URL:", url)
        response = requests.get(url)
        response.raise_for_status()
        
        # Save HTML content for debugging (optional)
        with open(temp_html_file, "w") as f:
            f.write(response.text)
        print("Saved HTML content for debugging")
        
        # Parse the HTML
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Find all tables
        tables = soup.find_all('table')
        print(f"Found {len(tables)} tables on the page")
        
        # Find the table with ONI data
        data_table = None
        for i, table in enumerate(tables):
            # Print first row text of each table for debugging
            first_row = table.find('tr')
            if first_row:
                first_row_text = first_row.get_text().strip()
                print(f"Table {i} first row: {first_row_text[:50]}...")
                
                # Check if this looks like our ONI table
                if "DJF" in first_row_text and "JFM" in first_row_text and "Year" in first_row_text:
                    data_table = table
                    print(f"Found ONI table at index {i}")
                    break
        
        if not data_table:
            # Try alternate method - look for a table with years as first column
            for i, table in enumerate(tables):
                rows = table.find_all('tr')
                if len(rows) > 1:  # Need at least header + one data row
                    cells = rows[1].find_all('td')
                    if len(cells) > 0:
                        # Check if first cell has a year (1950-2025)
                        year_text = cells[0].get_text().strip()
                        print(f"Table {i}, first data cell: {year_text}")
                        try:
                            year = int(year_text)
                            if 1950 <= year <= 2025:
                                data_table = table
                                print(f"Found ONI table by year detection at index {i}")
                                break
                        except ValueError:
                            pass
        
        if not data_table:
            raise ValueError("Could not find the ONI data table on the webpage")
        
        # Extract the data
        rows = data_table.find_all('tr')
        print(f"Found {len(rows)} rows in the ONI table")
        
        # Define headers 
        headers = ['Year', 'DJF', 'JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA', 'JAS', 'ASO', 'SON', 'OND', 'NDJ']
        
        # Extract data from each row
        data = []
        for i, row in enumerate(rows):
            # Skip header row
            if i == 0:
                continue
                
            cols = row.find_all('td')
            
            # Skip completely empty rows
            if len(cols) == 0:
                print(f"Skipping row {i}: no columns")
                continue
            
            # Convert to proper format
            row_data = []
            
            # Extract year
            year_text = cols[0].get_text().strip()
            try:
                year = int(year_text)
                row_data.append(str(year))
            except ValueError:
                print(f"Skipping row {i}: invalid year format '{year_text}'")
                continue
            
            # Handle rows with incomplete data (like current year with only some months)
            if len(cols) < 13:
                print(f"Row {i} has only {len(cols)} columns, will pad with empty values")
                # Add available monthly values
                for j in range(1, len(cols)):
                    cell = cols[j].get_text().strip()
                    if cell:
                        try:
                            val = float(cell)
                            row_data.append(f"{val:.1f}")
                        except ValueError:
                            row_data.append("")
                    else:
                        row_data.append("")
                
                # Pad missing months with empty values
                for j in range(len(cols), 13):
                    row_data.append("")
            else:
                # Normal row with all 13 columns
                for j in range(1, 13):
                    cell = cols[j].get_text().strip()
                    if cell:
                        try:
                            val = float(cell)
                            row_data.append(f"{val:.1f}")
                        except ValueError:
                            row_data.append("")
                    else:
                        row_data.append("")
            
            data.append(row_data)
            if i < 5 or i > len(rows) - 7:  # Only print first 5 and last 6 rows
                print(f"Row {i}: {row_data}")
        
        # Create DataFrame
        df = pd.DataFrame(data, columns=headers)
        print(f"Created DataFrame with {len(df)} rows")
        
        return df
    
    finally:
        # Clean up temporary HTML file whether the function succeeds or fails
        if os.path.exists(temp_html_file):
            try:
                os.remove(temp_html_file)
                print(f"Removed temporary HTML file: {temp_html_file}")
            except Exception as e:
                print(f"Warning: Failed to remove temporary HTML file: {e}")

def format_dataframe_for_csv(df):
    """Format the dataframe to match the expected CSV format."""
    # Keep Year as simple integers without decimal points
    df['Year'] = df['Year'].astype(int)
    return df

def main():

    # Get the directory containing the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Go up one level to find project root (if needed)
    project_root = os.path.dirname(script_dir)

    # Create paths relative to project root
    csv_path = os.path.join(project_root, "ONI", "enso_oni_raw.csv")
    
    try:
        # Fetch new data
        print("Fetching latest ONI data...")
        new_data = fetch_oni_data()
        
        # Format the data to match the expected CSV format
        new_data = format_dataframe_for_csv(new_data)
        
        # Output sample of the new data for debugging
        print("\nSample of data to be written:")
        print(new_data.head())
        print(f"Total rows: {len(new_data)}")
        
        # Confirm the latest year data
        latest_year = new_data['Year'].astype(float).max()
        print(f"Latest year in data: {latest_year}")
        latest_year_data = new_data[new_data['Year'] == latest_year]
        print("Latest year data:")
        print(latest_year_data)
        
        # Save data to CSV
        print(f"\nOverwriting file: {csv_path}")
        new_data.to_csv(csv_path, index=False)
        
        print("CSV file updated successfully!")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()