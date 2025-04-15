#Script which checks CCVD_OUTPUTS and MINI_PPT and deletes files older than 7 days 
#!/usr/bin/env python3
import os
from pathlib import Path
import time
from datetime import datetime
import sys

#Defined absolute paths to directories
OUTPUTS_DIR = Path("/srv/shiny-server/PDKE/CCVD/CCVD_OUTPUTS")
MINI_PPT_DIR = Path("/srv/shiny-server/PDKE/CCVD/MINI_PPT")

# Age threshold (one week) in seconds
AGE_THRESHOLD = 7 * 86400

#Deletes files older than 7 days within given directory
def delete_old_files(directory):
    #Checks for given directory and logs error if not found
    if not directory.exists():
        print(f"Error: Directory not found: {directory}", file=sys.stderr)
        return

    #Returns current time (since epoch) in seconds
    current_time = time.time()
    
    # Walks through directory and its subdirectories
    for root, _, files in os.walk(directory):
        for filename in files:
            filepath = Path(root) / filename
            try:
                # Gets last modification time of file
                mtime = os.path.getmtime(filepath)
                age = current_time - mtime

                #Deletes file if last modification (or creation) older than 7 days
                if age > AGE_THRESHOLD:
                    try:
                        os.remove(filepath)
                        print(f"Deleted: {filepath}")
                    except OSError as e:
                        print(f"Error deleting {filepath}: {e}", file=sys.stderr)
            except OSError as e:
                print(f"Error accessing {filepath}: {e}", file=sys.stderr)
                continue

def main():
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Logs and calls function to delete files in given directory
    print(f"Deleting files older than 7 days in {OUTPUTS_DIR} on {current_time}")
    delete_old_files(OUTPUTS_DIR)
    
    print(f"Deleting files older than 7 days in {MINI_PPT_DIR} on {current_time}")
    delete_old_files(MINI_PPT_DIR)

if __name__ == "__main__":
    main()
