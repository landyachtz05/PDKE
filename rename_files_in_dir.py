import os

import os
import argparse

def rename_files(directory, old_prefix, new_prefix):
    """Renames files in a directory based on prefix replacement."""

    try:
        for filename in os.listdir(directory):
            if filename.startswith(old_prefix):
                new_filename = new_prefix + filename[len(old_prefix):]
                old_filepath = os.path.join(directory, filename)
                new_filepath = os.path.join(directory, new_filename)

                try:
                    os.rename(old_filepath, new_filepath)
                    print(f"Renamed '{filename}' to '{new_filename}'")
                except OSError as e:
                    print(f"Error renaming '{filename}': {e}")

    except FileNotFoundError:
        print(f"Directory '{directory}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename files based on prefix.")
    parser.add_argument("directory", help="The directory containing the files.")
    parser.add_argument("old_prefix", help="The prefix to replace.")
    parser.add_argument("new_prefix", help="The new prefix.")
    args = parser.parse_args()

    directory_path = args.directory
    old_prefix = args.old_prefix
    new_prefix = args.new_prefix

    rename_files(directory_path, old_prefix, new_prefix)

    # Example usage (you can hardcode the directory if you prefer):
    # python rename_files_in_dir.py "/path/to/your/directory"  "original filename start" "new filename start"
    # python rename_files_in_dir.py "/Users/jgeis/Work/PDKE/CCVD/CCVD_OUTPUTS/Honouliuli National Historic Site" "Honouliuli National Historic Site " "Honouliuli_" 
