"""Data preparation module for Lambda particle analysis.

This module handles reading CSV files and creating combined feather files
for efficient data loading.
"""

import pandas as pd
import glob


def concat_csvs_with_unique_events(files):
    """Load and concatenate CSV files with globally unique event IDs.
    
    Args:
        files: List of file paths to CSV files
        
    Returns:
        pd.DataFrame: Combined dataframe with unique event IDs
    """
    dfs = []
    offset = 0

    for file in files:
        df = pd.read_csv(file)
        df['event'] = df['event'] + offset
        offset = df['event'].max() + 1
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def create_combined_feather_file(input_pattern, output_file):
    """Read multiple CSV files and save as a single feather file.
    
    Args:
        input_pattern: Glob pattern for input CSV files
        output_file: Path for output feather file
        
    Returns:
        bool: True if successful, False otherwise
    """
    files = sorted(glob.glob(input_pattern))
    
    if len(files) == 0:
        print(f"No files found matching pattern: {input_pattern}")
        print("Check the path and file name pattern.")
        return False
    
    print(f"Found {len(files)} files to process")
    combined_df = concat_csvs_with_unique_events(files)
    combined_df.to_feather(output_file)
    print(f"Successfully created: {output_file}")
    print(f"Total rows: {len(combined_df)}")
    
    return True


def main():
    """Main execution for data preparation."""
    # Example usage
    input_pattern = 'data\\k_lambda_5x41_5000evt_*.mcpart_lambda.csv.zip'
    output_file = 'combined_all_lambda_data.feather'
    
    create_combined_feather_file(input_pattern, output_file)


if __name__ == "__main__":
    main()
