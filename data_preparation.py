"""Data preparation module for Lambda particle analysis.

This module handles reading CSV files and creating combined feather files
for efficient data loading.

Usage:
    python data_preparation.py file1.csv file2.csv -o output.feather
    python data_preparation.py data/*.csv -o combined.feather
    python data_preparation.py *5x41*.csv -o myoutput.feather

Real usage:
    python data_preparation.py /mnt/c/data/meson-structure/2025-08-root/*18x275*mcpart*zip -o /mnt/c/data/meson-structure/18x275_mcpart_500k.feather

"""

import argparse
import glob
import sys
import pandas as pd


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
        print(f"  Reading: {file}")
        df = pd.read_csv(file)
        df['event'] = df['event'] + offset
        offset = df['event'].max() + 1
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def create_combined_feather_file(files, output_file):
    """Read multiple CSV files and save as a single feather file.
    
    Args:
        files: List of file paths to CSV files
        output_file: Path for output feather file
        
    Returns:
        bool: True if successful, False otherwise
    """
    if len(files) == 0:
        print("Error: No input files provided")
        return False
    
    print(f"\nProcessing {len(files)} file(s)...")
    combined_df = concat_csvs_with_unique_events(files)
    
    print(f"\nSaving to: {output_file}")
    combined_df.to_feather(output_file)
    
    print(f"✓ Successfully created: {output_file}")
    print(f"  Total rows: {len(combined_df):,}")
    print(f"  Total columns: {len(combined_df.columns)}")
    
    return True


def expand_file_patterns(patterns):
    """Expand file patterns using glob.
    
    This is useful on Windows where the shell doesn't expand wildcards.
    
    Args:
        patterns: List of file patterns (may contain wildcards)
        
    Returns:
        list: Expanded list of file paths
    """
    files = []
    for pattern in patterns:
        expanded = glob.glob(pattern)
        if expanded:
            files.extend(expanded)
        else:
            # If no match, keep the original (might be a literal filename)
            files.append(pattern)
    return sorted(files)


def main():
    """Main execution for data preparation."""
    parser = argparse.ArgumentParser(
        description='Combine Lambda particle CSV files into feather format',
        epilog='''
                Examples:
                  %(prog)s file1.csv file2.csv -o output.feather
                  %(prog)s data/*.csv -o combined.feather
                  %(prog)s *5x41*.csv -o myoutput.feather
            ''')

    parser.add_argument('files', nargs='+', help='Input CSV file(s) to combine (wildcards supported)')
    parser.add_argument('-o', '--output', default='combined_all_lambda_data.feather',
                        help='Output feather file (default: combined_all_lambda_data.feather)')

    args = parser.parse_args()

    # Expand wildcards (important for Windows)
    input_files = expand_file_patterns(args.files)
    
    # Check if files exist
    if not input_files:
        print("Error: No files found matching the specified patterns")
        sys.exit(1)
    
    # Create combined feather file
    success = create_combined_feather_file(input_files, args.output)
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
