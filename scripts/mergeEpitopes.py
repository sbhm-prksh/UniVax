import pandas as pd


def merge_epitopes(file_paths):
    """Merges epitopes from multiple CSV files into a single DataFrame.

    Args:
        file_paths: A list of file paths to the CSV files.

    Returns:
        A pandas DataFrame containing the merged epitopes.
    """

    dfs = []
    for file_path in file_paths:
        try:
            df = pd.read_csv(file_path)

            # Determine cell type based on filename
            cell_type = file_path.split('/')[-1].split('_')[1]

            # Add cell type column
            df['Cell Type'] = cell_type

            # Handle missing columns for non-TH types
            if cell_type != 'TH(IFN)':
                for col in ['Result', 'Allele', 'IFM_Gamma']:
                    if col not in df.columns:
                        df[col] = 'NA'

            # Reorder columns for consistent output
            df = df[['Cell Type', 'Protein ID', 'Epitope', 'Allergen Test', 'Toxicity Test', 'Antigen Test', 'Score', 'Allele', 'Result', 'IFM_Gamma']]

        except FileNotFoundError:
            print(f"Error: File not found - {file_path}")
            continue

        dfs.append(df)

    # Concatenate DataFrames (handles potential missing columns)
    merged_df = pd.concat(dfs, ignore_index=True, sort=False)

    return merged_df


# Example usage
file_paths = [
    './outputs/filtered_bCell.csv',
    './outputs/filtered_TC(A1).csv',
    './outputs/filtered_TC(B58).csv',
    './outputs/filtered_TH(IFN).csv'
]

merged_epitopes = merge_epitopes(file_paths)

# Save the merged DataFrame to a CSV file
merged_epitopes.to_csv('./outputs/epitopes.csv', index=False)