# Written by Shubham Prakash on 26/11/2024
# Purpose: This script will filter out the sequence which succeed in all of the three tests of Allergen, Antigen(PENDING) and Signal P detection.
#Input: input_csv variable. I have named it output.csv. Feel free to change. Should contains csv after all three tests.
#Output: output_csv variable. I have named it filtered_sequence.csv. Feel free to change.
import pandas as pd

# Read the original CSV file
input_csv = 'outputs/converted.csv'  # INPUT CSV File. Update Accordingly.
  # INPUT CSV File. Change Accordingly.
df = pd.read_csv(input_csv)

# Apply the filter conditions
filtered_df = df[
    (df['Allergen Test'] == 'Non-Allergen') & 
    (df['Antigen Test'] == 'Pending') & #change this after proper implementation of antigen-test (vaxijen)
    (df['Signal P'] == 'Found')
]

# Save the filtered data to a new CSV file
output_csv = './outputs/filtered_sequence.csv'  # OUTPUT CSV File. Change Accordingly.
filtered_df.to_csv(output_csv, index=False)

print(f"Filtered data has been saved to {output_csv}")