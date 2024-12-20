# Written By Shubham Prakash on 23/10/24
# Purpose: To do the Signal P identification using Phobius
# This code will ask user to enter the number of sequence they want to proceed with. This functionality was not present in our S1_Algpred.py.
# Feel free to use any csv consisting of any number of sequences.
# Input: A CSV file with columns like- Protein ID, Sequence, Allergen Test etc.UPDATE inputFile variable accordingly. So copy this input from either S0_Fasta2CSV or S1_Algpred(recommended).
# Output: Will modify the Signal P column to Found, if SignalP is present, else not found. We need Found one.

import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
import sys
######################################################################
# Function to send POST request and get the result of Signal P presence.
######################################################################
def get_signalP_result(sequence, id, index):
    print("----------------------------------------------------")
    print(f"{index + 1}. Attempting: SignalP Detection for: {id}")
    payload = {
        "protseq": f">{id}\n{sequence}",
        "format": 'nog',
    }
    url = "https://phobius.sbc.su.se/cgi-bin/predict.pl"
    
    r = requests.post(url, data=payload)
    
    if r.status_code == 200:
        soup = BeautifulSoup(r.text, 'html.parser')
        try:
            text = soup.find("pre").get_text()
            if "signal" in text.lower():
                result = "Found"
            else:
                result = "Not Found"
        except IndexError:
            result = "ERROR:DATA EXTRACTION"
        print(f"{index + 1}. Success   : SignalP Detection for: {id}")
    else:
        result = "ERROR:SERVER"
        print(f"{index + 1}. Failed  : SignalP Detection for: {id}")
    return result


######################################################################
# DRIVER CODE STARTS
######################################################################
# Input CSV file name
def main():
    inputFile = 'outputs/converted.csv'  # INPUT CSV File. Update Accordingly.
    df = pd.read_csv(inputFile)

    # # Ask the user how many sequences they want to process
    # num_sequences = int(input(f"Total sequences available: {len(df)}\nHow many sequences do you want to process? "))

    # # Ensure the number of sequences to process is within the available range
    # num_sequences = min(num_sequences, len(df))
    start_time = time.time()

    # Process only the specified number of sequences
    for index, row in df.head(sequence_count).iterrows():
        sequence = row['Sequence']
        id = row['Protein ID']
        result = get_signalP_result(sequence, id, index)
        df.at[index, 'Signal P'] = result

    # Save the updated DataFrame to the same CSV
    df.to_csv(inputFile, index=False)

    print("----------------------------------------------------")
    print(f"_______________________________\n| COMPLETED: SignalP Detection |\n| (Check {inputFile}.csv) |")
    end_time = time.time()
    total_time = end_time - start_time
    print(f"|    Time Taken: {total_time:.2f} sec     |\n|_____________________________|")
if __name__ == "__main__":
    sequence_count = int(sys.argv[1])  # Sequence count passed from Flask
    # output_file = sys.argv[2]  # Path to the CSV file
    main()

