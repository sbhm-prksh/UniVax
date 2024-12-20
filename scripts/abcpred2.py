# Written by Shubham Prakash on 26/11/2024
# Purpose: B Cell Epitope detection after doing Allergen, Antigen and Signal P Detection using ABCPred. AND THEN WILL RUN ALLERGEN, ANTIGEN, TOXIN TEST FOR THE ACHIEVED EPITOPES
# Input: csv_file variable. I have named it filtered_sequence.csv. Feel free to change.
# Output: output_file variable. I have named it bCell.csv. Feel free to change. It will contain all the results of bcell epitopes with their antigen, allergen and toxicity test.
# Output: output_csv variable. I have named it filtered_bCell.csv. Feel free to change.
# Due to low quality code send by the server, we need to use LXML parser in beautifullsoup inplace of html parser beacause TD is not getting closed and that gonna cause problem in scapping the score. FUCK!
import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
import csv
import os
import math
######################################################################
# FUNCTION TO SEND POST REQUEST TO ABCPred
######################################################################
def get_bcell_epitof(sequence, id, index, d2, threshold):
    print("------------------------------------------------------")
    print(f"|   {index + 1}. Attempting: B Cell Epitope for: {id}  |")
    url = "https://webs.iiitd.edu.in/cgibin/abcpred/test1_main.pl"
    payload = {
        "SEQNAME": "",
        "SEQ": f"{sequence}",
        "Threshold": threshold,
        "window": 16,
        "filter": "on"
    }
    r = requests.post(url, data=payload)
    epitope_count = 0  # Initialize the count of epitopes
    if r.status_code == 200:
        soup = BeautifulSoup(r.text, 'lxml')
        try:
            i = 0
            for elements in soup.find_all('td', attrs={'width': '50%'}):
                if i == 0:  # Skip the first irrelevant row
                    i += 1
                    continue
                d2['Protein ID'].append(id)
                d2['Epitope'].append(elements.get_text())
                i += 1
                epitope_count += 1
            score_elements = soup.find_all('td', attrs={'width': '20%'}) # For score calculation
            for i in range(3, len(score_elements), 2):  # Start at 3, step by 2
                d2['Score'].append(score_elements[i].get_text(strip=True))
        except:
            d2['Protein ID'].append(id)
            d2['Epitope'].append("ERROR:Data Extraction")
            d2['Score'].append("NA")
        print(f"|   {index + 1}. Success   : B Cell Epitope for: {id}  |")
    else:
        print(f"|    {index + 1}. Failed  : B Cell Epitope for: {id}   |")
    return epitope_count  # Return the count of epitopes found

######################################################################
# FUNCTION TO PRINT COUNT OF EPITOPS FOR SAKE OF RE-ITERATION IF NEEDED
######################################################################
def print_summary(d3,threshold):
    print("------------------------------------------------------")
    print("\n\n------------------------------------------------------")
    print(f"|    Summary of Epitope Counts (Threshold: {threshold})    |")
    print("------------------------------------------------------")
    print(f"{'Protein ID':<26}{'Epitope Count':>26}")
    print("------------------------------------------------------")
    for i in range(len(d3['Protein ID'])):
        print(f"{d3['Protein ID'][i]:<26}{d3['Epitope Count'][i]:>26}")
    print("------------------------------------------------------")

######################################################################
# FUNCTION TO RUN ALLERGIC TEST FOR A EPITOPE
######################################################################
# def get_allergenicity_result(sequence,id,index):
#     print("------------------------------------------------------")
#     print(f"|   {index + 1}. Attempting: Allergen Test for: {sequence}  |")

#     payload = {
#         "name": "Job5",
#         "seq": f">Protein\n{sequence}",  # Pass the sequence here
#         "terminus": 4,
#         "svm_th": 0.3,
#     }
#     url = "https://webs.iiitd.edu.in/raghava/algpred2/batch_action.php"
    
#     r = requests.post(url, data=payload)
    
#     if r.status_code == 200:
#         soup = BeautifulSoup(r.text, 'html.parser')
#         try:
#             # Find the result in the table (you might need to adjust the index if it's different)
#             # result = soup.find_all("td")[6].get_text()  # Assuming the 7th <td> contains the result
#             result = soup.find_all("td")[6].get_text().strip()
#         except IndexError:
#             result = "ERROR:DATA EXTRACTION"
#         print(f"|   {index + 1}. Success   : Allergen Test for: {sequence}  |")
#     else:
#         result = "ERROR:SERVER"
#         print(f"|   {index + 1}. Failed   : Allergen Test for: {sequence}  |")
#     return result
def get_allergen_test(output_file):
    def process_batch(batch, batch_index, batch_fasta_path, url, retries=3, delay=2):
        """
        Processes a single batch by submitting it to the server and handling retries.
        
        Args:
            batch (pd.DataFrame): DataFrame containing the batch sequences.
            batch_index (int): Index of the current batch.
            batch_fasta_path (str): Path to the batch FASTA file.
            url (str): Server URL for submission.
            retries (int): Number of retry attempts for failed submissions.
            delay (int): Delay between retries in seconds.
        
        Returns:
            list: List of results for the batch.
        """
        payload = {
            "name": f"Job_Batch_{batch_index + 1}",
            "seq": open(batch_fasta_path, 'r').read(),
            "terminus": 4,
            "svm_th": 0.3,
        }
        
        for attempt in range(1, retries + 1):
            try:
                response = requests.post(url, data=payload)
                
                if response.status_code == 200:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    result_cells = soup.find_all("td")
                    
                    # Parse results for each sequence
                    results = []
                    for i in range(len(batch)):
                        try:
                            result = result_cells[6 + i * 6].get_text().strip()
                            results.append(result)
                        except IndexError:
                            results.append("Result Extraction Error")
                    
                    print(f"Batch {batch_index + 1} processed successfully on attempt {attempt}.")
                    return results
                
                else:
                    print(f"Batch {batch_index + 1}: Server returned status {response.status_code}. Retrying...")
            
            except Exception as e:
                print(f"Batch {batch_index + 1}: Error on attempt {attempt} - {e}")
            
            time.sleep(delay)  # Delay before retrying
        
        print(f"Batch {batch_index + 1} failed after {retries} attempts.")
        return ["Submission Failed"] * len(batch)


    def batch_allergen_test(output_file, total_sequences=None):
        """
        Process Allergen test in batches with improved error handling.
        
        Args:
            output_file (str): Path to input CSV file.
            total_sequences (int, optional): Number of sequences to process.
        """
        df = pd.read_csv(output_file)
        
        # Validate total sequences
        max_sequences = len(df)
        if total_sequences is None or total_sequences > max_sequences:
            total_sequences = max_sequences
        
        # Fixed batch size of 200
        batch_size = 200
        num_batches = math.ceil(total_sequences / batch_size)
        print(f"\nProcessing {total_sequences} sequences in {num_batches} batches of {batch_size} sequences each.")
        
        # Create temporary directory
        os.makedirs('temp', exist_ok=True)
        url = "https://webs.iiitd.edu.in/raghava/algpred2/batch_action.php"
        
        # for batch_start in range(0, total_sequences, batch_size):
        #     batch_index = batch_start // batch_size
        #     batch = df.iloc[batch_start:min(batch_start + batch_size, total_sequences)]
            
        #     # Create batch FASTA file
        #     batch_fasta_path = f'temp/batch_{batch_index + 1}.fa'
        #     with open(batch_fasta_path, 'w') as f:
        #         for _, row in batch.iterrows():
        #             f.write(f">Protein_{row['Protein ID']}\n{row['Epitope']}\n")
            
        epitope_counter = 1
        for batch_start in range(0, total_sequences, batch_size):
            batch_index = batch_start // batch_size
            batch = df.iloc[batch_start:min(batch_start + batch_size, total_sequences)]
            
            # Create batch FASTA file
            batch_fasta_path = f'temp/batch_{batch_index + 1}.fa'
            with open(batch_fasta_path, 'w') as f:
                for _, row in batch.iterrows():
                    # Assign a unique identifier to each epitope
                    f.write(f">Epitope_{epitope_counter}\n{row['Epitope']}\n")
                    epitope_counter += 1
            # Process the batch
            results = process_batch(batch, batch_index, batch_fasta_path, url)
            
            # Update DataFrame with results
            for i, row_index in enumerate(range(batch_start, min(batch_start + batch_size, total_sequences))):
                df.at[row_index, 'Allergen Test'] = results[i]
            
            # Remove batch FASTA file
            os.remove(batch_fasta_path)
        
        # Save updated results
        df.to_csv(output_file, index=False)
        print("All batches processed. Results saved.")
        
        # Remove temporary directory
        os.rmdir('temp')


    def main():
        # output_file = input("Enter the path to the output CSV file: ").strip()
        
        # Read total available sequences
        df = pd.read_csv(output_file)
        max_sequences = len(df)
        
        # # Ask user for number of sequences to process
        # while True:
        #     try:
        #         print(f"Total sequences available: {max_sequences}")
        #         user_input = input(f"How many EPITOPES do you want to process for ALLERGEN TEST? (1-{max_sequences}, or press Enter for all): ").strip()
                
        #         if user_input == '':
        #             total_sequences = None
        #             break
                
        #         total_sequences = int(user_input)
        #         if 1 <= total_sequences <= max_sequences:
        #             break
        #         else:
        #             print(f"Please enter a number between 1 and {max_sequences}.")
        #     except ValueError:
        #         print("Please enter a valid number.")
        
        # # Process sequences
        total_sequences=max_sequences
        batch_allergen_test(output_file, total_sequences)


    if __name__ == "__main__":
        main()
######################################################################
# FUNCTION TO RUN TOXICITY TEST FOR A EPITOPE
######################################################################
def get_toxicity(output_file,num_sequences):
    ######################################################################
    # FUNCTION TO SEND POST REQUEST TO TOXINPRED SERVER
    ######################################################################
    def toxicity_result(sequence,id,index):
        url = "https://webs.iiitd.edu.in/raghava/toxinpred/multiple_test.php"
        cookies = {
            "PHPSESSID": "75s8sth03frii1hosaer2f39j1",  # Replace with the actual session ID if required
        }
        headers = {
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36",
            "Referer": "https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php",
        }
        data = {
            "seq": sequence,
            "method": "1",
            "eval": "10",
            "thval": "0.0",
            "field[]": ["4", "7", "9", "11", "13"],
        }
        print("------------------------------------------------------")
        print(f"|   {index + 1}. Attempting: Toxicity Test for: {sequence}  |")
        with requests.Session() as session:
            # Step 1: Submit the POST request
            response = session.post(url, data=data, headers=headers, cookies=cookies)
            
            # Check the response
            if response.status_code == 200:
                # Step 2: Parse the meta-refresh URL
                soup = BeautifulSoup(response.text, "html.parser")
                meta_refresh = soup.find("meta", attrs={"http-equiv": "refresh"})
                
                if meta_refresh:
                    # Extract the URL from the content attribute
                    refresh_url = meta_refresh["content"].split("url=")[-1]
                    final_url = f"https://webs.iiitd.edu.in/raghava/toxinpred/{refresh_url}"
                    
                    # Step 3: Send a GET request to the redirected URL
                    final_response = session.get(final_url)
                    
                    if final_response.status_code == 200:
                        # print("Final Response Text:\n", final_response.text)
                        soup=BeautifulSoup(final_response.text,"html.parser")
                        result=soup.find_all("td", align="center")[3].get_text() #Our final result (index need to be checked routinely)
                        print(f"|   {index + 1}. Success   : Toxicity Test for: {sequence}  |")
                    else:
                        result = "ERROR:DATA EXTRACTION"
                        print(f"|   {index + 1}. Failed   : Toxicity Test for: {sequence}  |")

                else:
                    result = "ERROR:NO REDIRECTION PAGE FOUND"
                    print(f"|   {index + 1}. Failed   : Toxicity Test for: {sequence}  |")

            else:
                result = "ERROR: INITIAL REQUEST FAILED"
                print(f"|   {index + 1}. Failed   : Toxicity Test for: {sequence}  |")
            return result
        
    ######################################################################
    # DRIVER CODE STARTS
    ######################################################################
    input_file=output_file

    df=pd.read_csv(input_file)
    for index, row in df.head(num_sequences).iterrows():
        sequence = row['Epitope']
        id=row['Protein ID']
        result = toxicity_result(sequence,id,index) 
        # Add the 'Toxicity Test' column with the fetched result
        df.at[index, 'Toxicity Test'] = result
    df.to_csv(input_file, index=False)
######################################################################
# FUNCTION TO RUN ANTIGEN TEST FOR A EPITOPE (!!!PENDING!!!!)
######################################################################
def get_antigen(output_file,num_sequences):
    ######################################################################
    # FUNCTION TO SEND POST REQUEST TO VAXIJEN SERVER
    ######################################################################
    def antigenicity_result(sequence,id, index):
        return "Pending"
    
    ######################################################################
    # DRIVER CODE STARTS
    ######################################################################
    input_file=output_file
    df=pd.read_csv(input_file)
    for index, row in df.head(num_sequences).iterrows():
        sequence = row['Epitope']
        id = row['Protein ID']
        result = antigenicity_result(sequence, id, index)
        df.at[index, 'Antigen Test'] = result
    df.to_csv(input_file, index=False)

######################################################################
# MAIN CODE STARTS
######################################################################
start_time = time.time()
csv_file = "./outputs/filtered_sequence.csv"  # INPUT CSV File. Update Accordingly.
df = pd.read_csv(csv_file)
d2 = {'Protein ID': [], 'Epitope': [], 'Score': []}  # New dictionary to store epitopes.
d3 = {'Protein ID':[], 'Epitope Count':[]} # New Dictionary to store count of epitopes.

# Ensure the number of sequences to process is within the available range
num_sequences = len(df)
# num_sequences = int(input(f"|            Total sequences available: {len(df)}        |\n>How many sequences do you want to process? "))
num_sequences = min(num_sequences, len(df))

# Process only the specified number of sequences
while True:
    # threshold = float(input(f"> Enter threshold (default is 0.51): ") or 0.51)
    threshold = float( 0.9)
    for index, row in df.head(num_sequences).iterrows():
        sequence = row['Sequence']
        id = row['Protein ID']
        epitope_count = get_bcell_epitof(sequence, id, index, d2, threshold)
        d3['Protein ID'].append(id)
        d3['Epitope Count'].append(epitope_count)
    # print_summary(d3, threshold)
    # action = input("> Do you want to repeat the process? (y/n): ").strip().lower() or 'n'
    action = 'n'

    if action != 'y':
        break
    else: #if wants to repeat the process then empty the dictionary containing epitopes and their count
        for key in d2:
            d2[key] = []
        for key in d3:
            d3[key] = []

        

# Convert the updated d2 dictionary to df2 dataframe
df2 = pd.DataFrame(data=d2) #contains epitopes
output_file = "./outputs/bCell.csv"  # OUTPUT CSV File. Update Accordingly.
df2.to_csv(output_file, index=False)
######################################################################
# NOW NEW STEPS FOR THREE TESTS OF THE EPITOPES THAT WE GOT!
######################################################################
######################################################################
# ALLERGEN TEST
######################################################################
# for index, row in df.head(num_sequences).iterrows():
#     sequence = row['Epitope']  # Get the sequence
#     id=row['Protein ID']
#     result = get_allergenicity_result(sequence,id,index) 
    
#     # Add the 'Allergen Test' column with the fetched result
#     df.at[index, 'Allergen Test'] = result
# df.to_csv(output_file, index=False)
get_allergen_test(output_file)

######################################################################
# TOXICITY TEST
######################################################################
df = pd.read_csv(output_file)
# num_sequences = int(input(f"|            Total Epitopes available: {len(df)}        |\n>How many Epitopes do you want to process for TOXICITY AND ANTIGEN test? "))
num_sequences=len(df)
num_sequences = min(num_sequences, len(df))
get_toxicity(output_file,num_sequences)

######################################################################
# ANTIGEN TEST
######################################################################
get_antigen(output_file,num_sequences)

######################################################################
# FILTERING EPITOPES WHICH SUCCEED IN ALLERGEN, ANTIGEN AND TOXICITY TEST
######################################################################
# Read the original CSV file
df = pd.read_csv(output_file)

# Apply the filter conditions
filtered_df = df[
    (df['Allergen Test'] == 'Non-Allergen') & 
    (df['Toxicity Test']=='Non-Toxin') &
    (df['Antigen Test'] == 'Pending')#change this after proper implementation of antigen-test (vaxijen)
]

# Save the filtered data to a new CSV file
output_csv = './outputs/filtered_bCell.csv'  # OUTPUT CSV File. Change Accordingly.
filtered_df.to_csv(output_csv, index=False)


end_time = time.time()
total_time = end_time - start_time
print("------------------------------------------------------")
print(f"| Success: Detecting B Cell Epitopes.                |")
print(f"| Success: Allergen Test of Epitopes.                |")
print(f"| Success: Toxicity Test of Epitopes.                |")
print(f"| PENDING: Antigen Test of Epitopes.                 |")
print(f"| Check {output_file} & {output_csv}                            |")       
print(f"| Time Taken: {total_time:.2f} sec                              |")
print("======================================================")
