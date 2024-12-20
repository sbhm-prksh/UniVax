
import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
import csv
import os
import math

######################################################################
# FUNCTION TO SEND POST REQUEST TO NETMHC-II 4.0 SERVER
######################################################################
def get_th_epitope(sequence, id, index, d2):
    print("------------------------------------------------------")
    print(f"{index + 1}. Proccessing: {id}")

    url = "https://services.healthtech.dtu.dk/cgi-bin/webface2.fcgi"
    payload = {
    "configfile": "/var/www/services/services/NetMHCIIpan-4.0/webface.cf",
    "inp": "0",
    "SEQPASTE": sequence,
    "SEQSUB": "(binary)",
    "length": "15",
    "PEPPASTE": "",
    "PEPSUB": "(binary)",
    "termAcon": "on",
    "master": "1",
    "slaveB": [
        "DRB1_0101", "DRB1_0102", "DRB1_0103", "DRB1_0301", "DRB1_0305",
        "DRB1_0401", "DRB1_0402", "DRB1_0403", "DRB1_0404", "DRB1_0405",
        "DRB1_0408", "DRB1_0701", "DRB1_0801", "DRB1_0803", "DRB1_0901",
        "DRB1_1001", "DRB1_1101", "DRB1_1104", "DRB1_1201", "DRB1_1301"
    ],
    "allele": "DRB1_0101,DRB1_0102,DRB1_0103,DRB1_0301,DRB1_0305,DRB1_0401,DRB1_0402,DRB1_0403,DRB1_0404,DRB1_0405,DRB1_0408,DRB1_0701,DRB1_0801,DRB1_0803,DRB1_0901,DRB1_1001,DRB1_1101,DRB1_1104,DRB1_1201,DRB1_1301",
    "MHCSEQPASTEa": "",
    "MHCSEQSUBa": "(binary)",
    "MHCSEQPASTEb": "",
    "MHCSEQSUBb": "(binary)",
    "thrs": "1",
    "thrw": "5",
    "thrf": "10",
    "sort": "on"
    }
    
    # Submit the job
    response = requests.post(url, data=payload, allow_redirects=False)
    if response.status_code != 302:
        print("Error: Job submission failed.")
        d2['Protein ID'].append(id)
        d2['Epitope'].append("ERROR: Submission Failed")
        d2['Allele'].append("NA")
        d2['Score'].append("NA")
        d2['Result'].append("NA")
        return

    # Extract the redirect URL (Job Status Page)
    base_url = "https://services.healthtech.dtu.dk"
    job_status_url = base_url + response.headers['Location']
    print(f"> Job submitted. Checking status...")

    # Poll the status page until the job is complete
    while True:
        time.sleep(2)  # Wait before polling
        status_response = requests.get(job_status_url)
        soup = BeautifulSoup(status_response.text, 'lxml')
        pre_tag = soup.find('pre')  # Look for the <pre> tag containing results
        # Check for completion markers
        if pre_tag and "Number of weak binders" in pre_tag.text:
            print(f"> Job finished. Extracting results...")
            # Extract epitope data from the <pre> content
            lines = pre_tag.text.splitlines()
            for line in lines:
                if "<=SB" in line:  # Identify lines with epitopes
                    parts = line.split()
                    epitope = parts[2]  # Assuming 3th column is Epitope
                    allele=parts[1] # Assuming 2th column is Epitope
                    score = float(parts[-4])  # Assuming forth last column is COMB score
                    d2['Protein ID'].append(id)
                    d2['Epitope'].append(epitope)
                    d2['Allele'].append(allele)
                    d2['Score'].append(score)
                    d2['Result'].append("Strong Binder")
            print(f"{index + 1}. Successfull: {id}")
            return
        
        elif "Jobid not provided" in soup.text:
            print(f"Job {id} encountered an error.")
            d2['Protein ID'].append(id)
            d2['Epitope'].append("ERROR: Job Failed")
            d2['Score'].append("NA")
            return
        else:
            print(f"> Still processing...")

######################################################################
# FUNCTION TO RUN ALLERGIC TEST FOR A EPITOPE
######################################################################
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
        total_sequences=max_sequences
        # Process sequences
        batch_allergen_test(output_file, total_sequences)


    if __name__ == "__main__":
        main()
                  
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
    # DRIVER CODE STARTS FOR TOXICITY RESULT
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
    # DRIVER CODE STARTS FOR ANTIGEN TEST
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
# FUNCTION TO RUN IFN GAMMA TEST FOR A EPITOPE
######################################################################
def get_ifn(output_file):
    def batch_ifn_gamma_test(output_file, total_sequences=None):
        """
        Process IFN Gamma test in batches with user-specified sequence count.
        
        Args:
            output_file (str): Path to input CSV file
            total_sequences (int, optional): Number of sequences to process
        """
        # Read input file
        df = pd.read_csv(output_file)
        # total_sequences=int(total_sequences)
        # Validate total sequences
        max_sequences = len(df)
        if total_sequences is None or (total_sequences) > max_sequences:
            total_sequences = max_sequences
        
        # Fixed batch size of 10
        batch_size = 10
        
        # Calculate number of batches
        num_batches = math.ceil(total_sequences / batch_size)
        print(f"\nProcessing {total_sequences} sequences in {num_batches} batches of 10 sequences.")
        
        # Create temporary directory
        os.makedirs('temp', exist_ok=True)
        
        # Process in batches
        for batch_start in range(0, total_sequences, batch_size):
            batch = df.iloc[batch_start:min(batch_start+batch_size, total_sequences)]
            
            # Create batch FASTA file
            batch_fasta_path = f'temp/batch_{batch_start//batch_size + 1}.fa'
            with open(batch_fasta_path, 'w') as f:
                for _, row in batch.iterrows():
                    f.write(f">Protein_{row['Protein ID']}\n{row['Epitope']}\n")
            
            # Batch submission
            payload = {
                "job_name": f"Batch_{batch_start//batch_size + 1}",
                "sequence": open(batch_fasta_path, 'r').read(),
                "method": "hybrid",
                "model": "main",
            }
            
            url = "https://webs.iiitd.edu.in/raghava/ifnepitope/run_submit-old.php"
            
            try:
                # Submit job
                response = requests.post(url, data=payload, allow_redirects=False)
                soup = BeautifulSoup(response.text, "lxml")
                meta_refresh = soup.find("meta", attrs={"http-equiv": "refresh"})
                
                if not meta_refresh:
                    print(f"Batch {batch_start//batch_size + 1}: Submission Failed")
                    continue
                
                content = meta_refresh.get("content", "")
                refresh_url = f"https://webs.iiitd.edu.in{content.split('url=')[-1]}"
                result_url = refresh_url.replace("refresh_page", "result_page")
                
                # Poll for results
                while True:
                    time.sleep(5)
                    response = requests.get(result_url)
                    soup = BeautifulSoup(response.text, "lxml")
                    
                    if "Epitope Name" in soup.text:
                        # Extract results for each sequence in batch
                        result_cells = soup.find_all("td")
                        for i, row_index in enumerate(range(batch_start, min(batch_start+batch_size, total_sequences))):
                            try:
                                # Every 6th cell contains IFN Gamma result
                                result = result_cells[4 + i*6].get_text().strip()
                                df.at[row_index, 'IFM_Gamma'] = result
                            except IndexError:
                                df.at[row_index, 'IFM_Gamma'] = "Result Extraction Error"
                        
                        print(f"Batch {batch_start//batch_size + 1} processed successfully")
                        break
                    else:
                        print(f"Batch {batch_start//batch_size + 1}: Still processing...")
            
            except Exception as e:
                print(f"Error processing batch {batch_start//batch_size + 1}: {e}")
        
        # Save updated results
        df.to_csv(output_file, index=False)
        
        # Clean up temporary files
        for file in os.listdir('temp'):
            os.remove(os.path.join('temp', file))
        os.rmdir('temp')

    def main():
        # Input file
        # Read total available sequences
        df = pd.read_csv(output_file)
        max_sequences = len(df)
        
    
        total_sequences=max_sequences
        batch_ifn_gamma_test(output_file, total_sequences)

    if __name__ == "__main__":
        main()

######################################################################
# DRIVER CODE STARTS
######################################################################
csv_file = "./outputs/filtered_sequence.csv"  # INPUT CSV File. Update Accordingly.
df = pd.read_csv(csv_file)
d2 = {'Protein ID': [],'Allele':[], 'Epitope': [], 'Score': [], 'Result':[]}  # New dictionary to store epitopes.
# print("|        1 Sequence Processing takes ~2min        |")
# num_sequences = int(input(f"|            Total sequences available: {len(df)}        |\n> How many sequences do you want to process? "))
num_sequences=1 #shubham
# Ensure the number of sequences to process is within the available range
num_sequences = min(num_sequences, len(df))
start_time = time.time()

# Process only the specified number of sequences
for index, row in df.head(num_sequences).iterrows():
    sequence = row['Sequence']
    id = row['Protein ID']
    get_th_epitope(sequence, id, index, d2)


# Convert the updated d2 dictionary to df2 dataframe
df2 = pd.DataFrame(data=d2) #contains epitopes
# df3 = pd.DataFrame(data=d3) #contains epitopes count

output_file = "./outputs/TH.csv"  # OUTPUT CSV File. Update Accordingly.
# output_file2="summary.csv" # OUTPUT2 CSV File. Update Accordingly.
df2.to_csv(output_file, index=False)
# df3.to_csv(output_file2, index=False)
end_time = time.time()
total_time = end_time - start_time
print("------------------------------------------------------")
print(f"| Success: TH Epitope Prediction.                    |")
print(f"| Check {output_file}                                    |")       
print(f"| Time Taken: {total_time:.2f} sec                          |")
print("======================================================")
get_allergen_test(output_file)
######################################################################
# TOXICITY TEST
######################################################################
df = pd.read_csv(output_file)
# num_sequences = int(input(f"|            Total Epitopes available: {len(df)}        |\n>How many Epitopes do you want to process for allergen and toxicity test? "))
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
output_csv = './outputs/filtered_TH.csv'  # OUTPUT CSV File. Change Accordingly.
filtered_df.to_csv(output_csv, index=False)
######################################################################
# IFN GAMMA TEST
######################################################################
get_ifn(output_csv)
######################################################################
# FILTERING EPITOPES WHICH SUCCEED IN ALLERGEN, ANTIGEN AND TOXICITY TEST
######################################################################
df = pd.read_csv(output_csv)
# Apply the filter conditions
filtered_df = df[
    (df['Allergen Test'] == 'Non-Allergen') & 
    (df['Toxicity Test']=='Non-Toxin') &
    (df['IFM_Gamma'] == 'POSITIVE') &
    (df['Antigen Test'] == 'Pending')#change this after proper implementation of antigen-test (vaxijen)
]
# Save the filtered data to a new CSV file
output_csv = './outputs/filtered_TH(IFN).csv'  # OUTPUT CSV File. Change Accordingly.
filtered_df.to_csv(output_csv, index=False)
end_time = time.time()
total_time = end_time - start_time
print("------------------------------------------------------")
print(f"| Success: Detecting TH Cell Epitopes.               |")
print(f"| Success: Allergen Test of Epitopes.                |")
print(f"| Success: Toxicity Test of Epitopes.                |")
print(f"| PENDING: Antigen Test of Epitopes.                 |")
print(f"| Check {output_file} & {output_csv}             |")       
print(f"| Time Taken: {total_time:.2f} sec                              |")
print("======================================================")
