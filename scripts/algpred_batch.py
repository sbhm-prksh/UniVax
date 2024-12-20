import pandas as pd
import os
import math
import time
import requests
from bs4 import BeautifulSoup
import sys
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
    
    for batch_start in range(0, total_sequences, batch_size):
        batch_index = batch_start // batch_size
        batch = df.iloc[batch_start:min(batch_start + batch_size, total_sequences)]
        
        # Create batch FASTA file
        batch_fasta_path = f'temp/batch_{batch_index + 1}.fa'
        with open(batch_fasta_path, 'w') as f:
            for _, row in batch.iterrows():
                f.write(f">Protein_{row['Protein ID']}\n{row['Sequence']}\n")
        
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
    output_file = 'outputs/converted.csv'  # Path to your converted.csv file
    
    # # Read total available sequences
    # df = pd.read_csv(output_file)
    # max_sequences = len(df)
    
    # # Ask user for number of sequences to process
    # while True:
    #     try:
    #         print(f"Total sequences available: {max_sequences}")
    #         user_input = input(f"How many sequences do you want to process? (1-{max_sequences}, or press Enter for all): ").strip()
            
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
    
    # Process sequences
    batch_allergen_test(output_file, sequence_count)


if __name__ == "__main__":
    sequence_count = int(sys.argv[1])  # Sequence count passed from Flask
    # output_file = sys.argv[2]  # Path to the CSV file
    main()
