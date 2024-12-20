from flask import Flask, render_template, request, redirect, url_for, session, send_from_directory, abort
import subprocess
import os
import csv
import time
from collections import defaultdict

app = Flask(__name__)
app.secret_key = 'your_secret_key'

UPLOAD_FOLDER = 'uploads'
OUTPUT_FOLDER = 'outputs'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['OUTPUT_FOLDER'] = OUTPUT_FOLDER

# Ensure directories exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')  # Landing page

session=defaultdict(str)
@app.route('/step1', methods=['GET', 'POST'])
def upload():
    session['step'] = 1  # Mark the current step as Step 1--
    if request.method == 'POST':
        file = request.files['file']
        if file:
            filename = file.filename
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(file_path)
            session['file_path'] = file_path
        else:
            filename = 'default_file.fa'
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            session['file_path'] = file_path
            # Start the FinalF2C.py script (convert FA to CSV)
        try:
            subprocess.run(['python3', './scripts/finalF2C.py', file_path, os.path.join(app.config['OUTPUT_FOLDER'], 'converted.csv')], check=True)
            session['file_converted'] = True
            session['converted_filename'] = 'converted.csv'
            app.logger.info(session)
            return render_template('step1.html', file_uploaded=True, file_converted=True, converted_filename=session['converted_filename'])
        except subprocess.CalledProcessError:
            session['error_message'] = "Error while converting the file. Please try again."
            return render_template('step1.html', error_message=session['error_message'])
    return render_template('step1.html', file_uploaded=False)

@app.route('/step2', methods=['GET', 'POST'])
def step2():
    session['step'] = 2  # Update to Step 2
    print(session)
    if 'file_converted' in session:
        # Read the converted CSV to get the number of sequences available
        with open(os.path.join(app.config['OUTPUT_FOLDER'], session['converted_filename']), 'r') as csvfile:
            reader = csv.reader(csvfile)
            available_sequences = sum(1 for row in reader)-1  # Count rows in CSV file

        if request.method == 'POST':
            num_sequences = int(request.form['num_sequences'])
            session['num_sequences'] = num_sequences  # Save the number of sequences user selected

            # Start the allergenicity test in the background
            try:
                subprocess.run(['python3', './scripts/algpred_batch.py', str(num_sequences)], check=True)
                session['file_processed'] = True
                session['processed_filename'] = 'converted.csv'  # Example filename for processed results
            except subprocess.CalledProcessError:
                session['error_message'] = "Error during allergenicity test. Please try again."

            return render_template('step2.html', processing=False, file_processed=session.get('file_processed'), available_sequences=available_sequences, processed_filename=session.get('processed_filename'), error_message=session.get('error_message'))

        return render_template('step2.html', available_sequences=available_sequences, processing=False)

    return redirect(url_for('index'))

@app.route('/download/<filename>')
def download(filename):
    # Ensure the file exists in the output folder
    file_path = os.path.join(app.config['OUTPUT_FOLDER'], filename)
    
    if os.path.exists(file_path):
        return send_from_directory(app.config['OUTPUT_FOLDER'], filename, as_attachment=True)
    else:
        # If the file does not exist, show a 404 error
        abort(404, description="File not found.")

@app.route('/step3', methods=['GET', 'POST'])
def step3():
    session['step'] = 3  # Update to Step 3

    if 'file_processed' in session:
        # Read the converted CSV to get the number of sequences available
        with open(os.path.join(app.config['OUTPUT_FOLDER'], session['converted_filename']), 'r') as csvfile:
            reader = csv.reader(csvfile)
            available_sequences = sum(1 for row in reader) - 1  # Count rows in CSV file

        if request.method == 'POST':
            num_sequences = int(request.form['num_sequences'])
            session['num_sequences'] = num_sequences  # Save the number of sequences user selected

            # Start the Signal P detection in the background
            session['processing'] = True  # Set the processing state
            try:
                subprocess.run(['python3', './scripts/phobius.py', str(num_sequences)], check=True)
                session['file_processed'] = True
                session['processed_filename'] = 'converted.csv'  # Example filename for processed results
                session['processing'] = False  # Clear processing state after completion
            except subprocess.CalledProcessError:
                session['error_message'] = "Error during Signal P detection. Please try again."
                session['processing'] = False  # Clear processing state in case of an error

            return render_template('step3.html', processing=session.get('processing'), file_processed=session.get('file_processed'),
                                   available_sequences=available_sequences, processed_filename=session.get('processed_filename'),
                                   error_message=session.get('error_message'))

        return render_template('step3.html', available_sequences=available_sequences, processing=session.get('processing'))

    return redirect(url_for('index'))

@app.route('/step4', methods=['GET', 'POST'])
def step4():
    session['step'] = 4  # Update to Step 4
    if 'file_processed' in session:
        # Handle POST request when the user clicks "Start Filtration"
        if request.method == 'POST':
            # Start the filtration process
            session['processing'] = True  # Set the processing state
            try:
                subprocess.run(['python3', './scripts/filtering.py'], check=True)
                session['file_filtered'] = True
                session['file_processed'] = True
                session['filtered_filename'] = 'filtered_sequence.csv'  # Replace with actual output filename
                session['processing'] = False  # Clear processing state after completion
            except subprocess.CalledProcessError:
                session['error_message'] = "Error during filtration process. Please try again."
                session['processing'] = False  # Clear processing state in case of an error
        else:
            session['file_processed'] = False
        return render_template('step4.html', 
                               processing=session.get('processing'), 
                               file_processed=session.get('file_processed'), 
                               filtered_filename=session.get('filtered_filename'), 
                               error_message=session.get('error_message'))

    return redirect(url_for('index'))

@app.route('/step5', methods=['GET', 'POST'])
def step5():
    session['step'] = 5  # Update to Step 5

    if 'file_processed' in session:
        if request.method == 'POST':
            try:
                session['processing'] = True
                subprocess.run(['python3', './scripts/abcpred2.py'], check=True)  # Run your Python script
                session['bCell_filename'] = 'bCell.csv'  # Simulate result file names
                session['filtered_filename'] = 'filtered_bCell.csv'
                session['file_processed'] = True
            except subprocess.CalledProcessError as e:
                session['file_processed'] = False
                print(f"Error running abcpred2.py: {e}")
            finally:
                session['file_processed'] = True
                session['processing'] = False

            # return redirect(url_for('step5'))
        
        else:
            session['file_processed'] = False


        # Pass session values to the template
        return render_template(
            'step5.html',
            file_processed=session.get('file_processed', False),
            processing=session.get('processing', False),
            bCell_filename=session.get('bCell_filename', ''),
            filtered_filename=session.get('filtered_filename', '')
        )

    return redirect(url_for('index'))
@app.route('/step6', methods=['GET', 'POST'])
def step6():
    session['step'] = 6  # Update to Step 6

    if 'file_processed' in session:
        if request.method == 'POST':
            try:
                session['processing'] = True
                subprocess.run(['python3', './scripts/tca1_2.py'], check=True)  # Run your Python script
                session['tca1_filename'] = 'TC(A1).csv'  # Simulate result file names
                session['filtered_filename'] = 'filtered_TC(A1).csv'
                session['file_processed'] = True
            except subprocess.CalledProcessError as e:
                session['file_processed'] = False
                print(f"Error running tca1_2.py: {e}")
            finally:
                session['file_processed'] = True
                session['processing'] = False

            # return redirect(url_for('step5'))
        
        else:
            session['file_processed'] = False


        # Pass session values to the template
        return render_template(
            'step6.html',
            file_processed=session.get('file_processed', False),
            processing=session.get('processing', False),
            tca1_filename=session.get('tca1_filename', ''),
            filtered_filename=session.get('filtered_filename', '')
        )

    return redirect(url_for('index'))
@app.route('/step7', methods=['GET', 'POST'])
def step7():
    session['step'] = 7  # Update to Step 7

    if 'file_processed' in session:
        if request.method == 'POST':
            try:
                session['processing'] = True
                subprocess.run(['python3', './scripts/tcb58.py'], check=True)  # Run your Python script
                session['tcb58_filename'] = 'TC(B58).csv'  # Simulate result file names
                session['filtered_filename'] = 'filtered_TC(B58).csv'
                session['file_processed'] = True
            except subprocess.CalledProcessError as e:
                session['file_processed'] = False
                print(f"Error running tcb58.py: {e}")
            finally:
                session['file_processed'] = True
                session['processing'] = False

            # return redirect(url_for('step5'))
        
        else:
            session['file_processed'] = False


        # Pass session values to the template
        return render_template(
            'step7.html',
            file_processed=session.get('file_processed', False),
            processing=session.get('processing', False),
            tcb58_filename=session.get('tcb58_filename', ''),
            filtered_filename=session.get('filtered_filename', '')
        )

    return redirect(url_for('index'))
@app.route('/step8', methods=['GET', 'POST'])
def step8():
    session['step'] = 8  # Update to Step 7

    if 'file_processed' in session:
        if request.method == 'POST':
            try:
                session['processing'] = True
                subprocess.run(['python3', './scripts/th.py'], check=True)  # Run your Python script
                session['th_filename'] = 'filtered_TH.csv'  # Simulate result file names
                session['filtered_filename'] = 'filtered_TH(IFN).csv'
                session['file_processed'] = True
            except subprocess.CalledProcessError as e:
                session['file_processed'] = False
                print(f"Error running th.py: {e}")
            finally:
                session['file_processed'] = True
                session['processing'] = False

            # return redirect(url_for('step5'))
        
        else:
            session['file_processed'] = False


        # Pass session values to the template
        return render_template(
            'step8.html',
            file_processed=session.get('file_processed', False),
            processing=session.get('processing', False),
            th_filename=session.get('th_filename', ''),
            filtered_filename=session.get('filtered_filename', '')
        )

    return redirect(url_for('index'))
@app.route('/step9', methods=['GET', 'POST'])
def step9():
    session['step'] = 9  # Update to Step 9
    if 'file_processed' in session:
        # Handle POST request when the user clicks "Start Filtration"
        if request.method == 'POST':
            # Start the filtration process
            session['processing'] = True  # Set the processing state
            try:
                subprocess.run(['python3', './scripts/mergeEpitopes.py'], check=True)
                session['file_filtered'] = True
                session['file_processed'] = True
                session['filtered_filename'] = 'epitopes.csv'  # Replace with actual output filename
                session['processing'] = False  # Clear processing state after completion
            except subprocess.CalledProcessError:
                session['error_message'] = "Error during filtration process. Please try again."
                session['processing'] = False  # Clear processing state in case of an error
        else:
            session['file_processed'] = False
        return render_template('step9.html', 
                               processing=session.get('processing'), 
                               file_processed=session.get('file_processed'), 
                               filtered_filename=session.get('filtered_filename'), 
                               error_message=session.get('error_message'))

    return redirect(url_for('index'))
@app.route('/last')
def last():
    steps_files = {
            "Step1: Fasta to CSV": "converted.csv",
            "Step2: Allergenicity Test": "converted.csv",
            "Step3: Signal P Detection": "converted.csv",
            "Step4: Filtrating Sequences": "filtered_sequence.csv",
            "Step5: B Cell Epitope Detection": "filtered_bCell.csv",
            "Step6: TC(A1) Epitope Detection": "filtered_TC(A1).csv",
            "Step7: TC(B58) Epitope Detection": "filtered_TC(B58).csv",
            "Step8: TH Epitope Detection": "filtered_TH(IFN).csv",
            "Step9: Merging Eptiopes": "epitopes.csv",
        }
    final_filename = "final_merged_output.csv"

    return render_template(
            'last.html',
            steps_files=steps_files,
            final_filename=final_filename
        )
if __name__ == '__main__':
    app.run(debug=True)
