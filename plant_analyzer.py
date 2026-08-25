import numpy as np
import pandas as pd
import argparse
from scipy.signal import savgol_filter
import joblib  #to load pre-trained model matrix

def plant_spectra_analyzer(raw_csv_path, trained_model_path):
    """
    Inputs: 
      raw_csv_path: Path to a CSV with columns ['wavelength', 'reflectance']
      trained_model_path: Path to a saved, pre-trained sklearn PLSRegression model
    Output: 
      A dictionary of predicted phytochemical percentages
    """

    df = pd.read_csv(raw_csv_path)
    df = df[(df['wavelength'] >= 120) & (df['wavelength'] <= 800)]
    raw_spectra = df['reflectance'].values.reshape(1, -1)
    
    snv_spectra = (raw_spectra - np.mean(raw_spectra)) / np.std(raw_spectra) #preprocessing: Standard Normal Variate (SNV) to remove light scattering
    window_len = 2*int(0.5+0.1*len(df))+1
    processed_spectra = savgol_filter(snv_spectra, window_length=window_len, polyorder=min(window_len-1,3), deriv=1, axis=1) # 1st derivative preprocessing (window_length must be odd, polyorder must be less than window_length)

    plsr_model = joblib.load(trained_model_path)

    predictions = np.dot(processed_spectra, plsr_model["coefficients"])[0]
    #predictions = plsr_model.predict(processed_spectra)[0] #predicted compound percentages
    target_compounds = ['Total_Chlorophyll', 'Total_Carotenoids', 'Anthocyanins', 'Total_Phenolics'] #target compounds match the order the model was originally trained on
    
    results = {compound: round(float(percentage), 3) for compound, percentage in zip(target_compounds, predictions)}
    return results

parser = argparse.ArgumentParser()
parser.add_argument("-file_paths", nargs=2)
parser.add_argument("-download_type", nargs=1, choices=["csv", "txt", "json"])

if __name__ == "__main__": #so that plant_analyzer.py doesn't execute by itself but only when called by a subprocess in main.py
    args = parser.parse_args()
    try:
        if args.download_type:
            print("[this should download the output of plant_analyzer]") #temporary print statement; we will actually download it here for the user as a csv/json/txt
        else:
            print(plant_spectra_analyzer(args.file_paths[0], args.file_paths[1]))
    except Exception as e:
        print("plant_analyzer had an error running: ",type(e).__name__)
        print("Details:")
        print(e)

#argparse.add_argument("raw_csv_path", type=str, help="Path to the raw reflectance CSV file")

#along with this csv, pair the calibration pkl file. Example execution:

# composition = airtight_spectra_analyzer("leaf_scan_001.csv", "trained_leaf_matrix.pkl")
# print(composition)