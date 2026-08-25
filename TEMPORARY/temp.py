import pickle
import numpy as np

# 9 wavelengths x 4 target compounds coefficients matrix
# Each column represents a target compound's regression weights across the 9 spectrum bands
mock_coefficients = [
    # Chl,   Car,   Anth,  Phen
    [1.020, -0.015, 0.0003, 0.005],  # Wavelength band 1
    [0.980, -0.012, 0.0002, 0.004],  # Wavelength band 2
    [0.850, -0.010, 0.0001, 0.003],  # Wavelength band 3
    [0.720, -0.008, 0.0001, 0.002],  # Wavelength band 4
    [0.510, -0.005, 0.0000, 0.001],  # Wavelength band 5
    [0.340,  0.002, 0.0001, 0.001],  # Wavelength band 6
    [0.120,  0.015, 0.0003, 0.002],  # Wavelength band 7
    [0.050,  0.022, 0.0005, 0.003],  # Wavelength band 8
    [0.010,  0.035, 0.0008, 0.006]   # Wavelength band 9
]

calibration_data = {
    "instrument": "SpectroPlant-V1",
    "calibration_date": "2026-08-21",
    "coefficients": mock_coefficients,
    "reference_white": "standard_ptfe",
    "valid_range_nm": [120, 800] # Aligned range
}

with open("$HOME/Downloads/Botanist/TEMPORARY/plant_calibration3.pkl", "wb") as f:
    pickle.dump(calibration_data, f)

print("Successfully generated 9x4 calibration matrix!")
