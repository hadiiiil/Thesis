import os
import shutil

# Define the paths
input_folder = 'BCT/uncompressed_dicom'
output_folder = 'BCT/uncompressed_nrrd'

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Iterate through each subfolder in the input folder
for subdir in os.listdir(input_folder):
    subfolder_path = os.path.join(input_folder, subdir)
    
    # Check if it's a directory
    if os.path.isdir(subfolder_path):
        # Find the NRRD file in the subfolder
        nrrd_file = os.path.join(subfolder_path, f'{subdir}.nrrd')
        
        # Check if the NRRD file exists
        if os.path.exists(nrrd_file):
            # Create a corresponding subfolder in the output folder
            output_subfolder = os.path.join(output_folder, subdir)
            os.makedirs(output_subfolder, exist_ok=True)
            
            # Move the NRRD file to the output subfolder
            shutil.move(nrrd_file, os.path.join(output_subfolder, f'{subdir}.nrrd'))

print("Extraction complete.")

# python thesis/move-nrrd.py