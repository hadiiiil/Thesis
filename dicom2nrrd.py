import os
import re
import sys
import SimpleITK as sitk
import pandas as pd

def natural_sort_key(string):
    """Key function for natural sorting"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string)]

def convert_dicom_series_to_nrrd(dicom_folder, voxel_spacing):

    print(f'Reading the dicom files of {dicom_folder}')
    # Read the DICOM series
    reader = sitk.ImageSeriesReader()
    dicom_series = reader.GetGDCMSeriesFileNames(dicom_folder)
    reader.SetFileNames(dicom_series)
    image = reader.Execute()
    
    # Set the voxel spacing
    image.SetSpacing(voxel_spacing)
    
    # Save the NRRD file
    # nrrd_file = os.path.splitext(dicom_folder)[0] + '.nrrd'
    nrrd_file = os.path.join(dicom_folder, os.path.splitext(os.path.basename(dicom_folder))[0] + '.nrrd')
    print(f'Writing the nrrd {nrrd_file}')
    sitk.WriteImage(image, nrrd_file)

def main(main_folder, csv_file):
    # Read the CSV file containing voxel spacing information
    voxel_spacing_df = pd.read_csv(csv_file)

    # Get sorted list of subfolders
    subfolders = sorted([subfolder for subfolder in os.listdir(main_folder) if os.path.isdir(os.path.join(main_folder, subfolder))], key=natural_sort_key)

    # Loop through each subfolder and convert DICOM series to NRRD
    for subfolder in subfolders:
        subfolder_path = os.path.join(main_folder, subfolder)
        if subfolder.startswith('UncompressedBreast'):
            # Extract study number from the folder name
            study_number = int(subfolder.replace('UncompressedBreast', ''))
            
            # Find the corresponding row in the CSV file
            voxel_spacing_row = voxel_spacing_df[voxel_spacing_df['Breast n.'] == study_number]
            
            if not voxel_spacing_row.empty:
                # Extract voxel spacing information
                coronal_pixel_pitch = voxel_spacing_row['Coronal pixel pitch (mm)'].values[0]
                slice_thickness = voxel_spacing_row['slicethickness (mm)'].values[0]
                
                # Create voxel spacing tuple
                voxel_spacing = (coronal_pixel_pitch, coronal_pixel_pitch, slice_thickness)
                
                # Convert DICOM series to NRRD using the extracted voxel spacing
                convert_dicom_series_to_nrrd(subfolder_path, voxel_spacing)
            else:
                print(f"No voxel spacing information found for study {study_number}. Skipping.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <main_folder_path> <csv_file_path>")
        sys.exit(1)
    
    main_folder = sys.argv[1]
    csv_file = sys.argv[2]
    main(main_folder, csv_file)

# python dicom2nrrd.py BCT/uncompressed_dicom BCT/uncompressed_dicom/Breast_Metadata.csv