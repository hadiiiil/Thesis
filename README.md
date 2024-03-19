#Create Conda library for installing the dependencies that are in the requirements file:
      conda create --name myenv python=3.10.9
      pip install -r requirements.txt
      
#For installing pygalmesh:
      conda install cgal eigen
      conda install pygalmesh
      
#First, We convert the Dicom series into NRRD image, run the dicom2nrrd.py by writing in the terminal:
      python Thesis/dicom2nrrd.py BCT/uncompressed_dicom BCT/uncompressed_dicom/Breast_Metadata.csv

#Second, We resample the images into isotropic voxel spacing 0.273, run the resampling.py by writing in the terminal:
      python Thesis/resampling.py BCT/uncompressed_nrrd 0.273

#Third, We generate the volume mesh, run the extract-model.py by writing in the terminal:
      
