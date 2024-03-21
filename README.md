# Create conda environment for installing the dependencies that are in the requirements file:
      conda create --name myenv python=3.10.9
      pip install -r requirements.txt
      
# For installing pygalmesh:
      conda install cgal eigen
      conda install pygalmesh
# For installing Tetgen python library interface from source at GitHub: https://github.com/pyvista/tetgen
      git clone https://github.com/pyvista/tetgen
      cd tetgen
      pip install .
      
# First, We convert the Dicom series into NRRD image, run the dicom2nrrd.py by writing in the terminal:
      python Thesis/dicom2nrrd.py BCT/uncompressed_dicom BCT/uncompressed_dicom/Breast_Metadata.csv
      
# Then move the nrrd images to a new folder uncompressed_nrrd by running the following command:
      python Thesis/move-nrrd.py

# Second, We resample the images into isotropic voxel spacing 0.273, run the resampling.py by writing in the terminal:
      python Thesis/resampling.py BCT/uncompressed_nrrd 0.273

# Third, We generate the volume mesh, run the gen-mesh.py by writing in the terminal:
      python Thesis/gen-mesh.py BCT/uncompressed_nrrd/UncompressedBreast1/UncompressedBreast1_resampled.nrrd -o BCT/uncompressed_nrrd/UncompressedBreast1/volmesh1.vtk
