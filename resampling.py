import itk
import sys
import re
import os

def natural_sort_key(string):
    """Key function for natural sorting"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string)]


def resample_to_isotropic(main_folder, isotropic_spacing):

    # Get sorted list of subfolders
    subfolders = sorted([subfolder for subfolder in os.listdir(main_folder) if os.path.isdir(os.path.join(main_folder, subfolder))], key=natural_sort_key)

    # Loop through each subfolder and convert NRRD images
    for subfolder in subfolders:
        subfolder_path = os.path.join(main_folder, subfolder)

        # Load the input image
        PixelType = itk.US
        Dimension = 3
        ImageType = itk.Image[PixelType, Dimension]
        input_image = itk.imread(os.path.join(subfolder_path, f'{subfolder}.nrrd'), PixelType)
        
        # Define the desired isotropic spacing
        desired_spacing = (isotropic_spacing, isotropic_spacing, isotropic_spacing)

        # Define the new origin and direction for the image
        origin = [0.0, 0.0, 0.0]

        # Create an identity direction matrix
        identity_matrix = itk.Matrix[itk.D, 3, 3]()
        identity_matrix.SetIdentity()


        # Get the original image's spacing and size
        original_spacing = input_image.GetSpacing()
        original_size = input_image.GetLargestPossibleRegion().GetSize()

        # Calculate the new size of the resampled image
        new_size = [int(round(original_size[i] * (original_spacing[i] / desired_spacing[i]))) for i in range(3)]

        # Set up the resampling filter
        ResampleImageFilterType = itk.ResampleImageFilter[ImageType, ImageType]
        resample_filter = ResampleImageFilterType.New()
        resample_filter.SetInput(input_image)
        resample_filter.SetOutputSpacing(desired_spacing)
        resample_filter.SetSize(new_size)
        resample_filter.SetOutputOrigin(origin)
        resample_filter.SetOutputDirection(identity_matrix)

        # Use the nearest neighbor interpolator
        interpolator = itk.NearestNeighborInterpolateImageFunction[ImageType, itk.D].New()
        resample_filter.SetInterpolator(interpolator)

        # Set the boundary condition to 0
        resample_filter.SetDefaultPixelValue(0)
        # Update the filter to perform the resampling
        resample_filter.Update()

        # Save the resampled image
        resampled_file = os.path.join(subfolder_path, f'{subfolder}_resampled.nrrd')
        print(f'Writing the resampled nrrd {resampled_file}')

        itk.imwrite(resample_filter.GetOutput(), resampled_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python resample_isotropic.py <main_folder_path> <isotropic_spacing_mm>")
        sys.exit(1)

    main_folder = sys.argv[1]
    isotropic_spacing = float(sys.argv[2])

    resample_to_isotropic(main_folder, isotropic_spacing)
    
# python thesis/resampling.py BCT/uncompressed_nrrd 0.1908