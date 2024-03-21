import itk
import vtk
import argparse
import sys
import numpy as np
import pygalmesh
import pyvista
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkIOGeometry import vtkSTLWriter
import subprocess


class Arguments:
    def __init__(self):
        self.value_outside = 0
        self.facet_angle = 30.0
        self.facet_size = 3.0
        self.facet_distance = 2.0
        self.cell_radius_edge_ratio = 3.0
        self.cell_size = 2.0
        self.bc_thickness = 2.0
        self.pixel_spacing = [0.273, 0.273, 0.273]  # X, Y, Z size
        self.bc = "CT"


def main(input_filename, output_filename):

    # Initialize the arguments
    args = Arguments()

    # Print the values
    print(f"\nValue Outside: {args.value_outside}")
    print(f"Facet Angle: {args.facet_angle}")
    print(f"Facet Size: {args.facet_size}")
    print(f"Facet Distance: {args.facet_distance}")
    print(f"Cell Radius-Edge Ratio: {args.cell_radius_edge_ratio}")
    print(f"Cell Size: {args.cell_size}")
    print(f"Image Spacing: [{args.pixel_spacing[0]}, {args.pixel_spacing[1]}, {args.pixel_spacing[2]}]")
    print(f"Boundary Condition Type: {args.bc}")
    if args.bc == "CT":
        print(f"Boundary Condition Thickness: {args.bc_thickness}")
    
    if args.bc == "CT" and args.bc_thickness == 0:
        print("Check BC conditions. BC thickness cannot be 0 for CT images.")
        sys.exit(1)

    # Load your image
    input_image = itk.imread(input_filename, itk.US)

    # Get the dimensions of the image
    image_size = input_image.GetLargestPossibleRegion().GetSize()
    print(f"Original Image Dimensions: {image_size}")

    # Get the voxel size
    voxel_size = input_image.GetSpacing()
    print(f"Voxel Size: {voxel_size}")

    #Thresholding with ITK
    ThresholdType = itk.BinaryThresholdImageFilter.New(input_image)
    ThresholdType.SetOutsideValue(0)
    ThresholdType.SetInsideValue(1)
    ThresholdType.SetLowerThreshold(1)
    ThresholdType.SetUpperThreshold(255)  
    ThresholdType.Update()

    itk_image = ThresholdType.GetOutput()

    # Convert ITK image to a NumPy array for the mesh
    np_image = itk.GetArrayFromImage(itk_image).astype(np.uint8)

    # Reshape the array
    np_array_reshaped = np.transpose(np_image, (2, 1, 0))

    # create a volume mesh:
    mesh = pygalmesh.generate_from_array(
        np_array_reshaped,
        voxel_size=tuple(args.pixel_spacing),
        lloyd=False,
        odt=False,
        perturb=True,
        exude=True,
        max_edge_size_at_feature_edges=args.facet_size,
        min_facet_angle=args.facet_angle,
        max_radius_surface_delaunay_ball=0.0,  
        max_facet_distance=args.facet_distance,
        max_circumradius_edge_ratio=args.cell_radius_edge_ratio,
        max_cell_circumradius=args.cell_size,
        verbose=False
        )
    
    # reader = vtk.vtkUnstructuredGridReader()
    # reader.SetInputArray(mesh)
    # reader.Update()

    # Convert the mesh to a VTK object
    init_mesh = pyvista.wrap(mesh)

    # Convert the mesh to a PyVista mesh (UnstructuredGrid)
    mesh_pv = pyvista.wrap(mesh)

    # Extract the surface as PolyData
    surface = mesh_pv.extract_surface()
    print(output_filename)
    # Now, surface is a PyVista PolyData object. To write it as an STL file, you can use PyVista's built-in method:
    surface.save("output_filename.stl")

    # Define the first command
    command1 = ["tetgen", "-pqa15.0O0.80", "/home/hadeel/output_filename.stl"]

    # Execute the first command
    result1 = subprocess.run(command1, text=True, capture_output=True)

    # Check if the first command was successful
    if result1.returncode == 0:
        print("First command executed successfully!")
        print("Output:\n", result1.stdout)
    else:
        print("First command execution failed!")
        print("Error:\n", result1.stderr)

    # Define the second command
    command2 = ["tetgen", "-k", "/home/hadeel/output_filename.1.smesh"]

    # Execute the second command only if the first one was successful
    if result1.returncode == 0:
        result2 = subprocess.run(command2, text=True, capture_output=True)

        # Check if the second command was successful
        if result2.returncode == 0:
            print("Second command executed successfully!")
            print("Output:\n", result2.stdout)
        else:
            print("Second command execution failed!")
            print("Error:\n", result2.stderr)
    else:
        print("Skipping the second command because the first command failed.")


    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName("/home/hadeel/output_filename.2.vtk")
    reader.Update()

    # Clean the mesh
    cleaner = vtk.vtkCleanUnstructuredGrid()
    cleaner.SetInputData(reader.GetOutput())
    cleaner.Update()

    # final_points = cleaner.GetOutput().GetNumberOfPoints()
    # final_cells = cleaner.GetOutput().GetNumberOfCells()
    # print(f"Final number of points: {final_points}")
    # print(f"Final number of cells: {final_cells}")

    initial_points = init_mesh.GetNumberOfPoints()
    initial_cells = init_mesh.GetNumberOfCells()
    print(f"Initial number of points: {initial_points}")
    print(f"Initial number of cells: {initial_cells}")

    mesh = cleaner.GetOutput()

    # Iterate through cells and count non-tetrahedral cells
    non_tetra_count = 0
    for i in range(mesh.GetNumberOfCells()):
        cell = mesh.GetCell(i)
        if cell.GetCellType() != vtk.VTK_TETRA:
            non_tetra_count += 1

    print(f"Number of non-tetrahedral cells: {non_tetra_count}")

    # Iterate through cells and count tetrahedral cells
    tetra_count = 0
    for i in range(mesh.GetNumberOfCells()):
        cell = mesh.GetCell(i)
        if cell.GetCellType() == vtk.VTK_TETRA:  # Check if the cell is a tetrahedron
            tetra_count += 1

    print(f"Number of tetrahedral cells: {tetra_count}")

    # Write the original mesh
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(output_filename)
    writer.SetInputData(mesh)
    writer.Update()


    def labeling_elements(mesh, image_path, args):
        PixelType = itk.US
        image = itk.imread(image_path, PixelType)
        print("Starting labeling_elements function")

        # Create vtkCellCenters to find centers of cells in the unstructured grid
        centers = vtk.vtkCellCenters()
        centers.SetInputData(mesh)
        centers.VertexCellsOn()
        centers.Update()

        point_set = centers.GetOutput()

        # Prepare an array to store the material data
        data = vtk.vtkIntArray()
        data.SetName("materials")

        # ITK image to NumPy array for efficient access
        np_image = itk.array_from_image(image)

        np_image = np.transpose(np_image, (2, 1, 0))

        # Iterate over each cell center
        for i in range(point_set.GetNumberOfPoints()):
            pt = point_set.GetPoint(i)

            # Convert VTK point to ITK physical point
            itk_pt = itk.Point[itk.D, 3](pt)   # Double

            # Convert physical point to continuous index
            continuous_index = image.TransformPhysicalPointToContinuousIndex(itk_pt)

            # Round the continuous index to get the nearest neighbor
            index = tuple(int(round(idx)) for idx in continuous_index)

            # Safeguard against out-of-bounds indexing
            if all(0 <= idx < dim for idx, dim in zip(index, np_image.shape)):
                pixel_value = np_image[index]
                # Ensure pixel value is one of the expected labels
                assert pixel_value in [0, 1, 2, 3], f"Unexpected pixel value {pixel_value}"
                
                # Special handling: change pixel value 0 to 3
                if pixel_value == 0:
                    pixel_value = 3
                
                data.InsertNextValue(int(pixel_value))

        mesh.GetCellData().SetScalars(data)

        # Initialize a vtkIntArray for boundary conditions
        bcdata = vtk.vtkIntArray()
        bcdata.SetName("boundaryConditions")

        if args.bc == "VICTRE":
            print("Using VICTRE labeling for BC")
            # Convert cell data to point data
            cell2point = vtk.vtkCellDataToPointData()
            cell2point.SetInputData(mesh)
            cell2point.Update()

            # Apply boundary conditions based on converted point data
            for i in range(cell2point.GetOutput().GetPointData().GetArray(0).GetNumberOfTuples()):
                value = cell2point.GetOutput().GetPointData().GetArray(0).GetTuple1(i)
                bcdata.InsertNextValue(1 if value > 3 else 0)

        elif args.bc == "CT":
            print("Using CT labeling for BC")
            # Directly use points from the unstructured grid
            bcpoints = mesh.GetPoints()
            numberofpoints = bcpoints.GetNumberOfPoints()
            bbox = bcpoints.GetBounds()

            print(bbox)

            # Apply boundary conditions based on spatial location relative to the bounding box
            for i in range(numberofpoints):
                pt = bcpoints.GetPoint(i)
                # Apply BC if the point is within the defined thickness from the bounding box's minimum z-value
                if pt[2] < (bbox[4] + args.bc_thickness):
                    bcdata.InsertNextValue(1)
                else:
                    bcdata.InsertNextValue(0)

        # Attach the boundary conditions data to the point data of the unstructured grid
        mesh.GetPointData().SetScalars(bcdata)
        print("Completed labeling_elements function")

    labeling_elements(mesh, input_filename, args)

    writer.SetInputData(mesh)
    writer.Update()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a volume mesh from an input image.")
    
    # Add arguments
    parser.add_argument("input_filename", type=str, help="Path to the input image file")
    parser.add_argument("-o", "--output", type=str, help="Path to the output mesh file")

    args = parser.parse_args()
    
    main(args.input_filename, args.output)

#python thesis/gen-mesh.py BCT/uncompressed_nrrd/UncompressedBreast1/UncompressedBreast1_resampled.nrrd -o BCT/uncompressed_nrrd/UncompressedBreast1/volmesh1.vtk
