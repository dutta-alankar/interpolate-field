import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import vtk
from vtk.util import numpy_support as VN

# --- 1. Simulation Parameters and Placeholder Data Generation ---

# Set a fixed random seed for reproducibility
np.random.seed(12345)

# Box size
L = 100.0

# Number of particles (increased to 500,000)
n_baryons_total = 500000
n_cdms_total = 500000

# NOTE: The user described arrays like (n_baryons, nxb, nyb, nzb).
# For this script, we are generating a flat list of 3D points.
# If your data is on a grid, you will need to reshape it first.

print("Generating placeholder data...")
# Generate random 3D positions for baryons and CDMs within the box
baryon_points = np.random.rand(n_baryons_total, 3) * L
cdm_points = np.random.rand(n_cdms_total, 3) * L

# Generate random 3D velocities for the particles
vel_baryons = np.random.rand(n_baryons_total, 3)
vel_cdms = np.random.rand(n_cdms_total, 3)

print("Data generation complete.")

# --- 2. Interpolate Velocities ---

print("Interpolating CDM velocities onto baryon positions...")
vel_cdms_interp = griddata(cdm_points, vel_cdms, baryon_points, method='nearest')
print("Interpolation complete.")


# --- 3. Calculate Velocity Difference ---

velocity_difference = np.linalg.norm(vel_baryons - vel_cdms_interp, axis=1)


# --- 4. Find and Print Position of Maximum Difference ---

max_diff_index = np.argmax(velocity_difference)
max_diff_position = baryon_points[max_diff_index]
max_diff_value = velocity_difference[max_diff_index]

print(f"\nMaximum velocity difference is {max_diff_value:.4f}")
print(f"This occurs at position (x, y, z): {max_diff_position}\n")


# --- 5. Export Data to VTK File using vtk ---

print("Exporting data to VTK file using vtk...")

# Convert numpy arrays to VTK arrays
points_array = VN.numpy_to_vtk(baryon_points, deep=True)
vel_interp_array = VN.numpy_to_vtk(vel_cdms_interp, deep=True)
vel_interp_array.SetName("vel_cdms_interp")
vel_diff_array = VN.numpy_to_vtk(velocity_difference, deep=True)
vel_diff_array.SetName("velocity_difference")


# Create vtkPoints and add the points array
vtk_points = vtk.vtkPoints()
vtk_points.SetData(points_array)

# Create vtkCellArray to define the topology (vertices)
n_points = baryon_points.shape[0]
vtk_cells = vtk.vtkCellArray()
# Create a vertex for each point
for i in range(n_points):
    vtk_cells.InsertNextCell(1)
    vtk_cells.InsertCellPoint(i)


# Create an unstructured grid
ugrid = vtk.vtkUnstructuredGrid()
ugrid.SetPoints(vtk_points)
ugrid.SetCells(vtk.VTK_VERTEX, vtk_cells)

# Add the data arrays to the point data
ugrid.GetPointData().AddArray(vel_interp_array)
ugrid.GetPointData().AddArray(vel_diff_array)


# Write the unstructured grid to a .vtu file
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("velocity_data.vtu")
writer.SetInputData(ugrid)
writer.Write()

print("Data successfully exported to 'velocity_data.vtu'")


# --- 6. Create the Slice Plot ---

print("Creating slice plot...")
# Define the center of the box and the thickness of the slice
z_center = max_diff_position[2]
slice_thickness = L * 0.05 # Taking a 5% slice of the box height

# Select the baryon points that are within the central z-slice
slice_mask = np.abs(baryon_points[:, 2] - z_center) < (slice_thickness / 2.0)
slice_points = baryon_points[slice_mask]
slice_values = velocity_difference[slice_mask]

if len(slice_points) < 3:
    print("Error: Not enough points in the selected slice to create a plot.")
    print("Try increasing the number of particles or the slice_thickness.")
else:
    # Create a regular grid to interpolate the unstructured slice data onto
    grid_x, grid_y = np.mgrid[0:L:200j, 0:L:200j]

    # Interpolate the scattered data from the slice onto the regular grid
    grid_z = griddata(slice_points[:, :2], slice_values, (grid_x, grid_y), method='nearest')

    

    # --- Plotting ---

    # Check if there are any valid values to plot
    if np.all(np.isnan(grid_z)) or not np.any(grid_z > 0):
        print("Error: The interpolated grid contains no positive data to plot.")
        print("This can happen if the slice is empty or the data has no variation.")
    else:
        fig, ax = plt.subplots(figsize=(10, 8))

        # Use pcolormesh for the colormap.
        # We use LogNorm for the logarithmic scale.
        
        im = ax.pcolormesh(
            grid_x,
            grid_y,
            grid_z,
            cmap='viridis',
            
        )

        # Add a colorbar
        cbar = fig.colorbar(im)
        cbar.set_label('|vel_baryons - vel_cdms_interp|')

        # Add a red cross at the position of the maximum difference
        ax.plot(max_diff_position[0], max_diff_position[1], 'rX', markersize=20, markeredgewidth=4, alpha=0.6)

        # Set plot labels and title
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_title(f'Slice Plot of Velocity Difference at height \u2248 {z_center:.2f}')
        ax.set_aspect('equal')

        # Save the figure
        output_filename = 'velocity_difference_slice.png'
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved successfully as '{output_filename}'")
        # plt.show() # Uncomment to display the plot interactively