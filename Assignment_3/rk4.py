import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

def load_vti_file(filename):
    """Load a VTI file and return the VTK object."""
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def rk4_integration(vtk_field, seed_point, step_size, max_steps):
    """Perform RK4 integration from a seed point."""
    def get_vector_at_point(point):
        """Use VTKProbeFilter to interpolate the vector at a given point."""
        probe = vtk.vtkProbeFilter()
        probe.SetSourceData(vtk_field)
        pointSource = vtk.vtkPointSource()
        pointSource.SetNumberOfPoints(1)
        pointSource.SetCenter(point)
        probe.SetInputConnection(pointSource.GetOutputPort())
        probe.Update()
        data = probe.GetOutput().GetPointData().GetVectors()
        if data:
            return probe.GetOutput().GetPointData().GetVectors().GetTuple(0)
        else:
            return None
    
    def integrate(point, step_size):
        """Calculate the next point using RK4 method."""
        k1 = np.array(get_vector_at_point(point))
        if k1 is None: return None
        k2 = np.array(get_vector_at_point(point + 0.5 * step_size * k1))
        if k2 is None: return None
        k3 = np.array(get_vector_at_point(point + 0.5 * step_size * k2))
        if k3 is None: return None
        k4 = np.array(get_vector_at_point(point + step_size * k3))
        if k4 is None: return None
        return point + (step_size / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    
    def trace_streamline(direction):
        """Trace the streamline in a given direction."""
        points = []
        point = np.array(seed_point)
        for _ in range(max_steps):
            bounds = vtk_field.GetBounds()
            if not (bounds[0] <= point[0] <= bounds[1] and bounds[2] <= point[1] <= bounds[3] and bounds[4] <= point[2] <= bounds[5]):
                break
            next_point = integrate(point, direction * step_size)
            if next_point is None: break
            points.append(next_point)
            point = next_point
        return points
    
    forward_points = trace_streamline(1)
    backward_points = trace_streamline(-1)
    
    if backward_points:
        backward_points.reverse()
    streamline_points = np.vstack((backward_points, [seed_point], forward_points)) if backward_points and forward_points else []
    return streamline_points

def create_vtk_poly_data(points):
    """Create a VTKPolyData object from a list of points."""
    points_vtk = vtk.vtkPoints()
    for point in points:
        points_vtk.InsertNextPoint(point.tolist())
    
    poly_data = vtk.vtkPolyData()
    poly_data.SetPoints(points_vtk)
    
    lines = vtk.vtkCellArray()
    for i in range(len(points) - 1):
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, i)
        line.GetPointIds().SetId(1, i + 1)
        lines.InsertNextCell(line)
    
    poly_data.SetLines(lines)
    return poly_data

def save_poly_data(poly_data, filename):
    """Save a VTKPolyData object to a file."""
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(poly_data)
    writer.Write()

def main():
    # Load the dataset
    vtk_field = load_vti_file("tornado3d_vector.vti")
    
    # Get seed point from user input
    seed_x = float(input("Enter seed X coordinate: "))
    seed_y = float(input("Enter seed Y coordinate: "))
    seed_z = float(input("Enter seed Z coordinate: "))
    seed_point = [seed_x, seed_y, seed_z]
    
    # Perform RK4 integration to trace the streamline
    points = rk4_integration(vtk_field, seed_point, 0.05, 1000)
    if len(points) == 0:
        print("No streamline was generated. The seed point may be out of the vector field bounds or in a region with no vector data.")
        return
    
    # Create and save the VTKPolyData file
    poly_data = create_vtk_poly_data(points)
    save_poly_data(poly_data, "streamline.vtp")
    print("Streamline saved to streamline.vtp.")

if __name__ == "__main__":
    main()
