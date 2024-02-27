import vtk
from vtk import *

# Initialize a global vtkCellArray to store contour segments
global_contour_segments = vtkCellArray()


#loading the given Hurricane Simulation dataset
reader = vtkXMLImageDataReader()
reader.SetFileName('Isabel_2D.vti')         
reader.Update() 

#creating data object of the vtkImageData from the reader's output
data = reader.GetOutput()
#print(data, "\n")

# Extracting the vtkDataArray named 'Pressure' from the point data of the VTK dataset
dataArr= data.GetPointData().GetArray('Pressure')
#print(dataArr.GetTuple().GetCell(0))

# Finding Pressure Range
#x, y = data.GetScalarRange()

#Getting isoValue as input from user
# Specifying the possible range of isovalues for the given data set
isoValueMin = -1438
isoValueMax = 630

# Getting isoValue as input from user within the specified range
print(f"Enter the isovalue that you want to extract in the range ({isoValueMin} to {isoValueMax}): ")
isoVal = float(input())

# Check if the entered isoValue is within the specified range
while isoVal < isoValueMin or isoVal > isoValueMax:
    print(f"Invalid isovalue Entered. Please try to enter a value in the range ({isoValueMin} to {isoValueMax}): ")
    isoVal = float(input())

# Getting the number of cells in the uniform grid of the dataset
n= data.GetNumberOfCells()

#Converting vtk array to numpy array
import numpy as np
# Importing vtk_to_numpy for conversion
from vtk.util.numpy_support import vtk_to_numpy
# Getting the dimensions of the dataset (width, height, depth(not needed in this context,so used _))
width, height, _ = data .GetDimensions()
#Getting the VTK array (scalar values) associated with the points in the dataset
vtk_array = data.GetPointData().GetScalars()
# Determining the number of components in the VTK array
components = vtk_array.GetNumberOfComponents()
# Converting the VTK array to a NumPy array and reshaping based on dataset dimensions
arr = vtk_to_numpy(vtk_array).reshape(height, width, components)
# Printing an example element from the resulting NumPy array
#print(arr[0][1])


#Interpolation Function in order to find points on the active edges of the cell
def linear_interpolation(val1, val2, C, vertex1, vertex2):
   # Computing the slope of the interpolation line
   slope = (val1 - C) / (val1 - val2)
   # Interpolating x-coordinate based on the slope
   x = vertex1[0] + slope * (vertex2[0] - vertex1[0])
   # Interpolating y-coordinate based on the slope
   y = vertex1[1] + slope * (vertex2[1] - vertex1[1])
   # Returning the interpolated coordinates as a tuple (x, y)
   return (x, y)


#Iterating through the dataset
isoContour=[]
for i in range (arr.shape[0]-1):
  for j in range (arr.shape[1]-1):

    #Quering vertices in the counterclockwise direction
    v1=arr[i,j]
    v2=arr[i+1,j]
    v3=arr[i+1,j+1]
    v4=arr[i,j+1]
    #print(v1,v2,v3,v4)

    #Defining possible cases from 0 to 15
    case=0
    if(v1>isoVal):
      case+=1
    if(v2>isoVal):
      case+=2
    if(v3>isoVal):
      case+=4
    if(v4>isoVal):
      case+=8
      
    #Defining coordinates of cell
    x1=i
    y1=j
    x2=i+1
    y2=j
    x3=i+1
    y3=j+1
    x4=i
    y4=j+1

    p1=x1,y1
    p2=x2,y2
    p3=x3,y3
    p4=x4,y4

    #No active edges
    if(case==0 or case== 15):
       continue
      
    #Interpolating in cases of active edge found
    if(case == 1):
       isoContour.append(linear_interpolation(v1, v2, isoVal, p1, p2))
       isoContour.append(linear_interpolation(v1, v4, isoVal, p1, p4))
    elif(case==2):
       isoContour.append(linear_interpolation(v2, v3, isoVal, p2, p3))
       isoContour.append(linear_interpolation(v1, v2, isoVal, p1, p2))
    elif(case==3):
       isoContour.append(linear_interpolation(v2, v3, isoVal, p2, p3))
       isoContour.append(linear_interpolation(v1, v4, isoVal, p1, p4))
    elif(case==4):
       isoContour.append(linear_interpolation(v2, v3, isoVal, p2, p3))
       isoContour.append(linear_interpolation(v3, v4, isoVal, p3, p4))

    #Case when two line segments will be found
    elif(case==5):
       isoContour.append(linear_interpolation(v1, v4, isoVal, p1, p4))
       isoContour.append(linear_interpolation(v1, v2, isoVal, p1, p2))
       isoContour.append(linear_interpolation(v3, v2, isoVal, p3, p2))
       isoContour.append(linear_interpolation(v3, v4, isoVal, p3, p4))
    elif(case==6):
       isoContour.append(linear_interpolation(v3, v4, isoVal, p3, p4))
       isoContour.append(linear_interpolation(v1, v2, isoVal, p1, p2))
    elif(case==7):
       isoContour.append(linear_interpolation(v3, v4, isoVal, p3, p4))
       isoContour.append(linear_interpolation(v1, v4, isoVal, p1, p4))
    elif(case==8):
       isoContour.append(linear_interpolation(v3, v4, isoVal, p3, p4))
       isoContour.append(linear_interpolation(v1, v4, isoVal, p1, p4))
    elif(case==9):
       isoContour.append(linear_interpolation(v3, v4, isoVal, p3, p4))
       isoContour.append(linear_interpolation(v1, v2, isoVal, p1, p2))

    #Case when two line segments will be found
    elif(case==10):
       isoContour.append(linear_interpolation(v2, v3, isoVal, p2, p3))
       isoContour.append(linear_interpolation(v1, v2, isoVal, p1, p2))
       isoContour.append(linear_interpolation(v1, v4, isoVal, p1, p4))
       isoContour.append(linear_interpolation(v3, v4, isoVal, p3, p4))
    elif(case==11):
       isoContour.append(linear_interpolation(v2, v3, isoVal, p2, p3))
       isoContour.append(linear_interpolation(v3, v4, isoVal, p3, p4))
    elif(case==12):
       isoContour.append(linear_interpolation(v2, v3, isoVal, p2, p3))
       isoContour.append(linear_interpolation(v1, v4, isoVal, p1, p4))
    elif(case==13):
       isoContour.append(linear_interpolation(v2, v3, isoVal, p2, p3))
       isoContour.append(linear_interpolation(v1, v2, isoVal, p1, p2))
    elif(case==14):
       isoContour.append(linear_interpolation(v1, v4, isoVal, p1, p4))
       isoContour.append(linear_interpolation(v1, v2, isoVal, p1, p2))

#Converting list into numpy array
isoContour= np.array(isoContour)
#for i in isoContour:
  #print(i)
#print(isoContour)


#Creating PolyData object
polyData = vtkPolyData()

#Creating Points Array and setting points to the array
points = vtkPoints()
# Iterating through each vertex in isoContour array
for vertex in isoContour:
        # Extracting x and y coordinates from each vertex
        x, y= vertex
        # Inserting the point into the vtkPoints array with coordinates (y, x, 25)
        points.InsertNextPoint(y,x, 25)
# Setting the points in the PolyData object
polyData.SetPoints(points)


# Creating an empty vtkCellArray to store line segments
lines = vtkCellArray()
# Getting the total number of points in the isoContour array
n = len(isoContour)
# Iterating through isoContour array in steps of 2 (since each line segment has two points)
for i in range(0, n-1, 2):
    # Creating a new vtkLine object for each line segment
    line = vtk.vtkLine()
    # Setting the first point of the line segment using isoContour array indices
    line.GetPointIds().SetId(0, i)
    # Setting the second point of the line segment using isoContour array indices
    line.GetPointIds().SetId(1, i + 1)
    # Inserting the newly created line segment into the vtkCellArray
    lines.InsertNextCell(line)
# Setting the lines in the PolyData object
polyData.SetLines(lines)


#Writing the VTK XML PolyData file format
writer = vtkXMLPolyDataWriter()
#.vtp file named isoContour will be generated
writer.SetFileName("isoContour.vtp")
writer.SetInputData(polyData)
writer.Write()

print("The vtp file is generated!")

			
		
		
		
		


