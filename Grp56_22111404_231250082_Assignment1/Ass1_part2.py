import vtk
from vtk import *

#loading the given dataset
reader = vtkXMLImageDataReader()
reader.SetFileName('Isabel_3D.vti')         
reader.Update() 

#creating data object from the reader's output
data = reader.GetOutputPort()
#print(data, "\n")

# Creating a vtkSmartVolumeMapper instance for volume rendering
volumeMapper= vtkSmartVolumeMapper()
# Setting the blending mode to composite
volumeMapper.SetBlendModeToComposite()
# Setting the input connection for the volume mapper with the data object
volumeMapper.SetInputConnection(data)

# Creating a vtkColorTransferFunction instance for color transfer in volume rendering
ctfun = vtkColorTransferFunction()
# Adding an RGB point for color mapping (-4931.54, 0.0, 1.0, 1.0)
ctfun.AddRGBPoint(-4931.54, 0.0, 1.0, 1.0)
# Adding an RGB point for color mapping (-2508.95, 0.0, 0.0, 1.0)
ctfun.AddRGBPoint(-2508.95, 0.0, 0.0, 1.0)
# Adding an RGB point for color mapping (-1873.9, 0.0, 0.0, 0.5)
ctfun.AddRGBPoint(-1873.9, 0.0, 0.0, 0.5)
# Adding an RGB point for color mapping (-1027.16, 1.0, 0.0, 0.0)
ctfun.AddRGBPoint(-1027.16, 1.0, 0.0, 0.0)
# Adding an RGB point for color mapping (-298.031, 1.0, 0.4, 0.0)
ctfun.AddRGBPoint(-298.031, 1.0, 0.4, 0.0)
# Adding an RGB point for color mapping (2594.97, 1.0, 1.0, 0.0)
ctfun.AddRGBPoint(2594.97, 1.0, 1.0, 0.0)

#Creating instances of vtkPiecewiseFunction() and setting them with the values given in the table
ofun = vtkPiecewiseFunction()
ofun.AddPoint(-4931.54, 1.0)
ofun.AddPoint(101.815, 0.002)
ofun.AddPoint(2594.97, 0.0)

#vtkVolumeProperty to represent properties associated with volume rendering
vproperty= vtkVolumeProperty()
vproperty.SetColor(ctfun)
vproperty.SetScalarOpacity(ofun)

#Setting the type of interpolation
vproperty.SetInterpolationTypeToLinear()

#User input for Phong Shading
print("Whether you want to visualize the dataset with Phong Shading?(yes/no): ")
choice = input()
choice = choice.lower()

#turning on Phong shading while rendering
if(choice=="yes"):
 #setting the Ambient,Diffiuse,Spectacular coefficient values as given 
 vproperty.ShadeOn()
 vproperty.SetAmbient(0.5)
 vproperty.SetDiffuse(0.5)
 vproperty.SetSpecular(0.5)

#vtkVolume to represent volumetric entity in the rendering process
volume= vtkVolume()
volume.SetMapper(volumeMapper)
volume.SetProperty(vproperty)

#Using vtkOutlineFilter to add an outline to the volume rendered data
outline=vtkOutlineFilter()
outline.SetInputConnection(data)

#Creating Mapper
poly= vtkPolyDataMapper()
poly.SetInputConnection(outline.GetOutputPort())

#Creating Actor and setting mapper to it
actor= vtkActor()
actor.SetMapper(poly)
actor.GetProperty().SetLineWidth(3)
actor.GetProperty().SetColor(0,0,0)

#Creating renderer and giving actor to it
rend=vtkRenderer()
rend.AddActor(actor)
rend.AddVolume(volume)
rend.SetBackground(1,1,1)

#Creating renderWindow as per given size
window= vtkRenderWindow()
window.AddRenderer(rend)
window.SetSize(1000,1000)

#Creating interactor
interactor= vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)

#Initializing interactor for data visualization
interactor.Initialize()
interactor.Start()





