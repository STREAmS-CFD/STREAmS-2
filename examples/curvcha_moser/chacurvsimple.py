# script-version: 2.0
# Catalyst state generated using paraview version 5.12.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2792, 1662]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [-0.029533780595091486, 78.87294632899682, 2.094019763533848]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.029533780595091486, 78.87294632899682, 16.904395275599924]
renderView1.CameraFocalPoint = [-0.029533780595091486, 78.87294632899682, 2.094019763533848]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 6.79076046484115
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(2792, 1662)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Partitioned Dataset Reader'
grid = XMLPartitionedDatasetReader(registrationName='grid', FileName=['F:\\BUTTARE\\LI\\grid_000001.vtpd'])

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=grid)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [-0.029533780595091486, 78.87294632899682, 2.0798506919617012]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [-0.029533780595091486, 78.87294632899682, 2.0798506919617012]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'u'
uTF2D = GetTransferFunction2D('u')

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')
uLUT.TransferFunction2D = uTF2D
uLUT.RGBPoints = [-0.022668755185689195, 0.278431372549, 0.278431372549, 0.858823529412, 0.0071560402252809525, 0.0, 0.0, 0.360784313725, 0.036772270633377194, 0.0, 1.0, 1.0, 0.06680563104722125, 0.0, 0.501960784314, 0.0, 0.09642186145531748, 1.0, 1.0, 0.0, 0.12624665686628767, 1.0, 0.380392156863, 0.0, 0.1560714522772578, 0.419607843137, 0.0, 0.0, 0.18589624768822793, 0.878431372549, 0.301960784314, 0.301960784314]
uLUT.ColorSpace = 'RGB'
uLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'u']
slice1Display.LookupTable = uLUT
slice1Display.SelectTCoordArray = 'None'
slice1Display.SelectNormalArray = 'None'
slice1Display.SelectTangentArray = 'None'
slice1Display.OSPRayScaleArray = 'density'
slice1Display.OSPRayScaleFunction = 'Piecewise Function'
slice1Display.Assembly = 'Hierarchy'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 1.284532980215289
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.06422664901076444
slice1Display.SetScaleArray = ['POINTS', 'density']
slice1Display.ScaleTransferFunction = 'Piecewise Function'
slice1Display.OpacityArray = ['POINTS', 'density']
slice1Display.OpacityTransferFunction = 'Piecewise Function'
slice1Display.DataAxesGrid = 'Grid Axes Representation'
slice1Display.PolarAxes = 'Polar Axes Representation'
slice1Display.SelectInputVectors = [None, '']
slice1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [0.30348586096551183, 0.0, 0.5, 0.0, 1.2736184965892372, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.9835601349720835, 0.0, 0.5, 0.0, 1.001918218709672, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.9835601349720835, 0.0, 0.5, 0.0, 1.001918218709672, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.Orientation = 'Horizontal'
uLUTColorBar.WindowLocation = 'Any Location'
uLUTColorBar.Position = [0.3265633954154728, 0.7556919374247892]
uLUTColorBar.Title = 'u'
uLUTColorBar.ComponentTitle = ''
uLUTColorBar.ScalarBarLength = 0.3299999999999999

# set color bar visibility
uLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')
uPWF.Points = [-0.022668755185689195, 0.0, 0.5, 0.0, 0.18589624768822793, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = renderView1
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 0.0007716789055173805
animationScene1.StartTime = 0.0007716789055173805
animationScene1.EndTime = 1.0007716789055174

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'Time Step'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [2792, 1662]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'Time Step'
options.CatalystLiveTrigger = 'Time Step'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
