// Adaptor for getting Fortran simulation code into ParaView CoProcessor.
#ifdef SINGLE_PRECISION
#define REAL float
#define VTKREALARRAY vtkFloatArray
#else
#define REAL double
#define VTKREALARRAY vtkDoubleArray
#endif

// CoProcessor specific headers
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkRectilinearGrid.h"
#include "vtkUnsignedCharArray.h"
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkMultiProcessController.h>

// Fortran specific header
#include "vtkCPPythonAdaptorAPI.h"

// Josh specific headers
#include <iostream>
#include <algorithm>
//#include <vtkCPDataDescription.h>
//#include <vtkCPInputDataDescription.h>
//#include <vtkCPProcessor.h>
//#include <vtkCPPythonScriptPipeline.h>
//#include <vtkDoubleArray.h>
//#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
/////#include <vtkMPI.h> // never loaded before, but not used (perhaps)
#include <vtkNew.h>
//#include <vtkPointData.h>
#include <vtkPoints.h>
//#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
//#include <vtkSmartPointer.h>
#include <cuda_runtime.h>


// These will be called from the Fortran "glue" code"
// Completely dependent on data layout, structured vs. unstructured, etc.
// since VTK/ParaView uses different internal layouts for each.

// Creates the data container for the CoProcessor.
extern "C" void createcpimagedata(int* nxstart, int* nxend, int* nzstart, int* nzend,
 int* nx, int* nz,
 REAL* dx, REAL* dz,
 REAL* origin, char* zone, int* nrank)
{

  if (!vtkCPPythonAdaptorAPI::GetCoProcessorData())
  {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
  }
  if (*nrank==0)
    std::cout<<"Creating imageGrid channel for "<<zone<<"\n";

  // Name should be consistent between here, Fortran and Python client script.
  // Funzione che mi restituisce un oggetto.  mi definisce una cosa o, se ce l'ho già, me lo dà.
  // in questo caso io ce l'ho sempre
  vtkCPDataDescription* datadescrip = vtkCPPythonAdaptorAPI::GetCoProcessorData();
  datadescrip->AddInput(zone);

  // The simulation grid is a 3-dimensional topologically and geometrically
  // regular grid. In VTK/ParaView, this is considered an image data set.
  vtkSmartPointer<vtkImageData> grid = vtkSmartPointer<vtkImageData>::New();

  REAL zero=0.;
  REAL one =1.;
  grid->SetExtent(*nxstart - 1, *nxend - 1, 0, 0, *nzstart - 1, *nzend - 1);
  grid->SetSpacing(*dx, one, *dz);
  grid->SetOrigin( zero, *origin, zero);

  vtkCPInputDataDescription* inputdescrip = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(zone);
  inputdescrip->SetGrid(grid);
  inputdescrip->SetWholeExtent( 0, *nx - 1, 0, 0, 0, *nz - 1); 

}

// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// By hand name mangling for fortran.
extern "C" void addfieldtoimage(char* zone, REAL* scalars, char* name, int* nrank)
{
  vtkCPInputDataDescription* idd =
    vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(zone);
  if (!idd)
    return;
  if (*nrank==0)
    std::cout<<"Updating field "<<name<<" on imageGrid "<<zone<<"\n";

  vtkImageData* Image = vtkImageData::SafeDownCast(idd->GetGrid());
  if (!Image)
  {
    vtkGenericWarningMacro("image: No adaptor grid to attach field data to.");
    return;
  }

  // field name must match that in the fortran code.
  // old paraview: if (idd->IsFieldNeeded(name))
  if (idd->IsFieldNeeded(name, vtkDataObject::POINT) == true)
  {
    vtkSmartPointer<VTKREALARRAY> field = vtkSmartPointer<VTKREALARRAY>::New();
    field->SetNumberOfComponents(1);
    field->SetName(name);
    field->SetArray(scalars, Image->GetNumberOfPoints(), 1);
    Image->GetPointData()->AddArray(field);
  }

}

// Creates the data container for the CoProcessor.
extern "C" void createcprectilineardata(int* nxstart, int* nxend, int* nystart, int* nyend, int* nzstart, int* nzend,
 int* nxstartg, int* nystartg, int* nzstartg,
 int* nxendg, int* nyendg, int* nzendg,
 REAL* x, REAL* y, REAL* z,
 int* nrank, int* nproc,
 char* zone)
{

  if (!vtkCPPythonAdaptorAPI::GetCoProcessorData())
  {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
  }
  if (*nrank==0)
    std::cout<<"Creating rectilinearGrid channel for "<<zone<<"\n";

  // Name should be consistent between here, Fortran and Python client script.
  // Funzione che mi restituisce un oggetto.  mi definisce una cosa o, se ce l'ho già, me lo dà.
  // in questo caso io ce l'ho sempre
  vtkCPDataDescription* datadescrip = vtkCPPythonAdaptorAPI::GetCoProcessorData();
  datadescrip->AddInput(zone);







  // The simulation grid is a 3-dimensional topologically and geometrically
  // regular grid, not evenly spaced. In VTK/ParaView, this is considered a rectilinear data set.
  vtkSmartPointer<vtkRectilinearGrid> grid = vtkSmartPointer<vtkRectilinearGrid>::New();

  // Set extents of subdomain
  grid->SetExtent(*nxstart, *nxend, *nystart, *nyend, *nzstart, *nzend);

  vtkSmartPointer<VTKREALARRAY> xCoords = vtkSmartPointer<VTKREALARRAY>::New();
  xCoords->SetNumberOfComponents(1);
  int sizex = *nxend-*nxstart+1;
  xCoords->SetArray (x, sizex, 1);

  vtkSmartPointer<VTKREALARRAY> yCoords = vtkSmartPointer<VTKREALARRAY>::New();
  yCoords->SetNumberOfComponents(1);
  int sizey = *nyend-*nystart+1;
  yCoords->SetArray (y, sizey, 1);

  vtkSmartPointer<VTKREALARRAY> zCoords = vtkSmartPointer<VTKREALARRAY>::New();
  zCoords->SetNumberOfComponents(1);
  int sizez = *nzend-*nzstart+1;
  zCoords->SetArray (z, sizez, 1);

  grid->SetXCoordinates(xCoords);
  grid->SetYCoordinates(yCoords);
  grid->SetZCoordinates(zCoords);

  vtkCPInputDataDescription* inputdescrip = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(zone);

  inputdescrip->SetGrid(grid);
  inputdescrip->SetWholeExtent( *nxstartg, *nxendg, *nystartg, *nyendg, *nzstartg, *nzendg); 

//  // Mark ghost cells
//  vtkNew<vtkUnsignedCharArray> ghostCells;
//  ghostCells->SetNumberOfTuples(grid->GetNumberOfCells());
//  ghostCells->SetName(vtkDataSetAttributes::GhostArrayName());
//  ghostCells->Fill(0);
//  grid->GetCellData()->AddArray(ghostCells);
//
//  vtkIdType cellId;
//  // imax boundary
//  for (int k=*nzstart-1; k<*nzend; k++)
//  {
//    for (int j=*nystart-1; j<*nyend; j++)
//    {
//       int i=nxend-1;
//       cellId = i*ny*
//       ghostCells->SetValue(cellId,vtkDataSetAttributes::DUPLICATECELL);
//    }
//  }

  //vtkMultiBlockDataSet* MBGrid = vtkMultiBlockDataSet::New();
  //vtkSmartPointer<vtkMultiPieceDataSet> multiPiece = vtkSmartPointer<vtkMultiPieceDataSet>::New();
  //multiPiece->SetNumberOfPieces(*nproc);
  //multiPiece->SetPiece(*nrank, grid);
  //MBGrid->SetNumberOfBlocks(1);
  //MBGrid->SetBlock(0, multiPiece.GetPointer());

  //inputdescrip->SetGrid(MBGrid);

}

// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// By hand name mangling for fortran.
extern "C" void addfieldtorectilinear(char* zone, REAL* scalars, char* name, int* nrank)
{
  vtkCPInputDataDescription* idd =
    vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(zone);
  if (!idd)
    return;
  if (*nrank==0)
    std::cout<<"Updating field "<<name<<" on rectilinearGrid "<<zone<<"\n";




  vtkRectilinearGrid* RGrid = vtkRectilinearGrid::SafeDownCast(idd->GetGrid());
  if (!RGrid)
  {
    vtkGenericWarningMacro("RGrid: No adaptor grid to attach field data to.");
    return;
  }

  // field name must match that in the fortran code.
  // old paraview: if (idd->IsFieldNeeded(name))
  if (idd->IsFieldNeeded(name, vtkDataObject::POINT) == true)
  {
    vtkSmartPointer<VTKREALARRAY> field = vtkSmartPointer<VTKREALARRAY>::New();
    field->SetNumberOfComponents(1);
    field->SetName(name);
    field->SetArray(scalars, RGrid->GetNumberOfPoints(), 1);
    RGrid->GetPointData()->AddArray(field);
  }







//  // blocco madre
//  vtkMultiBlockDataSet* multiBlock = vtkMultiBlockDataSet::SafeDownCast(idd->GetGrid());
//  vtkMultiPieceDataSet* multiPiece = vtkMultiPieceDataSet::SafeDownCast(multiBlock->GetBlock(0));

//  if (!multiBlock)
//  {
//    vtkGenericWarningMacro("multiBlock: No adaptor grid to attach field data to.");
//    return;
//  }
//
//  if (!multiPiece)
//  {
//    vtkGenericWarningMacro("multiPiece: No adaptor grid to attach field data to.");
//    return;
//  }

//  // field name must match that in the fortran code.
//  // old paraview: if (idd->IsFieldNeeded(name))
//  if (idd->IsFieldNeeded(name, vtkDataObject::POINT) == true)
//  {
////    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(multiPiece->GetPiece(*nrank));
//    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(idd->GetGrid());
//    vtkSmartPointer<VTKREALARRAY> field = vtkSmartPointer<VTKREALARRAY>::New();
//    field->SetNumberOfComponents(1);
//    field->SetName(name);
//    field->SetArray(scalars, dataSet->GetNumberOfPoints(), 1);
//    dataSet->GetPointData()->AddArray(field);
//  }

}









// Creates the data container for the CoProcessor.
extern "C" void createcpstructureddata(int* nxp2, int* nyp2, int* nzp2,
 int* nxstart, int* nystart, int* nzstart, int* nxend, int* nyend, int* nzend,
 int* nxmaxp2, int* nymaxp2, int* nzmaxp2,
 REAL* x, REAL* y, REAL* z,
 int* nrank, int* nproc,
 char* zone,
 REAL* xyzc)
{

  if (!vtkCPPythonAdaptorAPI::GetCoProcessorData())
  {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
  }
  if (*nrank==0)
    std::cout<<"Creating structuredGrid channel for "<<zone<<"\n";

  // Name should be consistent between here, Fortran and Python client script.
  // Funzione che mi restituisce un oggetto.  mi definisce una cosa o, se ce l'ho già, me lo dà.
  // in questo caso io ce l'ho sempre

  //vtkNew<vtkCPDataDescription> datadescrip; // = static_cast<vtkNew<vtkCPDataDescription>>( vtkCPPythonAdaptorAPI::GetCoProcessorData() );
  //datadescrip->AddInput(zone);
  //datadescrip->SetTimeData(time, timeStep);

  vtkCPDataDescription* datadescrip = vtkCPPythonAdaptorAPI::GetCoProcessorData();
  datadescrip->AddInput(zone);

  // The simulation grid is a 3-dimensional topologically and geometrically
  // regular grid, not evenly spaced. In VTK/ParaView, this is considered a rectilinear data set.
  // But Paraview crashes when using VTKm contour on rectilinear grid, so here we build
  // the structured version
  vtkStructuredGrid* grid = vtkStructuredGrid::New();
  //vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();

  // JOSH /////////////////////////
  float* point = (float*) xyzc;
  for (int k = 0; k < *nzp2; ++k) {
    for (int j = 0; j < *nyp2; ++j) {
      for (int i = 0; i < *nxp2; ++i) {
        point[0] = x[i];
        point[1] = y[j];
        point[2] = z[k];
	point += 3;
      }
    }
  }
  vtkNew<vtkFloatArray> pointArray;
  pointArray->SetNumberOfComponents(3);
  pointArray->SetArray((float*) xyzc, static_cast<vtkIdType>((*nxp2) * (*nyp2) * (*nzp2) * 3), 1);
  vtkNew<vtkPoints> points;
  points->SetData(pointArray.GetPointer());

  grid->SetPoints(points);
  //grid->SetExtent(0,*nxp2-1,0,*nyp2-1,0,*nzp2-1);
  grid->SetExtent(*nxstart, *nxend, *nystart, *nyend, *nzstart, *nzend);

  datadescrip->GetInputDescriptionByName(zone)->SetGrid(grid);
  datadescrip->GetInputDescriptionByName(zone)->SetWholeExtent( 0, *nxmaxp2-1, 1, *nymaxp2, 0, *nzmaxp2-1); 

  //vtkCPInputDataDescription* inputdescrip = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(zone);
  //inputdescrip->SetGrid(grid);
  //inputdescrip->SetWholeExtent( 0, *nxmaxp2-1, 1, *nymaxp2, 0, *nzmaxp2-1); 

}









// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// By hand name mangling for fortran.
extern "C" void addfieldtostructured(char* zone, REAL* scalars, char* name, int* nrank)
{
  vtkCPInputDataDescription* idd =
    vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName(zone);
  if (!idd)
    return;
  if (*nrank==0)
    std::cout<<"Updating field "<<name<<" on structuredGrid "<<zone<<"\n";




  vtkStructuredGrid* SGrid = vtkStructuredGrid::SafeDownCast(idd->GetGrid());
  if (!SGrid)
  {
    vtkGenericWarningMacro("SGrid: No adaptor grid to attach field data to.");
    return;
  }

  // field name must match that in the fortran code.
  // old paraview: if (idd->IsFieldNeeded(name))
  if (idd->IsFieldNeeded(name, vtkDataObject::POINT) == true)
  {
    vtkSmartPointer<VTKREALARRAY> field = vtkSmartPointer<VTKREALARRAY>::New();
    field->SetName(name);
    field->SetNumberOfComponents(1);
    //field->SetArray(scalars, static_cast<vtkIdType>(*ntot), 1); // Josh style
    field->SetArray(scalars, SGrid->GetNumberOfPoints(), 1); // example style
    SGrid->GetPointData()->AddArray(field.GetPointer());
  }

}

