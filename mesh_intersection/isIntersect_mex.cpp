#include "mex.h"

#include <igl/readOBJ.h>
#include <igl/copyleft/boolean/mesh_boolean.h>
#include <igl/matlab/prepare_lhs.h>
#include <Eigen/Core>


void mexFunction(
     int          nlhs,
     mxArray      *plhs[],
     int          nrhs,
     const mxArray *prhs[]
     )
{
  using namespace Eigen;
  /* Check for proper number of arguments */

  if (nrhs != 2) 
  {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
        "readOBJ requires 1 input arguments, the path of the file to open");
  }

  // Read the file path
  char* file_path = mxArrayToString(prhs[0]);

  MatrixXd VA, VB, VC;
  MatrixXi FA, FB, FC;

  // Read the first mesh
  if(!igl::readOBJ(file_path,VA,FA))
  {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:fileio", "igl::readOBJ failed.");
  }
  
  // Read the second mesh
  if(!igl::readOBJ(file_path,VB,FB))
  {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:fileio", "igl::readOBJ failed.");
  }
  
  // Check for intersection
  igl::mesh_boolean(VA,FA,VB,FB,MESH_BOOLEAN_TYPE_UNION,VC,FC);
  // Return the matrices to matlab
  switch(nlhs)
  {
    case 2:
      igl::matlab::prepare_lhs_index(FC,plhs+1);
    case 1:
      igl::matlab::prepare_lhs_double(VC,plhs);
    default: break;
  }

  return;
}
