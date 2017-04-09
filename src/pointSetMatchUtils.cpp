

//
#include <pointSetMatchUtils.h>

// Lev-Mar
#include <levmar.h>  //

// STL
#include <ctime>     // time
#include <cmath>     // exp
#include <iostream>  // cout, cerr, endl
#include <string>    // string
#include <vector>    // vector
#include <cassert>   // assert

// Utils
#include <debugUtils.h>

// Maya
#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <maya/MMatrix.h>
#include <maya/MString.h>
#include <maya/MTransformationMatrix.h>


//// Get the distance between two 3D points.
//double distance(double *a, double *b)
//{
//  double dx = (a[0] - b[0]);
//  double dy = (a[1] - b[1]);
//  double dz = (a[2] - b[2]);
//  return sqrt((dx * dx) + (dy * dy) + (dz * dz));
//}
//
//
//void curveFunc(double *p, double *x, int m, int n, void *data)
//{
//  PointCloudData *userData = (PointCloudData *) data;
//
//  double params[9] = {
//      userData->locked[0] ? userData->initial[0] : p[0],
//      userData->locked[1] ? userData->initial[1] : p[1],
//      userData->locked[2] ? userData->initial[2] : p[2],
//      userData->locked[3] ? userData->initial[3] : p[3],
//      userData->locked[4] ? userData->initial[4] : p[4],
//      userData->locked[5] ? userData->initial[5] : p[5],
//      userData->locked[6] ? userData->initial[6] : p[6],
//      userData->locked[7] ? userData->initial[7] : p[7],
//      userData->locked[8] ? userData->initial[8] : p[8],
//  };
//
//  register int i;
//  MVector translateVec(params[0], params[1], params[2]);
//  double *rotate = &params[3];
//  double scale[3] = {params[6], params[7], params[8]};
//  if (userData->uniformScale)
//  {
//    scale[0] = params[6];
//    scale[1] = params[6];
//    scale[2] = params[6];
//  }
//
//  MTransformationMatrix tfmMatrix;
//  MTransformationMatrix::RotationOrder order;
//  order = MTransformationMatrix::RotationOrder::kXYZ;
//  tfmMatrix.setTranslation(translateVec, MSpace::kWorld);
//  tfmMatrix.setRotation(rotate, order, MSpace::kWorld);
//  tfmMatrix.setScale(&scale[0], MSpace::kWorld);
//  MMatrix matrix = tfmMatrix.asMatrix();
//
//  for (i = 0; i < n; ++i)
//  {
////    double weight = userData->weights[i];
//    double srcPoint[3] = {userData->refPoints[i * 3 + 0],
//                          userData->refPoints[i * 3 + 1],
//                          userData->refPoints[i * 3 + 2]};
//    double movPoint[3] = {userData->movPoints[i * 3 + 0],
//                          userData->movPoints[i * 3 + 1],
//                          userData->movPoints[i * 3 + 2]};
//
//    MPoint point(movPoint[0], movPoint[1], movPoint[2]);
//    point *= matrix;
//
//    movPoint[0] = point.x;
//    movPoint[1] = point.y;
//    movPoint[2] = point.z;
//
//    double diff = distance(srcPoint, movPoint);
//    // diff *= weight;
//    x[i] = 0.5 * (diff * diff);
//  }
//}
//
//
//void setVectorFromArray(double *inValues, int valuesNum,
//                        std::vector<double> &outVector)
//{
//  outVector.resize((unsigned long) (valuesNum * 3));
//  std::vector<double>::iterator it = outVector.begin();
//  for (int i = 0; it != outVector.end();)
//  {
//    *it = inValues[i];
//    it++;
//    i++;
//  }
//  return;
//}
//
//
//bool solveTransform(int iterMax,
//                    bool uniformScale,
//                    std::vector<double> &refPoints,
//                    std::vector<double> &movPoints,
//                    std::vector<double> &weights,
//                    std::vector<unsigned int> &lockedParams,
//                    std::vector<double> &initialTransform,
//                    std::vector<double> &outTransform)
//{
//  INFO("refPoints.size(): " << refPoints.size());
//  INFO("movPoints.size(): " << movPoints.size());
//  INFO("weights.size(): " << weights.size());
//  assert(refPoints.size() == movPoints.size());
//  assert((refPoints.size() / 3) == weights.size());
//  assert((lockedParams.size()) == outTransform.size());
//  outTransform.resize(9);
//
//  register int i, j;
//  int ret;
//
//  double opts[LM_OPTS_SZ];
//  double info[LM_INFO_SZ];
//  const int m = 9; // number of unknown parameters
//  double params[m];
//  int n = (int) weights.size(); // number of measurement errors
//
//  // Options
//  opts[0] = LM_INIT_MU;
//  opts[1] = 1E-15;
//  opts[2] = 1E-15;
//  opts[3] = 1E-20;
//  opts[4] = LM_DIFF_DELTA; // Relevant only if the Jacobian is approximated using finite differences (which we are not using)
//
//  struct PointCloudData userData;
//  userData.numPoints = n;
//  userData.refPoints = &refPoints[0];
//  userData.movPoints = &movPoints[0];
////  userData.weights = &weights[0];
//  userData.initial = &initialTransform[0];
//  userData.locked = &lockedParams[0];
////  userData.uniformScale = true;
//
//  // Initial Parameters
//  INFO("Initial Parameters: ");
//  for (i = 0; i < m; ++i)
//  {
//    params[i] = initialTransform[i];
//    INFO("-> " << params[i]);
//  }
//  INFO("");
//
//  // Allocate a memory block for both 'work' and 'covar', so that
//  // the block is close together in physical memory.
//  double *work, *covar;
//  work = (double *) malloc((LM_DIF_WORKSZ(m, n) + m * m) * sizeof(double));
//  if (!work)
//  {
//    ERR("Memory allocation request failed in myLibrary()");
//    return false;
//  }
//  covar = work + LM_DIF_WORKSZ(m, n);
//
//  // no Jacobian, caller allocates work memory, covariance estimated
//  ret = dlevmar_dif(
//
//      // Function to call (input only)
//      // Function must be of the structure:
//      //   func(double *params, double *x, int m, int n, void *data)
//      curveFunc,
//
//      // Parameters (input and output)
//      // Should be filled with initialTransform estimate, will be filled
//      // with aOutput parameters
//      params,
//
//      // Measurement Vector (input only)
//      // NULL implies a zero vector
//      NULL,
//
//      // Parameter Vector Dimension (input only)
//      // (i.e. #unknowns)
//      m,
//
//      // Measurement Vector Dimension (input only)
//      n,
//
//      // Maximum Number of Iterations (input only)
//      iterMax,
//
//      // Minimisation options (input only)
//      // opts[0] = tau      (scale factor for initialTransform mu)
//      // opts[1] = epsilon1 (stopping threshold for ||J^T e||_inf)
//      // opts[2] = epsilon2 (stopping threshold for ||Dp||_2)
//      // opts[3] = epsilon3 (stopping threshold for ||e||_2)
//      // opts[4] = delta    (step used in difference approximation to the Jacobian)
//      //
//      // If \delta<0, the Jacobian is approximated with central differences
//      // which are more accurate (but slower!) compared to the forward
//      // differences employed by default.
//      // Set to NULL for defaults to be used.
//      opts,
//
//      // Output Information (output only)
//      // information regarding the minimization.
//      // info[0] = ||e||_2 at initialTransform params.
//      // info[1-4] = (all computed at estimated params)
//      //  [
//      //   ||e||_2,
//      //   ||J^T e||_inf,
//      //   ||Dp||_2,
//      //   \mu/max[J^T J]_ii
//      //  ]
//      // info[5] = number of iterations,
//      // info[6] = reason for terminating:
//      //   1 - stopped by small gradient J^T e
//      //   2 - stopped by small Dp
//      //   3 - stopped by iterMax
//      //   4 - singular matrix. Restart from current params with increased \mu
//      //   5 - no further error reduction is possible. Restart with increased mu
//      //   6 - stopped by small ||e||_2
//      //   7 - stopped by invalid (i.e. NaN or Inf) "func" refPoints; a user error
//      // info[7] = number of function evaluations
//      // info[8] = number of Jacobian evaluations
//      // info[9] = number linear systems solved (number of attempts for reducing error)
//      //
//      // Set to NULL if don't care
//      info,
//
//      // Working Data (input only)
//      // working memory, allocated internally if NULL. If !=NULL, it is assumed to
//      // point to a memory chunk at least LM_DIF_WORKSZ(m, n)*sizeof(double) bytes
//      // long
//      work,
//
//      // Covariance matrix (output only)
//      // Covariance matrix corresponding to LS solution; Assumed to point to a mxm matrix.
//      // Set to NULL if not needed.
//      covar,
//
//      // Custom Data for 'func' (input only)
//      // pointer to possibly needed additional data, passed uninterpreted to func.
//      // Set to NULL if not needed
//      (void *) &userData);
//
//  INFO("Covariance of the fit:");
//  for (i = 0; i < m; ++i)
//  {
//    for (j = 0; j < m; ++j)
//    {
//      INFO(covar[i * m + j]);
//    }
//    INFO("");
//  }
//  INFO("");
//
//  free(work);
//
//  INFO("Results:");
//  INFO("Levenberg-Marquardt returned " << ret << " in " << (int) info[5]
//                                       << " iterations");
//
//  int reasonNum = (int) info[6];
//  INFO("Reason: " << reasons[reasonNum]);
//  INFO("Reason number: " << info[6]);
//  INFO("");
//
//  INFO("Solved Parameters:");
//  for (i = 0; i < m; ++i)
//  {
//    INFO("-> " << params[i]);
//    outTransform[i] = params[i];
//  }
//  INFO("");
//  const double degreesToRadians = 57.29577951308232;
//  INFO("tx=" << outTransform[0]);
//  INFO("ty=" << outTransform[1]);
//  INFO("tz=" << outTransform[2]);
//  INFO("rx=" << outTransform[3] * degreesToRadians);
//  INFO("ry=" << outTransform[4] * degreesToRadians);
//  INFO("rz=" << outTransform[5] * degreesToRadians);
//  INFO("sx=" << outTransform[6]);
//  INFO("sy=" << outTransform[7]);
//  INFO("sz=" << outTransform[8]);
//  INFO("");
//
//  INFO(std::endl << std::endl << "Solve Information:");
//  INFO("Initial Error: " << info[0]);
//
//  INFO("Overall Error: " << info[1]);
//  INFO("J^T Error: " << info[2]);
//  INFO("Dp Error: " << info[3]);
//  INFO("Max Error: " << info[4]);
//
//  INFO("Iterations: " << info[5]);
//  INFO("Termination Reason: " << reasons[reasonNum]);
//  INFO("Function Evaluations: " << info[7]);
//  INFO("Jacobian Evaluations: " << info[8]);
//  INFO("Attempts for reducing error: " << info[9]);
//  return true;
//}
//
