
#ifndef MAYA_POINT_SET_MATCH_UTILS_H
#define MAYA_POINT_SET_MATCH_UTILS_H

// Lev-Mar
#include <levmar.h>  //

// STL
#include <ctime>     // time
#include <cmath>     // exp
#include <iostream>  // cout, cerr, endl
#include <string>    // string
#include <vector>    // vector
#include <cassert>   // assert
//
#include <math.h>
#include <vector>
#include <string>

// Utils
#include <debugUtils.h>

// Maya
#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <maya/MMatrix.h>
#include <maya/MString.h>
#include <maya/MTransformationMatrix.h>


// Lev-Mar Termination Reasons:
const std::string reasons[8] = {
    // reason 0
    "no reason, should not get here",

    // reason 1
    "stopped by small gradient J^T e",

    // reason 2
    "stopped by small Dp",

    // reason 3
    "stopped by itmax",

    // reason 4
    "singular matrix. Restart from current p with increased \\mu",

    // reason 5
    "no further error reduction is possible. Restart with increased mu",

    // reason 6
    "stopped by small ||e||_2",

    // reason 7
    "stopped by invalid (i.e. NaN or Inf) \"func\" refPoints (user error)",
};


struct PointCloudData
{
  // Number of points
  int numPoints;

  // Source Points to match to
  double *refPoints;

  // Movable points
  double *movPoints;

//  // List of weights for each pair of points
//  double *weights;

  // Initial transform parameters
  double *initial;

  // Locked parameters
  unsigned int *locked;
};


// Get the distance between two 3D points.
inline
double distance(double *a, double *b)
{
  double dx = (a[0] - b[0]);
  double dy = (a[1] - b[1]);
  double dz = (a[2] - b[2]);
  return sqrt((dx * dx) + (dy * dy) + (dz * dz));
}


inline
void curveFunc(double *p, double *x, int m, int n, void *data)
{
  register int i, j;

  PointCloudData *userData = (PointCloudData *) data;

  double params[9] = {
      0, 0, 0,
      0, 0, 0,
      1, 1, 1,
  };

  j = 0;
  for (i=0; i < 9; ++i)
  {
    if (userData->locked[i])
    {
      params[i] = userData->initial[i];
    }
    else
    {
      params[i] = p[j];
      ++j;
    }
  }

  const double degreesToRadians = 57.29577951308232;
  double rotate[3] = {params[3] / degreesToRadians,
                      params[4] / degreesToRadians,
                      params[5] / degreesToRadians};

  MVector translateVec(params[0], params[1], params[2]);
  MTransformationMatrix tfmMatrix;
  MTransformationMatrix::RotationOrder order;
  order = MTransformationMatrix::RotationOrder::kXYZ;
  tfmMatrix.setTranslation(translateVec, MSpace::kWorld);
  tfmMatrix.setRotation(&rotate[0], order);
  tfmMatrix.setScale(&params[6], MSpace::kWorld);
  MMatrix matrix = tfmMatrix.asMatrix().inverse();

//  std::cerr << n << "=";
  for (i = 0; i < (n / 4); ++i)
  {
//    INFO("i=" << i);
//    double weight = userData->weights[i];
    double srcPoint[3] = {userData->refPoints[i * 3 + 0],
                          userData->refPoints[i * 3 + 1],
                          userData->refPoints[i * 3 + 2]};
    double movPoint[3] = {userData->movPoints[i * 3 + 0],
                          userData->movPoints[i * 3 + 1],
                          userData->movPoints[i * 3 + 2]};
//    INFO("srcX=" << srcPoint[0] << ","
//                 << srcPoint[1] << ","
//                 << srcPoint[2]);
//
//    INFO("movX=" << movPoint[0] << ","
//                 << movPoint[1] << ","
//                 << movPoint[2]);

    MPoint point(movPoint[0], movPoint[1], movPoint[2]);
    point *= matrix;

    movPoint[0] = point.x;
    movPoint[1] = point.y;
    movPoint[2] = point.z;

//    INFO("aftX=" << movPoint[0] << ","
//                 << movPoint[1] << ","
//                 << movPoint[2]);

    double diff = distance(srcPoint, movPoint);
//    INFO("diff=" << diff);
    // diff *= weight;
    double diffX = fabs(movPoint[0] - srcPoint[0]);
    double diffY = fabs(movPoint[1] - srcPoint[1]);
    double diffZ = fabs(movPoint[2] - srcPoint[2]);
    x[i * 4 + 0] = 0.5 * (diff * diff);
    x[i * 4 + 1] = 0.5 * (diffX * diffX);
    x[i * 4 + 2] = 0.5 * (diffY * diffY);
    x[i * 4 + 3] = 0.5 * (diffZ * diffZ);
  }
//  std::cerr << std::endl << std::endl;
}


inline
void setVectorFromArray(double *inValues, int valuesNum,
                        std::vector<double> &outVector)
{
  outVector.resize((unsigned long) (valuesNum * 3));
  std::vector<double>::iterator it = outVector.begin();
  for (int i = 0; it != outVector.end();)
  {
    *it = inValues[i];
    it++;
    i++;
  }
  return;
}


inline
bool solveTransform(int iterMax,
                    bool uniformScale,
                    std::vector<double> &refPoints,
                    std::vector<double> &movPoints,
                    std::vector<double> &weights,
                    std::vector<unsigned int> &lockedParams,
                    std::vector<double> &initialTransform,
                    std::vector<double> &outTransform,
                    double &outError)
{
//  INFO("refPoints.size(): " << refPoints.size());
//  INFO("movPoints.size(): " << movPoints.size());
//  INFO("weights.size(): " << weights.size());
  assert(refPoints.size() == movPoints.size());
  assert((refPoints.size() / 3) == weights.size());
  assert((lockedParams.size()) == outTransform.size());

  register int i, j;
  int ret;

  int unlockedNum = 0;
  for(i = 0; i < lockedParams.size(); ++i)
  {
    if(lockedParams[i] == 0)
    {
      ++unlockedNum;
    }
  }
  INFO("Unlocked Parameters: " << unlockedNum);

  double opts[LM_OPTS_SZ];
  double info[LM_INFO_SZ];
  const int m = unlockedNum; // number of unknown parameters
  double params[m];
  int n = (int) weights.size() * 4; // number of measurement errors

  // Options
  opts[0] = LM_INIT_MU * 100.0;
  opts[1] = 1E-15;
  opts[2] = 1E-15;
  opts[3] = 1E-20;
  opts[4] = -LM_DIFF_DELTA * 10.0;

  struct PointCloudData userData;
  userData.numPoints = n;
  userData.refPoints = &refPoints[0];
  userData.movPoints = &movPoints[0];
//  userData.weights = &weights[0];
  userData.initial = &initialTransform[0];
  userData.locked = &lockedParams[0];
//  userData.uniformScale = uniformScale;

  // Set Initial parameters
  j = 0;
  for (i = 0; i < 9; ++i)
  {
    if (!lockedParams[i])
    {
      params[j] = initialTransform[i];
      ++j;
    }
  }

  // Initial Parameters
  INFO("Initial Parameters: ");
  for (i = 0; i < m; ++i)
  {
    INFO("-> " << params[i]);
  }
  INFO("");

  // Allocate a memory block for both 'work' and 'covar', so that
  // the block is close together in physical memory.
  double *work, *covar;
  work = (double *) malloc((LM_DIF_WORKSZ(m, n) + m * m) * sizeof(double));
  if (!work)
  {
    ERR("Memory allocation request failed in myLibrary()");
    return false;
  }
  covar = work + LM_DIF_WORKSZ(m, n);

  // no Jacobian, caller allocates work memory, covariance estimated
  ret = dlevmar_dif(

      // Function to call (input only)
      // Function must be of the structure:
      //   func(double *params, double *x, int m, int n, void *data)
      curveFunc,

      // Parameters (input and output)
      // Should be filled with initial estimate, will be filled
      // with output parameters
      params,

      // Measurement Vector (input only)
      // NULL implies a zero vector
      NULL,

      // Parameter Vector Dimension (input only)
      // (i.e. #unknowns)
      m,

      // Measurement Vector Dimension (input only)
      n,

      // Maximum Number of Iterations (input only)
      iterMax,

      // Minimisation options (input only)
      // opts[0] = tau      (scale factor for initialTransform mu)
      // opts[1] = epsilon1 (stopping threshold for ||J^T e||_inf)
      // opts[2] = epsilon2 (stopping threshold for ||Dp||_2)
      // opts[3] = epsilon3 (stopping threshold for ||e||_2)
      // opts[4] = delta    (step used in difference approximation to the Jacobian)
      //
      // If \delta<0, the Jacobian is approximated with central differences
      // which are more accurate (but slower!) compared to the forward
      // differences employed by default.
      // Set to NULL for defaults to be used.
      opts,

      // Output Information (output only)
      // information regarding the minimization.
      // info[0] = ||e||_2 at initialTransform params.
      // info[1-4] = (all computed at estimated params)
      //  [
      //   ||e||_2,
      //   ||J^T e||_inf,
      //   ||Dp||_2,
      //   \mu/max[J^T J]_ii
      //  ]
      // info[5] = number of iterations,
      // info[6] = reason for terminating:
      //   1 - stopped by small gradient J^T e
      //   2 - stopped by small Dp
      //   3 - stopped by iterMax
      //   4 - singular matrix. Restart from current params with increased \mu
      //   5 - no further error reduction is possible. Restart with increased mu
      //   6 - stopped by small ||e||_2
      //   7 - stopped by invalid (i.e. NaN or Inf) "func" refPoints; a user error
      // info[7] = number of function evaluations
      // info[8] = number of Jacobian evaluations
      // info[9] = number linear systems solved (number of attempts for reducing error)
      //
      // Set to NULL if don't care
      info,

      // Working Data (input only)
      // working memory, allocated internally if NULL. If !=NULL, it is assumed to
      // point to a memory chunk at least LM_DIF_WORKSZ(m, n)*sizeof(double) bytes
      // long
      work,

      // Covariance matrix (output only)
      // Covariance matrix corresponding to LS solution; Assumed to point to a mxm matrix.
      // Set to NULL if not needed.
      covar,

      // Custom Data for 'func' (input only)
      // pointer to possibly needed additional data, passed uninterpreted to func.
      // Set to NULL if not needed
      (void *) &userData);

  INFO("Covariance of the fit:");
  for (i = 0; i < m; ++i)
  {
    for (j = 0; j < m; ++j)
    {
      INFO(covar[i * m + j]);
    }
    INFO("");
  }
  INFO("");

  free(work);

  INFO("Results:");
  INFO("Levenberg-Marquardt returned " << ret << " in " << (int) info[5]
                                       << " iterations");

  int reasonNum = (int) info[6];
  INFO("Reason: " << reasons[reasonNum]);
  INFO("Reason number: " << info[6]);
  INFO("");

  INFO("Solved Parameters:");
  for (i = 0; i < m; ++i)
  {
    INFO("-> " << params[i]);
  }
  INFO("");

  j = 0;
  for (i = 0; i < 9; ++i)
  {
    if (lockedParams[i])
    {
      outTransform[i] = initialTransform[i];
    }
    else
    {
      outTransform[i] = params[j];
      ++j;
    }
  }

  const double degreesToRadians = 57.29577951308232;
  INFO("tx=" << outTransform[0]);
  INFO("ty=" << outTransform[1]);
  INFO("tz=" << outTransform[2]);
  INFO("rx=" << outTransform[3]);
  INFO("ry=" << outTransform[4]);
  INFO("rz=" << outTransform[5]);
  INFO("sx=" << outTransform[6]);
  INFO("sy=" << outTransform[7]);
  INFO("sz=" << outTransform[8]);
  INFO("");

  INFO(std::endl << std::endl << "Solve Information:");
  INFO("Initial Error: " << info[0]);

  INFO("Overall Error: " << info[1]);
  INFO("J^T Error: " << info[2]);
  INFO("Dp Error: " << info[3]);
  INFO("Max Error: " << info[4]);

  INFO("Iterations: " << info[5]);
  INFO("Termination Reason: " << reasons[reasonNum]);
  INFO("Function Evaluations: " << info[7]);
  INFO("Jacobian Evaluations: " << info[8]);
  INFO("Attempts for reducing error: " << info[9]);

  outError = info[1];

  if (ret == -1)
  {
    return false;
  }
  return true;
}


#endif // MAYA_POINT_SET_MATCH_UTILS_H
