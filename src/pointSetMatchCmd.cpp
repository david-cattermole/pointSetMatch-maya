
#include <pointSetMatchCmd.h>
#include <pointSetMatchUtils.h>

#include <debugUtils.h>

#include <testDataOrigin.h>
#include <testData2.h>

PointSetMatchCmd::~PointSetMatchCmd()
{}

void *PointSetMatchCmd::creator()
{
  return new PointSetMatchCmd();
}

// Set keyframes to move selected object in a spiral
MStatus PointSetMatchCmd::doIt(const MArgList &)
{
  MStatus status;

  std::vector<double> movPoints;
  std::vector<double> refPoints;
  std::vector<double> weights;
  std::vector<unsigned int> lockedParams;
  std::vector<double> initialTransform;
  std::vector<double> outTransform;

  INFO("testData2_num: " << testData2_num);
  INFO("testDataOrigin1_num: " << testDataOrigin1_num);
  setVectorFromArray(&testData2_values[0], testData2_num, movPoints);
  setVectorFromArray(&testDataOrigin1_values[0], testDataOrigin1_num, refPoints);

  int iterMax = 10000;
  bool uniformScale = true;
  weights.resize((unsigned long) testData2_num, 0.0);
  outTransform.resize(9);
  lockedParams.resize(9, 0);

  // Initial transforms
  initialTransform.resize(9);

  // Translate
  initialTransform[0] = 0.0;
  initialTransform[1] = 0.0;
  initialTransform[2] = 0.0;
  lockedParams[0] = 0;
  lockedParams[1] = 0;
  lockedParams[2] = 0;

  // Rotate
  initialTransform[3] = 0.0;
  initialTransform[4] = 0.0;
  initialTransform[5] = 0.0;
  lockedParams[3] = 0;
  lockedParams[4] = 0;
  lockedParams[5] = 0;

  // Scale
  initialTransform[6] = 1.0;
  initialTransform[7] = 1.0;
  initialTransform[8] = 1.0;
  lockedParams[6] = 0;
  lockedParams[7] = 0;
  lockedParams[8] = 0;

  double outError = -1.0;
  solveTransform(iterMax, uniformScale,
                 refPoints, movPoints,
                 weights, lockedParams, initialTransform,
                 outTransform, outError);

  return MS::kSuccess;
}



