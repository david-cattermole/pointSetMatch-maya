
#ifndef POINT_SET_MATCH_CMD_H
#define POINT_SET_MATCH_CMD_H

//#include <math.h>
#include <maya/MIOStream.h>

#include <maya/MPxCommand.h>

#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <maya/MMatrix.h>
#include <maya/MDagPath.h>
#include <maya/MFnAnimCurve.h>
#include <maya/MFnDagNode.h>
#include <maya/MString.h>


// Match Objects Command
class PointSetMatchCmd : public MPxCommand
{
public:

  PointSetMatchCmd()
  {};

  virtual ~PointSetMatchCmd();

  virtual MStatus doIt(const MArgList &args);

  static void *creator();
};

#endif // POINT_SET_MATCH_CMD_H
