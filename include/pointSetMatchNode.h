
#ifndef MAYA_POINT_SET_MATCH_NODE_H
#define MAYA_POINT_SET_MATCH_NODE_H

// STL
#include <string>
#include <cmath>

// Maya
#include <maya/MIOStream.h>
#include <maya/MPxNode.h>
#include <maya/MString.h>
#include <maya/MTypeId.h>
#include <maya/MPlug.h>


#define POINT_SET_MATCH_NODE_ID  0x40001

class PointSetMatchNode : public MPxNode
{
public:
  PointSetMatchNode();

  virtual ~PointSetMatchNode();

  virtual MStatus compute(const MPlug &plug, MDataBlock &dataBlock);

  static void *creator();

  static MStatus initialize();

public:

  static MObject aIterations;
  static MObject aRotateOrder;
  static MObject aUniformScale;

  static MObject aInputPoints;
  static MObject aInputPointsReferencePoint;
  static MObject aInputPointsReferencePointX;
  static MObject aInputPointsReferencePointY;
  static MObject aInputPointsReferencePointZ;
  static MObject aInputPointsMovePoint;
  static MObject aInputPointsMovePointX;
  static MObject aInputPointsMovePointY;
  static MObject aInputPointsMovePointZ;
  static MObject aInputPointsWeight;

  static MObject aInitialTranslateX;
  static MObject aInitialTranslateY;
  static MObject aInitialTranslateZ;
  static MObject aInitialRotateX;
  static MObject aInitialRotateY;
  static MObject aInitialRotateZ;
  static MObject aInitialScaleX;
  static MObject aInitialScaleY;
  static MObject aInitialScaleZ;

  static MObject aLockedTranslateX;
  static MObject aLockedTranslateY;
  static MObject aLockedTranslateZ;
  static MObject aLockedRotateX;
  static MObject aLockedRotateY;
  static MObject aLockedRotateZ;
  static MObject aLockedScaleX;
  static MObject aLockedScaleY;
  static MObject aLockedScaleZ;

  static MObject aOutputTranslateX;
  static MObject aOutputTranslateY;
  static MObject aOutputTranslateZ;
  static MObject aOutputRotateX;
  static MObject aOutputRotateY;
  static MObject aOutputRotateZ;
  static MObject aOutputScaleX;
  static MObject aOutputScaleY;
  static MObject aOutputScaleZ;
  static MObject aOutputError;

  static MTypeId id;
};


//inline bool
//equivalent(double a, double b  )
//{
//  return fabs( a - b ) < .000001 ;
//}

#endif // POINT_SET_MATCH_NODE_H