

#include <pointSetMatchNode.h>
#include <pointSetMatchUtils.h>

#include <maya/MString.h>
#include <maya/MVector.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>

// Attributes
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnCompoundAttribute.h>

// Debug
#include <utilities/debugUtils.h>

MTypeId PointSetMatchNode::id(POINT_SET_MATCH_NODE_ID);

MObject PointSetMatchNode::aIterations;
MObject PointSetMatchNode::aRotateOrder;
MObject PointSetMatchNode::aUniformScale;

MObject PointSetMatchNode::aInputPoints;
MObject PointSetMatchNode::aInputPointsReferencePoint;
MObject PointSetMatchNode::aInputPointsReferencePointX;
MObject PointSetMatchNode::aInputPointsReferencePointY;
MObject PointSetMatchNode::aInputPointsReferencePointZ;
MObject PointSetMatchNode::aInputPointsMovePoint;
MObject PointSetMatchNode::aInputPointsMovePointX;
MObject PointSetMatchNode::aInputPointsMovePointY;
MObject PointSetMatchNode::aInputPointsMovePointZ;
MObject PointSetMatchNode::aInputPointsWeight;

MObject PointSetMatchNode::aInitialTranslateX;
MObject PointSetMatchNode::aInitialTranslateY;
MObject PointSetMatchNode::aInitialTranslateZ;
MObject PointSetMatchNode::aInitialRotateX;
MObject PointSetMatchNode::aInitialRotateY;
MObject PointSetMatchNode::aInitialRotateZ;
MObject PointSetMatchNode::aInitialScaleX;
MObject PointSetMatchNode::aInitialScaleY;
MObject PointSetMatchNode::aInitialScaleZ;

MObject PointSetMatchNode::aLockedTranslateX;
MObject PointSetMatchNode::aLockedTranslateY;
MObject PointSetMatchNode::aLockedTranslateZ;
MObject PointSetMatchNode::aLockedRotateX;
MObject PointSetMatchNode::aLockedRotateY;
MObject PointSetMatchNode::aLockedRotateZ;
MObject PointSetMatchNode::aLockedScaleX;
MObject PointSetMatchNode::aLockedScaleY;
MObject PointSetMatchNode::aLockedScaleZ;

MObject PointSetMatchNode::aOutputTranslateX;
MObject PointSetMatchNode::aOutputTranslateY;
MObject PointSetMatchNode::aOutputTranslateZ;
MObject PointSetMatchNode::aOutputRotateX;
MObject PointSetMatchNode::aOutputRotateY;
MObject PointSetMatchNode::aOutputRotateZ;
MObject PointSetMatchNode::aOutputScaleX;
MObject PointSetMatchNode::aOutputScaleY;
MObject PointSetMatchNode::aOutputScaleZ;
MObject PointSetMatchNode::aOutputError;

PointSetMatchNode::PointSetMatchNode()
{}

PointSetMatchNode::~PointSetMatchNode()
{}

MStatus PointSetMatchNode::compute(const MPlug &plug, MDataBlock &dataBlock)
{
  MStatus status;

  if ((plug == aOutputTranslateX) ||
      (plug == aOutputTranslateY) ||
      (plug == aOutputTranslateZ) ||
      (plug == aOutputRotateX) ||
      (plug == aOutputRotateY) ||
      (plug == aOutputRotateZ) ||
      (plug == aOutputScaleX) ||
      (plug == aOutputScaleY) ||
      (plug == aOutputScaleZ) ||
      (plug == aOutputError))
  {
    std::vector<double> movPoints;
    std::vector<double> refPoints;
    std::vector<double> weights;

    //
    MArrayDataHandle inputPointsArray = dataBlock.inputArrayValue(aInputPoints);
    unsigned int inputPointsNum = inputPointsArray.elementCount();
    INFO("inputPointsNum: " << inputPointsNum);
    for ( unsigned int i = 0; i < inputPointsNum; i++ )
    {
      inputPointsArray.jumpToElement(i);
      MDataHandle inputPointElement = inputPointsArray.inputValue();

      // inputPointsArray.
      MDataHandle weightHandle = inputPointElement.child(aInputPointsWeight);
      double weight = weightHandle.asDouble();
      INFO("weight=" << weight);
      weights.push_back(weight);

      MDataHandle refHandle = inputPointElement.child(aInputPointsReferencePoint);
      double *refPoint = refHandle.asDouble3();
      // MVector refPoint = refHandle.asFloatVector();
      refPoints.push_back(refPoint[0]);
      refPoints.push_back(refPoint[1]);
      refPoints.push_back(refPoint[2]);

      MDataHandle movHandle = inputPointElement.child(aInputPointsMovePoint);
      double *movPoint = movHandle.asDouble3();
      // MVector movPoint = movHandle.asVector();
      INFO("movPointX=" << movPoint[0]);
      INFO("movPointY=" << movPoint[1]);
      INFO("movPointZ=" << movPoint[2]);
      movPoints.push_back(movPoint[0]);
      movPoints.push_back(movPoint[1]);
      movPoints.push_back(movPoint[2]);

      // inputPointsArray.next();
    }

    if (!inputPointsNum)
    {
      return MS::kUnknownParameter;
    }

    std::vector<unsigned int> lockedParams(9);
    lockedParams[0] = (unsigned int) dataBlock.inputValue(aLockedTranslateX).asBool();
    lockedParams[1] = (unsigned int) dataBlock.inputValue(aLockedTranslateY).asBool();
    lockedParams[2] = (unsigned int) dataBlock.inputValue(aLockedTranslateZ).asBool();
    lockedParams[3] = (unsigned int) dataBlock.inputValue(aLockedRotateX).asBool();
    lockedParams[4] = (unsigned int) dataBlock.inputValue(aLockedRotateY).asBool();
    lockedParams[5] = (unsigned int) dataBlock.inputValue(aLockedRotateZ).asBool();
    lockedParams[6] = (unsigned int) dataBlock.inputValue(aLockedScaleX).asBool();
    lockedParams[7] = (unsigned int) dataBlock.inputValue(aLockedScaleY).asBool();
    lockedParams[8] = (unsigned int) dataBlock.inputValue(aLockedScaleZ).asBool();

    std::vector<double> initialTransform(9);
    initialTransform[0] = (double) dataBlock.inputValue(aInitialTranslateX).asDouble();
    initialTransform[1] = (double) dataBlock.inputValue(aInitialTranslateY).asDouble();
    initialTransform[2] = (double) dataBlock.inputValue(aInitialTranslateZ).asDouble();
    initialTransform[3] = (double) dataBlock.inputValue(aInitialRotateX).asDouble();
    initialTransform[4] = (double) dataBlock.inputValue(aInitialRotateY).asDouble();
    initialTransform[5] = (double) dataBlock.inputValue(aInitialRotateZ).asDouble();
    initialTransform[6] = (double) dataBlock.inputValue(aInitialScaleX).asDouble();
    initialTransform[7] = (double) dataBlock.inputValue(aInitialScaleY).asDouble();
    initialTransform[8] = (double) dataBlock.inputValue(aInitialScaleZ).asDouble();

    int iterMax = dataBlock.inputValue(aIterations).asInt();
    bool uniformScale = dataBlock.inputValue(aUniformScale).asBool();
    // TODO: Rotation Order too.

    double outError = -1.0;
    std::vector<double> outTransform(9);
    bool ret = solveTransform(iterMax, uniformScale,
                              refPoints, movPoints,
                              weights, lockedParams, initialTransform,
                              outTransform,
                              outError);
    if (ret == false)
    {
      WRN("Solver returned false!");
      return MS::kUnknownParameter;
    }

    MDataHandle outTx = dataBlock.outputValue(aOutputTranslateX);
    MDataHandle outTy = dataBlock.outputValue(aOutputTranslateY);
    MDataHandle outTz = dataBlock.outputValue(aOutputTranslateZ);

    MDataHandle outRx = dataBlock.outputValue(aOutputRotateX);
    MDataHandle outRy = dataBlock.outputValue(aOutputRotateY);
    MDataHandle outRz = dataBlock.outputValue(aOutputRotateZ);

    MDataHandle outSx = dataBlock.outputValue(aOutputScaleX);
    MDataHandle outSy = dataBlock.outputValue(aOutputScaleY);
    MDataHandle outSz = dataBlock.outputValue(aOutputScaleZ);

    outTx.setDouble(outTransform[0]);
    outTy.setDouble(outTransform[1]);
    outTz.setDouble(outTransform[2]);
    outTx.setClean();
    outTy.setClean();
    outTz.setClean();

    outRx.setDouble(outTransform[3]);
    outRy.setDouble(outTransform[4]);
    outRz.setDouble(outTransform[5]);
    outRx.setClean();
    outRy.setClean();
    outRz.setClean();

    if (!uniformScale)
    {
      outSx.setDouble(outTransform[6]);
      outSy.setDouble(outTransform[7]);
      outSz.setDouble(outTransform[8]);
    }
    else
    {
      outSx.setDouble(outTransform[6]);
      outSy.setDouble(outTransform[6]);
      outSz.setDouble(outTransform[6]);
    }
    outSx.setClean();
    outSy.setClean();
    outSz.setClean();

    // TODO: Output errors for each point.
    MDataHandle outErr = dataBlock.outputValue(aOutputError);
    outErr.setDouble(outError);
    outErr.setClean();
  }
  else
  {
    return MS::kUnknownParameter;
  }

  return MS::kSuccess;
}

void *PointSetMatchNode::creator()
{
  return new PointSetMatchNode();
}

MStatus PointSetMatchNode::initialize()
{
  MFnNumericAttribute nAttr;
  MStatus stat;

  // Iterations
  aIterations = nAttr.create("iterations", "itr", MFnNumericData::kInt, 1000);
  stat = nAttr.setStorable(true);
  CHECK_MSTATUS_AND_RETURN_IT(stat);
  stat = addAttribute(aIterations);
  CHECK_MSTATUS_AND_RETURN_IT(stat);

  // Rotate Order
  aRotateOrder = nAttr.create("rotateOrder", "rord", MFnNumericData::kInt, 0);
  stat = nAttr.setStorable(true);
  CHECK_MSTATUS_AND_RETURN_IT(stat);
  stat = addAttribute(aRotateOrder);
  CHECK_MSTATUS_AND_RETURN_IT(stat);

  // Uniform Scale
  aUniformScale = nAttr.create("uniformScale", "unfscl", MFnNumericData::kBoolean, 1);
  stat = nAttr.setStorable(true);
  CHECK_MSTATUS_AND_RETURN_IT(stat);
  stat = addAttribute(aUniformScale);
  CHECK_MSTATUS_AND_RETURN_IT(stat);

  {
    aInputPointsReferencePointX = nAttr.create("refPointX", "refx",
                                               MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInputPointsReferencePointY = nAttr.create("refPointY", "refy",
                                               MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInputPointsReferencePointZ = nAttr.create("refPointZ", "refz",
                                               MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    MFnCompoundAttribute refPoint;
    aInputPointsReferencePoint = refPoint.create("refPoint", "ref", &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    stat = refPoint.addChild(aInputPointsReferencePointX);
    stat = refPoint.addChild(aInputPointsReferencePointY);
    stat = refPoint.addChild(aInputPointsReferencePointZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    stat = addAttribute(aInputPointsReferencePoint);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInputPointsMovePointX = nAttr.create("movPointX", "movx",
                                               MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInputPointsMovePointY = nAttr.create("movPointY", "movy",
                                               MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInputPointsMovePointZ = nAttr.create("movPointZ", "movz",
                                               MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    MFnCompoundAttribute movPoint;
    aInputPointsMovePoint = movPoint.create("movPoint", "mov", &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    stat = movPoint.addChild(aInputPointsMovePointX);
    stat = movPoint.addChild(aInputPointsMovePointY);
    stat = movPoint.addChild(aInputPointsMovePointZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    stat = addAttribute(aInputPointsMovePoint);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInputPointsWeight = nAttr.create("weight", "wgt",
                                      MFnNumericData::kDouble, 1.0);
    stat = nAttr.setMin(0.0);
    stat = nAttr.setMax(1.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    MFnCompoundAttribute inPoints;
    aInputPoints = inPoints.create("inputPoints", "inpts", &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    stat = inPoints.addChild(aInputPointsReferencePoint);
    stat = inPoints.addChild(aInputPointsMovePoint);
    stat = inPoints.addChild(aInputPointsWeight);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    stat = inPoints.setArray(true);
    stat = inPoints.setIndexMatters(false);
    stat = inPoints.setDisconnectBehavior(MFnAttribute::kDelete);
    stat = addAttribute(aInputPoints);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Initial Translation
  {
    aInitialTranslateX = nAttr.create("initialTranslateX", "itx",
                                      MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aInitialTranslateX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInitialTranslateY = nAttr.create("initialTranslateY", "ity",
                                      MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aInitialTranslateY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInitialTranslateZ = nAttr.create("initialTranslateZ", "itz",
                                      MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aInitialTranslateZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Initial Rotation
  {
    aInitialRotateX = nAttr.create("initialRotateX", "irx",
                                   MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aInitialRotateX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInitialRotateY = nAttr.create("initialRotateY", "iry",
                                   MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aInitialRotateY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInitialRotateZ = nAttr.create("initialRotateZ", "irz",
                                   MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aInitialRotateZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Initial Scale
  {
    aInitialScaleX = nAttr.create("initialScaleX", "isx",
                                  MFnNumericData::kDouble, 1.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aInitialScaleX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInitialScaleY = nAttr.create("initialScaleY", "isy",
                                  MFnNumericData::kDouble, 1.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aInitialScaleY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aInitialScaleZ = nAttr.create("initialScaleZ", "isz",
                                  MFnNumericData::kDouble, 1.0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aInitialScaleZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Locked Translation
  {
    aLockedTranslateX = nAttr.create("lockedTranslateX", "ltx",
                                     MFnNumericData::kBoolean, 0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aLockedTranslateX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aLockedTranslateY = nAttr.create("lockedTranslateY", "lty",
                                     MFnNumericData::kBoolean, 0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aLockedTranslateY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aLockedTranslateZ = nAttr.create("lockedTranslateZ", "ltz",
                                     MFnNumericData::kBoolean, 0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aLockedTranslateZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Locked Rotation
  {
    aLockedRotateX = nAttr.create("lockedRotateX", "lrx",
                                  MFnNumericData::kBoolean, 0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aLockedRotateX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aLockedRotateY = nAttr.create("lockedRotateY", "lry",
                                  MFnNumericData::kBoolean, 0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aLockedRotateY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aLockedRotateZ = nAttr.create("lockedRotateZ", "lrz",
                                  MFnNumericData::kBoolean, 0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aLockedRotateZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Locked Scale
  {
    aLockedScaleX = nAttr.create("lockedScaleX", "lsx",
                                 MFnNumericData::kBoolean, 0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aLockedScaleX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aLockedScaleY = nAttr.create("lockedScaleY", "lsy",
                                 MFnNumericData::kBoolean, 0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aLockedScaleY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aLockedScaleZ = nAttr.create("lockedScaleZ", "lsz",
                                 MFnNumericData::kBoolean, 0);
    stat = nAttr.setStorable(true);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aLockedScaleZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Output Translation
  {
    aOutputTranslateX = nAttr.create("outputTranslateX", "otx",
                                     MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    stat = nAttr.setWritable(false);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aOutputTranslateX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aOutputTranslateY = nAttr.create("outputTranslateY", "oty",
                                     MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    stat = nAttr.setWritable(false);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aOutputTranslateY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aOutputTranslateZ = nAttr.create("outputTranslateZ", "otz",
                                     MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    stat = nAttr.setWritable(false);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aOutputTranslateZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Output Rotate
  {
    aOutputRotateX = nAttr.create("outputRotateX", "orx",
                                  MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    stat = nAttr.setWritable(false);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aOutputRotateX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aOutputRotateY = nAttr.create("outputRotateY", "ory",
                                  MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    stat = nAttr.setWritable(false);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aOutputRotateY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aOutputRotateZ = nAttr.create("outputRotateZ", "orz",
                                  MFnNumericData::kDouble, 0.0);
    stat = nAttr.setStorable(true);
    stat = nAttr.setWritable(false);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aOutputRotateZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Output Scale
  {
    aOutputScaleX = nAttr.create("outputScaleX", "osx",
                                 MFnNumericData::kDouble, 1.0);
    stat = nAttr.setStorable(true);
    stat = nAttr.setWritable(false);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aOutputScaleX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aOutputScaleY = nAttr.create("outputScaleY", "osy",
                                 MFnNumericData::kDouble, 1.0);
    stat = nAttr.setStorable(true);
    stat = nAttr.setWritable(false);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aOutputScaleY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    aOutputScaleZ = nAttr.create("outputScaleZ", "osz",
                                 MFnNumericData::kDouble, 1.0);
    stat = nAttr.setStorable(true);
    stat = nAttr.setWritable(false);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(aOutputScaleZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  // Output Error
  aOutputError = nAttr.create("outputError", "oer",
                              MFnNumericData::kDouble, -1.0);
  stat = nAttr.setStorable(true);
  stat = nAttr.setWritable(false);
  CHECK_MSTATUS_AND_RETURN_IT(stat);
  stat = addAttribute(aOutputError);
  CHECK_MSTATUS_AND_RETURN_IT(stat);

  // Attribute Affects
  {
    stat = attributeAffects(aInputPoints, aOutputTranslateX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = attributeAffects(aInputPoints, aOutputTranslateY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = attributeAffects(aInputPoints, aOutputTranslateZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = attributeAffects(aInputPoints, aOutputRotateX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = attributeAffects(aInputPoints, aOutputRotateY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = attributeAffects(aInputPoints, aOutputRotateZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = attributeAffects(aInputPoints, aOutputScaleX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = attributeAffects(aInputPoints, aOutputScaleY);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = attributeAffects(aInputPoints, aOutputScaleZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = attributeAffects(aInputPoints, aOutputError);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    stat = attributeAffects(aInputPointsWeight, aOutputTranslateX);
    stat = attributeAffects(aInputPointsWeight, aOutputTranslateY);
    stat = attributeAffects(aInputPointsWeight, aOutputTranslateZ);
    stat = attributeAffects(aInputPointsWeight, aOutputRotateX);
    stat = attributeAffects(aInputPointsWeight, aOutputRotateY);
    stat = attributeAffects(aInputPointsWeight, aOutputRotateZ);
    stat = attributeAffects(aInputPointsWeight, aOutputScaleX);
    stat = attributeAffects(aInputPointsWeight, aOutputScaleY);
    stat = attributeAffects(aInputPointsWeight, aOutputScaleZ);
    stat = attributeAffects(aInputPointsWeight, aOutputError);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  {
    {
      stat = attributeAffects(aInitialTranslateX, aOutputTranslateX);
      stat = attributeAffects(aInitialTranslateX, aOutputTranslateY);
      stat = attributeAffects(aInitialTranslateX, aOutputTranslateZ);
      stat = attributeAffects(aInitialTranslateX, aOutputRotateX);
      stat = attributeAffects(aInitialTranslateX, aOutputRotateY);
      stat = attributeAffects(aInitialTranslateX, aOutputRotateZ);
      stat = attributeAffects(aInitialTranslateX, aOutputScaleX);
      stat = attributeAffects(aInitialTranslateX, aOutputScaleY);
      stat = attributeAffects(aInitialTranslateX, aOutputScaleZ);
      stat = attributeAffects(aInitialTranslateX, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aInitialTranslateY, aOutputTranslateX);
      stat = attributeAffects(aInitialTranslateY, aOutputTranslateY);
      stat = attributeAffects(aInitialTranslateY, aOutputTranslateZ);
      stat = attributeAffects(aInitialTranslateY, aOutputRotateX);
      stat = attributeAffects(aInitialTranslateY, aOutputRotateY);
      stat = attributeAffects(aInitialTranslateY, aOutputRotateZ);
      stat = attributeAffects(aInitialTranslateY, aOutputScaleX);
      stat = attributeAffects(aInitialTranslateY, aOutputScaleY);
      stat = attributeAffects(aInitialTranslateY, aOutputScaleZ);
      stat = attributeAffects(aInitialTranslateY, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aInitialTranslateZ, aOutputTranslateX);
      stat = attributeAffects(aInitialTranslateZ, aOutputTranslateY);
      stat = attributeAffects(aInitialTranslateZ, aOutputTranslateZ);
      stat = attributeAffects(aInitialTranslateZ, aOutputRotateX);
      stat = attributeAffects(aInitialTranslateZ, aOutputRotateY);
      stat = attributeAffects(aInitialTranslateZ, aOutputRotateZ);
      stat = attributeAffects(aInitialTranslateZ, aOutputScaleX);
      stat = attributeAffects(aInitialTranslateZ, aOutputScaleY);
      stat = attributeAffects(aInitialTranslateZ, aOutputScaleZ);
      stat = attributeAffects(aInitialTranslateZ, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);
    }

    {
      stat = attributeAffects(aInitialRotateX, aOutputTranslateX);
      stat = attributeAffects(aInitialRotateX, aOutputTranslateY);
      stat = attributeAffects(aInitialRotateX, aOutputTranslateZ);
      stat = attributeAffects(aInitialRotateX, aOutputRotateX);
      stat = attributeAffects(aInitialRotateX, aOutputRotateY);
      stat = attributeAffects(aInitialRotateX, aOutputRotateZ);
      stat = attributeAffects(aInitialRotateX, aOutputScaleX);
      stat = attributeAffects(aInitialRotateX, aOutputScaleY);
      stat = attributeAffects(aInitialRotateX, aOutputScaleZ);
      stat = attributeAffects(aInitialRotateX, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aInitialRotateY, aOutputTranslateX);
      stat = attributeAffects(aInitialRotateY, aOutputTranslateY);
      stat = attributeAffects(aInitialRotateY, aOutputTranslateZ);
      stat = attributeAffects(aInitialRotateY, aOutputRotateX);
      stat = attributeAffects(aInitialRotateY, aOutputRotateY);
      stat = attributeAffects(aInitialRotateY, aOutputRotateZ);
      stat = attributeAffects(aInitialRotateY, aOutputScaleX);
      stat = attributeAffects(aInitialRotateY, aOutputScaleY);
      stat = attributeAffects(aInitialRotateY, aOutputScaleZ);
      stat = attributeAffects(aInitialRotateY, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aInitialRotateZ, aOutputTranslateX);
      stat = attributeAffects(aInitialRotateZ, aOutputTranslateY);
      stat = attributeAffects(aInitialRotateZ, aOutputTranslateZ);
      stat = attributeAffects(aInitialRotateZ, aOutputRotateX);
      stat = attributeAffects(aInitialRotateZ, aOutputRotateY);
      stat = attributeAffects(aInitialRotateZ, aOutputRotateZ);
      stat = attributeAffects(aInitialRotateZ, aOutputScaleX);
      stat = attributeAffects(aInitialRotateZ, aOutputScaleY);
      stat = attributeAffects(aInitialRotateZ, aOutputScaleZ);
      stat = attributeAffects(aInitialRotateZ, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);
    }

    {
      stat = attributeAffects(aInitialScaleX, aOutputTranslateX);
      stat = attributeAffects(aInitialScaleX, aOutputTranslateY);
      stat = attributeAffects(aInitialScaleX, aOutputTranslateZ);
      stat = attributeAffects(aInitialScaleX, aOutputRotateX);
      stat = attributeAffects(aInitialScaleX, aOutputRotateY);
      stat = attributeAffects(aInitialScaleX, aOutputRotateZ);
      stat = attributeAffects(aInitialScaleX, aOutputScaleX);
      stat = attributeAffects(aInitialScaleX, aOutputScaleY);
      stat = attributeAffects(aInitialScaleX, aOutputScaleZ);
      stat = attributeAffects(aInitialScaleX, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aInitialScaleY, aOutputTranslateX);
      stat = attributeAffects(aInitialScaleY, aOutputTranslateY);
      stat = attributeAffects(aInitialScaleY, aOutputTranslateZ);
      stat = attributeAffects(aInitialScaleY, aOutputRotateX);
      stat = attributeAffects(aInitialScaleY, aOutputRotateY);
      stat = attributeAffects(aInitialScaleY, aOutputRotateZ);
      stat = attributeAffects(aInitialScaleY, aOutputScaleX);
      stat = attributeAffects(aInitialScaleY, aOutputScaleY);
      stat = attributeAffects(aInitialScaleY, aOutputScaleZ);
      stat = attributeAffects(aInitialScaleY, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aInitialScaleZ, aOutputTranslateX);
      stat = attributeAffects(aInitialScaleZ, aOutputTranslateY);
      stat = attributeAffects(aInitialScaleZ, aOutputTranslateZ);
      stat = attributeAffects(aInitialScaleZ, aOutputRotateX);
      stat = attributeAffects(aInitialScaleZ, aOutputRotateY);
      stat = attributeAffects(aInitialScaleZ, aOutputRotateZ);
      stat = attributeAffects(aInitialScaleZ, aOutputScaleX);
      stat = attributeAffects(aInitialScaleZ, aOutputScaleY);
      stat = attributeAffects(aInitialScaleZ, aOutputScaleZ);
      stat = attributeAffects(aInitialScaleZ, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);
    }
  }

  {
    {
      stat = attributeAffects(aLockedTranslateX, aOutputTranslateX);
      stat = attributeAffects(aLockedTranslateX, aOutputTranslateY);
      stat = attributeAffects(aLockedTranslateX, aOutputTranslateZ);
      stat = attributeAffects(aLockedTranslateX, aOutputRotateX);
      stat = attributeAffects(aLockedTranslateX, aOutputRotateY);
      stat = attributeAffects(aLockedTranslateX, aOutputRotateZ);
      stat = attributeAffects(aLockedTranslateX, aOutputScaleX);
      stat = attributeAffects(aLockedTranslateX, aOutputScaleY);
      stat = attributeAffects(aLockedTranslateX, aOutputScaleZ);
      stat = attributeAffects(aLockedTranslateX, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aLockedTranslateY, aOutputTranslateX);
      stat = attributeAffects(aLockedTranslateY, aOutputTranslateY);
      stat = attributeAffects(aLockedTranslateY, aOutputTranslateZ);
      stat = attributeAffects(aLockedTranslateY, aOutputRotateX);
      stat = attributeAffects(aLockedTranslateY, aOutputRotateY);
      stat = attributeAffects(aLockedTranslateY, aOutputRotateZ);
      stat = attributeAffects(aLockedTranslateY, aOutputScaleX);
      stat = attributeAffects(aLockedTranslateY, aOutputScaleY);
      stat = attributeAffects(aLockedTranslateY, aOutputScaleZ);
      stat = attributeAffects(aLockedTranslateY, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aLockedTranslateZ, aOutputTranslateX);
      stat = attributeAffects(aLockedTranslateZ, aOutputTranslateY);
      stat = attributeAffects(aLockedTranslateZ, aOutputTranslateZ);
      stat = attributeAffects(aLockedTranslateZ, aOutputRotateX);
      stat = attributeAffects(aLockedTranslateZ, aOutputRotateY);
      stat = attributeAffects(aLockedTranslateZ, aOutputRotateZ);
      stat = attributeAffects(aLockedTranslateZ, aOutputScaleX);
      stat = attributeAffects(aLockedTranslateZ, aOutputScaleY);
      stat = attributeAffects(aLockedTranslateZ, aOutputScaleZ);
      stat = attributeAffects(aLockedTranslateZ, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);
    }

    {
      stat = attributeAffects(aLockedRotateX, aOutputTranslateX);
      stat = attributeAffects(aLockedRotateX, aOutputTranslateY);
      stat = attributeAffects(aLockedRotateX, aOutputTranslateZ);
      stat = attributeAffects(aLockedRotateX, aOutputRotateX);
      stat = attributeAffects(aLockedRotateX, aOutputRotateY);
      stat = attributeAffects(aLockedRotateX, aOutputRotateZ);
      stat = attributeAffects(aLockedRotateX, aOutputScaleX);
      stat = attributeAffects(aLockedRotateX, aOutputScaleY);
      stat = attributeAffects(aLockedRotateX, aOutputScaleZ);
      stat = attributeAffects(aLockedRotateX, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aLockedRotateY, aOutputTranslateX);
      stat = attributeAffects(aLockedRotateY, aOutputTranslateY);
      stat = attributeAffects(aLockedRotateY, aOutputTranslateZ);
      stat = attributeAffects(aLockedRotateY, aOutputRotateX);
      stat = attributeAffects(aLockedRotateY, aOutputRotateY);
      stat = attributeAffects(aLockedRotateY, aOutputRotateZ);
      stat = attributeAffects(aLockedRotateY, aOutputScaleX);
      stat = attributeAffects(aLockedRotateY, aOutputScaleY);
      stat = attributeAffects(aLockedRotateY, aOutputScaleZ);
      stat = attributeAffects(aLockedRotateY, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aLockedRotateZ, aOutputTranslateX);
      stat = attributeAffects(aLockedRotateZ, aOutputTranslateY);
      stat = attributeAffects(aLockedRotateZ, aOutputTranslateZ);
      stat = attributeAffects(aLockedRotateZ, aOutputRotateX);
      stat = attributeAffects(aLockedRotateZ, aOutputRotateY);
      stat = attributeAffects(aLockedRotateZ, aOutputRotateZ);
      stat = attributeAffects(aLockedRotateZ, aOutputScaleX);
      stat = attributeAffects(aLockedRotateZ, aOutputScaleY);
      stat = attributeAffects(aLockedRotateZ, aOutputScaleZ);
      stat = attributeAffects(aLockedRotateZ, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);
    }

    {
      stat = attributeAffects(aLockedScaleX, aOutputTranslateX);
      stat = attributeAffects(aLockedScaleX, aOutputTranslateY);
      stat = attributeAffects(aLockedScaleX, aOutputTranslateZ);
      stat = attributeAffects(aLockedScaleX, aOutputRotateX);
      stat = attributeAffects(aLockedScaleX, aOutputRotateY);
      stat = attributeAffects(aLockedScaleX, aOutputRotateZ);
      stat = attributeAffects(aLockedScaleX, aOutputScaleX);
      stat = attributeAffects(aLockedScaleX, aOutputScaleY);
      stat = attributeAffects(aLockedScaleX, aOutputScaleZ);
      stat = attributeAffects(aLockedScaleX, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aLockedScaleY, aOutputTranslateX);
      stat = attributeAffects(aLockedScaleY, aOutputTranslateY);
      stat = attributeAffects(aLockedScaleY, aOutputTranslateZ);
      stat = attributeAffects(aLockedScaleY, aOutputRotateX);
      stat = attributeAffects(aLockedScaleY, aOutputRotateY);
      stat = attributeAffects(aLockedScaleY, aOutputRotateZ);
      stat = attributeAffects(aLockedScaleY, aOutputScaleX);
      stat = attributeAffects(aLockedScaleY, aOutputScaleY);
      stat = attributeAffects(aLockedScaleY, aOutputScaleZ);
      stat = attributeAffects(aLockedScaleY, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);

      stat = attributeAffects(aLockedScaleZ, aOutputTranslateX);
      stat = attributeAffects(aLockedScaleZ, aOutputTranslateY);
      stat = attributeAffects(aLockedScaleZ, aOutputTranslateZ);
      stat = attributeAffects(aLockedScaleZ, aOutputRotateX);
      stat = attributeAffects(aLockedScaleZ, aOutputRotateY);
      stat = attributeAffects(aLockedScaleZ, aOutputRotateZ);
      stat = attributeAffects(aLockedScaleZ, aOutputScaleX);
      stat = attributeAffects(aLockedScaleZ, aOutputScaleY);
      stat = attributeAffects(aLockedScaleZ, aOutputScaleZ);
      stat = attributeAffects(aLockedScaleZ, aOutputError);
      CHECK_MSTATUS_AND_RETURN_IT(stat);
    }
  }

  {
    stat = attributeAffects(aRotateOrder, aOutputTranslateX);
    stat = attributeAffects(aRotateOrder, aOutputTranslateY);
    stat = attributeAffects(aRotateOrder, aOutputTranslateZ);
    stat = attributeAffects(aRotateOrder, aOutputRotateX);
    stat = attributeAffects(aRotateOrder, aOutputRotateY);
    stat = attributeAffects(aRotateOrder, aOutputRotateZ);
    stat = attributeAffects(aRotateOrder, aOutputScaleX);
    stat = attributeAffects(aRotateOrder, aOutputScaleY);
    stat = attributeAffects(aRotateOrder, aOutputScaleZ);
    stat = attributeAffects(aRotateOrder, aOutputError);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    stat = attributeAffects(aIterations, aOutputTranslateX);
    stat = attributeAffects(aIterations, aOutputTranslateY);
    stat = attributeAffects(aIterations, aOutputTranslateZ);
    stat = attributeAffects(aIterations, aOutputRotateX);
    stat = attributeAffects(aIterations, aOutputRotateY);
    stat = attributeAffects(aIterations, aOutputRotateZ);
    stat = attributeAffects(aIterations, aOutputScaleX);
    stat = attributeAffects(aIterations, aOutputScaleY);
    stat = attributeAffects(aIterations, aOutputScaleZ);
    stat = attributeAffects(aIterations, aOutputError);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
  }

  return MS::kSuccess;
}
