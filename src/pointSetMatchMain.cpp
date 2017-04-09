
//
#include <pointSetMatchNode.h>
// #include <pointSetMatchCmd.h>

#include <maya/MGlobal.h>
#include <maya/MFnPlugin.h>

MStatus initializePlugin(MObject obj)
{
  MStatus status;
  MFnPlugin plugin(obj, PLUGIN_COMPANY, "1.0", "Any");

  status = plugin.registerNode("pointSetMatch",
                               PointSetMatchNode::id,
                               PointSetMatchNode::creator,
                               PointSetMatchNode::initialize);
  CHECK_MSTATUS_AND_RETURN_IT(status);

//  status = plugin.registerCommand("pointSetMatch",
//                                  PointSetMatchCmd::creator);
//  CHECK_MSTATUS_AND_RETURN_IT(status);

  return status;
}

MStatus uninitializePlugin(MObject obj)
{
  MStatus status;
  MFnPlugin plugin(obj);

  status = plugin.deregisterNode(PointSetMatchNode::id);
  CHECK_MSTATUS_AND_RETURN_IT(status);

//  status = plugin.deregisterCommand("pointSetMatch");
//  CHECK_MSTATUS_AND_RETURN_IT(status);

  return status;
}
