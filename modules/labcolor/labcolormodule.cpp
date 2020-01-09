/*********************************************************************
 *  Author  : Tino Weinkauf
 *  Init    : Thursday, February 01, 2018 - 20:34:51
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labcolor/labcolormodule.h>
#include <modules/labcolor/colormixing.h>
#include <modules/labcolor/colorinterpolation.h>

namespace inviwo {

using namespace kth;

LabColorModule::LabColorModule(InviwoApplication* app) : InviwoModule(app, "LabColor") {
    // Add a directory to the search path of the Shadermanager
    // ShaderManager::getPtr()->addShaderSearchPath(getPath(ModulePath::GLSL));

    // Register objects that can be shared with the rest of inviwo here:

    // Processors
    registerProcessor<ColorMixing>();
    registerProcessor<ColorInterpolation>();
}

}  // namespace inviwo
