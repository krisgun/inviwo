/*********************************************************************
 *  Author  : Tino Weinkauf
 *  Init    : Thursday, February 01, 2018 - 20:34:51
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labsubdivision/labsubdivisionmodule.h>
#include <modules/labsubdivision/chaikin.h>

namespace inviwo {

using namespace kth;

LabSubdivisionModule::LabSubdivisionModule(InviwoApplication* app)
    : InviwoModule(app, "LabSubdivision") {
    // Add a directory to the search path of the Shadermanager
    // ShaderManager::getPtr()->addShaderSearchPath(getPath(ModulePath::GLSL));

    // Register objects that can be shared with the rest of inviwo here:

    // Processors
    registerProcessor<Chaikin>();
}

}  // namespace inviwo
