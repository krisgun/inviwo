/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 30, 2017 - 13:12:27
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labraycasting/labraycastingmodule.h>
#include <labraycasting/compositor.h>

namespace inviwo
{

LabRaycastingModule::LabRaycastingModule(InviwoApplication* app) : InviwoModule(app, "DD2320Lab3")
    
{
	registerProcessor<Compositor>();
}

} // namespace
