/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 30, 2017 - 13:12:27
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <labraycasting/labraycastingmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo
{

class IVW_MODULE_LABRAYCASTING_API LabRaycastingModule : public InviwoModule
{
public:
    LabRaycastingModule(InviwoApplication* app);
    virtual ~LabRaycastingModule() = default;
};

} // namespace
