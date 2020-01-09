/*********************************************************************
 *  Author  : Tino Weinkauf
 *  Init    : Thursday, February 01, 2018 - 20:34:51
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <modules/labcolor/labcolormoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo {

class IVW_MODULE_LABCOLOR_API LabColorModule : public InviwoModule {
public:
    LabColorModule(InviwoApplication* app);
    virtual ~LabColorModule() = default;
};

}  // namespace inviwo
