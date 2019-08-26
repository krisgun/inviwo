/*********************************************************************
 *  Author  : Tino Weinkauf
 *  Init    : Monday, August 27, 2018 - 20:34:51
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <modules/labsetup/labsetupmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo {

class IVW_MODULE_LABSETUP_API LabSetupModule : public InviwoModule {
public:
    LabSetupModule(InviwoApplication* app);
    virtual ~LabSetupModule() = default;
};

}  // namespace inviwo
