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

#include <modules/labsubdivision/labsubdivisionmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo {

class IVW_MODULE_LABSUBDIVISION_API LabSubdivisionModule : public InviwoModule {
public:
    LabSubdivisionModule(InviwoApplication* app);
    virtual ~LabSubdivisionModule() = default;
};

}  // namespace inviwo
