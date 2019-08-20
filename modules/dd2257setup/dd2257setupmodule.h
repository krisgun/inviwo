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

#include <modules/dd2257setup/dd2257setupmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo
{

class IVW_MODULE_DD2257SETUP_API DD2257SetupModule : public InviwoModule
{
public:
    DD2257SetupModule(InviwoApplication* app);
    virtual ~DD2257SetupModule() = default;
};

} // namespace
