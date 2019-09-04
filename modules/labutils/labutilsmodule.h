#pragma once

#include <labutils/labutilsmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo {

class IVW_MODULE_LABUTILS_API LabUtilsModule : public InviwoModule {
public:
    LabUtilsModule(InviwoApplication* app);
    virtual ~LabUtilsModule() = default;
};

}  // namespace inviwo
