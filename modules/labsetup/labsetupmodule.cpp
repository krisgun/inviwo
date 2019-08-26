/*********************************************************************
 *  Author  : Tino Weinkauf
 *  Init    : Monday, August 27, 2018 - 20:34:51
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labsetup/labsetupmodule.h>
#include <modules/labsetup/connectpoints.h>

namespace inviwo {

using namespace kth;

LabSetupModule::LabSetupModule(InviwoApplication* app) : InviwoModule(app, "LabSetup") {
    // Add a directory to the search path of the Shadermanager
    // ShaderManager::getPtr()->addShaderSearchPath(getPath(ModulePath::GLSL));

    // Register objects that can be shared with the rest of inviwo here:

    // Processors
    registerProcessor<ConnectPoints>();

    // Properties
    // registerProperty<LabSetupProperty>();

    // Readers and writes
    // registerDataReader(util::make_unique<LabSetupReader>());
    // registerDataWriter(util::make_unique<LabSetupWriter>());

    // Data converters
    // registerRepresentationConverter(util::make_unique<LabSetupDisk2RAMConverter>());

    // Ports
    // registerPort<LabSetupOutport>();
    // registerPort<LabSetupInport>();

    // PropertyWidgets
    // registerPropertyWidget<LabSetupPropertyWidget, LabSetupProperty>("Default");

    // Dialogs
    // registerDialog<LabSetupDialog>(LabSetupOutport);

    // Other varius things
    // registerCapabilities(util::make_unique<LabSetupCapabilities>());
    // registerSettings(util::make_unique<LabSetupSettings>());
    // registerMetaData(util::make_unique<LabSetupMetaData>());
    // registerPortInspector("LabSetupOutport", "path/workspace.inv");
    // registerProcessorWidget(std::string processorClassName, std::unique_ptr<ProcessorWidget>
    // processorWidget); registerDrawer(util::make_unique_ptr<LabSetupDrawer>());
}

}  // namespace inviwo
