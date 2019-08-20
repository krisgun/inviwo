/*********************************************************************
 *  Author  : Tino Weinkauf
 *  Init    : Monday, August 27, 2018 - 20:34:51
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <dd2257setup/dd2257setupmodule.h>
#include <modules/dd2257setup/connectpoints.h>

namespace inviwo
{

using namespace kth;

DD2257SetupModule::DD2257SetupModule(InviwoApplication* app) : InviwoModule(app, "DD2257Setup")
{
    // Add a directory to the search path of the Shadermanager
    // ShaderManager::getPtr()->addShaderSearchPath(getPath(ModulePath::GLSL));

    // Register objects that can be shared with the rest of inviwo here:
    
    // Processors
    registerProcessor<ConnectPoints>();
    
    // Properties
    // registerProperty<DD2257SetupProperty>();
    
    // Readers and writes
    // registerDataReader(util::make_unique<DD2257SetupReader>());
    // registerDataWriter(util::make_unique<DD2257SetupWriter>());
    
    // Data converters
    // registerRepresentationConverter(util::make_unique<DD2257SetupDisk2RAMConverter>());

    // Ports
    // registerPort<DD2257SetupOutport>();
    // registerPort<DD2257SetupInport>();

    // PropertyWidgets
    // registerPropertyWidget<DD2257SetupPropertyWidget, DD2257SetupProperty>("Default");
    
    // Dialogs
    // registerDialog<DD2257SetupDialog>(DD2257SetupOutport);
    
    // Other varius things
    // registerCapabilities(util::make_unique<DD2257SetupCapabilities>());
    // registerSettings(util::make_unique<DD2257SetupSettings>());
    // registerMetaData(util::make_unique<DD2257SetupMetaData>());   
    // registerPortInspector("DD2257SetupOutport", "path/workspace.inv");
    // registerProcessorWidget(std::string processorClassName, std::unique_ptr<ProcessorWidget> processorWidget);
    // registerDrawer(util::make_unique_ptr<DD2257SetupDrawer>());  
}

} // namespace
