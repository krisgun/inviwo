/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 30, 2017 - 13:13:06
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <labraycasting/labraycastingmoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/ports/imageport.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <inviwo/core/datastructures/image/layerram.h>
#include <inviwo/core/properties/optionproperty.h>
#include <inviwo/core/properties/transferfunctionproperty.h>

#include <inviwo/core/properties/ordinalproperty.h>


namespace inviwo
{

/** \docpage{org.inviwo.Compositor, Compositor}
    ![](org.inviwo.Compositor.png?classIdentifier=org.inviwo.Compositor)

    This processor provides different compositing schemes for raycasting in a 3D volume.
    An orthographic projection of the slices of the volume in the XY plane is used, so that
    the result is an image displaying the content of the colume in z-direction.

    ### Inports
        * __inVol__ Input Volume.

    ### Outports
        * __image__ Output Image.
*/


/** \class Compositor
    \brief Processor for Raycasting in a 3D volume

    @author Himangshu Saikia
*/
class IVW_MODULE_LABRAYCASTING_API Compositor : public Processor
{
    //Friends
    //Types
public:

    //Construction / Deconstruction
public:
    Compositor();
    virtual ~Compositor() = default;

    //Methods
public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

protected:
    ///Our main computation function
    virtual void process() override;

    vec4 getAverage(const VolumeRAM* vr, const size_t x, const size_t y);
    vec4 getMaximum(const VolumeRAM* vr, const size_t x, const size_t y);
    vec4 getFirstHit(const VolumeRAM* vr, const size_t x, const size_t y, double isovalue);
    vec4 accumulateFrontToBack(const VolumeRAM* vr, const size_t x, const size_t y);
    vec4 accumulateBackToFront(const VolumeRAM* vr, const size_t x, const size_t y);
    void renderImage(const VolumeRAM* vr, LayerRAM* lr);

    //Ports
public:
    // input volume
    VolumeInport volIn_;
    // output image
    ImageOutport image_;


    OptionPropertyInt method_;
    TransferFunctionProperty tf_;

    // iso value for the first hit
    FloatProperty firstHitIsoValue_;

    //Properties
public:

    //Attributes
private:
    size3_t dims_;
};

} // namespace
