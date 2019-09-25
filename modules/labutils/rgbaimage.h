#pragma once

#include <inviwo/core/datastructures/image/image.h>
#include <inviwo/core/datastructures/image/layerram.h>
#include <inviwo/core/datastructures/image/layerramprecision.h>
#include <labutils/labutilsmoduledefine.h>

namespace inviwo {

class IVW_MODULE_LABUTILS_API RGBAImage {

    // Typedefs
public:
    /** Integer pixel indices. */
    typedef size2_t IndexType;

    /** Fractional pixel indices */
    typedef dvec2 PositionType;

public:
    /** \brief Create an image from a given inviwo image. */
    static const RGBAImage createFromImage(std::shared_ptr<const inviwo::Image> image);

    RGBAImage(std::shared_ptr<inviwo::Image> image);

    /**
     * \brief Create new empty image with RGBA values in [0,255].
     *
     * @param size Number of pixels.
     */
    RGBAImage(const IndexType& size);

    /** Copy constructor */
    RGBAImage(const RGBAImage& other);

    /** Copy assignment */
    RGBAImage& operator=(const RGBAImage& other);

    /** Default destructor */
    ~RGBAImage() {
        if (ownsData_) delete data_;
    }

    void setPixel(IndexType idx, dvec4 color);

    void setPixelGrayScale(IndexType idx, double value);

    dvec4 sample(PositionType fracIdx) const;

    double sampleGrayScale(PositionType fracIdx) const;

    dvec4 readPixel(IndexType idx) const;

    double readPixelGrayScale(IndexType idx) const;

    // Data
protected:
    RGBAImage(std::shared_ptr<const inviwo::Image> image);

    LayerRAM* data_;
    bool ownsData_ = false;
    IndexType size_;
};

}  // namespace inviwo
