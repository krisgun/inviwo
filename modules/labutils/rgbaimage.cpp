#pragma once

#include <labutils/rgbaimage.h>

namespace inviwo {

const RGBAImage RGBAImage::createFromImage(std::shared_ptr<const inviwo::Image> image) {
    return RGBAImage(image);
}

RGBAImage::RGBAImage(std::shared_ptr<const inviwo::Image> image) : ownsData_(false) {
    size_ = image->getDimensions();
    IVW_ASSERT(image.get(), "No valid volume.");
    IVW_ASSERT(image->getColorLayer()->getRepresentation<LayerRAM>(),
               "No valid Layer RAM representation.");
    data_ = const_cast<LayerRAM*>(image->getColorLayer()->getRepresentation<LayerRAM>());
}

RGBAImage::RGBAImage(std::shared_ptr<inviwo::Image> image) : ownsData_(false) {
    size_ = image->getDimensions();
    IVW_ASSERT(image.get(), "No valid volume.");
    IVW_ASSERT(image->getColorLayer()->getEditableRepresentation<LayerRAM>(),
               "No valid Layer RAM representation.");
    data_ = image->getColorLayer()->getEditableRepresentation<LayerRAM>();
}

RGBAImage::RGBAImage(const IndexType& size) : ownsData_(true), size_(size) {
    data_ = new LayerRAMPrecision<glm::u8vec4>(size_);
}

RGBAImage::RGBAImage(const RGBAImage& other)
    : data_(other.data_->clone()), ownsData_(true), size_(other.size_) {}

RGBAImage& RGBAImage::operator=(const RGBAImage& other) {
    data_ = other.data_->clone();
    ownsData_ = true;
    size_ = other.size_;

    return *this;
}

void RGBAImage::setPixel(IndexType idx, dvec4 color) { data_->setFromDVec4(idx, color); }

void RGBAImage::setPixelGrayScale(IndexType idx, double value) {
    setPixel(idx, vec4(value, value, value, 255));
}

dvec4 RGBAImage::sample(PositionType fracIdx) const {

    IVW_ASSERT(fracIdx[0] < size_[0] && fracIdx[1] < size_[1],
               "RGBA image accessed outside of its bounds.");

    dvec2 locPix(fracIdx[0] - int(fracIdx[0]), fracIdx[1] - int(fracIdx[1]));

    auto f00 = data_->getAsDVec4(size2_t(fracIdx[0], fracIdx[1]));
    auto f10 = fracIdx[0] + 1 < size_[0] ? data_->getAsDVec4(size2_t(fracIdx[0] + 1, fracIdx[1]))
                                         : dvec4(0.0);
    auto f01 = fracIdx[1] + 1 < size_[1] ? data_->getAsDVec4(size2_t(fracIdx[0], fracIdx[1] + 1))
                                         : dvec4(0.0);
    auto f11 = fracIdx[0] + 1 < size_[0] && fracIdx[1] + 1 < size_[1]
                   ? data_->getAsDVec4(size2_t(fracIdx[0] + 1, fracIdx[1] + 1))
                   : dvec4(0.0);

    return (1 - locPix[1]) * ((1 - locPix[0]) * f00 + locPix[0] * f10) +
           locPix[1] * ((1 - locPix[0]) * f01 + locPix[0] * f11);
}

double RGBAImage::sampleGrayScale(PositionType fracIdx) const {
    dvec4 color = sample(fracIdx);
    return (color[0] + color[1] + color[2]) / 3;
}

dvec4 RGBAImage::readPixel(IndexType idx) const { return data_->getAsDVec4(idx); }

double RGBAImage::readPixelGrayScale(IndexType idx) const {
    dvec4 color = readPixel(idx);
    return (color[0] + color[1] + color[2]) / 3;
}

}  // namespace inviwo
