/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/volume/volumeram.h>
#include <lablic/licprocessor.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
    , propKernel_("kernel", "Kernel Size", 150, 1, 500)
    , propFastLIC_("fastLIC", "FastLIC", true)
    , propMean("mean", "Mean", 128, 0.0, 255)
    , propDeviation("deviation", "Standard Deviation", 25.5, 0.0, 255)
// TODO: Register additional properties
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties 
    // TODO: Register additional properties
    addProperty(propKernel_);
    addProperty(propFastLIC_);
    addProperty(propMean);
    addProperty(propDeviation);
}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    vectorFieldDims_ = vol->getDimensions();

    auto tex = noiseTexIn_.getData();
    const RGBAImage texture = RGBAImage::createFromImage(tex);
    texDims_ = tex->getDimensions();

    double value = texture.readPixelGrayScale(size2_t(0, 0));
    LogProcessorInfo("rand val: " << value);

    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 0.0));

    // Hint: Output an image showing which pixels you have visited for debugging
    std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));
    
    auto pixelValSum {0};
    auto pixelValSumSquared {0};
	auto numOfNonBlack {0};
    auto kernelSize {propKernel_ / 2};
    //Check FastLIC toggle
    if (propFastLIC_) {
        #pragma omp parallel
        #pragma omp for
        for (auto j = 0; j < texDims_.y; ++j) {
            for (auto i = 0; i < texDims_.x; ++i) {
                if(visited[i][j]) continue;

                std::vector<ivec2> forwardPixels = Integrator::integratePoints(vectorField, vec2(i, j), 1, kernelSize, texDims_);
                std::vector<ivec2> backwardPixels = Integrator::integratePoints(vectorField, vec2(i, j), -1, kernelSize, texDims_);
                forwardPixels.insert(forwardPixels.end(), backwardPixels.begin(), backwardPixels.end());

                //apply convolution kernel, fora genom forwardPixels och summera ihop deras v�rden i texturen och dela p� antalet

                double sum = 0;
                for (int k = 0; k < forwardPixels.size(); ++k) {
                    sum += texture.readPixelGrayScale(size2_t(forwardPixels[k][0], forwardPixels[k][1]));
                }

                if (forwardPixels.size() > 0) {
                    sum /= forwardPixels.size();
                }

                for (auto k = 0; k < forwardPixels.size(); ++k) {
                    auto prevPixVal {licImage.readPixelGrayScale(size2_t(forwardPixels[k][0], forwardPixels[k][1]))};
                    licImage.setPixelGrayScale(size2_t(forwardPixels[k][0], forwardPixels[k][1]), sum);
                    visited[forwardPixels[k][0]][forwardPixels[k][1]] = 1;
					pixelValSum -= prevPixVal;
					pixelValSum += sum;
					pixelValSumSquared -= prevPixVal * prevPixVal;
					pixelValSumSquared += sum * sum;

					if (sum != 0 && prevPixVal == 0) {
						numOfNonBlack++;
					}
                }
            }
        }
    }
    else {
        #pragma omp parallel
        #pragma omp for
        for (auto j = 0; j < texDims_.y; j++) {
            for (auto i = 0; i < texDims_.x; i++) {
                std::vector<ivec2> forwardPixels = Integrator::integratePoints(vectorField, vec2(i, j), 1, kernelSize, texDims_);
                std::vector<ivec2> backwardPixels = Integrator::integratePoints(vectorField, vec2(i, j), -1, kernelSize, texDims_);
                forwardPixels.insert(forwardPixels.end(), backwardPixels.begin(), backwardPixels.end());

                //apply convolution kernel, fora genom forwardPixels och summera ihop deras v�rden i texturen och dela p� antalet

                double sum = 0;
                for (int k = 0; k < forwardPixels.size(); k++) {
                    sum += texture.readPixelGrayScale(size2_t(forwardPixels[k][0], forwardPixels[k][1]));
                }
                if (forwardPixels.size() > 0) {
                    sum /= forwardPixels.size();
                }
                licImage.setPixelGrayScale(size2_t(i, j), sum);
				pixelValSum += sum;
				pixelValSumSquared += sum * sum;
				if (sum != 0) {
					numOfNonBlack++;
				}
            }
        }
    }
    
	auto totalPix = texDims_.y*texDims_.x;
	auto mean = pixelValSum / numOfNonBlack;
	auto stdDev = sqrt((pixelValSumSquared - (numOfNonBlack * (mean * mean))) / (numOfNonBlack- 1));
    
	auto f = propDeviation / stdDev;

	for (auto j = 0; j < texDims_.y; j++) {
		for (auto i = 0; i < texDims_.x; i++) {
			auto oldPix = licImage.readPixelGrayScale(size2_t(i,j));
			licImage.setPixelGrayScale(size2_t(i, j), propMean + f*(oldPix - mean));
		}
	}

    licOut_.setData(outImage);
}

}  // namespace inviwo
