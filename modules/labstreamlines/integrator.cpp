/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/integrator.h>

namespace inviwo {

// TODO: Implement a single integration step here

dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position, const float stepSize) {
	//x_i+1 = x_i + stepSize * v(x_i)
    dvec2 v_x_i = vectorField.interpolate(position);
    return position + dvec2(stepSize*v_x_i[0], stepSize*v_x_i[1]);
    //     Access the vector field with vectorField.interpolate(...)
}


dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position, const float stepSize) {
    dvec2 v1 = vectorField.interpolate(position);
    dvec2 v2 = vectorField.interpolate(position + dvec2((stepSize/2) * v1[0], (stepSize/2) * v1[1]));
    dvec2 v3 = vectorField.interpolate(position + dvec2((stepSize/2) * v2[0], (stepSize/2) * v2[1]));
    dvec2 v4 = vectorField.interpolate(position + dvec2(stepSize * v3[0], stepSize * v3[1]));
	
	double xCord = stepSize*((v1[0] + 2 * v2[0] + 2 * v3[0] + v4[0])/6);
	double yCord = stepSize*((v1[1] + 2 * v2[1] + 2 * v3[1] + v4[1])/6);
	
	return position + dvec2(xCord, yCord);
}
//translate 
/*
i både x och y-led
mappa 0 till bboxmin
mappa bboxmax till dims

vi måste rasteriza
				 
pixRatio = startPixel[x] / (dims[x])  
startPoint = pixRatio * (bboxmax[x] - bboxmin[x])
/   


*/

std::vector<ivec2> Integrator::integratePoints(const VectorField2& vectorField, vec2 startPixel, int direction, int steps, size2_t dims) {
    
    ivec2 dimensions = vectorField.getNumVerticesPerDim();
    dvec2 BBoxMin_ = vectorField.getBBoxMin();
    dvec2 BBoxMax_ = vectorField.getBBoxMax();
    
	double stepSize = direction * std::min((BBoxMax_[0] - BBoxMin_[0]) / dims[0], (BBoxMax_[1] - BBoxMin_[1]) / dims[1]);

	//std::cout << stepSize;

	std::vector<ivec2> pixels;
	std::vector<dvec2> ratio;

	vec2 pixRatio = {1 / dims[0], 1 / dims[1]};
	vec2 startPoint = {pixRatio[0] * startPixel[0] * (BBoxMax_[0] - BBoxMin_[0]) + BBoxMin_[0], pixRatio[1] * startPixel[1] * (BBoxMax_[1] - BBoxMin_[1]) + BBoxMin_[1]};

    dvec2 oldPoint = startPoint;
    //Add start point to vector
    //points.push_back(oldPoint);
	pixels.push_back({oldPoint[0]/pixRatio[0], oldPoint[1]/pixRatio[1]});

    for (int i = 0; i < steps; i++) {

        dvec2 fieldVec = vectorField.interpolate(startPoint);
        //Stop when encountering zeros in vector field 
        bool isVecFieldZero = (fieldVec[0] > -0.001 && fieldVec[0] < 0.001) && (fieldVec[1] > -0.001 && fieldVec[1] < 0.001);
        if (isVecFieldZero) break;

        //Get new point from Rung Kut
        startPoint = Integrator::RK4(vectorField, startPoint, stepSize);
		     
        //Stop at boundary of bbox
        if (!vectorField.isInside(startPoint)) {
            startPoint = vectorField.clampPositionToBBox(startPoint);
            //points.push_back(startPoint);
			pixels.push_back({ startPoint[0] / pixRatio[0], startPoint[1] / pixRatio[1] });
            break;
        }
        //Add new point to vector
        //points.push_back(startPoint);
		pixels.push_back({ startPoint[0] / pixRatio[0], startPoint[1] / pixRatio[1] });
        oldPoint = startPoint;
    }

    return pixels;
}

void Integrator::drawPoint(const dvec2& p, const vec4& color, IndexBufferRAM* indexBuffer,
                           std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(p[0], p[1], 0), vec3(0, 0, 1), vec3(p[0], p[1], 0), color});
}

// Alias for draw point
void Integrator::drawNextPointInPolyline(const dvec2& p, const vec4& color,
                                         IndexBufferRAM* indexBuffer,
                                         std::vector<BasicMesh::Vertex>& vertices) {
    Integrator::drawPoint(p, color, indexBuffer, vertices);
}

void Integrator::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                                 IndexBufferRAM* indexBuffer,
                                 std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

}  // namespace inviwo
