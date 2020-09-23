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
