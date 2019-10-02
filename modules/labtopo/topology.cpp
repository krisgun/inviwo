/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 **********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/integrator.h>
#include <labutils/scalarvectorfield.h>
#include <labtopo/topology.h>
#include <labtopo/utils/gradients.h>

namespace inviwo {

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 1, 0, 1),    // Saddle
    vec4(0, 0, 1, 1),    // AttractingNode
    vec4(1, 0, 0, 1),    // RepellingNode
    vec4(0.5, 0, 1, 1),  // AttractingFocus
    vec4(1, 0.5, 0, 1),  // RepellingFocus
    vec4(0, 1, 0, 1)     // Center
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",    // Class identifier
    "Vector Field Topology",  // Display name
    "KTH Lab",                // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const { return processorInfo_; }

Topology::Topology()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);

    // TODO: Register additional properties
    // addProperty(propertyName);
}

void Topology::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);

    // Initialize mesh, vertices and index buffers for seperatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    // Add a bounding box to the mesh (The 2D mesh renderer will automatically adapt to this)
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    auto indexBufferBBox = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    drawLineSegment(BBoxMin, vec2(BBoxMin[0], BBoxMax[1]), vec4(0, 0, 0, 1), indexBufferBBox.get(),
                    vertices);
    drawLineSegment(vec2(BBoxMin[0], BBoxMax[1]), BBoxMax, vec4(0, 0, 0, 1), indexBufferBBox.get(),
                    vertices);
    drawLineSegment(BBoxMax, vec2(BBoxMax[0], BBoxMin[1]), vec4(0, 0, 0, 1), indexBufferBBox.get(),
                    vertices);
    drawLineSegment(vec2(BBoxMax[0], BBoxMin[1]), BBoxMin, vec4(0, 0, 0, 1), indexBufferBBox.get(),
                    vertices);

    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type adjacency.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    size2_t dims = vectorField.getNumVerticesPerDim();

    // Looping through all values in the vector field.
    for (int j = 0; j < dims[1]; ++j) {
        for (int i = 0; i < dims[0]; ++i) {
            dvec2 vectorValue = vectorField.getValueAtVertex(size2_t(i, j));
            dvec2 pos = vectorField.getPositionAtVertex(size2_t(i, j));
            // Computing the jacobian at a position
            dmat2 jacobian = vectorField.derive(pos);
            // Doing the eigen analysis
            auto eigenResult = util::eigenAnalysis(jacobian);
            // The result of the eigen analysis has attributed eigenvaluesRe eigenvaluesIm and
            // eigenvectors
        }
    }

    // Accessing the colors
    vec4 colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                               IndexBufferRAM* indexBuffer,
                               std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

}  // namespace inviwo
