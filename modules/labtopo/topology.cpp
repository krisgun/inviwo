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
    vec4(1, 1, 0, 1),    // Saddle - Yellow
    vec4(1, 0, 0, 1),    // AttractingNode - Red
    vec4(0, 0, 1, 1),    // RepellingNode - Blue
    vec4(0.5, 0, 1, 1),  // AttractingFocus - Purple
    vec4(1, 0.5, 0, 1),  // RepellingFocus - Orange
    vec4(0, 1, 0, 1)     // Center - Green
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
    , meshBBoxOut("meshBBoxOut")
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

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

    // Add a bounding box to the mesh
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    // Initialize mesh, vertices and index buffers for seperatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type strip.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines,
    // ConnectivityType::Strip);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    std::vector<dvec2> criticalPoints {};
    //Set start corners to be the bounding box
    std::vector<dvec2> startCorners = {BBoxMin, vec2{BBoxMax[0], BBoxMin[1]}, BBoxMax, vec2{BBoxMin[0], BBoxMax[1]}};
    double epsilon = 0.01;
    //Extract all critical points
    findCriticalPoints(vectorField, criticalPoints, epsilon, startCorners);
    

    for (auto const& point : criticalPoints) {
        //LogProcessorInfo("Critical points: " << point);
        Integrator::drawPoint(point, vec4(1, 1, 1, 1), indexBufferPoints.get(), vertices);
    }

    size2_t dims = vectorField.getNumVerticesPerDim();
    // Looping through all values in the vector field.
    for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t i = 0; i < dims[0]; ++i) {
            dvec2 vectorValue = vectorField.getValueAtVertex(size2_t(i, j));
            dvec2 pos = vectorField.getPositionAtVertex(size2_t(i, j));
            // Computing the jacobian at a position
            dmat2 jacobian = vectorField.derive(pos);
            // Doing the eigen analysis
            auto eigenResult = util::eigenAnalysis(jacobian);
            // The result of the eigen analysis has attributed eigenvaluesRe eigenvaluesIm and
            // eigenvectors
            if ((i == 0 || i == dims[0]-1) && (j == 0 || j == dims[1]-1)) {
                LogProcessorInfo("index: ("<< i << ", " << j << "): pos: " << pos << ", val: " << vectorValue)
            }
        }
    }

    // Other helpful functions
    // dvec2 pos = vectorField.getPositionAtVertex(size2_t(i, j));
    // Computing the jacobian at a position
    // dmat2 jacobian = vectorField.derive(pos);
    // Doing the eigen analysis
    // auto eigenResult = util::eigenAnalysis(jacobian);
    // The result of the eigen analysis has attributed eigenvaluesRe eigenvaluesIm and
    // eigenvectors

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

/*
    3 *-------* 2
      |       |
      |       |
      |       |
    0 *-------* 1
*/

void Topology::findCriticalPoints(const VectorField2& vectorField, std::vector<dvec2>& criticalPoints, double epsilon, std::vector<dvec2>& corners) {
    std::vector<dvec2> cornerVectorValues{};
    
    double minPointDist = std::min((corners[1][0] - corners[0][0]), (corners[2][1] - corners[1][1]));

    for (int i = 0; i < corners.size(); ++i) {
        dvec2 fieldVec = vectorField.interpolate(corners[i]);
        cornerVectorValues.push_back(fieldVec);

        //Check if vector field is zero at the current corner
        bool isVecFieldZero = (fieldVec[0] > -0.001 && fieldVec[0] < 0.001) && (fieldVec[1] > -0.001 && fieldVec[1] < 0.001);
        if (isVecFieldZero) {
            if (std::find(criticalPoints.begin(), criticalPoints.end(), corners[i]) == criticalPoints.end()) {
                criticalPoints.push_back(corners[i]);
            }
            return;
        }
    }

    if (minPointDist < epsilon) {
        if (isCornersDifferentSigns(cornerVectorValues)) {
            dvec2 elem{ corners[0][0] + (corners[1][0] - corners[0][0]) / 2 , corners[0][1] + (corners[3][1] - corners[0][1]) / 2 };
            if (std::find(criticalPoints.begin(), criticalPoints.end(), (elem)) == criticalPoints.end()) {
                criticalPoints.push_back(elem);
            }
        }
        return;
    }

    //Get boundary points for inner quadrants
    dvec2 point01 {corners[0][0] + (corners[1][0] - corners[0][0])/2, corners[0][1]};
    dvec2 point12 {corners[1][0], corners[1][1] + (corners[2][1] - corners[1][1])/2};
    dvec2 point23 {corners[3][0] + (corners[2][0] - corners[3][0])/2, corners[2][1]};
    dvec2 point03 {corners[0][0], corners[1][1] + (corners[2][1] - corners[1][1]) / 2};
    dvec2 pointMid {corners[0][0] + (corners[1][0] - corners[0][0]) / 2, corners[1][1] + (corners[2][1] - corners[1][1]) / 2 };

    //Create new 4 sub-quadrants of the original one
    std::vector<dvec2> q0Corners {corners[0], point01, pointMid, point03};
    std::vector<dvec2> q1Corners {point01, corners[1], point12, pointMid};
    std::vector<dvec2> q2Corners {pointMid, point12, corners[2], point23};
    std::vector<dvec2> q3Corners {point03, pointMid, point23, corners[3]};

    //Recurse over the sub-quadrants
    findCriticalPoints(vectorField, criticalPoints, epsilon, q0Corners);
    findCriticalPoints(vectorField, criticalPoints, epsilon, q1Corners);
    findCriticalPoints(vectorField, criticalPoints, epsilon, q2Corners);
    findCriticalPoints(vectorField, criticalPoints, epsilon, q3Corners);


    //if (std::min(size[0],size[1]) < epsilon) {
        //nollkoll => om det finns l�gg till mitten av rutan
    //}
    //f�r in criticalpointslista, fyll p� om den hittar inom en ruta som �r 
    //kolla points i h�rnen -> om det inte kan finnas nolla returna, annars g�r rutan mindre i fyra delar
    //n�r vi kommer till threshhold och det kan finnas en nolla l�gg till punkten i mitten
}

int Topology::sign(double number) {
    if (number > 0) return 1;
    if (number < 0) return -0;
    return 0;
}

bool Topology::isCornersDifferentSigns(std::vector<dvec2>& corners) {
    int x0 = Topology::sign(corners[0][0]);
    int x1 = Topology::sign(corners[1][0]);
    int x2 = Topology::sign(corners[2][0]);
    int x3 = Topology::sign(corners[3][0]);

    bool differentSignX = !(x0 == x1 && x0 == x2 && x0 == x3);

    int y0 = Topology::sign(corners[0][1]);
    int y1 = Topology::sign(corners[1][1]);
    int y2 = Topology::sign(corners[2][1]);
    int y3 = Topology::sign(corners[3][1]);

    bool differentSignY = !(y0 == y1 && y0 == y2 && y0 == y3);

    return differentSignX && differentSignY;
}

bool Topology::isZeroWithinBox(std::vector<dvec2>& corners) {
    bool differentSignsX = ( 
        (   
            //Different X-signs between corner 0 and 1
            (corners[0][0] > 0 && corners[1][0] < 0) 
            || (corners[0][0] < 0 && corners[1][0] > 0)
        )
        ||  //or  
        (   //Different X-signs between corner 2 and 3
            (corners[2][0] > 0 && corners[3][0] < 0) 
            || (corners[2][0] < 0 && corners[3][0] > 0)
        ) 
    );

    bool differentSignsY = (
        (   
            //Different Y-signs between corner 0 and 3
            (corners[0][1] > 0 && corners[3][1] < 0) 
            || (corners[0][1] < 0 && corners[3][1] > 0)
        )
        || //or
        (   
            //Different Y-signs between corner 1 and 2
            (corners[1][1] > 0 && corners[2][1] < 0) 
            || (corners[1][1] < 0 && corners[2][1] > 0)
        )
    );

    return differentSignsX && differentSignsY;
}

}  // namespace inviwo
