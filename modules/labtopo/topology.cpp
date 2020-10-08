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
#include <labstreamlines/streamlineintegrator.h>

namespace inviwo {

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 1, 0, 1),    // Saddle - Yellow
    vec4(1, 0, 0, 1),    // AttractingNode - Red
    vec4(0, 0, 1, 1),    // RepellingNode - Blue
    vec4(0.5, 0, 1, 1),  // AttractingFocus - Purple
    vec4(1, 0.5, 0, 1),  // RepellingFocus - Orange
    vec4(0, 1, 0, 1),     // Center - Green
};

enum TypeCP {
    Saddle = 0,
    AttractingNode = 1,
    RepellingNode = 2,
    AttractingFocus = 3,
    RepellingFocus = 4,
    Center = 5
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
    , propDirectionField("directionField", "Integrate in the Direction Field")
    , propEigenScale("eigenScale", "Eigen Vector Scaling", 0.1, 0.01, 1.0)
    , propSeparatricesSteps("separatricesSteps", "Int. Steps", 50, 1, 100)
    , propSeparatricesStepSize("separatricesStepSize", "Step Size", 0.1, 0.01, 1,0)
    , propSeparatricesVelocity("separatricesVelocity", "Vel. Lim.", 0.0, 0.0, 1.0)
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
    addProperty(propDirectionField);
    addProperty(propEigenScale);
    addProperty(propSeparatricesSteps);
    addProperty(propSeparatricesStepSize);
    addProperty(propSeparatricesVelocity);
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
    std::vector<dvec2> startCorners = {};
    std::vector<dvec2> criticalPoints {};
    double epsilon = 0.0001;

    size2_t dims = vectorField.getNumVerticesPerDim();
    // Looping through all values in the vector field.
    for (size_t j = 0; j < dims[1]-1; ++j) {
        for (size_t i = 0; i < dims[0]-1; ++i) {
            dvec2 corner0 = vectorField.getPositionAtVertex(size2_t(i, j));
            dvec2 corner1 = vectorField.getPositionAtVertex(size2_t(i+1, j));
            dvec2 corner2 = vectorField.getPositionAtVertex(size2_t(i+1, j+1));
            dvec2 corner3 = vectorField.getPositionAtVertex(size2_t(i, j+1));
            startCorners = {corner0, corner1, corner2, corner3};
            findCriticalPoints(vectorField, criticalPoints, epsilon, startCorners, propDirectionField);
        }
    }

    std::vector<std::vector<dvec2>> coloredCritPoints = Topology::classifyCriticalPoints(criticalPoints, vectorField);

    for (auto i = 0; i < coloredCritPoints.size(); ++i) {
        for (auto j = 0; j < coloredCritPoints[i].size(); ++j) {
            int colorIndex = i;
            Integrator::drawPoint(coloredCritPoints[i][j], ColorsCP[colorIndex], indexBufferPoints.get(), vertices);
        }
    }

    //Compute separatrices
    drawSeparatrices(vectorField, coloredCritPoints[Saddle], indexBufferSeparatrices, vertices);


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

void Topology::drawSeparatrices(const VectorField2& vectorField, std::vector<dvec2>& saddlePoints, std::shared_ptr<inviwo::IndexBufferRAM>& indexBuffer, std::vector<BasicMesh::Vertex>& vertices) {
	
	for (const auto &saddlePoint : saddlePoints) {
		dmat2 jacobian{ vectorField.derive(saddlePoint) };
		auto eigenResult = util::eigenAnalysis(jacobian);
		auto eigenVectors = eigenResult.eigenvectors;
        dvec2 eigenVector1 = eigenVectors[0];
        dvec2 eigenVector2 = eigenVectors[1];

        const double eigenScale{ propEigenScale };

        std::vector<dvec2> startPoints {
            {saddlePoint + eigenScale*eigenVector1}, {saddlePoint - eigenScale*eigenVector1}, 
            {saddlePoint + eigenScale*eigenVector2}, {saddlePoint - eigenScale*eigenVector2}
        };

        for (auto i = 0; i < startPoints.size(); ++i) {
            double stepSize = propSeparatricesStepSize;
            int steps = propSeparatricesSteps;
            int currStep = 0;
            double velocityLimit = propSeparatricesVelocity;
            bool arnold = false;
            dvec2 oldPoint = startPoints[i];

            double eigenValue = eigenResult.eigenvaluesRe[0];

            if (i > 1) {
                eigenValue = eigenResult.eigenvaluesRe[1];
            }
            int direction = eigenValue > 0 ? 1 : -1;
            Topology::drawLineSegment(startPoints[i], saddlePoint, vec4(1, 1, 1, 1), indexBuffer.get(), vertices);
            while (!arnold &&  currStep < steps) {
                startPoints[i] = Integrator::RK4(vectorField, startPoints[i], stepSize, direction, propDirectionField, velocityLimit, arnold);
                //Draw line segments between points
                Topology::drawLineSegment(oldPoint, startPoints[i], vec4(1, 1, 1, 1), indexBuffer.get(), vertices);
                oldPoint = startPoints[i];
                ++currStep;
            }
        } 
	}
}


bool Topology::isDuplicate(std::vector<dvec2>& criticalPoints, dvec2 newCriticalPoint, double epsilon) {
	for (int i = 0; i < criticalPoints.size(); ++i) {
		if (abs(criticalPoints[i][0] - newCriticalPoint[0]) < epsilon * 10 && abs(criticalPoints[i][1] - newCriticalPoint[1]) < epsilon * 10) {
			return true;
		}
	}
	return false;
}

/*
    3 *-------* 2
      |       |
      |       |
      |       |
    0 *-------* 1
*/

void Topology::findCriticalPoints(const VectorField2& vectorField, std::vector<dvec2>& criticalPoints, double epsilon, std::vector<dvec2>& corners, bool normalizeVecField) {
    std::vector<dvec2> cornerVectorValues{};
    
    VectorField2 usedVectorField = vectorField;

    if (normalizeVecField) { usedVectorField = StreamlineIntegrator::normalizeVectorField(usedVectorField); }

    double minPointDist = std::min((corners[1][0] - corners[0][0]), (corners[2][1] - corners[1][1]));

    for (int i = 0; i < corners.size(); ++i) {
        dvec2 fieldVec = usedVectorField.interpolate(corners[i]);
        cornerVectorValues.push_back(fieldVec);

        
        //Check if vector field is zero at the current corner
        double fieldZeroTol = 0.001;
        bool isVecFieldZero = (abs(fieldVec[0]) < fieldZeroTol) && (abs(fieldVec[1]) < fieldZeroTol);
        if (isVecFieldZero) {
            if (std::find(criticalPoints.begin(), criticalPoints.end(), corners[i]) == criticalPoints.end()) {
                criticalPoints.push_back(corners[i]);
            }
            return;
        }
        
    }

    //Change-of-sign test 
    if (isCornersDifferentSigns(cornerVectorValues)) {
        //Check if size of corner-box is smaller than epsilon
        if (minPointDist < epsilon) {
            dvec2 boxMidPoint { corners[0][0] + (corners[1][0] - corners[0][0]) / 2 , corners[0][1] + (corners[3][1] - corners[0][1]) / 2 };
            if (std::find(criticalPoints.begin(), criticalPoints.end(), (boxMidPoint)) == criticalPoints.end()) {
                criticalPoints.push_back(boxMidPoint);
                return;
            }
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
        findCriticalPoints(vectorField, criticalPoints, epsilon, q0Corners, normalizeVecField);
        findCriticalPoints(vectorField, criticalPoints, epsilon, q1Corners, normalizeVecField);
        findCriticalPoints(vectorField, criticalPoints, epsilon, q2Corners, normalizeVecField);
        findCriticalPoints(vectorField, criticalPoints, epsilon, q3Corners, normalizeVecField);
    }
}

//Returns 1 if number is greater than 0
//Returns -1 if number is smaller than 0
//Returns 0 otherwise
int Topology::sign(double number) {
    if (number > 0) return 1;
    if (number < 0) return -0;
    return 0;
}

//Check if both components of each corner vector are not the same
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

std::vector<std::vector<dvec2>> Topology::classifyCriticalPoints(std::vector<dvec2> criticalPoints, const VectorField2& vectorField) {
     
    std::vector<std::vector<dvec2> > classifiedCritPoints(
        6, std::vector<dvec2> ());


    for (auto const &point : criticalPoints) {
        
        dmat2 jacobian {vectorField.derive(point)};
        auto eigenResult = util::eigenAnalysis(jacobian);
        dvec2 imEigen = eigenResult.eigenvaluesIm;
        dvec2 reEigen = eigenResult.eigenvaluesRe;

        double threshold = 0.0001;
        bool reWithinThresh = (abs(reEigen[0]) < threshold) && (abs(reEigen[1]) < threshold);

        //Saddle point
        if (((reEigen[0] < 0 && reEigen[1] > 0) || (reEigen[1] < 0 && reEigen[0] > 0)) && imEigen[0] == 0 && imEigen[1] == 0) {
            classifiedCritPoints[Saddle].push_back(point);
        }
        
        //Repelling Node
        else if (reEigen[0] > 0 && reEigen[1] > 0 && imEigen[0] == 0 && imEigen[1] == 0) {
            classifiedCritPoints[RepellingNode].push_back(point);
        }

        //Attracting Node
        else if (reEigen[0] < 0 && reEigen[1] < 0 && imEigen[0] == 0 && imEigen[1] == 0) {
            classifiedCritPoints[AttractingNode].push_back(point);
        }

        //Center
        else if (reWithinThresh && imEigen[0] == -imEigen[1] && -imEigen[1] != 0) {
            classifiedCritPoints[Center].push_back(point);
        }

        //Attracting Focus
        else if (reEigen[0] == reEigen[1] && reEigen[1] < 0 && imEigen[0] == -imEigen[1] && -imEigen[1] != 0) {
            classifiedCritPoints[AttractingFocus].push_back(point);
        }

        //Repelling Focus
        else if (reEigen[0] == reEigen[1] && reEigen[1] > 0 && imEigen[0] == -imEigen[1] && -imEigen[1] != 0) {
            classifiedCritPoints[RepellingFocus].push_back(point);
        }
    }
    return classifiedCritPoints;
}


}  // namespace inviwo
