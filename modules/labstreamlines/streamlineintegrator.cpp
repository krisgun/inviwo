/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/interaction/events/mouseevent.h>
#include <inviwo/core/util/utilities.h>
#include <labstreamlines/integrator.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming
// scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5f, 0.5f), vec2(-1.f), vec2(1.f), vec2(0.1))
    , propSeedMode("seedMode", "Seeds")
    , propNumStepsTaken("numstepstaken", "Number of actual steps", 0, 0, 100000)
    , mouseMoveStart("mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
                     MouseButton::Left, MouseState::Press | MouseState::Move)
    , propForward("forward", "Forward Integration", true)
    , propBackward("backward", "Backward Integration", true)
    , propDirectionField("directionField", "Integrate in the Direction Field")
    , propStepSize("stepSize", "Step Size", 0.1f, 0.001f, 1.0f)
    , propNumberOfSteps("numberOfSteps", "Steps", 100, 1, 10000)
    , propArcLength("arcLength", "Arc Length Limit", 10, 0, 100)
    , propVelocityLimit("velocityLimit", "Velocity Limit")
    , propNumberOfRandomSeeds("numberOfRandomSeeds", "Random Seeds", 0, 0, 100)
    , propUniformGridSeedsX("uniformGridSeedsX", "X Seeds", 6, 0, 25)
    , propUniformGridSeedsY("uniformGridSeedsY", "Y Seeds", 6, 0, 25)
				
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional),
// increment (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(inData);
    addPort(meshOut);
    addPort(meshBBoxOut);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
	
    addProperty(propSeedMode);
    addProperty(propStartPoint);
    addProperty(propNumStepsTaken);
    propNumStepsTaken.setReadOnly(true);
    propNumStepsTaken.setSemantics(PropertySemantics::Text);
    addProperty(mouseMoveStart);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(propForward);
    addProperty(propBackward);
    addProperty(propStepSize);
    addProperty(propDirectionField);
    addProperty(propNumberOfSteps);
    addProperty(propArcLength);
    addProperty(propVelocityLimit);
    addProperty(propNumberOfRandomSeeds);
    addProperty(propUniformGridSeedsX);
    addProperty(propUniformGridSeedsY);

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart, propNumStepsTaken);
            // util::hide(...)
        } else {
            util::hide(propStartPoint, mouseMoveStart, propNumStepsTaken);
            // util::show(...)
        }
    });
}

void StreamlineIntegrator::eventMoveStart(Event* event) {
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to bounding box range
    mousePos[0] *= static_cast<float>(BBoxMax_[0] - BBoxMin_[0]);
    mousePos[1] *= static_cast<float>(BBoxMax_[1] - BBoxMin_[1]);
    mousePos += static_cast<vec2>(BBoxMin_);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void StreamlineIntegrator::process() {
    srand(1.0);
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;

    // Make bounding box without vertex duplication, instead of line segments which duplicate
    // vertices, create line segments between each added points with connectivity type of the index
    // buffer
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin_[0], BBoxMax_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax_[0], BBoxMin_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

	auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

    if (propSeedMode.get() == 0) {
        vec2 startPoint = propStartPoint.get();
        // Draw start point
        Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

        // TODO: Create one stream line from the given start point
        vec4 color = vec4(1, 0, 0, 1);
        StreamlineIntegrator::createStreamLine(vectorField, startPoint, propForward, propBackward, propStepSize, 
            propDirectionField, propNumberOfSteps, propArcLength, propVelocityLimit, color, indexBufferLines, vertices);

        // TODO: Use the propNumStepsTaken property to show how many steps have actually been
        // integrated This could be different from the desired number of steps due to stopping
        // conditions (too slow, boundary, ...)
        propNumStepsTaken.set(0);

    } else {
		vec4 color = vec4(0.1, 0.1, 0.1, 1);
		vec2 startPoint;
        for (int i = 0; i < propNumberOfRandomSeeds; i++) {
            double f = (double)rand() / RAND_MAX;
            double randomX = BBoxMin_[0] + f * (BBoxMax_[0] - BBoxMin_[0]);
            f = (double)rand() / RAND_MAX;
            double randomY = BBoxMin_[1] + f * (BBoxMax_[1] - BBoxMin_[1]);
            startPoint = vec2(randomX, randomY);

            // Draw start point
            //Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

            // TODO: Create one stream line from the given start point
            
            StreamlineIntegrator::createStreamLine(
                vectorField, startPoint, propForward, propBackward, propStepSize,
                propDirectionField, propNumberOfSteps, propArcLength, propVelocityLimit, color,
                indexBufferLines, vertices);

            // TODO: Use the propNumStepsTaken property to show how many steps have actually been
            // integrated This could be different from the desired number of steps due to stopping
            // conditions (too slow, boundary, ...)
            propNumStepsTaken.set(0);
		}
		
		if (propUniformGridSeedsX > 1 && propUniformGridSeedsY > 1) {
			double xJump = (BBoxMax_[0] - BBoxMin_[0]) / (propUniformGridSeedsX - 1);
			double yJump = (BBoxMax_[1] - BBoxMin_[1]) / (propUniformGridSeedsY - 1);

			for (int i = 0; i < propUniformGridSeedsX; i++) {
				for (int j = 0; j < propUniformGridSeedsY; j++) {

					startPoint = vec2(BBoxMin_[0] + i * xJump, BBoxMin_[1] + j * yJump);
					//Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

					StreamlineIntegrator::createStreamLine(
						vectorField, startPoint, propForward, propBackward, propStepSize,
						propDirectionField, propNumberOfSteps, propArcLength, propVelocityLimit, color,
						indexBufferLines, vertices);
				}
			}
		}
        
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
    }

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}  // namespace inviwo

void StreamlineIntegrator::createStreamLine(VectorField2& vectorField, vec2& startPoint, bool forward, bool backward,
    const float stepSize, bool directionField, int steps, float arcLength, float velocityLimit, vec4 color, 
    std::shared_ptr<inviwo::IndexBufferRAM>& indexBuffer, std::vector<BasicMesh::Vertex>& vertices) {

    if (directionField) {
         vectorField = StreamlineIntegrator::normalizeVectorField(vectorField);
    }

    if (forward && !backward) {
        StreamlineIntegrator::integratePoints(vectorField, startPoint, stepSize, steps, velocityLimit, color, indexBuffer, vertices, arcLength);
    }
    
    else if (!forward && backward) {
        StreamlineIntegrator::integratePoints(vectorField, startPoint, -stepSize, steps, velocityLimit, color, indexBuffer, vertices, arcLength);
    }
    else if (forward && backward) {
        StreamlineIntegrator::integratePoints(vectorField, startPoint, stepSize, (steps/2)+steps%2, velocityLimit, color, indexBuffer, vertices, arcLength/2);
        StreamlineIntegrator::integratePoints(vectorField, startPoint, -stepSize, steps/2, velocityLimit, color, indexBuffer, vertices, arcLength/2);
    }
    
}

void StreamlineIntegrator::integratePoints(VectorField2& vectorField, vec2 startPoint,
    const float stepSize, int steps, float velocityLimit, vec4 color,
    std::shared_ptr<inviwo::IndexBufferRAM>& indexBuffer, std::vector<BasicMesh::Vertex>& vertices, float arcLength) {

    ivec2 dimensions = vectorField.getNumVerticesPerDim();

	float integratedArcLength = 0; 

    dvec2 oldPoint = startPoint;
    for (int i = 0; i < steps; i++) {
        
        dvec2 fieldVec = vectorField.interpolate(startPoint);
        //Stop when norm is smaller than velocity limit
       
        bool isLessThanVelocityLimit = StreamlineIntegrator::vectorNorm(fieldVec) < velocityLimit;

        //Stop when encountering zeros in vector field 
        bool isVecFieldZero = (fieldVec[0] > -0.001 && fieldVec[0] < 0.001) && (fieldVec[1] > -0.001 && fieldVec[1] < 0.001);
        
		//Stop at arc length
        bool isGreaterArcLengthOrZero = integratedArcLength > arcLength || arcLength == 0;  

		if (isLessThanVelocityLimit || isVecFieldZero || isGreaterArcLengthOrZero) break;
        
        startPoint = Integrator::RK4(vectorField, startPoint, stepSize);

		integratedArcLength += StreamlineIntegrator::vectorNorm(dvec2(startPoint[0]-oldPoint[0], startPoint[1]-oldPoint[1]));

        //Stop at boundary of bbox
        if (!vectorField.isInside(startPoint)) {
            startPoint = vectorField.clampPositionToBBox(startPoint);
			if (!((startPoint[0] == oldPoint[0] && (startPoint[0] == BBoxMin_[0] || startPoint[0] == BBoxMax_[0])) || 
				(startPoint[1] == oldPoint[1] && (startPoint[1] == BBoxMin_[1] || startPoint[1] == BBoxMax_[1])))) {
				Integrator::drawLineSegment(oldPoint, startPoint, color, indexBuffer.get(), vertices);
			}
            break;
        }

        //Draw line segments between points
        Integrator::drawLineSegment(oldPoint, startPoint, color, indexBuffer.get(), vertices);
        oldPoint = startPoint;
    }
}

VectorField2 StreamlineIntegrator::normalizeVectorField(VectorField2 vectorField) {
    ivec2 dimensions = vectorField.getNumVerticesPerDim();

    for (int i = 0; i < dimensions[0]; i++) {
        for (int j = 0; j < dimensions[1]; j++) {
            dvec2 vector = vectorField.getValueAtVertex({i, j});
            double length = StreamlineIntegrator::vectorNorm(vector);
            vectorField.setValueAtVertex({i, j}, dvec2(vector[0]/length,vector[1]/length));
        }
    }
    return vectorField;
}

double StreamlineIntegrator::vectorNorm(dvec2 vector) {
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2));
}

}  // namespace inviwo
