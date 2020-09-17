/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, Anke Friederici
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_{
    "org.inviwo.MarchingSquares",  // Class identifier
    "Marching Squares",            // Display name
    "KTH Lab",                     // Category
    CodeState::Experimental,       // Code state
    Tags::None,                    // Tags
};

const std::array<std::array<int, 4>, 16> lookupTable = {{
    {},
    {0, 3},
    {0, 1},
    {1, 3},
    {1, 2},
    {0, 1, 2, 3},
    {0, 2},
    {2, 3},
    {2, 3},
    {0, 2},
    {0, 1, 2, 3},
    {1, 2},
    {1, 3},
    {0, 1},
    {0, 3},
    {}}};

const ProcessorInfo MarchingSquares::getProcessorInfo() const { return processorInfo_; }

MarchingSquares::MarchingSquares()
    : Processor()
    , inData("volumeIn")
    , meshIsoOut("meshIsoOut")
    , meshGridOut("meshGridOut")
    , propShowGrid("showGrid", "Show Grid")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f),
                    vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput,
                    PropertySemantics::Color)
    , propDeciderType("deciderType", "Decider Type")
    , propRandomSeed("seed", "Random Seed", 0, 0, std::mt19937::max())
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f), vec4(1.0f),
                   vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData) {
    // Register ports
    addPort(inData);
    addPort(meshIsoOut);
    addPort(meshGridOut);

    // Register properties
    addProperty(propShowGrid);
    addProperty(propGridColor);

    addProperty(propDeciderType);
    propDeciderType.addOption("asymptotic", "Asymptotic", 0);
    propDeciderType.addOption("random", "Random", 1);

    addProperty(propRandomSeed);
    propRandomSeed.setSemantics(PropertySemantics::Text);

    addProperty(propMultiple);

    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);

    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clear();
    propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 1.0f, 0.0f, 1.0f));
    propIsoTransferFunc.get().add(1.0f, vec4(1.0f, 0.0f, 0.0f, 1.0f));
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propGridColor, propRandomSeed, propNumContours, propIsoTransferFunc);

    propDeciderType.onChange([this]() {
        if (propDeciderType.get() == 1) {
            util::show(propRandomSeed);
        } else {
            util::hide(propRandomSeed);
        }
    });

    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]() {
        if (propShowGrid.get()) {
            util::show(propGridColor);
        } else {
            util::hide(propGridColor);
        }
    });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]() {
        if (propMultiple.get() == 0) {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        } else {
            //util::hide(propIsoValue);
            //util::show(propIsoColor, propNumContours);

            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
            util::hide(propIsoValue, propIsoColor);
            util::show(propNumContours, propIsoTransferFunc);
        }
    });
}

void MarchingSquares::process() {
    if (!inData.hasData()) {
        return;
    }

    // Create a structured grid from the input volume
    auto vol = inData.getData();
    auto grid = ScalarField2::createFieldFromVolume(vol);

    // Extract the minimum and maximum value from the input data
    const double minValue = grid.getMinValue();
    const double maxValue = grid.getMaxValue();

    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    //LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue
                                                                 // << ".");
    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from


    // - bounding box {xmin, ymin} - {xmax, ymax}
    const dvec2 bBoxMin = grid.getBBoxMin();
    const dvec2 bBoxMax = grid.getBBoxMax();
    // - cell size {dx, dy}
    const dvec2 cellSize = grid.getCellSize();

	const ivec2 nVertPerDim = grid.getNumVerticesPerDim();  //[5, 5] for simplegrid

    // Values at the vertex positions can be accessed by the indices of the vertex
    // with index i ranging between [0, nx-1] and j in [0, ny-1]
    ivec2 ij = {0, 0};
    double valueAt00 = grid.getValueAtVertex(ij);

    for (int i = 0; i < nVertPerDim[0]; i++) {
        for (int j = 0; j < nVertPerDim[1]; j++) {
            double value = grid.getValueAtVertex({i, j});
            //LogProcessorInfo("Value at (" << i << "," << j << "): " << value);
        }
    }
    //LogProcessorInfo("Bounding box: min: " << bBoxMin << " max: " << bBoxMax);
    // Initialize the output: mesh and vertices for the grid and bounding box
    auto gridmesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> gridvertices;

    auto indexBufferBBox = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // bottomLeft to topLeft
    drawLineSegment(bBoxMin, vec2(bBoxMin[0], bBoxMax[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topLeft to topRight
    drawLineSegment(vec2(bBoxMin[0], bBoxMax[1]), bBoxMax, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topRight to bottomRight
    drawLineSegment(bBoxMax, vec2(bBoxMax[0], bBoxMin[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // bottomRight to bottomLeft
    drawLineSegment(vec2(bBoxMax[0], bBoxMin[1]), bBoxMin, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);

    // Set the random seed to the one selected in the interface
    randGenerator.seed(static_cast<std::mt19937::result_type>(propRandomSeed.get()));
    // You can create a random sample between min and max with
    float minRand = 0.0;
    float maxRand = 1.0;
    float rand = randomValue(minRand, maxRand);

    // Properties are accessed with propertyName.get()
    if (propShowGrid.get()) {

        // The function drawLineSegments creates two vertices at the specified positions,
        // that are placed into the Vertex vector defining our mesh.
        // An index buffer specifies which of those vertices should be grouped into to make up
        // lines/trianges/quads. Here two vertices make up a line segment.
        auto indexBufferGrid = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        // Grid lines X-axis
        for (double i = bBoxMin[0] + cellSize[0]; i < bBoxMax[0]; i += cellSize[0]) {
            vec2 v1 = vec2(i, bBoxMin[1]);
            vec2 v2 = vec2(i, bBoxMax[1]);
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);
        }

        // Grid lines Y-axis
        for (double i = bBoxMin[1] + cellSize[1]; i < bBoxMax[1]; i = i + cellSize[1]) {
            vec2 v1 = vec2(bBoxMin[0], i);
            vec2 v2 = vec2(bBoxMax[0], i);
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);
        }
    }

    // Set the created grid mesh as output
    gridmesh->addVertices(gridvertices);
    meshGridOut.setData(gridmesh);

    // TODO (Bonus) Gaussian filter
    // Our input is const (i.e. cannot be altered), but you need to compute smoothed data and write
    // it somewhere
    // Create an editable structured grid with ScalarField2 smoothedField =
    // ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin); Values can be set with
    // smoothedField.setValueAtVertex({0, 0}, 4.2);
    // and read again in the same way as before
    // smoothedField.getValueAtVertex(ij);
    // Initialize the output: mesh and vertices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    
    if (propMultiple.get() == 0) {
        //LogProcessorInfo("Hit kommer vi också");
        renderIsoline(propIsoValue, &propIsoColor.get(), &grid, mesh, &vertices);
    }
    else {
        
        // TODO: Draw the given number (propNumContours) of isolines between
        // the minimum and maximum value
		auto stepSize = (abs(maxValue - minValue))/(propNumContours + 1);

		//vec4 color

        //LogProcessorInfo("MinVal: " << minValue);

		for (float k = 1; k <= propNumContours; k++) {
			LogProcessorInfo("k: " << k);
			// propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
            vec4 color = propIsoTransferFunc.get().sample(k/(propNumContours+1));
            LogProcessorInfo("Color: " << color);
			renderIsoline(minValue + k*stepSize, &color, &grid, mesh, &vertices);
		}

		 LogProcessorInfo("MaxVal: " << maxValue);


        // TODO (Bonus): Use the transfer function property to assign a color
        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, that means
        // vec4 color = propIsoTransferFunc.get().sample(0.0f);
        // is the color for the minimum value in the data
        // vec4 color = propIsoTransferFunc.get().sample(1.0f);
        // is the color for the maximum value in the data
    }

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour
	
    mesh->addVertices(vertices);
    meshIsoOut.setData(mesh);
}

void MarchingSquares::renderIsoline(double isoValue, vec4* color, ScalarField2* grid,
                                    std::shared_ptr<inviwo::BasicMesh> mesh,
                                    std::vector<BasicMesh::Vertex>* vertices) {
    // Get the definition of our structured grid with
    // - number of vertices in each dimension {nx, ny}
    const ivec2 nVertPerDim = grid->getNumVerticesPerDim();  //[5, 5] for simplegrid
	
    // TODO: Draw a single isoline at the specified isovalue (propIsoValue)
    // and color it with the specified color (propIsoColor)

    // mark all vertices
    std::vector<std::vector<int>> binaryImage(nVertPerDim[0], std::vector<int>(nVertPerDim[1]));

    for (auto i = 0; i < nVertPerDim[0]; ++i) {
        for (auto j = 0; j < nVertPerDim[1]; ++j) {
            if (grid->getValueAtVertex({i, j}) >= isoValue) {
                binaryImage[i][j] = 1;
            } else {
                binaryImage[i][j] = 0;
            }
            // LogProcessorInfo("Binary (" << i << ", " << j << "): " << binaryImage[i][j])
        }
    }

    // Build binary index
    for (auto i = 0; i < nVertPerDim[0] - 1; i++) {
        for (auto j = 0; j < nVertPerDim[1] - 1; j++) {
            auto cIndex = 0;
            cIndex |= binaryImage[i][j];
            cIndex |= (binaryImage[i + 1][j] << 1);
            cIndex |= (binaryImage[i + 1][j + 1] << 2);
            cIndex |= (binaryImage[i][j + 1] << 3);

            vec3 p0 = {grid->getPositionAtVertex({i, j}), grid->getValueAtVertex({i, j})};
            vec3 p1 = {grid->getPositionAtVertex({i + 1, j}), grid->getValueAtVertex({i + 1, j})};
            vec3 p2 = {grid->getPositionAtVertex({i + 1, j + 1}),
                       grid->getValueAtVertex({i + 1, j + 1})};
            vec3 p3 = {grid->getPositionAtVertex({i, j + 1}), grid->getValueAtVertex({i, j + 1})};

            auto indexBufferIsoContour =
                mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

            switch (cIndex) {
                case (0):
                case (15):
                    break;
                case (1):
                case (14): {
                    double e0 =
                        inverseLinearInterpolation(isoValue, {p0[2], p1[2]}, {p0[0], p1[0]});
                    double e3 = inverseLinearInterpolation(isoValue, {p0[2], p3[2]},
                                                           {p0[1], p3[1]});  // y-dir

                    vec2 v1 = {e0, p0[1]};
                    vec2 v2 = {p3[0], e3};

                    drawLineSegment(v1, v2, *color, indexBufferIsoContour.get(),
                                    *vertices);
                    break;
                }
                case (2):
                case (13): {
                    double e0 =
                        inverseLinearInterpolation(isoValue, {p0[2], p1[2]}, {p0[0], p1[0]});
                    double e1 = inverseLinearInterpolation(isoValue, {p1[2], p2[2]},
                                                           {p1[1], p2[1]});  // y-dir

                    vec2 v1 = {e0, p0[1]};
                    vec2 v2 = {p1[0], e1};

                    drawLineSegment(v1, v2, *color, indexBufferIsoContour.get(),
                                    *vertices);

                    break;
                }
                case (3):
                case (12): {
                    double e1 = inverseLinearInterpolation(isoValue, {p1[2], p2[2]},
                                                           {p1[1], p2[1]});  // y-dir
                    double e3 = inverseLinearInterpolation(isoValue, {p0[2], p3[2]},
                                                           {p0[1], p3[1]});  // y-dir

                    vec2 v1 = {p1[0], e1};
                    vec2 v2 = {p3[0], e3};

                    drawLineSegment(v1, v2, *color, indexBufferIsoContour.get(),
                                    *vertices);

                    break;
                }
                case (4):
                case (11): {
                    double e1 =
                        inverseLinearInterpolation(isoValue, {p1[2], p2[2]}, {p1[1], p2[1]});
                    double e2 =
                        inverseLinearInterpolation(isoValue, {p3[2], p2[2]}, {p3[0], p2[0]});

                    vec2 v1 = {p1[0], e1};
                    vec2 v2 = {e2, p2[1]};

                    drawLineSegment(v1, v2, *color, indexBufferIsoContour.get(),
                                    *vertices);
                    break;
                }
                case (5):
                case (10): {

                    double e0 =
                        inverseLinearInterpolation(isoValue, {p0[2], p1[2]}, {p0[0], p1[0]});
                    double e1 = inverseLinearInterpolation(isoValue, {p1[2], p2[2]},
                                                           {p1[1], p2[1]});  // y-dir
                    double e2 =
                        inverseLinearInterpolation(isoValue, {p2[2], p3[2]}, {p2[0], p3[0]});
                    double e3 = inverseLinearInterpolation(isoValue, {p0[2], p3[2]},
                                                           {p0[1], p3[1]});  // y-dir

                    int caseFive = 1;

                    if (propDeciderType.get() == 0) {
                        if (e2 > e0) {
                            caseFive = 0;
                        }
                    } else {
                        caseFive = randomValue(0, 2);
                    }

                    if (caseFive == 1) {
                        // case 5
                        vec2 v1 = {e0, p0[1]};
                        vec2 v2 = {p1[0], e1};
                        vec2 v3 = {e2, p2[1]};
                        vec2 v4 = {p3[0], e3};

                        drawLineSegment(v1, v2, *color, indexBufferIsoContour.get(),
                                        *vertices);
                        drawLineSegment(v3, v4, *color, indexBufferIsoContour.get(),
                                        *vertices);

                    } else {
                        // case 10
                        vec2 v1 = {e0, p0[1]};
                        vec2 v2 = {p1[0], e1};
                        vec2 v3 = {e2, p2[1]};
                        vec2 v4 = {p3[0], e3};

                        drawLineSegment(v1, v4, *color, indexBufferIsoContour.get(),
                                        *vertices);
                        drawLineSegment(v2, v3, *color, indexBufferIsoContour.get(),
                                        *vertices);
                    }
                    break;
                }
                case (6):
                case (9): {
                    double e0 =
                        inverseLinearInterpolation(isoValue, {p0[2], p1[2]}, {p0[0], p1[0]});
                    double e2 =
                        inverseLinearInterpolation(isoValue, {p2[2], p3[2]}, {p2[0], p3[0]});

                    vec2 v1 = {e0, p0[1]};
                    vec2 v2 = {e2, p2[1]};

                    drawLineSegment(v1, v2, *color, indexBufferIsoContour.get(),
                                    *vertices);
                    break;
                }
                case (7):
                case (8): {
                    double e2 =
                        inverseLinearInterpolation(isoValue, {p2[2], p3[2]}, {p2[0], p3[0]});
                    double e3 = inverseLinearInterpolation(isoValue, {p0[2], p3[2]},
                                                           {p0[1], p3[1]});  // y-dir
                    vec2 v1 = {e2, p2[1]};
                    vec2 v2 = {p3[0], e3};
                    drawLineSegment(v1, v2, *color, indexBufferIsoContour.get(),
                                    *vertices);

                    break;
                }
                default:
                    break;
            }
        }
    }
	
}

double MarchingSquares::inverseLinearInterpolation(const double isoValue, vec2 z, vec2 p) {
    return p[0] + ((isoValue - z[0])/((z[1]-z[0])/(p[1]-p[0])));
}

float MarchingSquares::randomValue(const float min, const float max) const {
    return min + uniformReal(randGenerator) * (max - min);
}

void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

}  // namespace inviwo
