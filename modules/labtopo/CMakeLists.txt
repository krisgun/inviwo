#--------------------------------------------------------------------
# Inviwo LabTopo Module
ivw_module(LabTopo)

#--------------------------------------------------------------------
# Add header files
set(HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/integrator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/interpolator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/topology.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/gradients.h
    # tnt/jama headers for eigenvector computation
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/tnt/jama_eig.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/tnt/tnt_array1d.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/tnt/tnt_array2d.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/tnt/tnt_math_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/tnt/tnt_i_refvec.h
)

#--------------------------------------------------------------------
# Add source files
set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/integrator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/interpolator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/topology.cpp
)
ivw_group("Sources" ${SOURCE_FILES} ${HEADER_FILES})

#--------------------------------------------------------------------
# Create module
ivw_create_module(${SOURCE_FILES} ${HEADER_FILES} ${SHADER_FILES})

