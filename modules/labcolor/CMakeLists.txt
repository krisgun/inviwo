#--------------------------------------------------------------------
# Inviwo LabColor Module
ivw_module(LabColor)

#--------------------------------------------------------------------
# Add header files
set(HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/colorinterpolation.h
    ${CMAKE_CURRENT_SOURCE_DIR}/colormixing.h
    ${CMAKE_CURRENT_SOURCE_DIR}/colormixingbase.h
    ${CMAKE_CURRENT_SOURCE_DIR}/colorspace/src/ColorSpace.h
    ${CMAKE_CURRENT_SOURCE_DIR}/colorspace/src/Comparison.h
    ${CMAKE_CURRENT_SOURCE_DIR}/colorspace/src/Conversion.h
    ${CMAKE_CURRENT_SOURCE_DIR}/colorspace/src/Utils.h
)
#~ ivw_group("Header Files" ${HEADER_FILES})

#--------------------------------------------------------------------
# Add source files
set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/colorinterpolation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/colormixing.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/colormixingbase.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/colorspace/src/ColorSpace.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/colorspace/src/Comparison.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/colorspace/src/Conversion.cpp
)
ivw_group("Sources" ${SOURCE_FILES} ${HEADER_FILES})


#--------------------------------------------------------------------
# Add shaders
set(SHADER_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/glsl/labcolorprocessor.frag
)
ivw_group("Shader Files" ${SHADER_FILES})


#--------------------------------------------------------------------
# Add Unittests
set(TEST_FILES
	#${CMAKE_CURRENT_SOURCE_DIR}/tests/unittests/labcolor-test.cpp
)
ivw_add_unittest(${TEST_FILES})

#--------------------------------------------------------------------
# Create module
ivw_create_module(${SOURCE_FILES} ${HEADER_FILES} ${SHADER_FILES})

#--------------------------------------------------------------------
# Add shader directory to pack
# ivw_add_to_module_pack(${CMAKE_CURRENT_SOURCE_DIR}/glsl)
