# Define module and return the generated name as MODULE_NAME
ALLPIX_UNIQUE_MODULE(MODULE_NAME)

#First make the root library for the pixel class, then make the module library

# Find the RootNewMacros.cmake file and include it
GET_FILENAME_COMPONENT(ROOT_CMAKE_DIR ${ROOT_USE_FILE} DIRECTORY)
FIND_FILE(ROOT_MACROS_FILE
          NAMES RootNewMacros.cmake RootMacros.cmake
          HINTS ${ROOT_CMAKE_DIR} ${ROOT_CMAKE_DIR}/modules/ $ENV{ROOTSYS}/cmake/modules)

IF(NOT ROOT_MACROS_FILE)
    MESSAGE(WARNING "Cannot find ROOT macros, including generic file as fallback.")
    INCLUDE(${ROOT_USE_FILE})
ELSE()
    INCLUDE(${ROOT_MACROS_FILE})
ENDIF()

# Generate the ROOT dictionary
ROOT_GENERATE_DICTIONARY(CaloOutputWritersDictionary
    ${CMAKE_CURRENT_SOURCE_DIR}/CaloOutputWriterModule.hpp
    LINKDEF
    ${CMAKE_CURRENT_SOURCE_DIR}/Linkdef.h
    -inlineInputHeader
    -I${CMAKE_CURRENT_SOURCE_DIR}
    MODULE
    CaloOutputWriters
)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/CaloOutputWritersDictionary.cxx
  PROPERTIES  GENERATED TRUE COMPILE_FLAGS "-w" # or whatever ?
)


# Define the library adding the object file created above
ADD_LIBRARY(CaloOutputWriters SHARED
    ${CMAKE_CURRENT_BINARY_DIR}/CaloOutputWritersDictionary.cxx
)

# Link the standard dependencies
TARGET_LINK_LIBRARIES(CaloOutputWriters ${ALLPIX_DEPS_LIBRARIES})

# Specify install for the messages
INSTALL(TARGETS CaloOutputWriters
    COMPONENT modules
    EXPORT Allpix
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)

# Also install the dictionary objects
INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libCaloOutputWriters_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libCaloOutputWriters.rootmap
    COMPONENT modules
    DESTINATION lib)

# Add source files to library
ALLPIX_MODULE_SOURCES(${MODULE_NAME}
    CaloOutputWriterModule.cpp
    Cluster.cpp
    # ADD SOURCE FILES HERE...
)

TARGET_LINK_LIBRARIES(${MODULE_NAME} CaloOutputWriters ROOT::Tree)

# Provide standard install target
ALLPIX_MODULE_INSTALL(${MODULE_NAME})
