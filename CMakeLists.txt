cmake_minimum_required(VERSION 3.17)
project(actsTest)

include(GNUInstallDirs)

# soft installation path
set(SW $ENV{HOME}/mpd/sw/ubuntu2404_x86-64)
set(ACTS_SRC ${SW}/ACTS/latest)
set(Boost_DIR ${SW}/Boost/latest/lib/cmake/Boost-1.83.0)
set(Acts_DIR ${SW}/ACTS/latest/lib/cmake/Acts)
set(DD4hep_DIR ${SW}/DD4hep/latest/cmake)
#set(Geant4_DIR ${SW}/GEANT4/latest/lib/cmake/Geant4)
set(ROOT_DIR ${SW}/ROOT/latest/cmake)
set(nlohmann_json_DIR ${SW}/nlohmann_json/latest/share/cmake/nlohmann_json)
set(LCIO_DIR ${SW}/LCIO/latest/lib/cmake/LCIO)
set(SIO_DIR ${SW}/LCIO/latest/lib/cmake/SIO)
set(EDM4HEP_DIR ${SW}/EDM4hep/latest/lib/cmake/EDM4HEP)
set(podio_DIR ${SW}/podio/latest/lib/cmake/podio)
set(HEPMC3_DIR ${SW}/HepMC3/latest/share/HepMC3/cmake)

find_package(Geant4 REQUIRED vis_all ui_all)
include(${Geant4_USE_FILE})

find_package(Eigen3 REQUIRED HINTS "${SW}/Eigen3/latest")

find_package(ROOT REQUIRED)
if(ROOT_FOUND)
    message(STATUS "Found ROOT")
    message(STATUS "ROOT includes: ${ROOT_INCLUDE_DIRS}")
    message(STATUS "ROOT libs: ${ROOT_LIBRARY_DIR}, ${ROOT_LIBRARIES}")
endif()

find_package(Acts REQUIRED)

include(${ROOT_USE_FILE})

set(PROJECT_INCLUDE_DIRS
        $ENV{PWD}
        ${SW}/Boost/latest/include
        ${SW}/Eigen3/latest/include/eigen3
        ${SW}/ACTS/latest/include
        ${SW}/TBB/latest/include
        ${SW}/nlohmann_json/latest/include
        ${SW}/ROOT/latest/include
        ${SW}/GEANT4/latest/include
        ${HEPMC3_INCLUDE_DIR}
)

set(PROJECT_LINK_DIRS
        ${SW}/TBB/latest/lib
        ${SW}/ACTS/latest/lib
        ${SW}/GEANT4/latest/lib
        ${ROOT_LIBRARY_DIR}
)

include_directories(${PROJECT_INCLUDE_DIRS})

link_directories(${PROJECT_LINK_DIRS})

set(PROJECT_SOURCES
        MySpacePointMaker.cpp
        MyDigitizationAlgorithm.cpp
#        MyTrackFindingAlgorithm.cpp
        MyTrackWriter.cpp
        MyFtdGeo.cxx
        MyFtdDetector.cpp
        BaseTpcSectorGeo.cxx
)

set(PROJECT_HEADERS
        MySpacePointMaker.hpp
        MyDigitizationAlgorithm.hpp
#        MyTrackFindingAlgorithm.hpp
        MyTrackWriter.hpp
        MyFtdGeo.h
        MyFtdDetector.h
        BaseTpcSectorGeo.h
)

set(PROJECT_LINK_LIBS
        ActsCore
        ActsFatras
        ActsPluginGeant4
        ActsExamplesFramework
        ActsExamplesGenerators
        ActsExamplesDetectorTelescope
        ActsExamplesFatras
        ActsPluginFpeMonitoring
        ActsExamplesMagneticField
        ActsExamplesDigitization
        ActsExamplesIoJson
        ActsExamplesIoRoot
        ActsExamplesIoCsv
        ActsExamplesTrackFinding
        ActsExamplesTrackFitting
        ActsExamplesTruthTracking
        ActsExamplesIoHepMC3
        ActsExamplesAmbiguityResolution
        tbb
        ${Geant4_LIBRARIES}
        ${ROOT_LIBRARIES}
)

# add dynamic library
add_library(actsTestLib SHARED
        ${PROJECT_SOURCES}
        ${PROJECT_HEADERS}
)

target_include_directories(actsTestLib
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        PRIVATE
        ${PROJECT_INCLUDE_DIRS}
)

target_link_libraries(actsTestLib PUBLIC ${PROJECT_LINK_LIBS})

# main exec
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE actsTestLib)

set(CMAKE_INSTALL_LIBDIR lib)
set(CMAKE_INSTALL_INCLUDEDIR include)
set(CMAKE_INSTALL_PREFIX $ENV{PWD}/stage)

install(TARGETS actsTestLib
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES ${PROJECT_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
