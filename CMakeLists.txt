cmake_minimum_required(VERSION 3.17)
project(hw)

# soft installation path
set(SW $ENV{HOME}/mpd/sw/slc7_x86-64)
set(ACTS_SRC $ENV{HOME}/mpd/acts)
set(Boost_DIR ${SW}/boost/latest/lib/cmake/Boost-1.83.0)
find_package(Eigen3 REQUIRED HINTS "${SW}/Eigen3/latest")
find_package(Acts REQUIRED)

include_directories(
        ${SW}/TBB/latest/include
        ${ACTS_SRC}/Examples/Detectors/TelescopeDetector/include
        ${ACTS_SRC}/Examples/Detectors/MagneticField/include
        ${ACTS_SRC}/Examples/Algorithms/Generators
        ${ACTS_SRC}/Examples/Algorithms/Fatras/include
        ${ACTS_SRC}/Examples/Algorithms/Digitization/include
        ${ACTS_SRC}/Examples/Algorithms/TrackFinding/include
        ${ACTS_SRC}/Examples/Algorithms/TruthTracking/
        ${ACTS_SRC}/Examples/Run/Common/include
        ${ACTS_SRC}/Examples/Io/Json/include
        ${ACTS_SRC}/Examples/Io/Root/include
        ${ACTS_SRC}/Examples/Io/Csv/include
)
link_directories(
        ${SW}/TBB/latest/lib
        ${SW}/acts/latest/lib
        )


# main exec
add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC
        ActsCore
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
        ActsExamplesTruthTracking
        tbb
)

# copy run macros
set(SCRIPTS
        digi-smearing-config.json)

foreach (_script ${SCRIPTS})
    configure_file(${_script}
            ${PROJECT_BINARY_DIR}/${_script}
            COPYONLY)
endforeach ()