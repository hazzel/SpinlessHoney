# Start configuration
cmake_minimum_required(VERSION 2.8)
project(spinless-frg)

enable_language(CXX)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(GSL)

FIND_PACKAGE(Boost 1.40 COMPONENTS system filesystem REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11 -pipe -g ${DEFINES} -DPROFILER")
find_library(PROFILER NAMES profiler PATHS $ENV{HOME}/gperftools-2.5/install/lib)
include_directories($ENV{HOME}/gperftools-2.5/install/include)
message("Using profiler: ${PROFILER}")

# Create executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)
add_executable(spinless ${FKT_SRC} spinless.cpp)
if (PROFILER)
	target_link_libraries(spinless ${Boost_LIBRARIES} ${PROFILER} ${GSL_LIBRARIES})
else()
	target_link_libraries(spinless ${Boost_LIBRARIES} ${GSL_LIBRARIES})
endif()
