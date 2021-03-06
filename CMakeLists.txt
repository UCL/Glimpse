cmake_minimum_required(VERSION 3.0)

# Paths to external libraries
set(NICAEA_URL https://github.com/UCL/nicaea)
set(SPARSE2D_URL https://github.com/CosmoStat/Sparse2D.git)

# Adding customized cmake module
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(FindPkgConfig)
include(ExternalProject)

project(glimpse)

# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


# Include external libraries
find_package(GSL       REQUIRED)
find_package(CCFits    REQUIRED)
find_package(FFTW      REQUIRED)
find_package(NFFT      REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Boost     COMPONENTS program_options REQUIRED)

# Include external projects
## NICAEA
ExternalProject_Add(NICAEA
        PREFIX NICAEA
        GIT_REPOSITORY ${NICAEA_URL}
        GIT_TAG 03aa73cfaa34f38c67cf583dd0b476d8392d6337
        CMAKE_CACHE_ARGS "-DFFTW_ROOT:FILEPATH=${FFTW_ROOT}" "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
        set(NICAEA_LIBRARY_DIR ${CMAKE_BINARY_DIR}/NICAEA/src/NICAEA/lib)
        set(NICAEA_INCLUDE_DIR ${CMAKE_BINARY_DIR}/NICAEA/src/NICAEA/include)
        set(NICAEA_LIBRARIES -lnicaea)

## Sparse2D
ExternalProject_Add(SPARSE2D
         PREFIX SPARSE2D
         GIT_REPOSITORY ${SPARSE2D_URL}
         GIT_TAG v2.1.2
         CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/extern" "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
        set(SPARSE2D_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib)
        set(SPARSE2D_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include)
        set(SPARSE2D_LIBRARIES -lsparse2d -ltools -lsparse1d)

include_directories(${SPARSE2D_INCLUDE_DIR}
		    ${GSL_INCLUDE_DIR}
		    ${CCFITS_INCLUDE_DIR}
		    ${CFITSIO_INCLUDE_DIR}
		    ${Boost_INCLUDE_DIR}
		    ${NICAEA_INCLUDE_DIR}
		    ${NFFT_INCLUDE_DIR}
		    ${FFTW_INCLUDE_DIR}
		    ${ARMADILLO_INCLUDE_DIRS})

link_directories(${SPARSE2D_LIBRARY_DIR}
		    ${GSL_LIBRARY_DIR}
		    ${CCFITS_LIBRARY_DIR}
		    ${CFITSIO_LIBRARY_DIR}
		    ${Boost_LIBRARY_DIR}
		    ${NICAEA_LIBRARY_DIR}
		    ${NFFT_LIBRARY_DIR}
		    ${FFTW_LIBRARY_DIR})

set(GLIMPSE_SRC 
		src/glimpse_functions.cpp
		src/survey.cpp
		src/redshift_distribution.cpp
		src/field.cpp
		src/surface_reconstruction.cpp
		src/density_reconstruction.cpp
		src/starlet_2d.cpp
		src/wavelet_transform.cpp
		src/gpu_utils.c)

# VERSIONING from git
execute_process(
  COMMAND git describe --tags
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.h.in
  ${CMAKE_SOURCE_DIR}/src/version.h
)

find_package(CUDA)

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lomp -Xpreprocessor")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -fpermissive -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -DDEBUG_FITS")

if(${CUDA_FOUND})
    message("Compiling CUDA accelerated reconstruction code, with 3D support")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCUDA_ACC")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -O3")
    foreach(CAPABILITY ${CUDA_COMPUTE_CAPABILITIES})
      set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_${CAPABILITY},code=sm_${CAPABILITY}")
    endforeach()

    cuda_add_executable(glimpse src/glimpse.cpp ${GLIMPSE_SRC} src/spg.cu src/spg.cpp)
    cuda_add_cufft_to_target(glimpse)
else(${CUDA_FOUND})
    message("Compiling without CUDA acceleration")
    add_executable(glimpse src/glimpse.cpp ${GLIMPSE_SRC})
endif(${CUDA_FOUND})

add_dependencies(glimpse NICAEA SPARSE2D)

target_link_libraries(glimpse ${NICAEA_LIBRARIES}
                              ${SPARSE2D_LIBRARIES}
                              ${Boost_LIBRARIES}
                              ${GSL_LIBRARIES}
                              ${NFFT_LIBRARIES}
                              ${FFTW_LIBRARIES}
                              ${CCFITS_LIBRARY}
                              ${CFITSIO_LIBRARY}
                              ${ARMADILLO_LIBRARIES})

# Testing
option(BUILD_TESTING "Build the testing tree." OFF)
include(CTest)

if(BUILD_TESTING)
  find_package(Catch2 REQUIRED)
  enable_testing()

  set(TEST_WAVELET_TRANSFORM_SRC
    test/test_wavelet_transform.cpp
    src/wavelet_transform.cpp
    src/starlet_2d.cpp
    src/gpu_utils.c
    )
  if(${CUDA_FOUND})
    cuda_add_executable(test_wavelet_transform ${TEST_WAVELET_TRANSFORM_SRC})
    cuda_add_cufft_to_target(test_wavelet_transform)
  else(${CUDA_FOUND})
    add_executable(test_wavelet_transform ${TEST_WAVELET_TRANSFORM_SRC})
  endif(${CUDA_FOUND})
  target_include_directories(test_wavelet_transform PRIVATE /usr/local/include)
  target_link_libraries(test_wavelet_transform ${SPARSE2D_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${CCFITS_LIBRARY}
    ${CFITSIO_LIBRARY})
  add_dependencies(test_wavelet_transform SPARSE2D)

  add_test( test_wavelet_transform_test test_wavelet_transform )

  set(TEST_FIELD_SRC
    test/test_field.cpp
    src/survey.cpp
    src/redshift_distribution.cpp
    src/field.cpp
    # src/surface_reconstruction.cpp
    # src/density_reconstruction.cpp
    src/starlet_2d.cpp
    src/wavelet_transform.cpp
    src/gpu_utils.c
    )
  if(${CUDA_FOUND})
    cuda_add_executable(test_field ${TEST_FIELD_SRC})
    cuda_add_cufft_to_target(test_field)
  else(${CUDA_FOUND})
    add_executable(test_field ${TEST_FIELD_SRC})
  endif(${CUDA_FOUND})
  target_link_libraries(test_field ${NICAEA_LIBRARIES}
    ${SPARSE2D_LIBRARIES}
    ${Boost_LIBRARIES}
    ${GSL_LIBRARIES}
    ${NFFT_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${CCFITS_LIBRARY}
    ${CFITSIO_LIBRARY}
    ${ARMADILLO_LIBRARIES}
    Catch2::Catch2)
  add_dependencies(test_field NICAEA
    SPARSE2D)

  add_test( test_field_test test_field )
  set_tests_properties(test_field_test PROPERTIES WILL_FAIL TRUE)
endif(BUILD_TESTING)
