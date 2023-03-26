# This is an oclMersenneTwister project to show and test the usage of the FindOpenCL
# script.

cmake_minimum_required( VERSION 2.6 )
project( clMussels_Mu )

# This is only required for the script to work in the version control
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../..")

find_package( OpenCL REQUIRED )

include_directories( ${OPENCL_INCLUDE_DIRS} )

if( OPENCL_HAS_CPP_BINDINGS )
	message( "OpenCL has CPP bindings. Full include is: " ${OPENCL_INCLUDE_DIRS} )
else( OPENCL_HAS_CPP_BINDINGS )
	message( "No OpenCL CPP bindings found" )
endif( OPENCL_HAS_CPP_BINDINGS )

add_executable( clMussels_Multiscale 
		Main_Program.cpp
		Transfer_between_PDE_IBM.cpp
		Check_Utilities.cpp
		Initializations.cpp
		PDE_Computing_Kernel.cpp
		Save_data.cpp
		Device_Utilities.cpp
		IBM_program.cpp
		Initializations.cpp
		Profiler.cpp
		oclRadixSort/oclRadixSort.cpp
		oclRadixSort/oclCommon.cpp
		oclRadixSort/oclScan.cpp
)

###############################################################################
#  include folder locations for CUDA, CUDPP, and the project
###############################################################################

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	INCLUDE_DIRECTORIES(
		./oclMersenneTwister
		./oclRadixSort)
else()
	INCLUDE_DIRECTORIES(
		./oclMersenneTwister
		./oclRadixSort)
endif()

target_link_libraries( clMussels_Multiscale ${OPENCL_LIBRARIES} )

# Generating the data folder where the model stores its results
file(MAKE_DIRECTORY ./build/Data)

# Defining the "run" command (use as: make run) that runs the model
add_custom_target(run ./clMussels_Multiscale)

# Defining the “destroy: command (use as: make wipe) that clears the build folder
add_custom_target(destroy rm -rf ../build/*
                  COMMENT “Clears_the_build_directory”)
