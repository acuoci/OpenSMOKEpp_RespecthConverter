#--------------------------------------------------------------------
# Preliminary Settings
#--------------------------------------------------------------------

cmake_minimum_required (VERSION 3.16)

project (OpenSMOKEpp_RespecthConverter)

# --- Set module path in order to use custom CMake modules
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{EIGEN_LIBRARY_PATH}/cmake)
set (CMAKE_COLOR_MAKEFILE ON)
set (VERSION "2.0.0")
set (INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")

if( NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

message (STATUS "Compiling in mode: ${CMAKE_BUILD_TYPE}")

#--------------------------------------------------------------------
# Compulsory libraries
#--------------------------------------------------------------------

# --- Boost (Compulsory)
find_package (Boost REQUIRED COMPONENTS date_time filesystem program_options system regex timer chrono)
message (STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
message (STATUS "Boost_LIBRARIES    = ${Boost_LIBRARIES}")

# --- OpenSmoke++ (Compulsory)
message (STATUS "OpenSMOKEpp_RespecthConverter = ${PROJECT_SOURCE_DIR}")

set(COMPULSORY_LIBRARIES ${Boost_LIBRARIES})

set(COMPULSORY_INCLUDE ${Boost_INCLUDE_DIRS})

#--------------------------------------------------------------------
# Compilation flags and options
#--------------------------------------------------------------------

add_compile_options(-std=c++17 -fPIC -m64 -w)
add_link_options(-std=c++17 -fPIC -m64 -w)

#--------------------------------------------------------------------
# Include flags
#--------------------------------------------------------------------

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/src
                     ${COMPULSORY_INCLUDE})

#--------------------------------------------------------------------
# Library flags
#--------------------------------------------------------------------

LINK_LIBRARIES (${COMPULSORY_LIBRARIES}
                ${SEMI_COMPULSORY_LIBRARIES}
                ${OPTIONAL_LIBRARIES})

#--------------------------------------------------------------------
# Targets to compile
#--------------------------------------------------------------------

add_executable (Respecth2OSpp.sh 
  ${PROJECT_SOURCE_DIR}/src/OpenSMOKEpp_RespecthConverter.cpp 
  ${PROJECT_SOURCE_DIR}/src/Bibliography.cpp 
  ${PROJECT_SOURCE_DIR}/src/Composition.cpp 
  ${PROJECT_SOURCE_DIR}/src/Conversions.cpp 
  ${PROJECT_SOURCE_DIR}/src/DatabaseSpecies.cpp 
  ${PROJECT_SOURCE_DIR}/src/Respecth2OpenSMOKEpp_BurnerStabilizedFlameSpeciation.cpp 
  ${PROJECT_SOURCE_DIR}/src/Respecth2OpenSMOKEpp_ConcentrationTimeProfile.cpp 
  ${PROJECT_SOURCE_DIR}/src/Respecth2OpenSMOKEpp.cpp 
  ${PROJECT_SOURCE_DIR}/src/Respecth2OpenSMOKEpp_IgnitionDelay.cpp 
  ${PROJECT_SOURCE_DIR}/src/Respecth2OpenSMOKEpp_JetStirredReactor.cpp 
  ${PROJECT_SOURCE_DIR}/src/Respecth2OpenSMOKEpp_LaminarBurningVelocity.cpp 
  ${PROJECT_SOURCE_DIR}/src/Respecth2OpenSMOKEpp_OutletConcentration.cpp 
  ${PROJECT_SOURCE_DIR}/src/Utilities.cpp
)

#--------------------------------------------------------------------
# Installing options
#--------------------------------------------------------------------

install(TARGETS Respecth2OSpp.sh DESTINATION ${INSTALL_BIN_DIR})