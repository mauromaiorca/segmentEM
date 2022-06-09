#
# mauro Subtomogram Average Tools
#
#

PROJECT(mauroSubtomogram)
 cmake_minimum_required(VERSION 2.6)

SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread")

SET(BUILD_SHARED_LIBS OFF)
SET(ITK_DYNAMIC_LOADING OFF)

SET(SRC_DIR "./src/")

IF(UNIX AND NOT APPLE)
        SET(LIB_DIR "./libs/")
	SET(BIN_DIR "./bin/")
ENDIF(UNIX AND NOT APPLE)

IF(APPLE)
	SET(BIN_DIR "./bin_macos/")
        SET(LIB_DIR "./libs_mac/")
ENDIF(APPLE)

IF(WIN32)
	SET(BIN_DIR "./bin_win32/")
ENDIF(WIN32)

#SET( GCC_COVERAGE_COMPILE_FLAGS  "${GCC_COVERAGE_COMPILE_FLAGS} -O2" )
SET( GCC_COVERAGE_COMPILE_FLAGS  "${GCC_COVERAGE_COMPILE_FLAGS} -Wc++11-extensions" )
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

#find_package(VTK REQUIRED)
#include(${VTK_USE_FILE})
# find_package(LAPACK REQUIRED)
#find_package(LAPACK REQUIRED)






find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
include(${CMAKE_SOURCE_DIR}/cmake/FindFFTW.cmake)

#find_package(fftw3 REQUIRED)
#include(${fftw3})
#SET(FFTW_INCLUDE_DIR "/media/sf_virtualMachineSharedFolder/AmplitudeStructureOfImages/code_production/c++/bin_linux/fftw/include/")
#SET(FFTW_LIBRARIES /media/sf_virtualMachineSharedFolder/AmplitudeStructureOfImages/code_production/c++/bin_linux/fftw/lib/libfftw3.a)

INCLUDE_DIRECTORIES(
   ${SRC_DIR}
   ${FFTW_INCLUDE_DIR}
)



ADD_EXECUTABLE(${BIN_DIR}selectParticles ${SRC_DIR}selectParticles.cxx) #it was refineParticles
TARGET_LINK_LIBRARIES(${BIN_DIR}selectParticles ${FFTW_LIBRARIES}  ${ITK_LIBRARIES})

ADD_EXECUTABLE(${BIN_DIR}starProcess ${SRC_DIR}starProcess.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}starProcess ${FFTW_LIBRARIES}  ${ITK_LIBRARIES})

ADD_EXECUTABLE(${BIN_DIR}starProcess2 ${SRC_DIR}starProcess2.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}starProcess2 ${FFTW_LIBRARIES}  ${ITK_LIBRARIES})


ADD_EXECUTABLE(${BIN_DIR}starSimulation ${SRC_DIR}starSimulation.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}starSimulation ${FFTW_LIBRARIES}  ${ITK_LIBRARIES})


ADD_EXECUTABLE(${BIN_DIR}meanMinMax ${SRC_DIR}meanMinMax.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}meanMinMax ${FFTW_LIBRARIES}  ${ITK_LIBRARIES})



ADD_EXECUTABLE(${BIN_DIR}stackProcess ${SRC_DIR}stackProcess.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}stackProcess ${FFTW_LIBRARIES}  ${ITK_LIBRARIES})

ADD_EXECUTABLE(${BIN_DIR}mapProcess ${SRC_DIR}mapProcess.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}mapProcess ${FFTW_LIBRARIES}  ${ITK_LIBRARIES})

ADD_EXECUTABLE(${BIN_DIR}assessParticles ${SRC_DIR}assessParticles.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}assessParticles ${FFTW_LIBRARIES}  ${ITK_LIBRARIES} )

ADD_EXECUTABLE(${BIN_DIR}mrcTools ${SRC_DIR}mrcTools.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}mrcTools ${FFTW_LIBRARIES}  ${ITK_LIBRARIES} )


#ADD_EXECUTABLE(${BIN_DIR}morphologocalClosing ${SRC_DIR}morphologocalClosing.cxx)
#TARGET_LINK_LIBRARIES(${BIN_DIR}morphologocalClosing ${FFTW_LIBRARIES}  ${ITK_LIBRARIES} )


ADD_EXECUTABLE(${BIN_DIR}maskProcess ${SRC_DIR}maskProcess.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}maskProcess ${FFTW_LIBRARIES}  ${ITK_LIBRARIES} )


ADD_EXECUTABLE(${BIN_DIR}emTasks ${SRC_DIR}emTasks.cxx)
TARGET_LINK_LIBRARIES(${BIN_DIR}emTasks ${FFTW_LIBRARIES}  ${ITK_LIBRARIES} )

#ADD_EXECUTABLE(${BIN_DIR}testImage ${SRC_DIR}test.cxx)
#TARGET_LINK_LIBRARIES(${BIN_DIR}testImage ${FFTW_LIBRARIES}  ${ITK_LIBRARIES} )

#ADD_EXECUTABLE(${BIN_DIR}interlaced ${SRC_DIR}interlaced.cxx)
#TARGET_LINK_LIBRARIES(${BIN_DIR}interlaced ${FFTW_LIBRARIES}  ${ITK_LIBRARIES} )


################################################
#
#  Install into the system
#
################################################
IF(UNIX)
#	SET (CMAKE_INSTALL_PREFIX ${DESTDIR}) #in order to support the make install option
#	INSTALL(TARGETS mrc2mhd mhd2mrc RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
ENDIF(UNIX)
