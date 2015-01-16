# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#SET (CMAKE_CXX_COMPILER             "/usr/bin/clang++")
SET (CMAKE_CXX_COMPILER             "/usr/bin/g++")
set(BUILD_SHARED_LIBS OFF)

SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)

PROJECT(pssc CXX C)

FILE(GLOB pssc_src src/pssc_config.cpp
    "src/*.h"
    "src/*.cpp"
)

set(USE_STATIC_BOOST ON)

set(Boost_USE_STATIC_LIBS ${USE_STATIC_BOOST})
set(Boost_USE_MULTITHREADED ON)
SET (Boost_USE_STATIC_RUNTIME OFF)
unset(Boost_INCLUDE_DIR CACHE)
unset(Boost_LIBRARY_DIRS CACHE)



SET( CMAKE_MODULE_PATH ${cpppdbparser2_SOURCE_DIR} )

#SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -iquote src -g -Wall -Wno-multichar -pthread -O0 -std=c++11")


#SET ( Boost_USE_STATIC_LIBS ON )
#SET ( Boost_USE_MULTITHREADED ON )
#SET ( Boost_USE_STATIC_RUNTIME ON )
#SET (BOOST_ALL_DYN_LINK           OFF)   # NOT force dynamic linking for all libraries
#find_package( ZLIB REQUIRED /usr/lib/x86_64-linux-gnu/libbz2.a)
#find_library( ZLIB  /usr/lib/x86_64-linux-gnu/libbz2.a)
find_package( ZLIB REQUIRED )
find_package( BZip2 REQUIRED )  

add_library(ZLIB STATIC IMPORTED)
add_library(BZip2 STATIC IMPORTED)

  set_target_properties(ZLIB PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(ZLIB PROPERTIES LINK_SEARCH_END_STATIC 1)
 set_property(TARGET ZLIB PROPERTY IMPORTED_LOCATION /usr/lib/x86_64-linux-gnu/libbz2.a)
 
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/include )
#set(INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/include )
message(INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/include)
FIND_PACKAGE ( Boost COMPONENTS thread filesystem program_options iostreams system date_time )

FIND_PACKAGE ( Irrlicht )

if(USE_STATIC_BOOST)
   message(STATUS "Linking against boost static libraries")
else(USE_STATIC_BOOST)
   message(STATUS "Linking against boost dynamic libraries")
endif(USE_STATIC_BOOST)

if(Boost_FOUND)
  INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )
  link_directories(${Boost_LIBRARY_DIRS})
  message(" ${Boost_LIBRARIES}\n")
else()
  message(FATAL_ERROR "You need boost to build this, set the BOOST_ROOT or BOOSTROOT env variables, or pass them to cmake")
endif()


# for debugging/profiling:
#SET ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -O0 -pg -static-libstdc++ -Bstatic -static-libgcc")
#SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -iquote src -pg -Wall -Wno-multichar -pthread -g -O0 -std=c++11 -Bstatic ")

SET ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -static-libstdc++")
SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -iquote src -Wall -Wno-multichar -pthread -O3 -std=c++11 ")


message("CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}\n")
message("CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}\n")

ADD_EXECUTABLE ( pssc ${pssc_src} )
TARGET_LINK_LIBRARIES ( pssc ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})


