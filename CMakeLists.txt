cmake_minimum_required(VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (lz-rlbwt)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/internal)
include_directories(~/include) #SDSL headers are here
include_directories(${PROJECT_SOURCE_DIR}/extern/BWTIL)
include_directories(${PROJECT_SOURCE_DIR}/extern/BWTIL/extern/bitvector/include)

LINK_DIRECTORIES(~/lib) #SDSL lib are here

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

add_executable(debug debug.cpp)
TARGET_LINK_LIBRARIES(debug sdsl)
TARGET_LINK_LIBRARIES(debug divsufsort)
TARGET_LINK_LIBRARIES(debug divsufsort64)

add_executable(lz-rlbwt-build lz-rlbwt-build.cpp)
TARGET_LINK_LIBRARIES(lz-rlbwt-build sdsl)
TARGET_LINK_LIBRARIES(lz-rlbwt-build divsufsort)
TARGET_LINK_LIBRARIES(lz-rlbwt-build divsufsort64)

add_executable(lz-rlbwt-search lz-rlbwt-search.cpp)
TARGET_LINK_LIBRARIES(lz-rlbwt-search sdsl)
TARGET_LINK_LIBRARIES(lz-rlbwt-search divsufsort)
TARGET_LINK_LIBRARIES(lz-rlbwt-search divsufsort64)