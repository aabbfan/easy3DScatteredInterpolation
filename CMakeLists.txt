cmake_minimum_required(VERSION 3.8)
project(easy3DScatteredInterpolation VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE release)
if(ENABLE_DEBUG)
  add_definitions(-DDEBUG_OUTPUT)
endif()

find_package(OpenMP REQUIRED)

set(LIBNABO_BUILD_PYTHON false)
set(LIBNABO_BUILD_TESTS false)
set(LIBNABO_BUILD_EXAMPLES false)
add_subdirectory(./libnabo)

set(SOURCE_FILES 3DScatteredInterpolate.cpp
                 main.cpp
)

add_executable(main ${SOURCE_FILES})

target_link_libraries(main nabo)