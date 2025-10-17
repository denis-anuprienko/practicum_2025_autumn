cmake_minimum_required(VERSION 2.8)
project(inmost_project)

#set(INMOST_DIR %%%%%) # uncomment and set location manually if needed
find_package(inmost REQUIRED)

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

add_executable(diffusion_fem diffusion_fem.cpp)

target_link_libraries(diffusion_fem ${INMOST_LIBRARIES})
