cmake_minimum_required(VERSION 3.10)

# Set the project name
project(MyProject)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB INCLUDE_DIR ${PROJECT_SOURCE_DIR}/dependencies/include)
file(GLOB TEST_DIR ${PROJECT_SOURCE_DIR}/dependencies/Test)

include_directories(${INCLUDE_DIR} ${TEST_DIR} )

# Add source files
file(GLOB LIB_SOURCES "dependencies/lib/*.cpp")
file(GLOB SRC_SOURCES "src/*.cpp")
file(GLOB TEST_SOURCES "Test/*.cpp")

# Create an executable from the source files
add_executable(MyProject ${SRC_SOURCES} ${LIB_SOURCES} ${TEST_SOURCES})

# Link libraries (if you have any precompiled libraries to link, specify them here)
# target_link_libraries(MyProject <library_name>)