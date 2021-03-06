message("Checking cmake version...")                                                                # Printing message...
cmake_minimum_required(VERSION 3.0)                                                                 # Requiring cmake 3.0...

message("")                                                                                         # Printing message...
message("################################################################################")         # Printing message...
message("################################ C/C++ COMPILER ################################")         # Printing message...
message("################################################################################")         # Printing message...

message("Setting project name...")                                                                  # Printing message...
project (Cappuccino)                                                                                # Setting project name...
message("Project name = ${PROJECT_NAME}")                                                           # Printing message...

set(CMAKE_CXX_STANDARD 17)                                                                          # Setting C++17 standard...
set(CMAKE_CXX_STANDARD_REQUIRED ON)                                                                 # Setting complier standard requirement...
set(CMAKE_CXX_EXTENSIONS OFF)                                                                       # Disabling compiler-specific standards...

message("Setting build type...")                                                                    # Printing message...
set(CMAKE_BUILD_TYPE Release)                                                                       # Setting build type...
message("Build type = ${CMAKE_BUILD_TYPE}")                                                         # Printing message...

message("")                                                                                         # Printing message...
message("################################################################################")         # Printing message...
message("#################################### Spinor ####################################")         # Printing message...
message("################################################################################")         # Printing message...

if(APPLE)                                                                                           # Detecting APPLE...
  set(TARGET_1 "spinor")                                                                            # Setting executable name...
  set(DIRECTORY_1 "Code")                                                                           # Setting directory name...

  message("Adding source files for ${TARGET_1}...")                                                 # Printing message...
  aux_source_directory(${CMAKE_HOME_DIRECTORY}/${DIRECTORY_1}/src SRC_1)                            # Getting all Neutrino source files...
  set(SOURCES_1                                                                                     # Setting "SOURCES" variable...
    ${SRC_1})                                                                                       # All project source files.

  message("Adding build target as executable...")                                                   # Printing message...
  add_executable(${TARGET_1} ${SOURCES_1})                                                          # Adding executable...

  message("Adding include files...")                                                                # Printing message...
  set(INCLUDES_1                                                                                    # Setting "INCLUDES" variable...
    ${CMAKE_HOME_DIRECTORY}/${DIRECTORY_1}/include                                                  # Setting example include directory...
    ${NEUTRINO_PATH}/include)                                                                       # Setting Neutrino include directory...
  target_include_directories(                                                                       # Setting include directories...
    ${TARGET_1} PRIVATE                                                                             # Target name.
    ${INCLUDES_1})                                                                                  # All include directories.

  message("Adding linked libraries...")                                                             # Printing message...
  target_link_libraries(                                                                            # Setting other linked libraries...
    ${TARGET_1}                                                                                     # Target name.
    "-framework OpenGL"                                                                             # OpenGL library.
    "-framework OpenCL"                                                                             # OpenCL library.
    ${GLFW_PATH}/lib-macos/libglfw.3.dylib                                                          # GLFW library.
    "-lm"                                                                                           # "math" library.
    ${NEUTRINO_PATH}/lib/libnu.a)                                                                   # Neutrino library.
endif(APPLE)

if(UNIX AND NOT APPLE)                                                                              # Detecting LINUX...
  set(TARGET_1 "spinor")                                                                            # Setting executable name...
  set(DIRECTORY_1 "Code")                                                                           # Setting directory name...

  message("Adding source files for ${TARGET_1}...")                                                 # Printing message...
  aux_source_directory(${CMAKE_HOME_DIRECTORY}/${DIRECTORY_1}/src SRC_1)                            # Getting all Neutrino source files...
  set(SOURCES_1                                                                                     # Setting "SOURCES" variable...
    ${SRC_1})                                                                                       # All project source files.

  message("Adding build target as executable...")                                                   # Printing message...
  add_executable(${TARGET_1} ${SOURCES_1})                                                          # Adding executable...

  message("Adding include files...")                                                                # Printing message...
  set(INCLUDES                                                                                      # Setting "INCLUDES" variable...
    ${CMAKE_HOME_DIRECTORY}/include                                                                 # Example include directory.
    ${NEUTRINO_PATH}/include)                                                                       # Neutrino include directory.
  target_include_directories(                                                                       # Setting include directories...
    ${TARGET_1} PRIVATE                                                                             # Target name.
    ${INCLUDES_1})                                                                                  # All include directories.

  message("Adding linked libraries...")                                                             # Printing message...
  target_link_libraries(                                                                            # Setting other linked libraries...
    ${TARGET_1}                                                                                     # Target name.
    "-lOpenGL"                                                                                      # OpenGL library.
    "-lOpenCL"                                                                                      # OpenCL library.
    "-ldl"                                                                                          # "libdl" library.
    "-lglfw"                                                                                        # GLFW library.
    "-lm"                                                                                           # "math" library.
    ${NEUTRINO_PATH}/lib/libnu.a)                                                                   # "neutrino" library.
endif(UNIX AND NOT APPLE)

if(WIN32)                                                                                           # Detecting WINDOWS...
  set(TARGET_1 "spinor")                                                                            # Setting executable name...
  set(DIRECTORY_1 "Code")                                                                           # Setting directory name...

  message("Adding source files for ${TARGET_1}...")                                                 # Printing message...
  aux_source_directory(${CMAKE_HOME_DIRECTORY}/${DIRECTORY_1}/src SRC_1)                            # Getting all Neutrino source files...
  string(REPLACE "\\" "/" GLAD_PATH "${GLAD_PATH}")                                                 # Adjusting backslashes...
  string(REPLACE "\\" "/" GLFW_PATH "${GLFW_PATH}")                                                 # Adjusting backslashes...
  string(REPLACE "\\" "/" GMSH_PATH "${GMSH_PATH}")                                                 # Adjusting backslashes...
  string(REPLACE "\\" "/" CL_PATH "${CL_PATH}")                                                     # Adjusting backslashes...
  string(REPLACE "\\" "/" NEUTRINO_PATH "${NEUTRINO_PATH}")                                         # Adjusting backslashes...
  set(SOURCES_1                                                                                     # Setting "SOURCES" variable...
    ${SRC_1})                                                                                       # All project source files.

  message("Adding build target as executable...")                                                   # Printing message...
  add_executable(${TARGET_1} ${SOURCES_1})                                                          # Adding executable...

  message("Adding include files...")                                                                # Printing message...
  set(INCLUDES_1                                                                                    # Setting "INCLUDES" variable...
    ${CMAKE_HOME_DIRECTORY}/include                                                                 # Example include directory.
    ${NEUTRINO_PATH}/include)                                                                       # Neutrino include directory.
  target_include_directories(                                                                       # Setting include directories...
    ${TARGET_1} PRIVATE                                                                             # Target name.
    ${INCLUDES_1})                                                                                  # All include directories.

  message("Adding linked libraries...")                                                             # Printing message...
  target_link_libraries(                                                                            # Setting other linked libraries...
    ${TARGET_1}                                                                                     # Target name.
    ${CL_PATH}/lib/x64/OpenCL.lib                                                                   # OpenCL library.
    ${GLFW_PATH}/lib-vc2019/glfw3.lib                                                               # GLFW library.
    ${NEUTRINO_PATH}/lib/nu.lib)                                                                    # "neutrino" library.
endif(WIN32)

message("Setting build directory...")                                                               # Printing message...
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/build)                                   # Setting build directory...
message(${CMAKE_HOME_DIRECTORY}/build)                                                              # Printing message...

message("DONE!")                                                                                    # Printing message...

message("")                                                                                         # Printing message...
message("################################################################################")         # Printing message...
message("################################# INSTRUCTIONS #################################")         # Printing message...
message("################################################################################")         # Printing message...

message("Now go to the \"build\" directory created by this command and do one (or more)")           # Printing message...
message("of the following things:")                                                                 # Printing message...
message("")                                                                                         # Printing message...
message("1. Type: \"make\" or \"make all\" in order to build all example executables.")             # Printing message...
message("2. Type: \"make <example_name>\" in order to build only a specific example executable.")   # Printing message...
message("   The name of each example is the lowercase name of its corresponding directory:")        # Printing message...
message("   e.g. EXAMPLE = Sinusoid --> EXECUTABLE = sinusoid")                                     # Printing message...
message("        make sinusoid")                                                                    # Printing message...
message("3. Type: \"make doc\" in order to build the Doxygen documentation of the project.")        # Printing message...
message("")                                                                                         # Printing message...
message("################################################################################")         # Printing message...
message("############################# CONFIGURATION REPORT #############################")         # Printing message...
message("################################################################################")         # Printing message...
