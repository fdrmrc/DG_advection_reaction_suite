##
#  CMake script for your program:
##


# Set the name of the project and target:
SET(TARGET "DG_upwind")

# Declare all source files the target consists of:
SET(TARGET_SRC
  #source/base_problem.cc
  #source/base_block_problem.cc
  #source/stokes.cc
  source/DG_upwind.cc
  # include/DG_upwind.h
  main.cc
  # You can specify additional files here!
)

PROJECT(${TARGET} CXX)



# Declare all source files the target consists of. Here, this is only
# the one step-X.cc file, but as you expand your project you may wish
# to add other source files as well. If your project becomes much larger,
# you may want to either replace the following statement by something like
#  FILE(GLOB_RECURSE TARGET_SRC  "source/*.cc")
#  FILE(GLOB_RECURSE TARGET_INC  "include/*.h")
#  SET(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# or switch altogether to the large project CMakeLists.txt file discussed
# in the "CMake in user projects" page accessible from the "User info"
# page of the documentation.
#SET(TARGET_SRC
#  ${TARGET}.cc
#  )

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

FIND_PACKAGE(deal.II 9.3.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()


###############

PROJECT(gtest)

FIND_PACKAGE(deal.II 9.3.0 REQUIRED
HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})


DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Fix warning on mac
SET(CMAKE_MACOSX_RPATH 1)

# Executable for the exerciese
ADD_EXECUTABLE(poisson source/dg_upwind.cc main.cc)
DEAL_II_SETUP_TARGET(poisson)

# Library of the executable
ADD_LIBRARY(dg_upwind-lib SHARED $<TARGET_OBJECTS:poisson>)
DEAL_II_SETUP_TARGET(dg_upwind-lib)

# Tester executable
FIND_PACKAGE(GTest)
FILE(GLOB test_files tests/*cc)
ADD_EXECUTABLE(gtest ${test_files})
TARGET_LINK_LIBRARIES(gtest ${GTEST_LIBRARY} dg_upwind-lib)   
DEAL_II_SETUP_TARGET(gtest)

INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS} ./include/)

ENABLE_TESTING()
GTEST_DISCOVER_TESTS(gtest)
