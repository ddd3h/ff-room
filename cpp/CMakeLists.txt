cmake_minimum_required(VERSION 3.26...4.0)
project(ffroom_core CXX)

find_package(OpenMP)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release build flags
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
if(MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(-Wall -Wextra)
endif()

# Use the Python specified by -DPython3_EXECUTABLE (conda/venv-aware).
# Must be set before find_package(Python3) and before pybind11.
if(DEFINED Python3_EXECUTABLE)
  set(Python3_FIND_STRATEGY LOCATION)
endif()
set(PYBIND11_FINDPYTHON ON)
find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development.Module)
message(STATUS "Python3: ${Python3_EXECUTABLE}")
message(STATUS "Python3 include: ${Python3_INCLUDE_DIRS}")

# pybind11: try pip-installed first, fall back to FetchContent
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_VARIABLE pybind11_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE pybind11_pip_result
  ERROR_QUIET
)
if(pybind11_pip_result EQUAL 0 AND pybind11_DIR)
  find_package(pybind11 3.1 CONFIG QUIET HINTS "${pybind11_DIR}")
endif()
if(NOT pybind11_FOUND)
  find_package(pybind11 3.1 CONFIG QUIET)
endif()
if(NOT pybind11_FOUND)
  message(STATUS "pybind11 not found locally; using FetchContent")
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v3.1.0
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(pybind11)
endif()

set(CORE_SOURCES
    src/core/Grid.cpp
    src/core/FanModel.cpp
    src/core/PoissonSolver.cpp
    src/core/BoundaryManager.cpp
    src/core/FluidSolver.cpp
    src/core/LBMSolver.cpp
)

# Pybind11 extension module
pybind11_add_module(_ffroom_core
    ${CORE_SOURCES}
    src/bindings/module.cpp
)

target_include_directories(_ffroom_core PRIVATE src)

if(OpenMP_CXX_FOUND)
    target_link_libraries(_ffroom_core PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(_ffroom_core PRIVATE FFROOM_HAS_OPENMP)
    message(STATUS "OpenMP found — enabling parallelization")
else()
    message(STATUS "OpenMP not found — building without parallelization")
endif()

# Wheels install relative to their staging prefix; build.sh supplies the package
# directory as its prefix for standalone builds.
if(SKBUILD)
  set(FFROOM_INSTALL_DIR ff_room)
else()
  set(FFROOM_INSTALL_DIR .)
endif()
install(TARGETS _ffroom_core
        LIBRARY DESTINATION ${FFROOM_INSTALL_DIR} COMPONENT python
        RUNTIME DESTINATION ${FFROOM_INSTALL_DIR} COMPONENT python)

# Unit tests (optional, requires catch2)
option(BUILD_TESTS "Build C++ unit tests" OFF)
if(BUILD_TESTS)
  enable_testing()
  find_package(Catch2 3 QUIET)
  if(Catch2_FOUND)
    add_executable(ffroom_tests
        ${CORE_SOURCES}
        tests/test_grid.cpp
        tests/test_poisson.cpp
        tests/test_couette.cpp
    )
    target_include_directories(ffroom_tests PRIVATE src)
    target_link_libraries(ffroom_tests PRIVATE Catch2::Catch2WithMain)
    add_test(NAME ffroom_tests COMMAND ffroom_tests)
  else()
    message(WARNING "Catch2 not found; skipping C++ tests")
  endif()
endif()
