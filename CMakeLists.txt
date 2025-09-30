cmake_minimum_required(VERSION 3.16)
project(SPHinXsysMedicalEngineering VERSION 1.0.0 LANGUAGES CXX)

set(SPHinXsysMedicalEngineering_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# ------ include the SPHinXsys library ------
add_subdirectory(sphinxsys_lib)

# ------ include the source code ------
add_subdirectory(src)

enable_testing()
add_subdirectory(tests)