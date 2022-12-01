# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(sCVD)

if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(sCVD main.cpp)
target_link_libraries(sCVD Garfield::Garfield)
