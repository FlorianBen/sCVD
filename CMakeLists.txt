# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(sCVD)

if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(sCVD_alpha_trim main_alpha_trim.cpp)
target_link_libraries(sCVD_alpha_trim Garfield::Garfield)

add_executable(sCVD_alpha_srim main_alpha_srim.cpp)
target_link_libraries(sCVD_alpha_srim Garfield::Garfield)

add_executable(sCVD_mip main_mip.cpp)
target_link_libraries(sCVD_mip Garfield::Garfield)

add_executable(sCVD_gamma main_gamma.cpp)
target_link_libraries(sCVD_gamma Garfield::Garfield)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
