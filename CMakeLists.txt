# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

project (deep)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_subdirectory(helpers)

add_executable(gen_eep gen_eep.cpp spWf.cpp sigmaep_bound.cpp)
target_link_libraries(gen_eep ${ROOT_LIBRARIES} helpers)

#add_executable(test_cc1_cc2 test_cc1_cc2.cpp sigmaep_bound.cpp)
#target_link_libraries(test_cc1_cc2 ${ROOT_LIBRARIES} helpers)

#add_executable(digitize digitize.cpp sigmaep_bound.cpp)
#target_link_libraries(digitize ${ROOT_LIBRARIES} helpers)
