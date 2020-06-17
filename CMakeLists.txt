cmake_minimum_required(VERSION 2.8)

project(Gcare CXX)


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

project (Gcare)

add_executable(gcare_graph ./src/main.cc ./src/util.cc ./src/data_graph.cc ./src/query_graph.cc ./src/wander_join.cc ./src/cset.cc ./src/sumrdf.cc ./src/jsub.cc ./src/impr.cc)
set_target_properties(gcare_graph PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(gcare_graph PRIVATE ./include)
target_compile_options(gcare_graph PRIVATE -g -O2 -std=c++11 -Bstatic)
target_link_libraries(gcare_graph boost_program_options stdc++fs)

add_executable(gcare_relation ./src/main.cc ./src/util.cc ./src/ndvector.cc ./src/data_relations.cc ./src/query_relations.cc ./src/correlated_sampling.cc ./src/bound_sketch.cc)
set_target_properties(gcare_relation PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(gcare_relation PRIVATE ./include)
target_compile_options(gcare_relation PRIVATE -g -O2 -std=c++11 -Bstatic)
target_compile_definitions(gcare_relation PRIVATE -DRELATION)
target_link_libraries(gcare_relation boost_program_options stdc++fs)
