cmake_minimum_required(VERSION 3.10)
project(SMVacuumDecay LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(smvacuumdecay
  src/RGE.cpp
  src/EffectivePotential.cpp
  src/CanonicalBounce.cpp
  src/FubiniLipatov.cpp
)
target_include_directories(smvacuumdecay PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(benchmark_point apps/benchmark_point.cpp)
target_link_libraries(benchmark_point PRIVATE smvacuumdecay)

add_executable(generate_phase_diagram apps/generate_phase_diagram.cpp)
target_link_libraries(generate_phase_diagram PRIVATE smvacuumdecay)

add_executable(write_reference_tables apps/write_reference_tables.cpp)
target_link_libraries(write_reference_tables PRIVATE smvacuumdecay)

add_executable(convergence_scan apps/convergence_scan.cpp)
target_link_libraries(convergence_scan PRIVATE smvacuumdecay)

enable_testing()

# All tests run from the source root: the regression tests read the frozen
# tables under reference/ via relative paths.
set(REGRESSION_TESTS test_rge test_potential test_FL test_bounce test_scientific)
foreach(t ${REGRESSION_TESTS})
  add_executable(${t} tests/regression/${t}.cpp)
  target_link_libraries(${t} PRIVATE smvacuumdecay)
  add_test(NAME regression/${t} COMMAND ${t}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()

set(UNIT_TESTS test_numerics test_matching test_threshold)
foreach(t ${UNIT_TESTS})
  add_executable(${t} tests/unit/${t}.cpp)
  target_link_libraries(${t} PRIVATE smvacuumdecay)
  add_test(NAME unit/${t} COMMAND ${t}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()

set(PHYSICS_TESTS test_fubini_limit test_convergence test_boundary)
foreach(t ${PHYSICS_TESTS})
  add_executable(${t} tests/physics/${t}.cpp)
  target_link_libraries(${t} PRIVATE smvacuumdecay)
  add_test(NAME physics/${t} COMMAND ${t}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
