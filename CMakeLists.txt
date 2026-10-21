cmake_minimum_required(VERSION 3.20)
project(nuflow VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)          # ROOT 6.28 in the spack stack is built with C++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
add_compile_options(-Wall -Wextra -Wno-unused-parameter)

option(NUFLOW_WITH_HEPMC3 "Build the NuHepMC (HepMC3) output writer" ON)
option(NUFLOW_BUILD_TESTS "Build the Catch2 test suite" ON)

find_package(ROOT 6 REQUIRED COMPONENTS Core RIO Tree Hist)
find_package(yaml-cpp REQUIRED)
find_package(Threads REQUIRED)

add_library(nuflow SHARED
  core/Config.cc
  core/Nucleus.cc
  flux/FluxSource.cc
  flux/FluxReaders.cc
  initial_state/NuclearModel.cc
  initial_state/RelativisticFermiGas.cc
  initial_state/FreeNucleon.cc
  kinematics/FormFactors.cc
  kinematics/QECurrent.cc
  kinematics/QuasiElasticKinematics.cc
  kinematics/PionProductionCurrent.cc
  kinematics/ResonantPionKinematics.cc
  channel/CCQEChannel.cc
  channel/CCRESChannel.cc
  channel/ChannelSelector.cc
  fsi/TransparentNucleus.cc
  orchestration/ModelFactory.cc
  orchestration/XSecTableBuilder.cc
  orchestration/EnergySampler.cc
  orchestration/Generator.cc
  orchestration/RootTreeWriter.cc
)
target_include_directories(nuflow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(nuflow PUBLIC NUFLOW_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(nuflow PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist yaml-cpp::yaml-cpp Threads::Threads)

if(NUFLOW_WITH_HEPMC3)
  find_package(HepMC3 REQUIRED)
  target_sources(nuflow PRIVATE orchestration/NuHepMCWriter.cc)
  target_compile_definitions(nuflow PUBLIC NUFLOW_HAS_HEPMC3)
  target_include_directories(nuflow PUBLIC ${HEPMC3_INCLUDE_DIR})
  target_link_libraries(nuflow PUBLIC ${HEPMC3_LIBRARIES})
endif()

add_executable(nuflow-gen apps/nuflow-gen.cc)
target_link_libraries(nuflow-gen PRIVATE nuflow)
add_executable(nuflow-xsec apps/nuflow-xsec.cc)
target_link_libraries(nuflow-xsec PRIVATE nuflow)
add_executable(nuflow-qe-table apps/nuflow-qe-table.cc)
target_link_libraries(nuflow-qe-table PRIVATE nuflow)

if(NUFLOW_BUILD_TESTS)
  find_package(Catch2 3 REQUIRED)
  enable_testing()
  add_executable(nuflow-tests
    tests/test_formfactors.cc
    tests/test_qe_current.cc
    tests/test_nuclear.cc
    tests/test_flux.cc
    tests/test_generator.cc
    tests/test_pion_production.cc
  )
  target_link_libraries(nuflow-tests PRIVATE nuflow Catch2::Catch2WithMain)
  include(Catch)
  catch_discover_tests(nuflow-tests)
endif()
