cmake_minimum_required(VERSION 3.16)
project(OCRT C)
file(GLOB_RECURSE OCRT_SRC CONFIGURE_DEPENDS src/*.c)
add_executable(ocrt ${OCRT_SRC})
target_include_directories(ocrt PRIVATE src)
target_compile_options(ocrt PRIVATE -std=c11 -O3 -ffp-contract=fast -fassociative-math -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp)
target_link_libraries(ocrt m gomp)
enable_testing()
if(EXISTS ${CMAKE_SOURCE_DIR}/tests/test_fr631_mie_reader.c)
  add_executable(test_fr631_mie_reader
    tests/test_fr631_mie_reader.c
    src/shared/mie_io.c
    src/shared/io_utils.c
    src/shared/numerics.c)
  target_include_directories(test_fr631_mie_reader PRIVATE src)
  target_link_libraries(test_fr631_mie_reader m)
  add_test(NAME fr631_mie_reader COMMAND test_fr631_mie_reader ${CMAKE_SOURCE_DIR}/inputs/tsm_ahn/Red_clay_AHN.mie)
endif()
add_test(NAME ocrt_nadir_smoke COMMAND ocrt --surface ocean --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 --sza 0 --vza 0.001 --raa 0 --wind-speed 0 --wavelength 555 --pressure 0 --aod 0 --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 --n-mu-water 16 --n-layers 60 --m-max 2 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(ocrt_nadir_smoke PROPERTIES ENVIRONMENT "OCRT_ADVANCED=1")
