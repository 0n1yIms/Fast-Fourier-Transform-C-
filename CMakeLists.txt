cmake_minimum_required(VERSION 3.6.0)


set(Halide_DIR ./Halide/lib/cmake/Halide)
set(HalideHelpers_DIR ./Halide/lib/cmake/HalideHelpers)
find_package(Halide REQUIRED)


add_executable(main main.cpp)

target_include_directories(main PUBLIC
        ${Halide_INCLUDE}
        )
target_link_libraries(main
      Halide::Halide)