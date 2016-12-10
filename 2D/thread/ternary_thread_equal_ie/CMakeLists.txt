cmake_minimum_required (VERSION 2.8.12.2)

project (Ternary_thread_class)

#set(CMAKE_BUILD_TYPE Debug)
set(CXX "/gpfs/runtime/opt/gcc/5.2.0/bin/gcc")
set(CMAKE_CXX_FLAGS "-lfftw3_threads -lm -std=c++11 -O0  -g")
INCLUDE_DIRECTORIES(/gpfs/runtime/opt/fftw/3.3.4/*)
LINK_DIRECTORIES(/gpfs/runtime/opt/fftw/3.3.4/*)

# define library location
#set_target_properties(Ternary_thread PROPERTIES LINKER_LANGUAGE CXX)

# create the executable file
add_executable (ternary_thread ternary_thread.cpp)
