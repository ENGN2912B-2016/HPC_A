cmake_minimum_required (VERSION 2.8.12.2)

project (Ternary_serial_class)

#set(CMAKE_BUILD_TYPE Debug)
set(CXX "/gpfs/runtime/opt/gcc/5.2.0/bin/gcc")
set(CMAKE_CXX_FLAGS "-lfftw3 -lm -std=c++11 -O0  -g")
INCLUDE_DIRECTORIES(/gpfs/runtime/opt/fftw/3.3.4/*)
LINK_DIRECTORIES(/gpfs/runtime/opt/fftw/3.3.4/*)

# define library location
#set_target_properties(Ternary_serial PROPERTIES LINKER_LANGUAGE CXX)

# create the executable file
add_executable (ternary_serial ternary_serial.cpp)
#target_link_libraries (ternary_serial Ternary_Serial)
