cmake_minimum_required(VERSION 3.10)
project(zfp)

set(CMAKE_CXX_STANDARD 14)

#set(ZFP_DIR $ENV{ZFP_DIR})
#set(ZFP_LIB $ENV{ZFP_LIB})
set(ZFP_DIR /home/majidrp/Software/zfp-0.5.5)
set(ZFP_LIB ${ZFP_DIR}/build/lib/libzfp.so)

if(NOT DEFINED ZFP_DIR)
    message( FATAL_ERROR "ZFP_DIR is not set!" )
endif(NOT DEFINED ZFP_DIR)

if(NOT DEFINED ZFP_LIB)
    message( FATAL_ERROR "ZFP_LIB is not set!" )
endif(NOT DEFINED ZFP_LIB)

include_directories(
        ${ZFP_DIR}/array
        ${ZFP_DIR}/array/zfp
        ${ZFP_DIR}/include
)

add_executable(zfp main.cpp)
target_link_libraries(zfp ${ZFP_LIB})