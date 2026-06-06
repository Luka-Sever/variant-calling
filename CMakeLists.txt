# verzija CMAKE-a i potrebni jezični procesori
cmake_minimum_required(VERSION 3.14)
project(VariantCaller C CXX)

# standard C++-a
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# instaliranje potrebnog libraryja
include(FetchContent)
FetchContent_Declare(
    htslib_download
    GIT_REPOSITORY https://github.com/samtools/htslib.git
    GIT_TAG 1.19.1
)

FetchContent_MakeAvailable(htslib_download)
execute_process(
    COMMAND ./configure --disable-bz2 --disable-lzma --disable-libcurl --with-libdeflate=no
    WORKING_DIRECTORY ${htslib_download_SOURCE_DIR}
)
execute_process(
    COMMAND make -j4
    WORKING_DIRECTORY ${htslib_download_SOURCE_DIR}
)


# početna lokacija
add_executable(variant_caller
    src/main.cpp
)

# includamo direktorij
target_include_directories(variant_caller PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${htslib_download_SOURCE_DIR}
)

# linkamo direktorij
target_link_libraries(variant_caller PRIVATE 
    ${htslib_download_SOURCE_DIR}/libhts.a
    z
    bz2
    lzma
    curl
    pthread
)