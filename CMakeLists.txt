"
    Uputa za pokretanje:
    u korijenskom direktoriju pokrenuti:

    mkdir build && cd build
    cmake ..
    cmake --build .
"

// verzija CMAKE-a i potrebni jezični procesori
cmake_minimum_required(VERSION 3.14)
project(VariantCaller C CXX)

// standard C++-a
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

// instaliranje potrebnog libraryja
include(FetchContent)
FetchContent_MakeAvailable(htslib)

// početna lokacija
add_executable(variant_caller
    main.cpp
)

// linkanje
target_link_libraries(variant_caller PRIVATE hts)
target_include_directories(variant_caller PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})