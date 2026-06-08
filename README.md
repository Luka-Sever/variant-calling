# Varaint Caller

This repository has C++ source code for implementation of variant calling over single-nucleotide polymorphism (SNP). It is a CLI tool that is used on Linux based systems.
This project is part of subject on Faculty of Electrical Engineering and Computing, [Bioninformatics 1](https://www.fer.unizg.hr/en/course/enbio1)

## Features
 - Time profiling: each phase has time it took for each phasd to complete.
 - Memory profiling: for memory profiling, use `/usr/bin/time`  with parameter `-v` in front of main program.

## Prerequisits
For this project you need:
 - Cmake, => 3.14
 - g++ compiler, which supports C++17 standard
 - all required libraries, run:
```
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    cmake \
    libhts-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev
```

Optionally, you should have minimap2 installed if you want to create SAM files

## Installing
Clone this repo and run following commands:
```
mkdir build && cd build  
cmake ..
cmake --build .
```
To run a program, inside created build directory, run it as:
```
./variant_caller <path to FASTA file> <path to SAM file> <path to save>
```

## Example
```
./variant_caller ../data/lambda.fasta ../data/alignment.sam ../data/output.csv 
```
Inside data folder, you will find your output. <br>
If you want to get memory usage, run:
```
/usr/bin/time -v ./variant_caller ../data/lambda.fasta ../data/alignment.sam ../data/output.csv 
```

# Contirbutors
Luka Sever <br>
Ivan Dundović
