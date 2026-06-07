#ifndef PARSER_HPP
#define PARSER_HPP
#include <iostream>
#include <fstream>
#include <string>

/*
    TODO:
    1. Create a function that reads a FASTA file and stores the genome as a single string so it can be used for comparison :)

    COMPLETED: YES
*/

/*
    Function takes a file path, skips headers (">") and returns a concatenated string.
    Written by: Ivan Dundović
*/
std::string read_fasta(std::string path){
    std::ifstream file(path);
    if(!file.is_open()) throw std::runtime_error("Cannot open the file: " + path);
    std::string genom, line;
    while (std::getline(file, line)) if(!line.empty() && line[0] != '>') genom += line;

    return genom;
}

#endif