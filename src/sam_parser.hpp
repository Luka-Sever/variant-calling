#ifndef SAM_PARSER_HPP
#define SAM_PARSER_HPP

#include <iostream>
#include <string>
#include <htslib/sam.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include "types.hpp"

/**
 * We convert the SAM file into a BAM file to reduce memory usage.
 * Written by: Luka Sever
 */
void convert_sam_to_bam(const char* inputFile, const char* outputFile) {
    /*
        We attempt to read the input file and save it to the output file.
    */
    samFile* file = sam_open(inputFile, "r");
    if (!file) std::cout << "Cannot open the file.\n";

    sam_hdr_t *header = sam_hdr_read(file);
    if (!header) std::cout << "Cannot read the header.\n";

    samFile *out = sam_open(outputFile, "wb");
    if (!out) std::cout << "Cannot make the BAM file header.\n";

    if (sam_hdr_write(out, header) < 0) 
        std::cout << "Cannot write in the BAM file:.\n";


    // Empty BAM1 structure
    bam1_t *aln = bam_init1();

    // Reading the file
    while (sam_read1(file, header, aln) >= 0) {
        if (sam_write1(out, header, aln) < 0) {
            std::cerr << "Error writing record.\n";
            break;
        }
    }

    /*
        Memory cleanup
    */

    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(file);
    sam_close(out);

    std::cout << "Conversion completed successfully!" << std::endl;
}
/*
    Helper function: reverse_complement swapping base pairs
    Written by: Ivan Dundović
*/
std::string reverse_complement(std::string seq) {
    std::map<char,char> comp = {{'A','T'},{'T','A'},{'C','G'},{'G','C'}};
    std::string temp;
    for (char c : seq) temp.push_back(comp[c]);
    std::reverse(temp.begin(), temp.end());
    return temp;
}
/*
    SAM parser
    Written by: Ivan Dundović
*/
std::vector<SamRecord> parse_sam(std::string filename){
    std::vector<SamRecord> records;
    std::ifstream file(filename);
    if(!file.is_open()) throw std::runtime_error("Cannot open the file:  " + filename);
    std::string line;
    while (getline(file, line)){
        if(line[0] == '@') continue;
        // Splitting lines into columns
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;
        while (ss >> token){
            columns.push_back(token);
        }
        // Saving selected columns into SAM records
        SamRecord record;
        record.flag = std::stoi(columns[1]);
        record.pos = std::stoi(columns[3]) - 1; // SAM is 1-based, so we subtract 1 to convert to 0-based indexing
        record.cigar = columns[5];
        record.seq = columns[9];

        if(record.flag & 4) continue;

        if(record.cigar == "*") continue;

        if(record.flag & 16) record.seq = reverse_complement(record.seq);
        records.push_back(record);
    }
    return records;
}
#endif