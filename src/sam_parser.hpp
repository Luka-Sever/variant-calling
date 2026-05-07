#ifndef SAM_PARSER_HPP
#define SAM_PARSER_HPP

#include <iostream>
#include <string>
#include <htslib/sam.h>
#include <types.hpp>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

/**
 * Pretvaramo sam datoteku u bam datoteku da samnjimo
 * potrošenu memoriju
 */

void convert_sam_to_bam(const char* inputFile, const char* outputFile) {
    /*
    
        Pokušavamo čitati input file i spremati
        u output file
    
    */
    samFile* file = sam_open(inputFile, "r");
    if (!file) cout << "Ne mogu otvoriti datoteku.\n";

    sam_hdr_t *header = sam_hdr_read(file);
    if (!header) cout << "Ne mogu pročitati header.\n";

    samFile *out = sam_open(outputFile, "wb");
    if (!out) cout << "Ne mogu stvoriti bam datoteku header.\n";

    if (sam_hdr_write(out, header) < 0) 
        cout << "Ne mogu pisati u bam datoteku.\n";


    // Prazna bam1 struktura
    bam1_t *aln = bam_init1();

    // Čitanje datoteke
    while (sam_read1(file, header, aln) >= 0) {
        if (sam_write1(out, header, aln) < 0) {
            std::cerr << "Greška pri pisanju zapisa.\n";
            break;
        }
    }

    /*
        Čišćenje memorije
    */

    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(file);
    sam_close(out);

    std::cout << "Konverzija uspješno završena!" << std::endl;
}
/*
    pomoćna funkcija : reverse_comploment
*/
std::string reverse_complement(std::string seq) {
    std::map<char,char> comp = {{'A','T'},{'T','A'},{'C','G'},{'G','C'}};
    std::string temp;
    for (char c : seq) temp.push_back(comp[c]);
    std::reverse(temp.begin(), temp.end());
    return temp;
}
/*
    parser sam-a, ako bude potrebno možemo dodavati još vrijednosti koje treba parsirati
*/
std::vector<SamRecord> parse_sam(std::string filename){
    std::vector<SamRecord> records;
    std::ifstream file(filename);
    if(!file.is_open()) throw std::runtime_error("Ne mogu otvoriti fajl: " + filename);
    std::string line;
    while (getline(file, line)){
        if(line[0] == '@') continue;

        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;
        while (ss >> token){
            columns.push_back(token);
        }
        SamRecord record;
        record.flag = std::stoi(columns[1]);
        record.pos = std::stoi(columns[3]) - 1; //  SAM je 1-based, pa zato -1 da imamo 0-based
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