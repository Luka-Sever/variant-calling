#ifndef VARIANT_CALLER_HPP
#define VARIANT_CALLER_HPP

#include <iostream>
#include <vector>
#include "types.hpp"

// vraća indeks slova
int base_to_index(char base) {
    switch (base)
    {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return 4;
    }
}

// funkcija vraća bazu od indeksa
char index_to_base(int index) {
    switch(index) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

typedef struct BasePileup {
    int count[5] = {0}; // A, C, G, T ili Unkown
} BasePileup;

// glavna funkcija određivanja mutacija, vraća
// vektor objekata mutacija
std::vector<Mutation> variant_caller(std::vector<SamRecord> records, std::string reference) {
    std::vector<Mutation> mutations;

    if (reference.empty()) return mutations;

    // struktura u kojoj se broje reference
    std::vector<BasePileup> pileup(reference.length());

    // olakšano čitanje CIGAR stringa
    std::regex patrn("[0-9]+(M|I|D|S|H)");

    for (const auto &record : records) {
        // ne zanimaju nas nepreslikani redovi
        if (record.flag & 4) continue;
        if (record.cigar == "*") continue;

        int ref_pos = record.pos;
        int read_pos = 0;

        std::string exp = record.cigar;
        std::smatch match;
        
        while(std::regex_search(exp, match, patrn)) {
            std::string m = match.str();
            char operation = m[m.length() - 1];
            int repetition = std::stoi(m.substr(0, m.length() - 1));
            exp = match.suffix();
        }

    }

    return mutations;

}

#endif