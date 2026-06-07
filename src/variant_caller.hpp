#ifndef VARIANT_CALLER_HPP
#define VARIANT_CALLER_HPP

#include <iostream>
#include <vector>
#include "types.hpp"

// Returns the index of the character
// Written by: Luka Sever
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

// Function returns the base at the given index
// Written by: Luka Sever
char index_to_base(int index) {
    switch(index) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

// Base pileup structure
typedef struct BasePileup {
    int count[5] = {0}; // A, C, G, T or Unknown
    int insertion_count = 0;
    int deletion_count = 0;
} BasePileup;

// Main function for determining mutations, returns
// a vector of mutation objects
// Written by: Luka Sever
std::vector<Mutation> variant_caller(std::vector<SamRecord> records, std::string reference) {
    std::vector<Mutation> mutations;

    // Just in case
    if (reference.empty()) return mutations;

    // Structure used to count references
    std::vector<BasePileup> pileup(reference.length());

    // Simplifies reading the CIGAR string
    std::regex patrn("[0-9]+(M|I|D|S|H)");
    for (const auto &record : records) {

        // We are not interested in unmapped reads
        if (record.flag & 4) continue;
        if (record.cigar == "*") continue;

        // We are interested in the position in the reference genome
        // and the current read position
        int ref_pos = record.pos;
        int read_pos = 0;

        // We follow the CIGAR string to read the changes in the alignment
        std::string exp = record.cigar;
        std::smatch match;
        
        // While we haven't finished processing all reads with the same CIGAR
        while(std::regex_search(exp, match, patrn)) {
            // The result we obtained
            std::string m = match.str();
            // Which operation
            char operation = m[m.length() - 1];
            // How many times it is repeated
            int repetition = std::stoi(m.substr(0, m.length() - 1));
            exp = match.suffix();

            // For specific operations 
            switch (operation) {
                // match
                case 'M':
                    for (int i = 0; i < repetition; i++) {
                        if (read_pos >= record.seq.length() || ref_pos >= reference.length()) break;

                        // How many bases match the reference
                        if (ref_pos >= 0 && ref_pos < (int)reference.length() &&
                            read_pos >= 0 && read_pos < (int)record.seq.length()) {

                                char base = record.seq[read_pos];
                                int baseIndex = base_to_index(base);

                                pileup[ref_pos].count[baseIndex]++;

                            }
                            ref_pos++;
                            read_pos++;
                        }
                        break;
                
                // Insertion
                case 'I':
                    if (ref_pos >= 0 && ref_pos < (int)reference.length()) {
                        pileup[ref_pos].insertion_count++;
                    }
                    read_pos += repetition;
                    break;
                    
                // Deletion
                case 'D':
                    if (ref_pos >= 0 && ref_pos < (int)reference.length())  {
                        pileup[ref_pos].deletion_count++;
                    }
                    ref_pos += repetition;
                    break;
                    
                // Substitution
                case 'S':
                    read_pos += repetition;
                    break;
                
                default:
                    break;
                }
                    
            }
        }

        // Main read parameters,
        // depth indicates how many reads we have,
        // and frequency indicates how much is needed
        // to conclude that a mutation has occurred
        // (set to 50%)
        const int MIN_DEPTH = 1;
        const float MIN_FREQ = 0.5f;


        for (int pos = 0; pos < (int)reference.length(); pos++) {
            int ref_base_index = base_to_index(reference[pos]);
            
            // Size of the substitution
            int substitute_depth = 0;
            for (int i = 0; i < 5; i++) {
                substitute_depth += pileup[pos].count[i];
            }
            int total_depth = substitute_depth + pileup[pos].insertion_count + pileup[pos].deletion_count;
            
            if (total_depth < MIN_DEPTH) continue;

            int max_index = 0;
            int max_count = pileup[pos].count[0];
            for (int i = 1; i < 5; i++) {
                if (pileup[pos].count[i] > max_count) {
                    max_count = pileup[pos].count[i];
                    max_index = i;
                }
            }
            //
            // Decision: whether it is an insertion, deletion, or substitution
            //
            if (pileup[pos].insertion_count > max_index && pileup[pos].insertion_count > pileup[pos].deletion_count) {
                float freq = (float)pileup[pos].insertion_count / total_depth;
                if (freq >= MIN_FREQ) {
                    if (freq >= MIN_FREQ) {
                        Mutation m;
                        m.pos = pos;
                        m.count = pileup[pos].insertion_count;
                        m.alt_base = index_to_base(max_index);
                        m.type = MutationType::INSERTION;
                        mutations.push_back(m);
                        continue;
                    }
                }
            }

            if (pileup[pos].deletion_count > max_index && pileup[pos].deletion_count > pileup[pos].insertion_count) {
                float freq = (float)pileup[pos].deletion_count / total_depth;
                if (freq >= MIN_FREQ) {
                    if (freq >= MIN_FREQ) {
                        Mutation m;
                        m.pos = pos;
                        m.count = pileup[pos].deletion_count;
                        m.alt_base = 'D';
                        m.type = MutationType::DELETION;
                        mutations.push_back(m);
                        continue;
                    }
                }
            }

            if (max_index != ref_base_index && max_index != 4) {
                float freq = (float)max_count / total_depth;
                if (freq >= MIN_FREQ) {
                    Mutation m;
                    m.pos = pos;
                    m.alt_base = index_to_base(max_index);
                    m.count = max_count;
                    m.type = MutationType::SUBSTITUTION;
                    mutations.push_back(m);
                }
            }
        }


    return mutations;

}

#endif