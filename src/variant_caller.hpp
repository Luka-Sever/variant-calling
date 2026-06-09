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
    int match_count = 0;
    // Stores the inserted base sequence for reporting purposes
    std::string inserted_bases;
} BasePileup;

// Main function for determining mutations, returns
// a vector of mutation objects
// Written by: Luka Sever
std::vector<Mutation> variant_caller(const std::vector<Variant>& allVariants, std::string& reference) {
    std::vector<Mutation> mutations;

    // Just in case
    if (reference.empty()) return mutations;

    // Structure used to count references
    std::vector<BasePileup> pileup(reference.length());

    // Simplifies reading the CIGAR string
    for (const auto &variant : allVariants) {

        // Out of bounds safety
        if (variant.pos < 0 || variant.pos >= (int)reference.length()) continue;

        if (variant.type == 'X') {
            // Mismatch: increment count for alternative base
            if (!variant.bases.empty()) {
                char alt_base = variant.bases[0];
                int index = base_to_index(alt_base);
                pileup[variant.pos].count[index]++;
            }
        } else if (variant.type == '=') {
            pileup[variant.pos].match_count++;
        } else if (variant.type == 'I') {
            // Insertion: store the inserted sequence for later reporting
            pileup[variant.pos].insertion_count++;
            if (!variant.bases.empty()) {
                pileup[variant.pos].inserted_bases = variant.bases;
            }
        } else if (variant.type == 'D') {
            // Deletion
            pileup[variant.pos].deletion_count++;
        }
    }
        // Main read parameters,
        // depth indicates how many reads we have,
        // and frequency indicates how much is needed
        // to conclude that a mutation has occurred
        const int MIN_DEPTH = 13;
        const float MIN_FREQ = 0.4f;

        for (int pos = 0; pos < (int)reference.length(); pos++) {
            int ref_base_index = base_to_index(reference[pos]);
            
            // Size of the substitution
            int substitute_depth = 0;
            for (int i = 0; i < 5; i++) {
                substitute_depth += pileup[pos].count[i];
            }
            int total_depth = substitute_depth + pileup[pos].insertion_count + pileup[pos].deletion_count + pileup[pos].match_count;
            
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
            if (pileup[pos].insertion_count > max_count && pileup[pos].insertion_count > pileup[pos].deletion_count) {
                float freq = (float)pileup[pos].insertion_count / total_depth;
                if (freq >= MIN_FREQ) {
                    Mutation m;
                    m.pos = pos;
                    m.count = pileup[pos].insertion_count;
                    m.alt_base = pileup[pos].inserted_bases.empty() ? 'I' : pileup[pos].inserted_bases[0];
                    m.type = MutationType::INSERTION;
                    mutations.push_back(m);
                    continue;
                }
            }

            if (pileup[pos].deletion_count > max_count && pileup[pos].deletion_count > pileup[pos].insertion_count) {
                float freq = (float)pileup[pos].deletion_count / total_depth;
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