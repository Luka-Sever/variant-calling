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
    int insertion_count = 0;
    int deletion_count = 0;
} BasePileup;

// glavna funkcija određivanja mutacija, vraća
// vektor objekata mutacija
// Napisao: Luka Sever
std::vector<Mutation> variant_caller(std::vector<SamRecord> records, std::string reference) {
    std::vector<Mutation> mutations;

    // za svaki slučaj
    if (reference.empty()) return mutations;

    // struktura u kojoj se broje reference
    std::vector<BasePileup> pileup(reference.length());

    // olakšano čitanje CIGAR stringa
    std::regex patrn("[0-9]+(M|I|D|S|H)");
    for (const auto &record : records) {

        // ne zanimaju nas nepreslikani redovi
        if (record.flag & 4) continue;
        if (record.cigar == "*") continue;

        // zanima nas pozicija u referentnom genomu
        // i pozicija koju trenutno čitamo
        int ref_pos = record.pos;
        int read_pos = 0;

        // pratimo CIGAR string, čitat ćemo kakve su promjene u očitanju
        std::string exp = record.cigar;
        std::smatch match;
        
        // dok nismo do kraja pretražili koliko očitanja ima isti CIGAR
        while(std::regex_search(exp, match, patrn)) {
            // rezultat koji smo dobili
            std::string m = match.str();
            // koja operacija
            char operation = m[m.length() - 1];
            // koliko puta se ponavlja
            int repetition = std::stoi(m.substr(0, m.length() - 1));
            exp = match.suffix();

            // za određene operacije 
            switch (operation) {
                // match
                case 'M':
                    for (int i = 0; i < repetition; i++) {
                        // koliko se baza poklapa s referencom
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
                
                // umetanje
                case 'I':
                    if (ref_pos >= 0 && ref_pos < (int)reference.length()) {
                        pileup[ref_pos].insertion_count++;
                    }
                    read_pos += repetition;
                    break;
                    
                // brisanje
                case 'D':
                    if (ref_pos >= 0 && ref_pos < (int)reference.length())  {
                        pileup[ref_pos].deletion_count++;
                    }
                    ref_pos += repetition;
                    break;
                    
                // substitucija
                case 'S':
                    read_pos += repetition;
                    break;
                
                default:
                    break;
                }
                    
            }
        }

        // glavni parametri očitanja,
        // depth nam govori koliko očitanja imamo,
        // a frekvencija govori koliko nam treba
        // da zaključimo da se dogodila mutacija
        // (postavljeno na 50%)
        const int MIN_DEPTH = 1;
        const float MIN_FREQ = 0.5f;


        for (int pos = 0; pos < (int)reference.length(); pos++) {
            int ref_base_index = base_to_index(reference[pos]);
            
            // koliko je velika substitucija
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
            // Odluka: je li umetanje, brisanje ili substitucija
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