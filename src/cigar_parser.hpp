#ifndef CIGAR_PARSER_HPP
#define CIGAR_PARSER_HPP

#include <iostream>
#include <vector>
#include <regex>
#include <string>
#include "types.hpp"

/*
    TODO: create a function that parses the CIGAR string :)

    DONE: Yes!
*/
/*
    Description: we define a regex pattern; in our case, the CIGAR string consists of a SEQUENCE of digits
    and a SINGLE letter (M, I, D, S, H). We set initial positions (ref_pos, read_pos),
    and for each recognized match for the given regex and CIGAR string, e.g. 89M, we do the following:

    1. Define the operation (always at the last position, e.g. 'M') and the number of repetitions of that operation.
    2. If the operation is M:
            it can be a MISMATCH or MATCH, so we check all cases
            if it is a MISMATCH we add a Variant ('X', ref_pos, s) to the vector
            we increment both ref_pos and read_pos

    3. If the operation is I:
            means bases exist in the sequence but not in the reference
            we add a Variant ('I', ref_pos, list of bases) to the vector
            we increase read_pos by repetition

    4. If the operation is D:
            means bases do not exist in the sequence but exist in the reference
            we add a Variant ('D', ref_pos, "-") to the vector | *repetition
            we increase ref_pos by one | *repetition

    5. If the operation is S:
            bases exist but are not mapped, so we skip them

    6. If the operation is H:
            ignore
    
   Written by: Ivan Dundović
*/
std::vector<Variant>parse_cigar(SamRecord record, std::string reference){
    std::vector<Variant> variants;
    // Defining the regex
    std::regex patrn("[0-9]+(M|I|D|S|H)");
    std::string exp = record.cigar;
    std::smatch match;
    int ref_pos = record.pos;
    int read_pos = 0;
    // Finding the CIGAR string using regex
    while (std::regex_search(exp, match, patrn))
    {
        std::string m = match.str();
        char operation = m[m.length() - 1];
        int repetition = std::stoi(m.substr(0, m.length() - 1));
        exp = match.suffix();
        // Operations are described in detail above before the function definition starts
        if(operation == 'M'){
            int i = 0;
            while (i < repetition)
            {
                if (read_pos >= record.seq.length() || ref_pos >= reference.length()) break;
                if(record.seq[read_pos] != reference[ref_pos]){
                    std::string s(1, record.seq[read_pos]);
                    Variant v = {'X', ref_pos, s};
                    variants.push_back(v);
                }
                i++;
                ref_pos++;
                read_pos++;
            }    
        }else if (operation == 'I'){
            if (read_pos + repetition <= record.seq.length()) {
                std::string bases = record.seq.substr(read_pos, repetition);
                Variant v = {'I', ref_pos, bases};
                variants.push_back(v);
            }
            read_pos += repetition;

        }else if(operation == 'D'){
            for(int i =0; i < repetition; i++){
                Variant v = {'D', ref_pos, "-"};
                variants.push_back(v);
            }
            ref_pos += repetition;

        }else if(operation == 'S') read_pos += repetition;
        else if(operation == 'H') continue;
    }
    return variants;
}

#endif