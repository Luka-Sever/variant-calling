#ifndef CIGAR_PARSER_HPP
#define CIGAR_PARSER_HPP

#include <iostream>
#include <vector>
#include <regex>
#include <string>
#include "types.hpp"

/*
    TO DOO: napraviti funkciju koja parsira CIGAR string :)

    DOVRŠENO: Da!

*/

/*
    opis: definiramo jedan regex niz u našem slučaju CIGAR string se sastoji od NIZA znamenki i JEDNOG slova (M, I, D, S, H).
    Postavimo početne pozicije (ref_pos, read_pos) te za svaki prepoznat match za zadani regex i cigar string npr. 89M radimo sljedeće:
    1. definiramo operaciju (uvijek na zadnjem mijestu m) i broj ponavljanja zadane operacije (repetition)

    2. ako je opearcija M:
            može biti MISMATCH ili MATCH pa provjeravamo sve
            ako je MISMATCH dodajemo jedan Variant ('X', ref_pos, s) u vektor 
            povećavamo i ref_pos i read_pos
    3. ako je opearcija I:
            znaći baze postoje u seq ali ne i u referenci
            dodajemo Variant ('I', ref_pos, niz baza) u vektor
            povećavamo read_pos za repetition
    4. ako je operacija D:
            znaći baze ne postoje u seq ali postoje u referenci
            dodajemo Variant ('D', ref_pos, "-") u vektor | *repetition
            povećavamo ref_pos za jedan  | *repetition
    5. ako je operacija S:
            baze postoje ali nisu mapirane pa preskaćemo
    6. ako je operacija H:
            ignoriraj
*/
std::vector<Variant>parse_cigar(SamRecord record, std::string reference){
    std::vector<Variant> variants;
    std::regex patrn("[0-9]+(M|I|D|S|H)");
    std::string exp = record.cigar;
    std::smatch match;
    int ref_pos = record.pos;
    int read_pos = 0;
    
    while (std::regex_search(exp, match, patrn))
    {
        std::string m = match.str();
        char operation = m[m.length() - 1];
        int repetition = std::stoi(m.substr(0, m.length() - 1));
        exp = match.suffix();
        if(operation == 'M'){
            int i = 0;
            while (i < repetition)
            {
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
            std::string bases = record.seq.substr(read_pos, repetition);
            Variant v = {'I', ref_pos, bases};
            read_pos += repetition;
            variants.push_back(v);

        }else if(operation == 'D'){
            for(int i =0; i < repetition; i++){
                Variant v = {'D', ref_pos, "-"};
                ref_pos++;
                variants.push_back(v);
            }

        }else if(operation == 'S') read_pos += repetition;
        else if(operation == 'H') continue;
    }
    return variants;
}

#endif