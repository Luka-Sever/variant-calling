#ifndef TYPES_HPP
#define TYPES_HPP
#include <string>

/*
    TO DOO:
    1. Dodati sve tipove što će nam trebati. :)

    DOVRŠENO: MOŽDA
*/
/*
    SamRecord struktura
    Detaljniji opis:
    Sadrži varijablu flag tipa int koja predstavlja vrijednost zastavice u .sam fajlu.
    (Značenje zastavica: 0 -> read je pravilno mapiran, 4 -> read nije mapiran, 16 -> read je u reversu)
    Zatim imamo int pos koja predstavlja poziciju na referenci, te dva stringa cigar u kojem se nalazi cigar string te
    seq u kojoj se nalazi sekvenca.
*/

struct SamRecord
{
    int flag; // za provjeru je li mapiran, reversed
    int pos; // pozicija na referenci
    std::string cigar; // cigar string tipa "5M2I3M1D4M"
    std::string seq; // procitana sekvenca
};

/*
    Variant struktura
*/
struct Variant
{
    char type;
    int pos;    // pozicija na referenci
    std::string bases; //koja baza/e 
};

struct Mutation
{
    char type;
    int pos;
    std::string base; 
};


#endif