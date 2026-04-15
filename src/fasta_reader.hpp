#ifndef PARSER_HPP
#define PARSER_HPP
#include <iostream>
#include <fstream>
#include <string>

/*
    TO DOO : 
    1. Napraviti funkciju koja cita fasta file te sprema genom kao jedan string kako bi ga mogli koristiti za usporedbu :)

    DOVRŠENO : DA
*/
/*
    Detaljniji opis funkcije:
    Funkcija prima putanju do fajl (path). Kreiramo varijablu file tipa infstream (input file stream). (linija 25)
    Zatim provjeravamo jesmo li uspijeli otvoriti fajl ako da nastavljamo inače bacamo grešku. (linija 26)
    Kad smo prošli if petlju nastavljamo tako što stvaramo string genom i line. Genom koristimo kako bi u njega spremili
    cijeli genom, dok line koristimo kako bi u nju stavljali liniju fajla koje ćemo dodavati u genom. (linija 27)
    Pomću getline(file, line) uzimamo liniju po liniju i to radimo u while petlji sve dok nismo ostali bez linija. (linija 28)
    Koristimo if provjeru ako fajl zapocinje s < to znaci da je header i ne treba nam u genomu pa ga ne dodajemo. (linija 28)
    I na kraju veaćamo genom. (linija 30)
*/

std::string read_fasta(std::string path){
    std::ifstream file(path);
    if(!file.is_open()) throw std::runtime_error("Ne mogu otvoriti fajl: " + path);
    std::string genom, line;
    while (std::getline(file, line)) if(!line.empty() && line[0] != '<') genom += line;

    return genom;
}

#endif