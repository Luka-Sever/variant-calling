#ifndef TYPES_HPP
#define TYPES_HPP
#include <string>

/*
    TODO:
    1. Add all the types we will need. :)

    DONE: YES
*/
// Enumeration of mutation types
// Written by: Luka Sever
enum MutationType {
    SUBSTITUTION,
    INSERTION,
    DELETION
};

// Structure representing a SAM record
// Written by: Ivan Dundović
struct SamRecord
{
    int flag; // indicates whether the read is mapped, reversed, etc.
    int pos; // position on the reference
    std::string cigar; // CIGAR string (e.g. "5M2I3M1D4M")
    std::string seq; // read sequence
};

/*
    Variant structure
    Written by: Ivan Dunodvić
*/
struct Variant
{
    char type;
    int pos;    // position on the reference
    std::string bases; // which base(s)
};

// Structure representing a mutation
// Written by: Luka Sever
struct Mutation
{
    int pos;
    int count;
    char alt_base;
    MutationType type;
};

#endif