#ifndef VARIANT_CALLER_HPP
#define VARIANT_CALLER_HPP

#include <iostream>
#include <vector>
#include "types.hpp"
/*
    TO DOO: napisati funkciju variant_caller :)

    Dovršeno: NE
*/

std::vector<Mutation> variant_caller(std::vector<Variant> allVariants, std::vector<SamRecord> records, std::string reference);

#endif