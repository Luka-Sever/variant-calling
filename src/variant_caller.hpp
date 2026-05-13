#ifndef VARIANT_CALLER_HPP
#define VARIANT_CALLER_HPP

#include <iostream>
#include <vector>
#include "types.hpp"
/*
    TO DOO: napisati funkciju variant_caller :)

    Dovršeno: NE
*/

typedef struct BaseInfo {
    int count[5]; // A, C, G, T ili Unkown
    std::vector<float> qualities;
} BaseInfo;

std::vector<Mutation> variant_caller(std::vector<SamRecord> records, std::string reference) {

    std::vector<Mutation> variants;
    for (auto &record : records) {

        if (record.flag & 4) continue;

    }

    return variants;

}

#endif