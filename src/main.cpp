#include <iostream>
#include "fasta_reader.hpp"
#include "sam_parser.hpp"
#include "cigar_parser.hpp"
#include "variant_caller.hpp"
#include "types.hpp"
#include <vector>

int main(){
    std::string ref_path = "../data/lambda.fasta";
    std::string sam_path = "../data/aligment.sam";

    std::string reference = read_fasta(ref_path);

    std::vector<SamRecord> records = parse_sam(sam_path);

    std::vector<Variant> allVariants;
    for(auto record : records){
        std::vector<Variant> v = parse_cigar(record, reference);
        allVariants.insert(allVariants.end(), v.begin(), v.end());
    }

    // TO DOOO :)

    std::vector<Mutation> mutations = variant_caller(all_variants, records, reference);
}