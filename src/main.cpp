#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "fasta_reader.hpp"
#include "sam_parser.hpp"
#include "cigar_parser.hpp"
#include "variant_caller.hpp"
#include "types.hpp"

using namespace std::string_literals;

int main(int argc, char* argv[]){
    // početni argumenti
    if (argc < 4) {
        std::cout << "Greška! Uporaba: " << argv[0] << " <putanja FASTA datoteke> <putanja SAM očitanja> <putanja za spremanje>\n";
        exit(1); 
    }
    std::string ref_path(argv[1]);
    std::string sam_path(argv[2]);

    std::ofstream output(argv[3]);

    std::string reference = read_fasta(ref_path);

    std::vector<SamRecord> records = parse_sam(sam_path);

    std::vector<Variant> allVariants;
    for(auto record : records){
        std::vector<Variant> v = parse_cigar(record, reference);
        allVariants.insert(allVariants.end(), v.begin(), v.end());
    }

    // traženje mutacije i određivanje pogreške
    std::vector<Mutation> mutations = variant_caller(records, reference);

    // psianje u CSV formatu
    std::string mutations_type[3] = {"X", "I", "D"};
    for (const auto& mutation : mutations) {
        std::string line;

        // stvaranje CSV reda
        line.append(mutations_type[mutation.type] + ",");
        line.append(std::to_string(mutation.pos));
        
        // određivanje vrste mutacije
        switch (mutation.type) {
            case MutationType::SUBSTITUTION:
                line.append(","s + mutation.alt_base);
                break;

            case MutationType::INSERTION:
                line.append(","s + mutation.alt_base);
                break;
            
            case MutationType::DELETION:
                line.append(",-");
                break;
        }
        
        output << line << "\n";

    }
    output.close();
}