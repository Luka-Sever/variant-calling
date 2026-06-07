#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>
#include "fasta_reader.hpp"
#include "sam_parser.hpp"
#include "cigar_parser.hpp"
#include "variant_caller.hpp"
#include "types.hpp"

using namespace std::string_literals;

int main(int argc, char* argv[]){
    // Argument check
    if (argc < 4) {
        std::cout << "Error! Usage: " << argv[0] << " <path to FASTA file> <path to SAM file> <path to save>\n";
        exit(1); 
    }

    // reading arguments
    std::string ref_path(argv[1]);
    std::string sam_path(argv[2]);
    std::ofstream output(argv[3]);

    // start measuring time
    auto total_time_start = std::chrono::high_resolution_clock::now();

    // reading FASTA
    std::cout << "[1/4] Reading FASTA and SAM files..." << std::endl;
    auto time_start1 = std::chrono::high_resolution_clock::now();
    std::string reference = read_fasta(ref_path);
    std::vector<SamRecord> records = parse_sam(sam_path);
    auto time_end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_time1 = time_end1 - time_start1;
    std::cout << "Completed in: " << read_time1.count() << "s" << std::endl;

    // reading all variants
    std::cout << "[2/4] Getting variants from CIGAR string..." << std::endl;
    auto time_start2 = std::chrono::high_resolution_clock::now();
    std::vector<Variant> allVariants;
    for(auto record : records){
        std::vector<Variant> v = parse_cigar(record, reference);
        allVariants.insert(allVariants.end(), v.begin(), v.end());
    }
    auto time_end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_time2 = time_end2 - time_start2;
    std::cout << "Completed in: " << read_time2.count() << "s" << std::endl;

    // main funcftion
    std::cout << "[3/4] Using variant caller..." << std::endl;
    auto time_start3 = std::chrono::high_resolution_clock::now();
    std::vector<Mutation> mutations = variant_caller(allVariants, reference);
    auto time_end3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_time3 = time_end3 - time_start3;
    std::cout << "Completed in: " << read_time3.count() << "s" << std::endl;

    // write in CSV 
    // Written by: Luka Sever
    std::cout << "[4/4] Writing to CSV..." << std::endl;
    auto time_start4 = std::chrono::high_resolution_clock::now();
    std::string mutations_type[3] = {"X", "I", "D"};
    for (const auto& mutation : mutations) {
        std::string line;

        // creating CSV row
        line.append(mutations_type[mutation.type] + ",");
        line.append(std::to_string(mutation.pos));
        
        // resolving mutation type
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
    auto time_end4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_time4 = time_end4 - time_start4;
    std::cout << "Completed in: " << read_time4.count() << "s" << std::endl;
    auto total_time_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time = total_time_end - total_time_start;
    std::cout << "Total time: " << total_time.count() << "s" << std::endl;

    return 0;

}