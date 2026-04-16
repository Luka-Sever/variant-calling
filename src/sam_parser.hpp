#include <iostream>
#include <string>
#include <htslib/sam.h>

using namespace std;

/**
 * Pretvaramo sam datoteku u bam datoteku da samnjimo
 * potrošenu memoriju
 */

void convert_sam_to_bam(const char* inputFile, const char* outputFile) {
    /*
    
        Pokušavamo čitati input file i spremati
        u output file
    
    */
    samFile* file = sam_open(inputFile, "r");
    if (!file) cout << "Ne mogu otvoriti datoteku.\n";

    sam_hdr_t *header = sam_hdr_read(file);
    if (!header) cout << "Ne mogu pročitati header.\n";

    samFile *out = sam_open(outputFile, "wb");
    if (!out) cout << "Ne mogu stvoriti bam datoteku header.\n";

    if (sam_hdr_write(out, header) < 0) 
        cout << "Ne mogu pisati u bam datoteku.\n";


    // Prazna bam1 struktura
    bam1_t *aln = bam_init1();

    // Čitanje datoteke
    while (sam_read1(file, header, aln) >= 0) {
        if (sam_write1(out, header, aln) < 0) {
            std::cerr << "Greška pri pisanju zapisa.\n";
            break;
        }
    }

    /*
        Čišćenje memorije
    */

    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(file);
    sam_close(out);

    std::cout << "Konverzija uspješno završena!" << std::endl;
}