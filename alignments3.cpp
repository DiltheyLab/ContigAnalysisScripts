#include <iostream>
#include <fstream>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>


using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 2){
        std::cout << "Wrong number of inputs, please provide a path to the input fasta file" << std::endl;
        return 1;  // Invalid number of arguments.
    }

    typedef String<Dna> TSequence;                             // sequence type
    typedef StringSet<DnaString> TStringSet;                    // container for strings
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;   // dependent string set
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;       // alignment graph

    SeqFileIn seqFileIn(toCString(argv[1]));
    //readRecord(id, seq, seqFileIn);
    

    StringSet<CharString> ids;
    Dna5String seq;
    TStringSet all_bridges;

    // Reads all remaining records.
    readRecords(ids, all_bridges, seqFileIn, 4);

    TSequence seq1 = "GTCTTCAAATCCTTCTAGTGAATTTTGAGGTTTTTTTTTT";
    
    int nr_of_seqs = 4;

    TStringSet sequences;
    TAlignGraph alignG; 
    int score;
    
    for (int i = 0; i < nr_of_seqs; ++i){
        //std::cout << all_bridges[i] << std::endl;
        clear(sequences);
        appendValue(sequences, seq1);
        appendValue(sequences, all_bridges[i]);
        clear(alignG);
        alignG = TAlignGraph(sequences);
        score = globalAlignment(alignG, Score<int, Simple>(3, -3, -3), AlignConfig<true, false, true, false>(), LinearGaps());
        std::cout << "Score: " << score << std::endl;
        std::cout << alignG << std::endl;
        int tpos = (int) getLastCoveredPosition(alignG,0);
        int nid = -1;
        int npos = -1;
        getProjectedPosition(alignG,0, tpos-1, nid, npos);
        std::cout << "last covered: " << (int) 0 << "\t" << (int) tpos << std::endl;
        std::cout << "last covered: " << (int) nid << "\t" << (int) npos << std::endl;

    }



    /* this could be a good check
    int score = globalAlignment(alignG1, Score<int, Simple>(1, -1, -1), AlignConfig<true, true, true, true>(), LinearGaps());
    std::cout << "Score: " << score << std::endl;
    std::cout << alignG1 << std::endl;
    */

    /* this doesnt work
    score = globalAlignment(alignG2, Score<int, Simple>(1, -1, -1), AlignConfig<false, true, false, true>(), LinearGaps());
    std::cout << "Score: " << score << std::endl;
    std::cout << alignG2 << std::endl;
    */

    /* this works
    score = globalAlignment(alignG3, Score<int, Simple>(1, -1, -1), AlignConfig<true, false, true, false>(), LinearGaps());
    std::cout << "Score: " << score << std::endl;
    std::cout << alignG3 << std::endl;
    */

    return 0;
}

