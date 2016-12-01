# KCNFIN
RNA-Transcript Based Genome Scaffolding Tool 





KCNFIN is a transcript based genome scaffolding tool written in perl. It requires perl 5.15 or greater, Moose, and graph.pm.   
Input: psl file, genome assembly fasta file, transcript fasta file.
Output: new DNA assembly.


Usage:  --scaffold  --psl <string>  --rna <string>  --dna <string>  [options] 
    
            --scaffold        Takes both DNA and RNA fasta file and regular psl file describing the alignments between the sequences in the given fasta files  
            --gap             Takes same input as -scaffold option and will attempt to fill gaps and update alignment information prior to scaffolding          
            --dna             DNA fasta file
            --rna             RNA fasta file 
            --psl             Standard BLAT psl desribing alignments between sequences located in DNA and RNA fasta files
            -h|--help         prints this message

  [options]:

           [--pscore <Num> scaffolding path coverage score used to filter paths, ranges from 0 to 100 ( default 50 )]
           [--nscore <Num> alignment coverage score used to filter alignments from psl, ranges from 0 to 100 ( default 1.6 )]   
           [--pid <Num>   percentage of identity used to filter psl alignments ( default 98.1 )]
           [--dcount <Int> maximum number of DNA sequences that each RNA sequence can be aligned to ( default 50 )]
           [--match <Int> filter psl alignments with matched number of bases below this number (default 100)]  
           [--intron <Int> maximum intron size that can be represented by an edge ( default 100000 )]  
           [--overlap <Int> maximum base pairs that two DNA sequences can be overlapping on RNA before overlap is considered significant ( default 5 )] 

    




