RepeatMatcher
=============

Simple analysis tool to manually annotate RepeatModeler sequences.

** Under development, please contact me before using it.**

MOTIVATION
After running RepeatModeler you probably will have a bunch of sequences to be manually annotated. RepeatMatcher is designed to help in this hard task. We provide a simple pipeline to generate all the required files and a graphical user interface to manipulate and annotate the sequences. The GUI is intented to show all the information to make the decisions easy for each sequence.

THE PIPELINE
We started with the RepeatModeler output file (fasta), the steps are:
  1. Masking all the low complexity regions with RepeatMasker and discard sequences with filters for a minimal composition.
  2. All sequences will be compared with itself using Cross_match.
  3. All sequences are compared in nucleotide space with known repeats consensi with Cross_match.
  4. All sequences are compared in peptide space with known repeats consensi with Blastx.
  5. All sequences are compared in peptide space withthe peptides in NCBI NR database with Blastx.
  6. All sequences are folded with RNAVienna and the principal structure is plotted with R:R4RNA.

To run the pipeline, first you need to resolve all the dependencies:
  - RepeatMasker
  - Blast
  - Cross_match
  - RNAVienna 
  - R 
  - R:R4RNA
  - Databases: NCBI NR, known repeats (in nucleotides and peptides).
  - Perl:Tk (for the GUI)

Then, check if it's required to adjust the parameters in RepeatMatcher.conf.

Launching the pipeline:
    perl RepeatMatcher.pl  -s sequences.fa -k knownrepeats.fa -o NewAnnotation

for a complete list of options:
    perl RepeatMatcher.pl --help
 
THE GUI
After you ran the pipeline, you have all the files required for RepeatMatcherGUI.pl. The GUI uses a text file to keep all the parameters and annotations made (the LOG file). The first lines of that file record the files used, the folowing lines records all changes for each sequence.

Launching the GUI for the first time:
    perl RepeatMatcherGUI.pl -o OUT -i SEQS -s SELF -a ALIGN -b BLAST -n BLAST -f FOLD -l LOG

Launching the GUI with a generated LOG:
    perl RepeatMatcherGUI.pl -r LOG

Juan Caballero, Institute for Systems Biology @ 2012
