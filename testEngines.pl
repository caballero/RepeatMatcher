#!/usr/bin/perl

use strict;
use warnings;
use lib './lib';
use NCBIBlastSearchEngine;
use WUBlastSearchEngine;
use CrossmatchSearchEngine;
use SearchEngineI;
use SearchResultCollection;
use Data::Dumper;
use ArrayListIterator;

my $Engine = WUBlastSearchEngine->new( pathToEngine=>"/usr/local/wublast/blastn" );
$Engine->setMatrix( "/home/asmit/Matrices/nt/wumatrix" );
$Engine->setQuery( "./gator_annotation/rep" );
$Engine->setSubject( "./gator_annotation/allMis0.fa" );
my $searchResults = $Engine->search();

my $hits = $searchResults->size();
print "WU-Blast found $hits hits\n";
for ( my $i = 0 ; $i < $hits; $i++ ) {
    my $qName  = $searchResults->get( $i )->queryName;
    my $qStart = $searchResults->get( $i )->queryStart;
    my $qEnd   = $searchResults->get( $i )->queryEnd;
    my $hName  = $searchResults->get( $i )->subjName;
    my $hStart = $searchResults->get( $i )->subjStart;
    my $hEnd   = $searchResults->get( $i )->subjEnd;
    my $dir    = $searchResults->get( $i )->orientation;
    my $score  = $searchResults->get( $i )->score;
    print join "\t", $qName, $qStart, $qEnd, $hName, $hStart, $hEnd, $dir, "$score\n";
}

exit 1;
my $Engine2 = NCBIBlastSearchEngine->new( pathToEngine=>"/usr/local/rmblast/bin/rmblastn" );
$Engine2->setMatrix( "/home/asmit/Matrices/simple.matrix" );
$Engine2->setQuery( "./gator_annotation/rep" );
$Engine2->setSubject( "./gator_annotation/allMis0.fa" );
my $searchResults2 = $Engine2->search();

print "NCBI-Blast\n";
while (my $res2 = $searchResults2->next()) {
    print Dumper $res2;
    print "=" x 80, "\n";
}
