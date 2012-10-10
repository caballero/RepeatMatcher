#!/usr/bin/perl

use strict;
use warnings;
use lib './lib';
use NCBIBlastSearchEngine;
use WUBlastSearchEngine;
use CrossmatchSearchEngine;
use SearchEngineI;
use SearchResultCollection;

my $Engine = WUBlastSearchEngine->new( pathToEngine=>"/usr/local/wublast/blastn" );
$Engine->setMatrix( "./Matrices/wublast/nt/wumatrix" );
$Engine->setQuery( "repeat.fa" );
$Engine->setSubject( "allMis0.fa" );
my $searchResults = $Engine->search();

open WU, ">tests.wublast" or die;
my $hits = $searchResults->size();
print WU "WU-Blast found $hits hits\n";
for ( my $i = 0 ; $i < $hits; $i++ ) {
    my $qName  = $searchResults->get( $i )->getQueryName;
    my $qStart = $searchResults->get( $i )->getQueryStart;
    my $qEnd   = $searchResults->get( $i )->getQueryEnd;
    my $hName  = $searchResults->get( $i )->getSubjName;
    my $hStart = $searchResults->get( $i )->getSubjStart;
    my $hEnd   = $searchResults->get( $i )->getSubjEnd;
    my $dir    = $searchResults->get( $i )->getOrientation;
    my $score  = $searchResults->get( $i )->getScore;
    print WU join "\t", $qName, $qStart, $qEnd, $hName, $hStart, $hEnd, $dir, "$score\n";
}
close WU;

my $Engine2 = NCBIBlastSearchEngine->new( pathToEngine=>"/usr/local/rmblast/bin/rmblastn" );
$Engine2->setMatrix( "./Matrices/ncbi/nt/simple.matrix" );
$Engine2->setQuery( "repeat.fa" );
$Engine2->setSubject( "allMis0.fa" );
my $searchResults2 = $Engine2->search();

open RM, ">tests.rmblast" or die;
my $hits2 = $searchResults2->size();
print RM "RM-Blast found $hits2 hits\n";
for ( my $i = 0 ; $i < $hits2; $i++ ) {
    my $qName  = $searchResults2->get( $i )->getQueryName;
    my $qStart = $searchResults2->get( $i )->getQueryStart;
    my $qEnd   = $searchResults2->get( $i )->getQueryEnd;
    my $hName  = $searchResults2->get( $i )->getSubjName;
    my $hStart = $searchResults2->get( $i )->getSubjStart;
    my $hEnd   = $searchResults2->get( $i )->getSubjEnd;
    my $dir    = $searchResults2->get( $i )->getOrientation;
    my $score  = $searchResults2->get( $i )->getScore;
    print RM join "\t", $qName, $qStart, $qEnd, $hName, $hStart, $hEnd, $dir, "$score\n";
}
close RM;

