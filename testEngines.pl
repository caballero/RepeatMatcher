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
$Engine->setMatrix( "/home/asmit/Matrices/nt/wumatrix" );
$Engine->setQuery( "./gator_annotation/rep" );
$Engine->setSubject( "./gator_annotation/allMis0.fa" );
my $searchResults = $Engine->search();

my $hits = $searchResults->size();
print "WU-Blast found $hits hits\n";
for ( my $i = 0 ; $i < $hits; $i++ ) {
    my $qName  = $searchResults->get( $i )->getQueryName;
    my $qStart = $searchResults->get( $i )->getQueryStart;
    my $qEnd   = $searchResults->get( $i )->getQueryEnd;
    my $hName  = $searchResults->get( $i )->getSubjName;
    my $hStart = $searchResults->get( $i )->getSubjStart;
    my $hEnd   = $searchResults->get( $i )->getSubjEnd;
    my $dir    = $searchResults->get( $i )->getOrientation;
    my $score  = $searchResults->get( $i )->getScore;
    print join "\t", $qName, $qStart, $qEnd, $hName, $hStart, $hEnd, $dir, "$score\n";
}

my $Engine2 = NCBIBlastSearchEngine->new( pathToEngine=>"/usr/local/rmblast/bin/rmblastn" );
$Engine2->setMatrix( "/home/asmit/Matrices/simple.matrix" );
$Engine2->setQuery( "./gator_annotation/rep" );
$Engine2->setSubject( "./gator_annotation/allMis0.fa" );
my $searchResults2 = $Engine2->search();

my $hits2 = $searchResults2->size();
print "RM-Blast found $hits2 hits\n";
for ( my $i = 0 ; $i < $hits2; $i++ ) {
    my $qName  = $searchResults2->get( $i )->getQueryName;
    my $qStart = $searchResults2->get( $i )->getQueryStart;
    my $qEnd   = $searchResults2->get( $i )->getQueryEnd;
    my $hName  = $searchResults2->get( $i )->getSubjName;
    my $hStart = $searchResults2->get( $i )->getSubjStart;
    my $hEnd   = $searchResults2->get( $i )->getSubjEnd;
    my $dir    = $searchResults2->get( $i )->getOrientation;
    my $score  = $searchResults2->get( $i )->getScore;
    print join "\t", $qName, $qStart, $qEnd, $hName, $hStart, $hEnd, $dir, "$score\n";
}

