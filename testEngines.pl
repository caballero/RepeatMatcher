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

while (my $res = $searchResults->next()) {
    print Dumper $res;
    print "=" x 80, "\n";
}
