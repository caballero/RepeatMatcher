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

my $NCBIEngine = NCBIBlastSearchEngine->new( pathToEngine=>"/usr/local/rmblast/bin/rmblastn" );
$NCBIEngine->setMatrix( "/home/asmit/Matrices/simple.matrix" );
$NCBIEngine->setQuery( "./gator_annotation/rep" );
$NCBIEngine->setSubject( "./gator_annotation/allMis0.fa" );
my $searchResults = $NCBIEngine->search();

while (my $res = $searchResults->next()) {
    print Dumper $res;
    print "=" x 80, "\n";
}
