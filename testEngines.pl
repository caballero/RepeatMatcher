#!/usr/bin/perl

use strict;
use warnings;
use lib './lib';
use NCBIBlastSearchEngine;
use WUBlastSearchEngine;
use CrossmatchSearchEngine;
use SearchEngineI;
use SearchResultCollection;

my $NCBIEngine = NCBIBlastSearchEngine->new(pathToEngine=>"/usr/local/rmblast/bin/rmblastn" );
$NCBIEngine->setMatrix( "/users/asmit/Matrices/simple.matrix" );
$NCBIEngine->setQuery( "rep" );
$NCBIEngine->setSubject( "allMis0.fa" );
my $searchResults = $NCBIEngine->search();

foreach my $obj (@$searchResults) {
    print "$obj\n";
}
