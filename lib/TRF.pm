#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) TRF.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A wrapper class for the Tandem Repeat Finder ( TRF )
##      program written by Gary Benson, Department of Biomathematical
##      Sciences Mount Sinai School of Medicine.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2004 Developed by
#* Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
# bless( {
#          'matchWeight' => 2,
#          'mismatchPenalty' => 7,
#          'delta' => 7,
#          'pm' => 80,
#          'pi' => 10,
#          'minScore' => 50,
#          'maxPeriod' => 500,
#          'maskedOutputFile' => '/users/bob/seq.fasta.masked',
#          'pathToEngine' => '/usr/local/bin/trf',
#          'workDir' => '/tmp',
#          'sequenceFile' => '/users/bob/seq.fasta'
#        },
#      'TRF' );
#
###############################################################################
# ChangeLog
#
#     $Log: TRF.pm,v $
#     Revision 1.30  2011/04/26 22:41:44  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
# To Do:
#
#

=head1 NAME

TRF

=head1 SYNOPSIS

use TRF

  my $trfEngine = TRF->new( 
                       pathToEngine=>"/usr/local/bin/trf",
                       workDir => "/tmp" );

  $trfEngine->setMatchWeight( 2 );
  $trfEngine->setMismatchPenalty( 7 );
  $trfEngine->setDelta( 7 );
  $trfEngine->setPm( 80 );
  $trfEngine->setPi( 10 );
  $trfEngine->setMinScore( 50 );
  $trfEngine->setMaxPeriod( 500 );
  $trfEngine->setMaskedOutputFile( "/users/bob/seq.fasta.masked" );
  $trfEngine->setSequenceFile( "/users/bob/seq.fasta" );
  my $results = $trfEngine->search();

Usage: 

=head1 DESCRIPTION

An encapsulation class for the trf search engine.  It is capable
for running searches and returning the results as a perl array of 
TRFResult objects.

The trf parser built into this object captures several types of 
data from the trf output files.  


  Indices     Period  Copy   Cons    Perc    Perc    Score  A  C  G  T Entropy Consensus
              Size    Number Size    Matches Indels                     (0-2)
  ---------------------------------------------------------------------------------------------
  11376 11412 16      2.3    16      95      0       65    40  8 18 32 1.80    GTAAAGATTTGCACAT
  ...

  Indices      = Indices of the repeat relative to the start of the sequence.
  Period Size  = Period size of the repeat.
  Copy Number  = Number of copies aligned with the consensus pattern.
  Cons Size    = Size of consensus pattern (may differ slightly from the period size).
  Perc Matches = Percent of matches between adjacent copies overall.
  Perc Indels  = Percent of indels between adjacent copies overall.
  Score        = Alignment score.
  A, C, G, T   =  Percent composition for each of the four nucleotides.
  Entropy      = Entropy measure based on percent composition.
  Consensus    = Consensus sequence


Limitations:

  TRF can only accept one sequence in the input file. The input sequence is 
  limited to 5MB. 

=head1 SEE ALSO

=over 4

TRFResult

=back

=head1 COPYRIGHT

Copyright 2005 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package TRF;
use strict;
use POSIX qw(:sys_wait_h);
use Data::Dumper;
use File::Basename;
use File::Spec;
use Carp;
use FileHandle;
use IPC::Open3;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter SearchEngineI);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "TRF";
my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  croak $CLASS
      . "::new: Missing path to search engine!\n\n"
      . "use \$searchEngine = $CLASS->new( pathToEngine=>\"/usr/local/"
      . "bin/trf\")\n"
      if ( not defined $nameValuePairs{'pathToEngine'} );

  croak $CLASS
      . "::new: Missing a working directory!\n\n"
      . "use \$searchEngine = $CLASS->new( workDir=>\"/tmp\");\n"
      if ( not defined $nameValuePairs{'workDir'} );

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  $this->setPathToEngine( $nameValuePairs{'pathToEngine'} );
  $this->setWorkDir( $nameValuePairs{'workDir'} );

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setPathToEngine()

  Use: my $value    = getPathToEngine( );
  Use: my $oldValue = setPathToEngine( $value );

  Get/Set the fully qualified path to the search engine
  binary file.

=cut

##-------------------------------------------------------------------------##
sub getPathToEngine {
  my $this = shift;

  return $this->{'pathToEngine'};
}

sub setPathToEngine {
  my $this  = shift;
  my $value = shift;

  croak $CLASS. "::setPathToEngine( $value ): Program does not exist!"
      if ( not -x $value || `which $value` );

  my $oldValue = $this->{'pathToEngine'};
  $this->{'pathToEngine'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setMatchWeight()

  Use: my $value    = getMatchWeight( );
  Use: my $oldValue = setMatchWeight( $value );

  Get/Set the MatchWeight attribute.  

=cut

##---------------------------------------------------------------------##
sub getMatchWeight {
  my $this = shift;

  return $this->{'matchWeight'};
}

sub setMatchWeight {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'matchWeight'};
  $this->{'matchWeight'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setMismatchPenalty()

  Use: my $value    = getMismatchPenalty( );
  Use: my $oldValue = setMismatchPenalty( $value );

  Get/Set the MismatchPenalty attribute.  

=cut

##---------------------------------------------------------------------##
sub getMismatchPenalty {
  my $this = shift;

  return $this->{'mismatchPenalty'};
}

sub setMismatchPenalty {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'mismatchPenalty'};
  $this->{'mismatchPenalty'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setDelta()

  Use: my $value    = getDelta( );
  Use: my $oldValue = setDelta( $value );

  Get/Set the Delta attribute.  

=cut

##---------------------------------------------------------------------##
sub getDelta {
  my $this = shift;

  return $this->{'delta'};
}

sub setDelta {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'delta'};
  $this->{'delta'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setPm()

  Use: my $value    = getPm( );
  Use: my $oldValue = setPm( $value );

  Get/Set the Pm attribute.  

=cut

##---------------------------------------------------------------------##
sub getPm {
  my $this = shift;

  return $this->{'pm'};
}

sub setPm {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'pm'};
  $this->{'pm'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setPi()

  Use: my $value    = getPi( );
  Use: my $oldValue = setPi( $value );

  Get/Set the Pi attribute.  

=cut

##---------------------------------------------------------------------##
sub getPi {
  my $this = shift;

  return $this->{'pi'};
}

sub setPi {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'pi'};
  $this->{'pi'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setMinScore()

  Use: my $value    = getMinScore( );
  Use: my $oldValue = setMinScore( $value );

  Get/Set the MinScore attribute.  

=cut

##---------------------------------------------------------------------##
sub getMinScore {
  my $this = shift;

  return $this->{'minScore'};
}

sub setMinScore {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'minScore'};
  $this->{'minScore'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setMaxPeriod()

  Use: my $value    = getMaxPeriod( );
  Use: my $oldValue = setMaxPeriod( $value );

  Get/Set the MaxPeriod attribute.  

=cut

##---------------------------------------------------------------------##
sub getMaxPeriod {
  my $this = shift;

  return $this->{'maxPeriod'};
}

sub setMaxPeriod {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'maxPeriod'};
  $this->{'maxPeriod'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setMaskedOutputFile()

  Use: my $value    = getMaskedOutputFile( );
  Use: my $oldValue = setMaskedOutputFile( $value );

  Get/Set the MaskedOutputFile attribute.  

=cut

##---------------------------------------------------------------------##
sub getMaskedOutputFile {
  my $this = shift;

  return $this->{'maskedOutputFile'};
}

sub setMaskedOutputFile {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'maskedOutputFile'};
  $this->{'maskedOutputFile'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setWorkDir()

  Use: my $value    = getWorkDir( );
  Use: my $oldValue = setWorkDir( $value );

  Get/Set the WorkDir attribute.  

=cut

##---------------------------------------------------------------------##
sub getWorkDir {
  my $this = shift;

  return $this->{'workDir'};
}

sub setWorkDir {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'workDir'};
  $this->{'workDir'} = $value;

  return $oldValue;
}

##---------------------------------------------------------------------##

=head2 get_setSequenceFile()

  Use: my $value    = getSequenceFile( );
  Use: my $oldValue = setSequenceFile( $value );

  Get/Set the SequenceFile attribute.  

=cut

##---------------------------------------------------------------------##
sub getSequenceFile {
  my $this = shift;

  return $this->{'sequenceFile'};
}

sub setSequenceFile {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'sequenceFile'};
  $this->{'sequenceFile'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 search()

  Use: my ( $resultCode, $trfArrayRef ) = search( );
 or
  Use: my ( $resultCode, $trfArrayRef )
                          = search( minScore=> 300,
                                    ...
                                  );

  Run the search and return a perl array containing TRFResult.pm
  objects.

=cut

##-------------------------------------------------------------------------##
sub search {
  my $this           = shift;
  my %nameValuePairs = @_;

  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak( $CLASS . "::search: Instance variable $name doesn't exist." );
      }
      $this->$method( $value );
    }
  }

  # Test if engine is available
  my $engine = $this->getPathToEngine();
  if ( !defined $engine || !-f "$engine" ) {
    croak $CLASS
        . "::search: The path to the search engine is undefined or\n"
        . "is set incorrectly: $engine\n";
  }

  # Generate parameter line
  my $parameters = " ";
  my $value;
  my $outputPrefix = "";

  if ( defined( $value = $this->getSequenceFile() ) ) {
    if ( -f $value ) {
      $value = File::Spec->rel2abs( $value );
      $parameters .= "$value ";
      $outputPrefix = basename( $value );
    }
    else {
      croak $CLASS
          . "::search: Error...sequence file ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error sequence file undefined!\n";
  }

  if ( defined( $value = $this->getMatchWeight() ) ) {
    $parameters   .= "$value ";
    $outputPrefix .= ".$value";
  }
  else {
    $parameters   .= "2 ";
    $outputPrefix .= ".2";
  }

  if ( defined( $value = $this->getMismatchPenalty() ) ) {
    $parameters   .= "$value ";
    $outputPrefix .= ".$value";
  }
  else {
    $parameters   .= "7 ";
    $outputPrefix .= ".7";
  }

  if ( defined( $value = $this->getDelta() ) ) {
    $parameters   .= "$value ";
    $outputPrefix .= ".$value";
  }
  else {
    $parameters   .= "7 ";
    $outputPrefix .= ".7";
  }

  if ( defined( $value = $this->getPm() ) ) {
    $parameters   .= "$value ";
    $outputPrefix .= ".$value";
  }
  else {
    $parameters   .= "80 ";
    $outputPrefix .= ".80";
  }

  if ( defined( $value = $this->getPi() ) ) {
    $parameters   .= "$value ";
    $outputPrefix .= ".$value";
  }
  else {
    $parameters   .= "10 ";
    $outputPrefix .= ".10";
  }

  if ( defined( $value = $this->getMinScore() ) ) {
    $parameters   .= "$value ";
    $outputPrefix .= ".$value";
  }
  else {
    $parameters   .= "50 ";
    $outputPrefix .= ".50";
  }

  if ( defined( $value = $this->getMaxPeriod() ) ) {
    $parameters   .= "$value ";
    $outputPrefix .= ".$value";
  }
  else {
    $parameters   .= "500 ";
    $outputPrefix .= ".500";
  }

  if ( defined $this->getMaskedOutputFile() ) {
    $parameters .= "-m ";
  }

  $parameters .= "-d ";

  my $workDir = $this->getWorkDir();

  my $POUTPUT = new FileHandle;
  my $pid;

  $pid = open( $POUTPUT, "cd $workDir; $engine $parameters 2>/dev/null |" );

  if ( $DEBUG ) {
    my $outFile;
    do {
      $outFile = "$workDir/trfResults-" . time() . ".out";
    } while ( -f $outFile );
    open OUT, ">$outFile";
    while ( <$POUTPUT> ) {
      print OUT $_;
    }
    close OUT;
  }
  close $POUTPUT;
  my $trfResultsRef =
      parseOutput( searchOutput => "$workDir/$outputPrefix.dat" );

  unless ( $DEBUG ) {
    unlink( "$workDir/$outputPrefix.dat" )
        if ( -e "$workDir/$outputPrefix.dat" );
    unlink( "$workDir/$outputPrefix.1.html" )
        if ( -e "$workDir/$outputPrefix.1.html" );
    unlink( "$workDir/$outputPrefix.1.txt.html" )
        if ( -e "$workDir/$outputPrefix.1.txt.html" );
  }

  my $resultCode = $?;

  return ( $resultCode, $trfResultsRef );
}

##---------------------------------------------------------------------##

=head1 CLASS METHODS

=cut

##---------------------------------------------------------------------##

##---------------------------------------------------------------------##

=head2 parseOutput()

  Use: my $trfResultArrayRef = parseOutput(
                                     searchOutput => $filename|$FH,
                                          );

  Parse the result of a search and return a perl array containing
  TRFResult objects.

=cut

##---------------------------------------------------------------------##
sub parseOutput {
  my %nameValueParams = @_;

  croak $CLASS. "::parseOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  my $TRFFILE;
  if ( ref( $nameValueParams{'searchOutput'} ) !~ /GLOB|FileHandle/ ) {
    print $CLASS
        . "::parseOutput() Opening file "
        . $nameValueParams{'searchOutput'} . "\n"
        if ( $DEBUG );
    open $TRFFILE, $nameValueParams{'searchOutput'}
        or die $CLASS
        . "::parseOutput: Unable to open "
        . "results file: $nameValueParams{'searchOutput'} : $!";
  }
  else {
    $TRFFILE = $nameValueParams{'searchOutput'};
  }

  my @trfFieldKeys = qw( start end period copyNumber consSize
      percMatches percIndels score percA percC percG percT
      entropy consensus );

  my @results = ();
  while ( <$TRFFILE> ) {

    # Look for a header line (e.g.):
    if ( /^\s*\d+\s+\d+/ ) {
      print "Found header line: $_\n" if ( $DEBUG );
      my @dataArray = split;    # Create an array out of this line
           # Version 3.21 includes the subject sequence of the match.
           # For now lets toss this.
      pop @dataArray if ( $#dataArray == 14 );
      if ( $#dataArray == 13 ) {

        my %nVPairs = ();
        map {
          my $datum = shift( @dataArray );
          $nVPairs{$_} = $datum;
        } @trfFieldKeys;

        push @results, TRFResult->new( %nVPairs );

      }
    }

  }
  close $TRFFILE;

  return \@results;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _ucFirst( $string );
##
##   Uppercases the first character in a string and returns it.
##
##-------------------------------------------------------------------------##
sub _ucFirst {
  my $string = shift;

  if ( defined $string && $string ne "" ) {
    substr( $string, 0, 1 ) = uc( substr( $string, 0, 1 ) );
  }
  return $string;
}

##-------------------------------------------------------------------------##
## Serialization & Debug Routines
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $string = toString([$this]);
##
##      $this         : Normally passed implicitly
##
##  Returns
##
##      Uses the Data::Dumper to create a printable reprentation
##      of a data structure.  In this case the object data itself.
##
##-------------------------------------------------------------------------##
sub toString {
  my $this = shift;
  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  return $data_dumper->Dump();
}

##-------------------------------------------------------------------------##
## Use: my serializeOUT( $filename );
##
##	  $filename	: A filename to be created
##
##  Returns
##
##	Uses the Data::Dumper module to save out the data
##	structure as a text file.  This text file can be
##	read back into an object of this type.
##
##-------------------------------------------------------------------------##
sub serializeOUT {
  my $this     = shift;
  my $fileName = shift;

  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  open OUT, ">$fileName";
  print OUT $data_dumper->Dump();
  close OUT;
}

##-------------------------------------------------------------------------##
## Use: my serializeIN( $filename );
##
##	$filename	: A filename containing a serialized object
##
##  Returns
##
##	Uses the Data::Dumper module to read in data
##	from a serialized PERL object or data structure.
##
##-------------------------------------------------------------------------##
sub serializeIN {
  my $this         = shift;
  my $fileName     = shift;
  my $fileContents = "";
  my $oldSep       = $/;
  undef $/;
  my $in;
  open $in, "$fileName";
  $fileContents = <$in>;
  $/            = $oldSep;
  close $in;
  return eval( $fileContents );
}

1;
