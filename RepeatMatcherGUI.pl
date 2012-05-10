#!/usr/bin/perl

=head1 NAME

RepeatMatcherGUI.pl

=head1 DESCRIPTION

Create a graphic user interface to annotate RepeatModeler consensi.

=head1 USAGE

    perl RepeatMatcherGUI [PARAMETERS]

    Parameter     Description
    -o --out      Output file (Fasta)
    -i --input    Input file (Fasta)
    -s --self     Self-comparison of input (cross_match output)
    -a --align    Alignments of input to a reference (cross_match output)
    -b --blastx   Blast output of the input to reference petides (blast output)
    -l --log      Write log file here (used to keep track of progress)
    -r --reload   Reload a previous project

    -v --verbose  Verbose mode
    -h --help     Print this screen
    --version     Print version
    
=head1 EXAMPLES

    1. Start a new project
    perl RepeatMatcherGUI.pl -o OUT -i SEQS -s SELF -a ALIGN -b BLASTX -l LOG

    2. Reload a started project
    perl RepeatMatcherGUI.pl -r LOG

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2012

=head1 CONTACT

jcaballero@systemsbiology.org

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with code.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Tk;

# Default parameters
my $help     = undef;         # Print help
my $verbose  = undef;         # Verbose mode
my $version  = undef;         # Version call flag
my $out      = undef;
my $in       = undef;
my $self     = undef;
my $align    = undef;
my $blastx   = undef;
my $log      = undef;
my $reload   = undef;


# Main variables
my $our_version = 0.1;        # Script version number
my %data;
my ($cnt, $id);

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'version'          => \$version,
    'o|out=s'          => \$out,
    'i|in=s'           => \$in,
    's|self=s'         => \$self,
    'a|align=s'        => \$align,
    'b|blastx=s'       => \$blastx,
    'l|log=s'          => \$log,
    'r|reload:s'       => \$reload
) or pod2usage(-verbose => 2);
printVersion() if (defined $version);    
pod2usage(-verbose => 2) if (defined $help);

reloadProject() 
if (defined $reload) {
    reloadProject($log);
}
else {
    startLog($log);
}


###################################
####   S U B R O U T I N E S   ####
###################################

sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

sub startLog {
    
}

sub reloadProject {
    
}

sub loadIn {
    warn "Loading sequences from $in\n" if (defined $verbose);
    open F, "$in" or die "cannot open file $in\n";
    my $id;
    while (<F>) {
        if (m/>(.+?)#/) {
            $id = $1;
            $data{$id}{'label'} .= $_;
        }
        else {
            $data{$id}{'seq'}   .= $_;
        }
    }
    close F;
}

sub loadAlign {
    warn "loading aligments in $align\n" if (defined $verbose);
    open A, "$align" or die "cannot open file $align\n";
    my $id = 'skip';
    while (<A>) {
       if (m/^\s*\d+.+(rnd-\d+_family-\d+)#/) {
           $id = $1;
       }
       $data{$id}{'align'} .= $_ if ($id ne 'skip');
    }
    close A;
}

sub loadSelf {
    warn "loading self-comparison in $self\n" if (defined $verbose);
    my $id;
    open S, "$self" or die "cannot open file $self\n";
    while (<S>) {
       if (m/^\s*\d+.+(rnd-\d+_family-\d+)#/) {
           $id = $1;
           $data{$id}{'self'} .= $_;
       }
    }
    close S;
}

sub loadBlastx {
    warn "loading blastx aligments in $blastx\n" if (defined $verbose);
    my $id;
    open B, "$blastx" or die "cannot open file $blastx\n";
    local $/ = "\nBLASTX";
    while (<B>) {
        m/Query= (rnd-\d+_family-\d+)#/;
        $id = $1;
        if (m/No hits found/) {
            $data{$id}{'blastx'} .= 'No hits found';
        }
        else {
            s/.+Searching\.+done\n\n//;
            $data{$id}{'blastx'} .= $_;
        }
    }
    close B;
}

sub revcomp {
    my $in  =  shift @_;
    my @ln  =  split (/\n/, $in);
    my $rc  =  shift @ln;
       $rc .=  "\n";
    my $sq  =  join "", @ln;
       $sq  =  reverse $sq;
       $sq  =~ tr/ACGTacgt/TGCAtgca/;
    while ($sq) {
        $rc .= substr ($sq, 0, 50), "\n";
        substr ($sq, 0, 50) = '';
    }
    return $rc;
}
