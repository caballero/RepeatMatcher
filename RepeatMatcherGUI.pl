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
use Tk::Hlist;

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
# %data => is the main struct
#       $data{$id}
#                  -> 'label'  : sequence fasta comment
#                  -> 'seq'    : nucleotide sequence
#                  -> 'self'   : sequence self-comparisons
#                  -> 'align'  : sequence alignments to reference nucleotides
#                  -> 'blastx' : sequence alignments to reference peptides
#                  -> 'status' : sequence analysis status
my %data;
my $box_width  = 100;
my $box_height = 30;
my $call_id;
my $delete;

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

if (defined $reload) {
    reloadProject($log);
}
else {
    startLog($log);
}

# load data from files in %data
loadIn();
loadSelf();
loadAlign();
loadBlastx();
my @ids = sort (keys %data);

# Create the GUI
my $mw = MainWindow -> new; # Main Window
$mw -> title('RepeatMatcherGUI');
my $data_frame      = $mw -> Frame();

# Edit frame
my $edit_frame      = $mw -> Frame();
my $label_entry     = $edit_frame -> Entry(-width => $box_width, -background => 'white');
my $update_button   = $edit_frame -> Button(-text => 'Update', -command => \&updateSeq);
my $revcomp_button  = $edit_frame -> Button(-text => 'ReverseSeq', -command => \&reverseSeq);
my $del_checkbox    = $edit_frame -> Checkbutton(-text => 'Delete', -variable => \$delete, -background => 'white');

# IDs list frame
my $id_hlist        = $mw -> Scrolled ('HList', 
                                        -itemtype   => 'text',
                                        -separator  => '/',
                                        -selectmode => 'single',
                                        -width      => 20, 
                                        -height     => 2 * $box_height,,
                                        -background => 'white',
                                        -browsecmd  => sub { 
                                                            $call_id = shift;
                                                            &callID();
                                                       }
                                        );
foreach my $id (@ids) {
    $id_hlist -> add($id, -text => $id);
}


# Sequence frame
my $seq_txt         = $data_frame -> Scrolled('Text', 
                                               -scrollbars => "osoe", 
                                               -width      => $box_width, 
                                               -height     => $box_height,
                                               -background => 'white');

# Align frame
my $align_txt       = $data_frame -> Scrolled('Text', 
                                               -scrollbars => "osoe", 
                                               -width      => $box_width, 
                                               -height     => $box_height,
                                               -background => 'white');

# Self frame
my $self_txt        = $data_frame -> Scrolled('Text', 
                                               -scrollbars => "osoe", 
                                               -width      => $box_width, 
                                               -height     => $box_height,
                                               -background => 'white');


# Blastx frame
my $blastx_txt      = $data_frame -> Scrolled('Text', 
                                               -scrollbars => "osoe", 
                                               -width      => $box_width, 
                                               -height     => $box_height,
                                               -background => 'white');

# Geometry managment
$label_entry    -> grid (-row => 1, -column => 1);
$update_button  -> grid (-row => 1, -column => 2);
$revcomp_button -> grid (-row => 1, -column => 3);
$del_checkbox   -> grid (-row => 1, -column => 4);
$edit_frame     -> grid (-row => 1, -column => 1, -columnspan => 4);

$id_hlist       -> grid (-row => 2, -column => 1);

$seq_txt        -> grid (-row => 1, -column => 2);
$align_txt      -> grid (-row => 1, -column => 3);
$self_txt       -> grid (-row => 2, -column => 2);
$blastx_txt     -> grid (-row => 2, -column => 3);
$data_frame     -> grid (-row => 2, -column => 2, -columnspan => 3, -rowspan => 2);

MainLoop();

###################################
####   S U B R O U T I N E S   ####
###################################

sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

sub startLog {
    my $file = shift @_;
    my $ver  = 0;
    $ver = 1 if (defined $verbose); 
    open LOG, ">$file" or die "cannot open $file\n";
    print LOG <<_LOG_   
# RepeatMatcherGUI log file
seq_file: $in
out_file: $out
self_file: $self
align_file: $align
blastx_file: $blastx
verbose_mode: $ver
_LOG_
;
}

sub reloadProject {
    my $file = shift @_;
    open LOG, "+>$file" or die "cannot open $file\n";
    while (<LOG>) {
        chomp;
        next if (m/^#/);
        if    (m/^seq_file: (.+)/)     { $in      = $1; }
        elsif (m/^out_file: (.+)/)     { $out     = $1; }
        elsif (m/^self_file: (.+)/)    { $self    = $1; }
        elsif (m/^align_file: (.+)/)   { $align   = $1; }
        elsif (m/^blastx_file: (.+)/)  { $blastx  = $1; }
        elsif (m/^verbose_mode: (.+)/) { $verbose = $1; }
        else {
            my ($id, $status) = split (/: /, $_);
            $data{$id}{'status'} = $status;
        }
    }
    $verbose = undef if ($verbose == 0);
}

sub loadIn {
    warn "Loading sequences from $in\n" if (defined $verbose);
    open F, "$in" or die "cannot open file $in\n";
    my $id;
    while (<F>) {
        if (m/>(.+?)#/) {
            $id = $1;
            s/>//;
            chomp;
            $data{$id}{'label'} .= $_;
            $data{$id}{'status'} = 'u' unless (defined $data{$id}{'status'});
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
            my @hit = split (/\n/, $_);
            pop @hit;
            while (1) {
                my $del = shift @hit;
                last if ($del =~ /Searching/);
            }
            shift @hit;
            shift @hit;
            shift @hit;
            $data{$id}{'blastx'} .= join ("\n", @hit);
        }
    }
    close B;
}

sub revcomp {
    my $rc  =  '';
    my $sq  =  shift @_;
       $sq  =~ s/\n+//;
       $sq  =  reverse $sq;
       $sq  =~ tr/ACGTacgt/TGCAtgca/;
    while ($sq) {
        $rc .= substr ($sq, 0, 50);
        $rc .= "\n";
        substr ($sq, 0, 50) = '';
    }
    return $rc;
}

sub callID {
    my $lab_     = $call_id;
    my $seq_     = 'No sequence';
    my $self_    = 'No matches';
    my $align_   = 'No matches';
    my $blastx_  = 'No matches';
    
    $lab_        = $data{$call_id}{'label'}  if (defined $data{$call_id}{'label'}); 
    $seq_        = $data{$call_id}{'seq'}    if (defined $data{$call_id}{'seq'});
    $self_       = $data{$call_id}{'self'}   if (defined $data{$call_id}{'self'});
    $align_      = $data{$call_id}{'align'}  if (defined $data{$call_id}{'align'});
    $blastx_     = $data{$call_id}{'blastx'} if (defined $data{$call_id}{'blastx'});
    
    $seq_txt    -> selectAll;
    $seq_txt    -> deleteSelected;
    $seq_txt    -> insert('end', ">$lab_\n$seq_");
    
    $align_txt  -> selectAll;
    $align_txt  -> deleteSelected;
    $align_txt  -> insert('end', $align_);
    
    $self_txt   -> selectAll;
    $self_txt   -> deleteSelected;
    $self_txt   -> insert('end', $self_);
    
    $blastx_txt -> selectAll;
    $blastx_txt -> deleteSelected;
    $blastx_txt -> insert('end', $blastx_);
 
    $label_entry -> configure(-text => $lab_);
}

sub reverseSeq {
    $data{$call_id}{'seq'} = revcomp($data{$call_id}{'seq'});
    callID();
}

sub updateSeq {
    
}
