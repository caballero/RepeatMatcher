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
#                  -> label    : sequence fasta comment
#                  -> seq      : nucleotide sequence
#                  -> self     : sequence self-comparisons
#                  -> align    : sequence alignments to reference nucleotides
#                  -> blastx   : sequence alignments to reference peptides
#                  -> delete   : delete flag
#                  -> reverse  : reverse flag
#                  -> newlabel : new (edited) label for sequence
my %data;
my $box_width  = 90;
my $box_height = 30;
my ($call_id, $delete, $reverse);

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
pod2usage(-verbose => 2) unless (defined $log or defined $reload);

if (defined $reload) {
    reloadProject($reload);
}
else {
    die "missing sequence file (-i)\n" unless (defined $in);
    die "missing self-comparison file (-s)\n" unless (defined $self);
    die "missing alignments file (-a)\n" unless (defined $align);
    die "missing blastx file (-b)\n" unless (defined $blastx);
    die "missing output file (-o)\n" unless (defined $out);
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

# Edit frame
my $edit_frame    = $mw         ->       Frame();
my $label_entry   = $edit_frame ->       Entry(
                                                -width      => int(1.5 * $box_width), 
                                                -background => 'white'
                                              );
my $update_button = $edit_frame ->      Button(
                                                -text       => 'Update', 
                                                -command    => \&updateSeq
                                              );
my $rev_checkbox  = $edit_frame -> Checkbutton(
                                                -text       => 'Reverse', 
                                                -variable   => \$reverse,
                                                -background => 'white',
                                                -command    => sub {
                                                               $data{$call_id}{'reverse'} = $reverse;
                                                               }
                                              );
my $del_checkbox  = $edit_frame -> Checkbutton(
                                                -text       => 'Delete',
                                                -variable   => \$delete, 
                                                -background => 'white',
                                                -command    => sub {
                                                               $data{$call_id}{'delete'} = $delete;
                                                               }
                                              );
my $fasta_button  = $edit_frame ->      Button(
                                                -text       => 'Export Fasta', 
                                                -command    => \&writeFasta
                                              );

# IDs list frame
my $id_hlist      = $mw         ->    Scrolled(
                                                'HList', 
                                                -itemtype   => 'text',
                                                -separator  => '/',
                                                -selectmode => 'single',
                                                -width      => 25, 
                                                -height     => int(2 * $box_height),
                                                -background => 'white',
                                                -browsecmd  => sub { 
                                                               $call_id = shift;
                                                               &callID();
                                                               }
                                              );
foreach my $id (@ids) {
    my $class = $data{$id}{'class'};
    $id_hlist -> add($id, -text => "$id#$class");
}


# Put all objects in a frame
my $data_frame    = $mw         ->       Frame();
# Sequence frame
my $seq_txt       = $data_frame ->    Scrolled(
                                                'Text', 
                                                -scrollbars => "osoe", 
                                                -width      => $box_width, 
                                                -height     => $box_height,
                                                -background => 'white'
                                              );

# Align frame
my $align_txt     = $data_frame ->    Scrolled(
                                                'Text', 
                                                -scrollbars => "osoe", 
                                                -width      => $box_width, 
                                                -height     => $box_height,
                                                -background => 'white'
                                              );

# Self frame
my $self_txt      = $data_frame ->    Scrolled(
                                                'Text', 
                                                -scrollbars => "osoe", 
                                                -width      => $box_width, 
                                                -height     => $box_height,
                                                -background => 'white'
                                              );


# Blastx frame
my $blastx_txt    = $data_frame ->    Scrolled(
                                                'Text', 
                                                -scrollbars => "osoe", 
                                                -width      => $box_width, 
                                                -height     => $box_height,
                                                -background => 'white'
                                              );

# Geometry
$label_entry      -> grid (-row => 1, -column => 1);
$rev_checkbox     -> grid (-row => 1, -column => 2);
$del_checkbox     -> grid (-row => 1, -column => 3);
$update_button    -> grid (-row => 1, -column => 4);
$fasta_button     -> grid (-row => 1, -column => 6);
$edit_frame       -> grid (-row => 1, -column => 1, -columnspan => 6);
$id_hlist         -> grid (-row => 2, -column => 1);
$seq_txt          -> grid (-row => 1, -column => 2);
$align_txt        -> grid (-row => 1, -column => 3);
$self_txt         -> grid (-row => 2, -column => 2);
$blastx_txt       -> grid (-row => 2, -column => 3);
$data_frame       -> grid (-row => 2, -column => 2, -columnspan => 3, -rowspan => 2);

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
    warn "creating LOG in $file\n" if (defined $verbose);
    open LOG, ">$file" or die "cannot open $file\n";
    print LOG <<_LOG_   
# RepeatMatcherGUI log file
seq_file: $in
out_file: $out
log_file: $log
self_file: $self
align_file: $align
blastx_file: $blastx
verbose_mode: $ver
_LOG_
;
close LOG;
}

sub reloadProject {
    my $file = shift @_;
    warn "reading info from LOG in $file\n" if (defined $verbose);
    open LOG, "$file" or die "cannot open $file\n";
    while (<LOG>) {
        chomp;
        next if (m/^#/);
        if    (m/^seq_file: (.+)/)     { $in      = $1; }
        elsif (m/^out_file: (.+)/)     { $out     = $1; }
        elsif (m/^log_file: (.+)/)     { $log     = $1; }
        elsif (m/^self_file: (.+)/)    { $self    = $1; }
        elsif (m/^align_file: (.+)/)   { $align   = $1; }
        elsif (m/^blastx_file: (.+)/)  { $blastx  = $1; }
        elsif (m/^verbose_mode: (.+)/) { $verbose = $1; }
        else {
            my ($id, $del, $rev, $new) = split (/\t/, $_);
            $data{$id}{'delete'}   = $del;
            $data{$id}{'reverse'}  = $rev;
            $data{$id}{'newlabel'} = $new;
        }
    }
    $verbose = undef if ($verbose == 0);
    close LOG;
}

sub loadIn {
    warn "Loading sequences from $in\n" if (defined $verbose);
    open FASTA, "$in" or die "cannot open file $in\n";
    my $id;
    my $class;
    while (<FASTA>) {
        if (m/>(.+?)#(.+?) /) {
            $id    = $1;
            $class = $2;
            s/>//;
            chomp;
            $data{$id}{'label'}    = $_;
            $data{$id}{'class'}    = $class;
            $data{$id}{'delete'}   = 0  unless (defined $data{$id}{'delete'});
            $data{$id}{'reverse'}  = 0  unless (defined $data{$id}{'reverse'});
            $data{$id}{'newlabel'} = $_ unless (defined $data{$id}{'newlabel'});
        }
        else {
            $data{$id}{'seq'}   .= $_;
        }
    }
    close FASTA;
}

sub loadAlign {
    warn "loading aligments in $align\n" if (defined $verbose);
    open ALIGN, "$align" or die "cannot open file $align\n";
    my $id = 'skip';
    while (<ALIGN>) {
       if (m/^\s*\d+.+(rnd-\d+_family-\d+)#/) {
           $id = $1;
       }
       $data{$id}{'align'} .= $_ if ($id ne 'skip');
    }
    close ALIGN;
}

sub loadSelf {
    warn "loading self-comparison in $self\n" if (defined $verbose);
    my $id;
    open SELF, "$self" or die "cannot open file $self\n";
    while (<SELF>) {
       if (m/^\s*\d+.+(rnd-\d+_family-\d+)#/) {
           $id = $1;
           $data{$id}{'self'} .= $_;
       }
    }
    close SELF;
}

sub loadBlastx {
    warn "loading blastx aligments in $blastx\n" if (defined $verbose);
    my $id;
    open BLAST, "$blastx" or die "cannot open file $blastx\n";
    local $/ = "\nBLASTX";
    while (<BLAST>) {
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
    close BLAST;
}

sub callID {
    my $lab_      = $data{$call_id}{'label'};
    my $seq_      = 'No sequence';
    my $self_     = 'No matches';
    my $align_    = 'No matches';
    my $blastx_   = 'No matches';
    
    $lab_         = $data{$call_id}{'newlabel'} if (defined $data{$call_id}{'newlabel'}); 
    $seq_         = $data{$call_id}{'seq'}      if (defined $data{$call_id}{'seq'});
    $self_        = $data{$call_id}{'self'}     if (defined $data{$call_id}{'self'});
    $align_       = $data{$call_id}{'align'}    if (defined $data{$call_id}{'align'});
    $blastx_      = $data{$call_id}{'blastx'}   if (defined $data{$call_id}{'blastx'});
    
    $seq_txt      -> selectAll;
    $seq_txt      -> deleteSelected;
    $seq_txt      -> insert('end', ">$lab_\n$seq_");
    
    $align_txt    -> selectAll;
    $align_txt    -> deleteSelected;
    $align_txt    -> insert('end', $align_);
    
    $self_txt     -> selectAll;
    $self_txt     -> deleteSelected;
    $self_txt     -> insert('end', $self_);
    
    $blastx_txt   -> selectAll;
    $blastx_txt   -> deleteSelected;
    $blastx_txt   -> insert('end', $blastx_);
 
    $label_entry  -> configure(-text => $lab_);
    
     
    if ($data{$call_id}{'reverse'} == 1) {
        $rev_checkbox -> select;
        $reverse = 1;
    }
    else {
        $reverse = 0;
    }
    
    if ($data{$call_id}{'delete'}  == 1) {
        $del_checkbox -> select;
        $delete = 1;
    }
    else {
        $delete = 0;
    }
}


sub updateSeq {
    my $rec = '';
    my $del = 0;
    my $rev = 0;
    my $lab = '';
    $del    = 1 if ($data{$call_id}{'delete'}  == 1);
    $rev    = 1 if ($data{$call_id}{'reverse'} == 1);
    
    my $new = $label_entry -> get;
    if ($new ne $data{$call_id}{'label'}) {
        $data{$call_id}{'newlabel'} = $new;
        $lab = $new;
    }
    
    open  LOG, ">>$log" or die "cannot open $log\n";
    print LOG "$call_id\t$del\t$rev\t$lab\n";
    warn "LOG> $call_id\t$del\t$rev\t$lab\n" if (defined $verbose);
    close LOG;
}

sub revcomp {
    my $rc  =  '';
    my $sq  =  shift @_;
       $sq  =~ s/\n+//g;
       $sq  =  reverse $sq;
       $sq  =~ tr/ACGTacgt/TGCAtgca/;
    while ($sq) {
        $rc .= substr ($sq, 0, 50);
        $rc .= "\n";
        substr ($sq, 0, 50) = '';
    }
    return $rc;
}

sub writeFasta {
    open OUT, ">$out" or die "cannot write $out";
    warn "writing sequences to $out\n" if (defined $verbose);
    foreach my $id (@ids) {
        next if ($data{$id}{'delete'} == 1);
        my $lab = $data{$id}{'label'};
        $lab = $data{$id}{'newlabel'} if (defined $data{$id}{'newlabel'});
        my $seq = $data{$id}{'seq'};
        $seq = revcomp($seq) if ($data{$id}{'reverse'} == 1);
        print OUT ">$lab\n$seq";
    }
    close OUT;
    warn "    done\n" if (defined $verbose);
}
