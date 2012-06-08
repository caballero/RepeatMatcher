#!/usr/bin/perl

=head1 NAME

RepeatMatcherGUI.pl

=head1 DESCRIPTION

Create a graphic user interface to annotate RepeatModeler consensi.

=head1 USAGE

    perl RepeatMatcherGUI [PARAMETERS]

    Parameter     Description
    -o --out      Output file with final sequences (Fasta)
    -x --exclude  Output file with excluded sequences (Fasta)
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
use Tk::ROText;
use Tk::Label;
use Tk::ItemStyle;

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
my $exclude  = undef;

# Main variables
my $our_version = 0.1;        # Script version number
# %data => is the main struct
# $data{$id}
#            |-> label    : sequence fasta comment
#            |-> seq      : nucleotide sequence
#            |-> self     : sequence self-comparisons
#            |-> align    : sequence alignments to reference nucleotides
#            |-> blastx   : sequence alignments to reference peptides
#            |-> delete   : delete flag
#            |-> reverse  : reverse flag
#            |-> newlabel : new (edited) label for sequence
#            |-> question : question mark flag
#            |-> list     : position in Hlist
#            |-> status   : flag for finished sequences
my %data;
my %classes;
my @classes;
my $box_width  = 100;
my $box_height = 30;
my ($call_id, $delete, $reverse, $question);

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
    'r|reload:s'       => \$reload,
    'x|exclude:s'      => \$exclude
) or pod2usage(-verbose => 2);
printVersion() if (defined $version);    
pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $log or defined $reload);

if (defined $reload) {
    reloadProject($reload);
}
else {
    die "missing sequence file (-i)\n"        unless (defined $in);
    die "missing self-comparison file (-s)\n" unless (defined $self);
    die "missing alignments file (-a)\n"      unless (defined $align);
    die "missing blastx file (-b)\n"          unless (defined $blastx);
    die "missing output file (-o)\n"          unless (defined $out);
    die "missing exclude file (-e)\n"         unless (defined $exclude);
    startLog($log);
}

# load data from files in %data
loadIn();
loadSelf();
loadAlign();
loadBlastx();
loadRepClasses();
my @ids = sort (keys %data);

# Create the GUI
my $mw = MainWindow -> new; # Main Window
$mw -> title('RepeatMatcherGUI');

# Edit frame
my $edit_frame  = $mw         ->      Frame();
my $id_label    = $edit_frame ->      Label(-width  => 30);
my $class_entry = $edit_frame ->      Entry(-width  => 20, -background => 'white');
my $rep_button  = $edit_frame -> Menubutton(-text   => 'Repeat Family',
                                            -relief => 'raised'
                                              );
my $rep_menu    = $rep_button ->       Menu();
my $dna_menu    = $rep_menu   ->    cascade(-label  => 'DNA');
my $dnahat_menu = $rep_menu   ->    cascade(-label  => 'DNA/hAT');
my $dnatcm_menu = $rep_menu   ->    cascade(-label  => 'DNA/TcMar');
my $line_menu   = $rep_menu   ->    cascade(-label  => 'LINE');
my $ltr_menu    = $rep_menu   ->    cascade(-label  => 'LTR');
my $sine_menu   = $rep_menu   ->    cascade(-label  => 'SINE');
my $sat_menu    = $rep_menu   ->    cascade(-label  => 'Satellite');
my $other_menu  = $rep_menu   ->    cascade(-label  => 'Other');
foreach my $rep (@classes) {
    if    ($rep =~ m#^DNA/hAT#) {
        $dnahat_menu -> command(-label => $rep, -command => \&confClass($rep));
    }
    elsif ($rep =~ m#^DNA/TcMar#) {
        $dnatcm_menu -> command(-label => $rep, -command => \&confClass($rep));
    }
    elsif ($rep =~ m/^DNA/) {
        $dna_menu    -> command(-label => $rep, -command => \&confClass($rep));
    }
    elsif ($rep =~ m/^LINE/) {
        $line_menu   -> command(-label => $rep, -command => \&confClass($rep));
    }
    elsif ($rep =~ m/^LTR/) {
        $ltr_menu    -> command(-label => $rep, -command => \&confClass($rep));
    }
    elsif ($rep =~ m/^SINE/) {
        $sine_menu   -> command(-label => $rep, -command => \&confClass($rep));
    }
    elsif ($rep =~ m/^Satellite/) {
        $sat_menu    -> command(-label => $rep, -command => \&confClass($rep));
    }
    else {
        $other_menu  -> command(-label => $rep, -command => \&confClass($rep));
    }
}
$rep_button -> configure(-menu => $rep_menu);

my $qm_checkbox   = $edit_frame -> Checkbutton(-text       => '?', 
                                               -variable   => \$question,
                                               -command    => sub {
                                                                   $data{$call_id}{'question'} = $question;
                                                              });

my $info_entry    = $edit_frame ->       Entry(-width      => 80, 
                                               -background => 'white');

my $update_button = $edit_frame ->      Button(-text       => 'Update', 
                                               -command => \&updateSeq);

my $rev_checkbox  = $edit_frame -> Checkbutton(-text       => 'Reverse', 
                                               -variable => \$reverse,
                                               -command    => sub {
                                                                   $data{$call_id}{'reverse'} = $reverse;
                                                               });

my $del_checkbox  = $edit_frame -> Checkbutton(-text       => 'Exclude', 
                                               -variable => \$delete, 
                                               -command    => sub {
                                                                   $data{$call_id}{'delete'} = $delete;
                                                               }
                                              );
my $fasta_button  = $edit_frame ->      Button(-text       => 'Export Fasta', 
                                               -command    => \&writeFasta
                                              );

# IDs list frame
my $id_hlist      = $mw         ->    Scrolled('HList', 
                                               -itemtype   => 'text',
                                               -separator  => '#',
                                               -selectmode => 'single',
                                               -width      => 25, 
                                               -height     => int(2 * $box_height) - 2,
                                               -background => 'white',
                                               -browsecmd  => sub {
                                                                   my $call = shift;
                                                                   my @call = split(/#/, $call);
                                                                   $call_id = pop @call;
                                                                   &callID();
                                                               });

my $done_style   = $id_hlist    ->   ItemStyle('text', 
                                                -foreground => 'blue',
                                                -background => 'white');

my $done2_style  = $id_hlist    ->   ItemStyle('text', 
                                               -foreground => 'red',
                                               -background => 'white');

my $undone_style = $id_hlist    ->   ItemStyle('text', 
                                               -foreground => 'black', 
                                               -background => 'white');

$id_hlist -> add("#",  # root node
                 -text => "#", 
                 -style => $undone_style); 

foreach my $class (sort keys %classes) {
    $id_hlist -> add("#$class", 
                     -text => $class, 
                     -style => $undone_style);
}

foreach my $id (@ids) {
    my $class  = $data{$id}{'class'};
    my $status = $data{$id}{'status'};
    if ($status == 1) {
        if ($data{$id}{'label'} ne $data{$id}{'newlabel'} or 
            $data{$id}{'reverse'} == 1                    or 
            $data{$id}{'delete'}  == 1) {
                $id_hlist -> add("#$class#$id", 
                                 -text => $id, 
                                 -style => $done2_style);
        }
        else {
                $id_hlist -> add("#$class#$id", 
                                 -text => $id, 
                                 -style => $done_style);
        }
    }
    else {
                $id_hlist -> add("#$class#$id", 
                                 -text => $id, 
                                 -style => $undone_style);
    }
}


# Put all objects in a frame
my $data_frame    = $mw         ->       Frame();
# Sequence frame
my $seq_txt       = $data_frame ->    Scrolled('ROText', 
                                               -scrollbars => "osoe", 
                                               -width      => $box_width, 
                                               -height     => $box_height,
                                               -background => 'white');

# Align frame
my $align_txt     = $data_frame ->    Scrolled('ROText', 
                                               -scrollbars => "osoe", 
                                               -width      => $box_width, 
                                               -height     => $box_height,
                                               -background => 'white');

# Self frame
my $self_txt      = $data_frame ->    Scrolled('ROText', 
                                               -scrollbars => "osoe", 
                                               -width      => $box_width, 
                                               -height     => $box_height,
                                               -background => 'white');


# Blastx frame
my $blastx_txt    = $data_frame ->    Scrolled('ROText', 
                                               -scrollbars => "osoe", 
                                               -width      => $box_width, 
                                               -height     => $box_height,
                                               -background => 'white');

# Geometry
$id_label      -> grid (-row => 1, -column => 1);
$class_entry   -> grid (-row => 1, -column => 2);
$rep_button    -> grid (-row => 1, -column => 3);
$qm_checkbox   -> grid (-row => 1, -column => 4);
$info_entry    -> grid (-row => 1, -column => 5);
$rev_checkbox  -> grid (-row => 1, -column => 6);
$del_checkbox  -> grid (-row => 1, -column => 7);
$update_button -> grid (-row => 1, -column => 8);
$fasta_button  -> grid (-row => 1, -column => 9);
$edit_frame    -> grid (-row => 1, -column => 1, -columnspan => 9);
$id_hlist      -> grid (-row => 2, -column => 1);
$seq_txt       -> grid (-row => 1, -column => 2);
$align_txt     -> grid (-row => 1, -column => 3);
$self_txt      -> grid (-row => 2, -column => 2);
$blastx_txt    -> grid (-row => 2, -column => 3);
$data_frame    -> grid (-row => 2, -column => 2, 
                        -columnspan => 3, -rowspan => 2);

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
exclude_file: $exclude
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
        elsif (m/^exclude_file: (.+)/) { $exclude = $1; }
        elsif (m/^verbose_mode: (.+)/) { $verbose = $1; }
        else {
            my ($id, $del, $rev, $new) = split (/\t/, $_);
            $data{$id}{'delete'}   = $del;
            $data{$id}{'reverse'}  = $rev;
            if ($new =~ m/.+?#(.+?) /) {
                my $class = $1;
                $data{$id}{'class'}    = $class;
                $data{$id}{'newlabel'} = $new;
                $data{$id}{'question'} = 1 if ($class =~ m/\?/);
                $classes{$class}       = 1;
            }
            $data{$id}{'status'}   = 1;
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
            unless (defined $data{$id}{$class}) {
                $data{$id}{'class'}    = $class;
                $classes{$class} = 1;
            }
                
            $data{$id}{'delete'}   = 0  unless (defined $data{$id}{'delete'});
            $data{$id}{'reverse'}  = 0  unless (defined $data{$id}{'reverse'});
            $data{$id}{'newlabel'} = $_ unless (defined $data{$id}{'newlabel'});
            $data{$id}{'question'} = 0  unless (defined $data{$id}{'question'});
            $data{$id}{'status'}   = 0  unless (defined $data{$id}{'status'});            
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
    return unless ($call_id =~ m/rnd-\d+_family-\d+/);
    my $lab_     = $data{$call_id}{'label'};
    my $seq_     = 'No sequence';
    my $self_    = 'No matches';
    my $align_   = 'No matches';
    my $blastx_  = 'No matches';
    my $class_   = 'No class';
    
    $class_      = $data{$call_id}{'class'}    if (defined $data{$call_id}{'class'});
    $lab_        = $data{$call_id}{'newlabel'} if (defined $data{$call_id}{'newlabel'}); 
    $seq_        = $data{$call_id}{'seq'}      if (defined $data{$call_id}{'seq'});
    $self_       = $data{$call_id}{'self'}     if (defined $data{$call_id}{'self'});
    $align_      = $data{$call_id}{'align'}    if (defined $data{$call_id}{'align'});
    $blastx_     = $data{$call_id}{'blastx'}   if (defined $data{$call_id}{'blastx'});
    
    $seq_txt     -> selectAll;
    $seq_txt     -> deleteSelected;
    $seq_txt     -> insert('end', ">$lab_\n$seq_");
    
    $align_txt   -> selectAll;
    $align_txt   -> deleteSelected;
    $align_txt   -> insert('end', $align_);
    
    $self_txt    -> selectAll;
    $self_txt    -> deleteSelected;
    $self_txt    -> insert('end', $self_);
    
    $blastx_txt  -> selectAll;
    $blastx_txt  -> deleteSelected;
    $blastx_txt  -> insert('end', $blastx_);

    $id_label    -> configure(-text => "Repeat: $call_id");
    $class_entry -> configure(-text => $class_);
    $lab_ =~ s/^.+? //;
    $info_entry  -> configure(-text => $lab_);
     
    if ($data{$call_id}{'reverse'} == 1) {
        $rev_checkbox -> select;
        $reverse = 1;
    }
    else {
        $rev_checkbox -> deselect;
        $reverse = 0;
    }
    
    if ($data{$call_id}{'delete'}  == 1) {
        $del_checkbox -> select;
        $delete = 1;
    }
    else {
        $del_checkbox -> deselect;
        $delete = 0;
    }
    
    if ($data{$call_id}{'question'}  == 1) {
        $qm_checkbox -> select;
        $question = 1;
    }
    else {
        $qm_checkbox -> deselect;
        $question = 0;
    }    
}

sub updateSeq {
    my $rec = '';
    my $del = 0;
    my $rev = 0;
    my $lab = '';
    $del    = 1 if ($data{$call_id}{'delete'}  == 1);
    $rev    = 1 if ($data{$call_id}{'reverse'} == 1);
    
    my $class = $data{$call_id}{'class'};
    
    my $new_class = $class_entry -> get;
    my $new_info  = $info_entry  -> get;
    $new_class .= '?' if ($question == 1 and $new_class !~ m/\?/);
    my $new = "$call_id#$new_class $new_info";
    
    if ($new ne $data{$call_id}{'label'}) {
        $data{$call_id}{'newlabel'} = $new;
        $lab = $new;
    }
    
    open  LOG, ">>$log" or die "cannot open $log\n";
    print LOG "$call_id\t$del\t$rev\t$lab\n";
    warn "LOG> $call_id\t$del\t$rev\t$lab\n" if (defined $verbose);
    close LOG;
    
    if ($lab =~ m/\w+/ or $del == 1 or $rev == 1) {
        $id_hlist -> itemConfigure("#$class#$call_id", 0, -style => $done2_style);
    }
    else {
        $id_hlist -> itemConfigure("#$class#$call_id", 0, -style => $done_style);
    }

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
    open OUT, ">$out" or die "cannot write $out\n";
    open BAD, ">$exclude" or die "cannot write $exclude\n";
    warn "writing final sequences to $out\n" if (defined $verbose);
    warn "writing exclude sequences to $exclude\n" if (defined $verbose);
    foreach my $id (@ids) {
        my $lab = $data{$id}{'label'};
        $lab = $data{$id}{'newlabel'} if (defined $data{$id}{'newlabel'});
        my $seq = $data{$id}{'seq'};
        $seq = revcomp($seq) if ($data{$id}{'reverse'} == 1);
        if ($data{$id}{'delete'} == 1) {
            print BAD ">$lab\n$seq";
        }
        else {
            print OUT ">$lab\n$seq";
        }
    }
    close OUT;
    warn "Done.\n" if (defined $verbose);
}

sub confClass {
    my $rep = shift @_;
    $class_entry -> configure(-text => $rep);
}
 
sub loadRepClasses {
    warn "loading repeat families\n" if (defined $verbose);
    @classes = qw#
ARTEFACT
DNA
DNA/Academ
DNA/CMC-Chapaev
DNA/CMC-Chapaev-3
DNA/CMC-EnSpm
DNA/CMC-Mirage
DNA/CMC-Transib
DNA/Chapaev
DNA/Crypton
DNA/Ginger
DNA/Harbinger
DNA/Kolobok
DNA/Kolobok-Hydra
DNA/Kolobok-T2
DNA/MULE-F
DNA/MULE-MuDR
DNA/MULE-NOF
DNA/Maverick
DNA/Merlin
DNA/Novosib
DNA/P
DNA/P-Fungi
DNA/PIF-Harbinger
DNA/PIF-ISL2EU
DNA/PiggyBac
DNA/Sola
DNA/TcMar
DNA/TcMar-Ant1
DNA/TcMar-Fot1
DNA/TcMar-Gizmo
DNA/TcMar-ISRm11
DNA/TcMar-Mariner
DNA/TcMar-Mogwai
DNA/TcMar-Pogo
DNA/TcMar-Sagan
DNA/TcMar-Stowaway
DNA/TcMar-Tc1
DNA/TcMar-Tc2
DNA/TcMar-Tc4
DNA/TcMar-Tigger
DNA/TcMar-m44
DNA/Transib
DNA/Zator
DNA/Zisupton
DNA/hAT
DNA/hAT-Ac
DNA/hAT-Blackjack
DNA/hAT-Charlie
DNA/hAT-Pegasus
DNA/hAT-Restless
DNA/hAT-Tag1
DNA/hAT-Tip100
DNA/hAT-Tol2
DNA/hAT-hAT1
DNA/hAT-hAT5
DNA/hAT-hATm
DNA/hAT-hATw
DNA/hAT-hATx
DNA/hAT-hobo
LINE
LINE/Ambal
LINE/CR1
LINE/CR1-Zenon
LINE/CRE
LINE/DRE
LINE/Dong-R4
LINE/Genie
LINE/I
LINE/Jockey
LINE/L1
LINE/L1-Tx1
LINE/L2
LINE/L2-Hydra
LINE/LOA
LINE/Odin
LINE/Penelope
LINE/Proto1
LINE/Proto2
LINE/R1
LINE/R2
LINE/R2-Hero
LINE/RTE
LINE/RTE-BovB
LINE/RTE-RTE
LINE/RTE-X
LINE/Rex-Babar
LINE/Tad1
LINE/Zorro
LINE/telomeric
LTR
LTR/Caulimovirus
LTR/Copia
LTR/Copia(Xen1)
LTR/DIRS
LTR/ERV
LTR/ERV-Foamy
LTR/ERV-Lenti
LTR/ERV1
LTR/ERVK
LTR/ERVL
LTR/ERVL-MaLR
LTR/Gypsy
LTR/Gypsy-Troyka
LTR/Ngaro
LTR/Pao
LTR/TATE
LTR/Viper
Low_complexity
Other
Other/Composite
Other/DNA_virus
Other/centromeric
Other/subtelomeric
RC/Helitron
RNA
Retroposon
SINE
SINE/5S
SINE/7SL
SINE/Alu
SINE/B2
SINE/B4
SINE/BovA
SINE/C
SINE/CORE
SINE/Deu
SINE/Dong-R4
SINE/I
SINE/ID
SINE/L1
SINE/L2
SINE/MIR
SINE/Mermaid
SINE/R1
SINE/R2
SINE/RTE
SINE/RTE-BovB
SINE/Salmon
SINE/Sauria
SINE/V
SINE/tRNA
SINE/tRNA-7SL
SINE/tRNA-CR1
SINE/tRNA-Glu
SINE/tRNA-L2
SINE/tRNA-Lys
SINE/tRNA-R2
SINE/tRNA-RTE
Satellite
Satellite/W-chromosome
Satellite/Y-chromosome
Satellite/acromeric
Satellite/centromeric
Satellite/macro
Satellite/subtelomeric
Satellite/telomeric
Segmental
Simple_repeat
Unknown
Unknown/Y-chromosome
Unknown/centromeric
rRNA
scRNA
snRNA
tRNA
#;
}
