#!/usr/bin/perl

=head1 NAME

RepeatMatcher.pl

=head1 DESCRIPTION

Automatic pipeline to annotate repeats from RepeatModeler.

WARNING: you must have RepeatMasker, R:, Blast and Cross_Match in your system.

=head1 USAGE

RepeatMatcher.pl [PARAMETERS]

     Parameters       Description                 Type
     -s --sequences   RepeatModeler sequences     Fasta
     -k --known       Known repeat annotation     Fasta
     -c --config      Configuration file          Text
     -o --out         Base name for files         
     -e --edit        Edit configuration
     -d --demo        Don't run the commands, show the commands instead
     
     Skipping some steps:
     --no-low         Don't mask low complexity repeats
     --no-self-comp   Don't self compare repeats
     --no-known-comp  Don't compare with the annotation
     --no-rep-blast   Don't blast to repeat peptides
     --no-nr-blast    Don't run blast to NR database
     --no-fold        Don't fold the sequence
     
     Other parameters:
     -h --help        Print this screen
     -v --verbose     Verbose mode
     --version        Print version information and exit
     
=head1 EXAMPLES

     RepeatMatcher.pl -s RepeatModeler.fa -k RepBase.fa -o NewAnnotation

=head1 AUTHOR

Juan Caballero, RepeatMaker.org - Institute for Systems Biology @ 2012

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

# Default parameters
my $help          = undef;
my $verbose       = undef;
my $version       = undef;
my $orig_seq      = undef;
my $out           = undef;
my $known_seq     = undef;
my $no_low        = undef;
my $no_self_comp  = undef;
my $no_known_comp = undef;
my $no_rep_blast  = undef;
my $no_nr_blast   = undef;
my $no_fold       = undef;
my $edit_conf     = undef;
my $conf_file     = 'RepeatMatcher.conf';
my $demo          = undef;

# Exec programs location/call, change it if they're not in your PATH
my $repeatmasker = 'RepeatMasker';
my $crossmatch   = 'cross_match';
my $blast        = 'blastall';
my $dnafold      = 'RNAfold';
my $fold2png     = 'plotR4RNA.sh';

# Main variables
my $our_version = 0.1;
my %conf;
my %seq;
my %aln;
my ($mask_seq, $self_comp, $known_comp, $rep_blast, $nr_blast, $fold);

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    's|sequences=s'    => \$orig_seq,
    'k|known:s'        => \$known_seq,
    'c|config:s'       => \$conf_file,
    'o|out:s'          => \$out,
    'e|edit'           => \$edit_conf,
    'no-low'           => \$no_low,
    'no-self-comp'     => \$no_self_comp,
    'no-known-comp'    => \$no_known_comp,
    'no-rep-blast'     => \$no_rep_blast,
    'no-nr-blast'      => \$no_nr_blast,
    'no-fold'          => \$no_fold,
    'd|demo'           => \$demo
) or pod2usage(-verbose => 2);
printVersion()           if  (defined $version);    
pod2usage(-verbose => 2) if  (defined $help);
pod2usage(-verbose => 2) if !(defined $orig_seq);

unless (defined $out) {
    $out = $orig_seq;
    $out =~ s/\.\w+?$//;
}

# Load Configuration
readConfig($conf_file);
editConfig() if (defined $edit_conf);

# STEP 1. Mask low complexity sequences
if (defined $no_low) {
    $mask_seq = $orig_seq;
}
else { # mask low complexity sequences, remove bad sequences
    $mask_seq = "$out.mask";
    maskLow($orig_seq, $mask_seq);
}

# STEP 2. Self-check-up
unless (defined $no_self_comp) {
    $self_comp = "$out.self";
    selfComp($mask_seq, $self_comp);
}   

# STEP 3. Sequence alignments
unless (defined $no_known_comp) {
    $known_comp = "$out.known";
    annComp($mask_seq, $known_seq, $known_comp);
}

# STEP 4. Blast to known repeats
unless (defined $no_rep_blast) {
    $rep_blast = "$out.repblast";
    blastxRep($mask_seq, $rep_blast);
}

# STEP 5. Blast to NR
unless (defined $no_nr_blast) {
    $nr_blast = "$out.nrblast";
    blastxNR($mask_seq, $nr_blast);
}

# STEP 6. Fold sequences
unless (defined $no_fold) {
    $fold = "$out.fold";
    foldSeq($mask_seq, $fold);
}

warn "RepeatMatcher => done\n" if (defined $verbose);

###################################
####   S U B R O U T I N E S   ####
###################################

sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

sub readConfig {
    my $in = shift @_;
    die  "Configuration file is missing\n"  if !(-e $in);
    warn "Reading configuration from $in\n" if  (defined $verbose);
    open CONF, "$in" or die "cannot open file $in\n";
    while (<CONF>) {
        next if (m/^#/);
        next if (m/^\n/);
        chomp;
        my ($var, $val) = split (/: /, $_);
        $conf{$var} = $val;
    }
    close CONF;
}

sub editConfig {
    my @par = qw/min_mask min_size crossmatch_self crossmatch_comp blastx_rep blastx_nr fold plot_w plot_h/;
    my $par = undef;
    my $res = undef;
    print "Manual configuration of parameters\n";
    foreach $par (@par) {
        print "$par: ", $conf{$par}, ", edit? [y/N] ";
        chomp($res = <>);
        if ($res =~ m/y/i) {
            print "$par: ";
            chomp($res  = <>);
            $conf{$par} = $res;
        }
    }
    print "Save new parameters? [y/N] ";
    chomp($res  = <>);
    if ($res =~ m/y/i) {
        print "File name: ";
        chomp($res  = <>);
        open CONF, ">$res" or die "cannot write $res\n";
        foreach $par (@par) {
            print CONF "$par: ", $conf{$par}, "\n";
        }
        close CONF;
    }
}

sub maskLow {
    my ($in, $msk) = @_;
    my $min_msk    = $conf{'min_mask'};
    my $min_len    = $conf{'min_size'};
    warn "Masking low complexity sequences in $in\n" if (defined $verbose);
    my $cmd = "$repeatmasker -low $in";
    warn "CMD: $cmd\n" if (defined $verbose);
    if (defined $demo) {
        print "$cmd\n";
        return;
    }
    system ($cmd);
    
    die "cannot run $cmd\n" unless (-e "$in.masked");
    warn "Filtering sequences in $in.masked\n" if (defined $verbose);
    open O, ">$msk" or die "cannot open file $msk\n";
    open F, "$in.masked" or die "cannot open file $in.masked\n";
    local $/ = "\n>"; # slurp mode for Fasta
    while (<F>) {
        s/>//g;
        my @seq  =  split (/\n/, $_);
        my $id   =  shift @seq; # remove header
        my $seq  =  join "", @seq;
           $seq  =~ s/[^ACGT]/N/g;
        my $nnum =  $seq =~ tr/N/N/;
        my $len  =  length $seq;
        my $freq =  100 * $nnum / $len;
        # print if sequence has >$min_len effective bases and is <$min_msk masked
        if ( $seq =~ m/[ACGT]{$min_len,}/ and $freq < $min_msk) {
            print O ">$_" 
        }
        else {
            warn "$id: rejected, low complexity\n" if (defined $verbose);
        }
    }
    close O;
    close F;
    
    # clean up
    my @suffix = qw/cat log out ref tbl masked/;
    foreach my $suffix (@suffix) {
        unlink "$in.$suffix";
    }
}

sub selfComp {
    my ($in, $out) = @_;
    warn "Running self-comparison for $in\n" if (defined $verbose);
    my $param = $conf{'crossmatch_self'};
    my $cmd = "$crossmatch $in $param > $out";
    warn "CMD: $cmd\n" if (defined $verbose);
    if (defined $demo) {
        print "$cmd\n";
        return;
    }
    system ($cmd);
    die "cannot run $cmd" if (-z $out);
}

sub annComp {
    my ($in, $ann, $out) = @_;
    warn "Running annotation comparison for $in\n" if (defined $verbose);
    my $param = $conf{'crossmatch_comp'};
    my $cmd = "$crossmatch $in $ann $param -alignments > $out";
    warn "CMD: $cmd\n" if (defined $verbose);
    if (defined $demo) {
        print "$cmd\n";
        return;
    }
    system ($cmd);
    die "cannot run $cmd" if (-z $out);
}

sub blastxRep {
    my ($in, $out) = @_;
    warn "Comparing $in to repeat proteins with BlastX\n" if (defined $verbose);
    my $param = $conf{'blastx_rep'};
    my $cmd = "$blast -p blastx -i $in -o $out $param";
    warn "CMD: $cmd\n" if (defined $verbose);
    if (defined $demo) {
        print "$cmd\n";
        return;
    }
    system ($cmd);
    die "cannot run $cmd\n" if (-z $out);
}

sub blastxNR {
    my ($in, $out) = @_;
    warn "Comparing $in to NR with BlastX\n" if (defined $verbose);
    my $param = $conf{'blastx_nr'};
    my $cmd = "$blast -p blastx -i $in -o $out $param";
    warn "CMD: $cmd\n" if (defined $verbose);
    if (defined $demo) {
        print "$cmd\n";
        return;
    }
    system ($cmd);
    die "cannot run $cmd\n" if (-z $out);
}

sub foldSeq {
    my ($in, $out) = @_;
    warn "Folding sequences in $in\n" if (defined $verbose);
    unless (-d $out) {
        mkdir $out or die "cannot create directory\n";
    }
    my $fold_param   = $conf{'fold'};
    my $ps2png_param = $conf{'ps2png'};
    local $/ = "\n>"; # Fasta slurp
    open F, "$in" or die "cannot open file\n";
    while (<F>) {
        s/>//g;
        my @seq = split (/\n/, $_);
        my $sid = shift @seq;
        $sid =~ s/#.*//g;
        my $seq = join "", @seq;
        my $cmd = "echo $seq | $dnafold $fold_param > $out/$sid.fold";

        warn "CMD: $cmd\n"  if (defined $verbose);
        if (defined $demo) {
            print "$cmd\n";
        }
        else {
            system ($cmd);
            die "cannot run $cmd\n" if (-z "$out/$sid.fold");
        }
    }
    close F;
    convertFold2PNG($out); 
}

sub convertFold2PNG {
    my $dir = shift;
    my $w   = $conf{'plot_w'}; $w ||= 600;
    my $h   = $conf{'plot_h'}; $h ||= 300;
    my $cmd = 'R --vanilla --quiet < R.scr';
    chdir $dir;
    opendir D, "." or die "cannot read $dir\n";
    while (my $f = readdir D) {
        next unless ($f =~ m/fold$/);
        my $png = $f;
        $png =~ s/fold/png/;
        my $rcmd = '';
        $rcmd   .= "library(R4RNA)\n";
        $rcmd   .= "rna <- readVienna(\"$f\")\n";
        $rcmd   .= "png(\"$png\", width=$w, height=$h)\n";
        $rcmd   .= "plotHelix(rna)\n";
        if (defined $demo) {
           print "$rcmd\n$cmd\n";
        }
        else { 
           open  R, ">R.scr" or die "cannot create R.scr\n";
           print R $rcmd;
           close R;
           system($cmd);
        }
    }
    unlink "R.scr";
}
