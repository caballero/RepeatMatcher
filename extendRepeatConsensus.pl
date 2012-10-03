#!/usr/bin/perl

=head1 NAME

extendRepeatConsensus.pl

=head1 DESCRIPTION

Iterative program to extend the borders in for a RepeatModeler consensus.

=head1 USAGE

    perl extendRepeatConsensus.pl [PARAM]

    Parameter       Description                                       Default
    -i --in         Consensus Fasta
    -g --genome     Genome Fasta
    -o --out        Output fasta
    -s --size       Step size per iteration                           [10]
    -x --engine     Alignment engine                                  [blastn]
    -e --evalue     Minimal evalue for matches                        [1e-20]
    -n --numseqs    Maximal number of sequences to try extending      [500]
    -m --minseqs    Minimal number of sequences to continue extending [3]
    -l --minlen     Minimal length of sequences                       [200]
    -d --div        Divergence level (14,18,20,25)                    [14]
    -x --maxn       Maximal number of no-bases in extension           [2]
    -w --win        Extension window                                  [100]
    --no3p          Don't extend to 3'
    --no5p          Don't extend to 5'
    -t --temp       Temporary file names                              [temp]
    
    -h --help       Print this screen and exit
    -v --verbose    Verbose mode on
    --version       Print version and exit

=head1 EXAMPLES

    perl extendRepeatConsensus.pl -i repeat.fa -g genome.fa -o new_repeat.fa

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

# Default parameters
my $help     =    undef;         # Print help
my $verbose  =    undef;         # Verbose mode
my $version  =    undef;         # Version call flag
my $in       =    undef;
my $genome   =    undef;
my $out      =    undef;
my $size     =       10;
my $engine   = 'blastn';
my $evalue   = 1 / 1e20;
my $numseq   =      500;
my $minseq   =        3;
my $minlen   =      200;
my $div      =       14;
my $maxn     =        2;
my $win      =      100;
my $no3p     =    undef;
my $no5p     =    undef;
my $temp     =   'temp';

# Main variables
my $our_version = 0.1;        # Script version number
my $blast       = '/usr/local/bin/blastall';
my $formatdb    = '/usr/local/blast/bin/formatdb';
my $cross_match = '/usr/local/bin/cross_match';
my $linup       = '/home/asmit/bin/Linup';
my $matrix_dir  = '/home/asmit/Matrices';
my $new         = '';
my %genome      = ();
my %genome_len  = ();

# Calling options
GetOptions(
    'h|help'            => \$help,
    'v|verbose'         => \$verbose,
    'i|in=s'            => \$in,
    'o|out=s'           => \$out,
    'g|genome=s'        => \$genome,
    'd|divergence:i'    => \$div,
    's|size:i'          => \$size,
    'e|evalue:s'        => \$evalue,
    'l|minlen:i'        => \$minlen,
    'x|engine:s'        => \$engine,
    't|temp:s'          => \$temp,
    'n|numseq:i'        => \$numseq,
    'm|minseq:i'        => \$minseq,
    'z|maxn:i'          => \$maxn,
    'w|win:i'           => \$win,
    'no3p'              => \$no3p,
    'no5p'              => \$no5p
) or pod2usage(-verbose => 2);
printVersion()           if  (defined $version);
pod2usage(-verbose => 2) if  (defined $help);
pod2usage(-verbose => 2) if !(defined $in);
pod2usage(-verbose => 2) if !(defined $out);
pod2usage(-verbose => 2) if !(defined $genome);

checkCmd();
checkIndex($engine, $genome);
my $cm_param    = checkDiv($div);
my ($lab, $rep) = readFasta($in);
my $iter        = 1;
loadGenome($genome);

###################################
####        M A I N            ####
###################################
while (1) {
    warn "extending repeat, iter: $iter\n" if (defined $verbose);
    $new = extendRepeat($rep);
    my $len_old = length $rep;
    my $len_new = length $new;
    last if ($len_old == $len_new);
    $iter++;
}
printFasta("$lab | extended", $new, $out);

###################################
####   S U B R O U T I N E S   ####
###################################
sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

sub readFasta {
    my $file = shift @_;
    warn "reading file $file\n" if (defined $verbose);
    my ($name, $seq);
    open F, "$file" or die "cannot open $file\n";
    while (<F>) {
        chomp;
        if (m/>/) {
            $name = $_;
        }
        else {
            $seq .= $_;
        }
    }
    close F;
    return $name, $seq;
}

sub printFasta {
    my ($head, $seq, $file) = @_;
    warn "writing file $file\n" if (defined $verbose);
    open  F, ">$file" or die "cannot write $file\n";
    print F "$head\n";
    while ($seq) {
        my $s = substr($seq, 0, 70);
        print F "$s\n";
        substr($seq, 0, 70) = '';
    }
    close F;
}

sub checkIndex {
    my ($engine, $genome) = @_;
    if ($engine eq 'blastn') {
        unless (-e "$genome.nhr" and -e "$genome.nin" and -e "$genome.nsq") {
            warn "missing indexes for $genome, generating them\n" if (defined $verbose);
            system ("$formatdb -i $genome -p F -o F");
        }
    }
}

sub checkCmd {
    die "blastall is missing\n"    if !(-e $blast);
    die "formatdb is missing\n"    if !(-e $formatdb);
    die "cross_match is missing\n" if !(-e $cross_match);
    die "Linup is missing\n"       if !(-e $linup);
}   

sub loadGenome {
    my ($file) = @_;
    warn "reading file $file\n" if (defined $verbose);
    open F, "$file" or die "cannot open $file\n";
    my $name = '';
    while (<F>) {
        chomp;
        if (m/^>/) {
            s/>//;
            s/\s+.*//;
            $name = $_;
        }
        else {
            $genome{$name} .= $_;
            $genome_len{$name} += length $_;
        }
    }
    close F;
}

sub extendRepeat {
    my ($rep)       = @_;
    my @left_seqs   = ();
    my @right_seqs  = ();
    my $left        = '';
    my $right       = '';
    my $cons        = '';
    my $null        = '';
    my $new         = $rep;
    my $ext         = 'Z' x $size;
    open  F, ">$temp.fa" or die "cannot write $temp.fa\n";
    print F  ">repeat\n$rep\n";
    close F;
    if ($engine eq 'blastn') {
        system ("$blast -p blastn -e $evalue -i $temp.fa -d $genome -m 8 -o $temp.out");
        @left_seqs  = parseBlastLeft ("$temp.out") unless (defined $no5p);
        @right_seqs = parseBlastRight("$temp.out", length $rep) unless (defined $no3p);
    }
    else { die "search engine $engine isn't supported\n"; }
      
    if (($#left_seqs + 1)  >= $minseq) {
        $cons  = createConsensus("$ext$rep", @left_seqs);
        $left  = substr($cons, 0, $size);
        $null  = $left =~ tr/N/N/;
        $left  = '' if ($null >= $maxn);
    }
    if (($#right_seqs + 1) >= $minseq) {
        $cons  = createConsensus("$rep$ext", @right_seqs);
        $right = substr($cons, (length $cons) - $size, $size);
        $null  = $right =~ tr/N/N/;
        $right = '' if ($null >= $maxn);
    }
    
    warn "extensions: left=$left, right=$right\n" if (defined $verbose);
    $new   = "$left$rep$right";
    return $new;
}

sub parseBlastLeft {
    my ($file) = @_;
    warn "parsing left matches from $file\n" if (defined $verbose);
    my @seqs   = ();
    open F, "$file" or die "cannot open $file\n";
    while (<F>) {
        chomp;
        my ($qry, $hit, $div, $len, $mis, $gap, $qini, $qend, $hini, $hend, $e, $sco) = split (/\t/, $_);
        next if ($qini > $win);
        next if ($len  < $minlen);
        if ($hini < $hend) { # hit f-f
            next if ($hini < $size);
            push @seqs, substr($genome{$hit}, $hini - $qini - 1, ($hend - $hini) + $size);
        }
        else {               # hit f-r
            next if ($hend < $size);
            push @seqs, revcomp(substr($genome{$hit}, $hini + $qini - 1, ($hini - $hend) + $size));
        }
        last if (($#seqs + 1 ) >= $numseq); 
    }
    close F;
    return @seqs;
}

sub parseBlastRight {
    my ($file, $len) = @_;
    warn "parsing right matches from $file\n" if (defined $verbose);
    my @seqs   = ();
    open F, "$file" or die "cannot open $file\n";
    while (<F>) {
        chomp;
        my ($qry, $hit, $div, $len, $mis, $gap, $qini, $qend, $hini, $hend, $e, $sco) = split (/\t/, $_);
        next if ($qend < ($len - $win));
        next if ($len  < $minlen);
        if ($hini < $hend) { # hit f-f
            next if (($hend + $size) > $genome_len{$hit}); 
            push @seqs, substr($genome{$hit}, $hend + ($len - $qend) - 1, ($hend - $hini) + $size);
        }
        else {               # hit f-r
            next if (($hini + $size) > $genome_len{$hit});
            push @seqs, revcomp(substr($genome{$hit}, $hend - 1, ($hini - $hend) + $size));
        }
        last if (($#seqs + 1 ) >= $numseq); 
    }
    close F;
    return @seqs;
}

sub createConsensus {
    my $rep = shift @_;
    open  R, ">$temp.rep.fa" or die "cannot write $temp.rep.fa\n";
    print R  ">rep0\n$rep\n";
    close R;
    
    open  F, ">$temp.repseq.fa" or die "cannot write $temp.repseq.fa\n";
    my $i = 1;
    while (my $seq = shift @_) {
        print F ">rep$i\n$seq\n";
        $i++;
    }
    close F;
    
    system "$cross_match $temp.repseq.fa $temp.rep.fa $cm_param -alignments > $temp.cm_out";
    
    system "$linup $temp.cm_out $matrix_dir/linupmatrix > $temp.ali";
    
    my $con = '';
    open A, "$temp.ali" or die "cannot open file $temp.ali\n";
    while (<A>) {
        chomp;
        next unless (m/^consensus/);
        s/consensus\s+\d+\s+//;
        s/-//g;
        $con .= $_;
    }
    close A;
    
    return $con;
}

sub checkDiv {
    my ($div) = @_;
    my $par   = '';
    if    ($div == 14) { $par = "-M $matrix_dir/14p41g.matrix -gap_init -33 -gap_ext -6 -minscore 7 -minmatch 200"; }
    elsif ($div == 18) { $par = "-M $matrix_dir/18p41g.matrix -gap_init -30 -gap_ext -6 -minscore 7 -minmatch 200"; }
    elsif ($div == 20) { $par = "-M $matrix_dir/20p41g.matrix -gap_init -28 -gap_ext -6 -minscore 7 -minmatch 200"; }
    elsif ($div == 25) { $par = "-M $matrix_dir/25p41g.matrix -gap_init -25 -gap_ext -5 -minscore 7 -minmatch 200"; }
    else  { die "Wrong divergence value, use [14,18,20,25]\n"; }
    
    warn "div=$div, cm_param=$par\n" if (defined $verbose);
    return $par;
}

sub revcomp{
    my ($s) = @_;
    my $r = reverse $s;
    $r =~ tr/ACGTacgt/TGCAtgca/;
    return $r;
}
