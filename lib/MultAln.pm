#!/usr/bin/perl

package MultAln;

#BEGIN {push @INC, '/home/asmit/bin')};

use CROSSMATCH ();
use ReadMatrix ();

$CGparam = 12;  # TG or CA match is 19, TG <-> CA mismatch -12
                #Previously set at 14. Seems to overestimate in very old elements
$TAparam = -5;  # TG or CA to TA is -4, so slightly worse than that
                # CG -> TA mismatch would have been -8

##########################################################################
#
#   A MultAln object has the following structure:
#
#   An array of hashes given by the following keys:
#      name
#      seq
#      start
#      end
#
##########################################################################

sub new {
    my $class = shift;
    my $object = {};
    bless $object, $class;
    $object;
}

sub length {
    my $object = shift;
    @_ ? $object->{length} = shift : $object->{length};
}

sub seq {
    my $object = shift;
    my $i = shift;
    @_ ? $object->{list}[$i]{seq} = shift : $object->{list}[$i]{seq};
}

sub start {
    my $object = shift;
    my $i = shift;
    @_ ? $object->{list}[$i]{start} = shift : $object->{list}[$i]{start};
}

sub end {
    my $object = shift;
    my $i = shift;
    @_ ? $object->{list}[$i]{end} = shift : $object->{list}[$i]{end};
}

sub seqStart {
    my $object = shift;
    my $i = shift;
    @_ ? $object->{list}[$i]{seqStart} = shift : $object->{list}[$i]{seqStart};
}

sub name {
    my $object = shift;
    my $i = shift;
    @_ ? $object->{list}[$i]{name} = shift : $object->{list}[$i]{name};
}

sub list {
    my $object = shift;
    return $object->{list};
}

sub numSequences {
    my $object = shift;
    return scalar(@{$object->list});
}

sub div {
    my $object = shift;
    my $n = shift;
    @_ ? $object->{list}[$n]{div} = shift : $object->{list}[$n]{div};
}

sub divergence {
    my $object = shift;
    my $consensus = shift;
    foreach $n (0..$object->numSequences - 1) {
        my $total = 0;
        my $change = 0;
        my $j = 0;
        foreach $i ($object->start($n)..$object->end($n)) {
            my $a = substr($consensus, $i, 1);
            my $b = substr($object->seq($n), $j, 1);
            $j++;
            next if ($a eq '-');
            next if ($b eq '-');
            $change++ if ($a ne $b);
            $total++;
        }
        $object->div($n, $change / $total);
    }
}

sub peek {
    my $object = shift;
    my ($seqNum, $pos) = @_;
    my $start = $object->start($seqNum);
    if ($pos < $start) { return ''; }
    elsif ($pos  > $object->end($seqNum)) { return ''; }
    else { 
        return substr($object->seq($seqNum), $pos - $start, 1);
    }
}

sub splice {
    my $object = shift;
    my $offset = shift;
    my $length = shift;
    if (@_) {
        splice @{$object->{list}}, $offset, $length, @_;
    }
    else {
        splice @{$object->{list}}, $offset, 1;
    }
}

sub align {
    my $object = shift;
    my $file = shift;

    my $xm = CROSSMATCH->new();
    $xm->Parse ( $file );
#
#   we require that all alignments have a common target sequence
#
    my $nh = $xm->numhits;
    my $tn = $xm->targetId(0);
    foreach $i (1..$nh - 1) {
	$it = ($xm->targetId($i));
	if ($it ne $tn) {
	    print STDERR "$it is not $tn\n";
	    die "crossmatch file requires a common sequence";
	}
  }

#
#   Form reference sequence
#

    my $tMin = 1000000;
    my $tMax = -1;
    foreach $i (0..$nh - 1) {
        my $ts = $xm->targetStart($i);
        $tMin = $ts if ($ts < $tMin);
        my $te = $xm->targetEnd($i);
        $tMax = $te if ($te > $tMax);
    }
    $refSeq = ' ' x ($tMax - $tMin + 1);
    foreach $i (0..$nh - 1) {
        my $seq = $xm->targetSeq($i);
        $seq =~ s/-//g;
        my $len = length $seq;
        my $ts = $xm->targetStart($i);
        substr($refSeq, $ts - $tMin, $len) = $seq;
    }

#
#   Form gap pattern arrays for both query and subject sequences:
#   The gap pattern array of sequence is the array @gap defined by:
#       $gap[i] = number of gaps between positions i and i+1;
#   Here @gapPattern[$i] is the gap pattern for $xm->targetSeq($i) and
#   @refGapPattern is the gap pattern for $refSeq
#
#   EXAMPLE:
#      gapped seq: ACGC--GCA---CGGTGC-CGT-C
#      sequence:   ACGCGCACGGTGCCGTC
#      gapPattern: 00020030000010010
#

    foreach $i (0..$nh - 1) {
        my @unGaps = split(/[^-]/, $xm->targetSeq($i));
        shift @unGaps;
        @gapPattern[$i] = [map ( length, @unGaps)];
        my $t = $xm->targetStart($i) - $tMin;
        unshift @{$gapPattern[$i]}, (0) x $t;
    }

    my $len = length($refSeq);
    foreach $j (0..$len - 1) {
        $refGapPattern[$j] = 0;
        foreach $i (0..$nh - 1) {
            $refGapPattern[$j] = $gapPattern[$i][$j] 
                if ($gapPattern[$i][$j] > $refGapPattern[$j]);
        }
    }

    my $seq = '';
    foreach $j (0..$len - 1) {
        $seq .= substr($refSeq, $j, 1);
        $seq .= '-' x $refGapPattern[$j];
    }

    $object->length(length($seq));      # $length of the alignment
    $object->seq(0, $seq);      # $seq is the first sequence of the alignment
    $object->start(0, 0);
    $object->end(0, ($object->length) - 1);
    $object->name(0, $xm->targetId(0));
    $object->seqStart(0, $tMin);

#
#   compute the start of each alignment relative to the start of the reference seq
#
#
#   $totalGaps[$j] = the total number of gaps prior to position $j 
#   in the reference sequence
#
    my @totalGaps;
    $totalGaps[0] = 0;
    foreach $j (1..$len - 1) {
        $totalGaps[$j] = $totalGaps[$j - 1] + $refGapPattern[$j - 1];
    }

    foreach $i (0..$nh - 1) {
        my $start = $xm->targetStart($i) - $tMin;
        $start += $totalGaps[$start];
        $object->start($i + 1, $start);
    }

    foreach $i (0..$nh - 1) {
        my $start = $object->start($i + 1);
        my $seq = '';
        my $len = length ( $xm->querySeq($i) );
        my $k = $xm->targetStart($i) - $tMin; # position in ungapped ref
        foreach $j (0..$len - 1) {
            my $n = substr($xm->querySeq($i), $j, 1);
            my $a = substr($xm->targetSeq($i), $j, 1);
            $seq .= $n;
            if ($a ne '-') {
                my $numgaps = $refGapPattern[$k] - $gapPattern[$i][$k];
                $seq .= '-' x $numgaps;
                $k++;
            }
        }
        $object->seq($i + 1, $seq);
        $len = length $seq;
        $object->end($i + 1, $start + $len - 1);
        $object->name($i + 1, $xm->queryName($i));
        $object->seqStart($i + 1, $xm->queryStart($i));
    }
}

sub read {
    my $object = shift;
    my $file = shift;

    open AF, $file || die "unable to open file $file";
    my $alignmentLength = 0;
    my @list = ();
    WHILELOOP: while (<AF>) {
        last if (/^consensus sequence:/);
        if ($_ !~ /^\s*$/) {
            chomp;
            my ($name, $start, $spaces, $seq, $end) = 
                /(\S+)\s*(\d+)\s{2}(\s*)(\S+)\s+(\d+)$/;
            my %x = (name => $name,
                     seqStart => $start,
                     spaces => length($spaces),
                     seq => $seq,
                     seqEnd => $end);
            my $i;
            for ($i = $#list; $i > -1; $i--) {
                my $item = $list[$i];
                if ($item->{name} eq $name && $item->{seqEnd} == $start - 1) {
                    $item->{seq} .= $seq;
                    $item->{seqEnd} = $end;
                    next WHILELOOP;
                }
            }
            $x{start} = $alignmentLength + $x{spaces};
            push @list, \%x;
        }
        else {
            foreach $item (@list) {
                $alignmentLength = length $item->{seq} 
                if (length($item->{seq}) > $alignmentLength);
            }
        }
    }
    close AF;

    $object->length($alignmentLength);
    foreach $i (0..$#list) {
        $object->name($i, $list[$i]{name});
        $object->start($i, $list[$i]{start});
        $object->seqStart($i, $list[$i]{seqStart});
        $object->seq($i, $list[$i]{seq});
        $object->end($i, $list[$i]{start} + length($list[$i]{seq}) - 1);
    }
}

sub sort {
#    my $object = shift;
    @{$_[0]->{list}} = sort byStart @{$_[0]->{list}};
}

sub Print {
    my $object = shift;
    my @list = @{$object->list};
#    @list = sort byStart @list;

    my $lineStart = 0;
    my $lineEnd = 99;
    while ($lineStart < $object->length) {
        foreach (@list) {
            next if ($_->{start} >= $lineEnd);
            next if ($_->{end} <= $lineStart);
            my $seq = '';
            if ($_->{start} > $lineStart) {
                $seq = ' ' x ($_->{start} - $lineStart);
            }
            my $seqStart = $lineStart - $_->{start};
            $seqStart = 0 if ($seqStart < 0);
            my $seqEnd = $lineEnd - $_->{start};
            $seqEnd = $_->{end} if ($seqEnd > $_->{end});
            $seq .= substr($_->{seq}, $seqStart, $seqEnd - $seqStart + 1);

            my $name = $_->{name};
            my $start;
            if ($seqStart == 0) {
                $start = $_->{seqStart};
            }
            else {
                my $priorSeq = substr($_->{seq}, 0, $seqStart);
                my $numLetters = ($priorSeq =~ tr/A-Z/A-Z/);
                $start = $_->{seqStart} + $numLetters;
            }
            $numLetters = ($seq =~ tr/A-Z/A-Z/);
            my $end = $start + $numLetters - 1;
            ($wname, $wstart, $wseq, $wend) = ($name, $start, $seq, $end);
            write;
            format STDOUT =
@<<<<<<<<<<<<<<<  @#####  @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  @#####
$wname, $wstart, $wseq, $wend
.
        }
        print "\n";
        $lineStart = $lineEnd + 1;
        $lineEnd = $lineStart + 99;
    }
}

sub byStart {
    $a->{start} <=> $b->{start};
}

sub profile {
    my $object = shift;

    my @list = @{$object->list};
    my @profile = ();
    foreach $seq (@list) {
        my $start = $seq->{start};
        my $i = 0;
        grep $profile[$start + $i++]{$_}++, split('', $seq->{seq});
    }


    return @profile;
}

sub consensus {
    my $object = shift;
    my $matrixFile = shift;
    my ($alphabet_r, $matrix_r) = ReadMatrix($matrixFile);
    foreach $n (@$alphabet_r) {
        $matrix_r->{$n, '-'} = $matrix_r->{'-', $n} = -6;
    }
    push @$alphabet_r, '-';
    $matrix_r->{'-', '-'} = 3;

    my @profile = $object->profile;
    my $consensus = '';
    my @cScore = ();
    foreach $i (0..$object->length - 1) {
        my $maxScore = -1000000;
        my $n = '';
        foreach $a (@$alphabet_r) {
            my $score = 0;
            foreach $b (keys %{$profile[$i]}) {
                $score += $profile[$i]{$b} * $matrix_r->{$a, $b};
            }
            if ($score > $maxScore) {
                $n = $a;
                $maxScore = $score;
            }
        }
        $consensus .= $n;
        push @cScore, $maxScore;
    }
#
#   go through the consensus and consider changing each dinucleotide
#   to a 'CG'
#
    FLOOP: foreach $i (0..length($consensus) - 2) {
        next if (substr($consensus, $i, 1) eq '-');
        my $CGscore = 0;
        my $dnScore = 0;
        my $dn1 = substr($consensus, $i, 1);
        my $k = $i + 1;
        while (substr($consensus, $k, 1) eq '-') {
            $k++;
            last FLOOP if ($k >= length($consensus));
        }
        my $dn2 = substr($consensus, $k, 1);
        foreach (@{$object->list}) {
            next if ($_->{start} > $i);
            next if ($_->{end} < $i);
            my $start = $_->{start};
            my $j = $i - $start;
            my $a = substr($_->{seq}, $j, 1);
            my $b = substr($_->{seq}, $k - $start, 1);
            my $c = $a . $b;
            $dnScore += ($matrix_r->{$dn1, $a}) + ($matrix_r->{$dn2, $b});
            if ($c eq 'CA' || $c eq 'TG') {
	      $CGscore += $CGparam;
	    } elsif ($c eq 'TA') {
	      $CGscore += $TAparam;
            } elsif ($c =~ /T[CT]/) {
              $CGscore += 2 + ($matrix_r->{G, $b});
	      # in other words; C->T transition scores +2
	    } elsif ($c =~ /[AG]A/) {
	      $CGscore += 2 + ($matrix_r->{C, $a});
            } else {
                $CGscore += ($matrix_r->{C, $a}) + ($matrix_r->{G, $b});
            }
        }
        if ($CGscore > $dnScore) {
            substr($consensus, $i, 1) = 'C';
            substr($consensus, $k, 1) = 'G';
        }
    }
    if (wantarray) {
        return ($consensus, \@cScore);
    }
    else {
        return $consensus;
    }
}

sub Mono {
    my $object = shift;
    $consensus = shift;
    my $nS = $object->numSequences;
    my %tF = ();
    my %tP =  ();
    my $numMN;
    my $i = 0;
    while (substr($consensus, $i, 1) eq '-') { $i++; }
    while ($i < length($consensus)) {
        my $a = substr($consensus, $i, 1);
        foreach $l (0..$nS - 1) {
            my $seq = $object->seq($l);
            my $start = $object->start($l);
            next if ($start > $i);
            my $end = $object->end($l);
            next if ($end < $j);
            my $b = substr($seq, $i - $start, 1);
            next if ($b eq '-');
            $tF{$a}{$b}++;
        }
        $i++;
        while (substr($consensus, $i, 1) eq '-') { $i++; }
    }
    foreach $a (A, C, G, T) {
        $numMN = 0;
        foreach $b (A, C, G, T) {
            $numMN += $tF{$a}{$b};
        }
        foreach $b (A, C, G, T) {
            $tP{$a}{$b} = $tF{$a}{$b} / $numMN;
        }
    }
    return (\%tF, \%tP);
}
1;
