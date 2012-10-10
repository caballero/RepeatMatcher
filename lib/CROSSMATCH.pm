#!/usr/bin/perl


package CROSSMATCH;
use Carp;
sub new {
    my $class = shift;
    my $object = [];
    bless $object, $class;
    $object;
}

#########################################################################
#
#   data structure
#
#   A crossmatch object consists of the following data:
#
#      list of hashes: each list element is parametrized by the following
#      keys:
#         score
#         mismatches
#         ins
#         del
#         queryName
#         queryStart
#         queryEnd
#         queryLength
#         targetId
#         targetType
#         targetStart
#         targetEnd
#         targetLength
#         direction
#         querySeq
#         targetSeq
#
#########################################################################

sub numhits {
    my $object = shift;
    $#{$object} + 1;
}

sub score {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{score} = shift : $object->[$i]{score};
}

sub mismatches {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{mismatches} = shift : $object->[$i]{mismatches};
}

sub ins {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{ins} = shift : $object->[$i]{ins};
}

sub del {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{del} = shift : $object->[$i]{del};
}

sub queryName {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{queryName} = shift : $object->[$i]{queryName};
}

sub queryStart {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{queryStart} = shift : $object->[$i]{queryStart};
}

sub queryEnd {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{queryEnd} = shift : $object->[$i]{queryEnd};
}

sub queryLength {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{queryLength} = shift : $object->[$i]{queryLength};
}

sub targetId {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{targetId} = shift : $object->[$i]{targetId};
}

sub targetName {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{targetName} = shift : $object->[$i]{targetName};
}

sub targetType {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{targetType} = shift : $object->[$i]{targetType};
}

sub targetStart {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{targetStart} = shift : $object->[$i]{targetStart};
}

sub targetEnd {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{targetEnd} = shift : $object->[$i]{targetEnd};
}

sub targetLength {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{targetLength} = shift : $object->[$i]{targetLength};
}

sub direction {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{direction} = shift : $object->[$i]{direction};
}

sub querySeq {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{querySeq} = shift : $object->[$i]{querySeq};
}

sub targetSeq {
    my $object = shift;
    my $i = shift;
    @_ ? $object->[$i]{targetSeq} = shift : $object->[$i]{targetSeq};
}

sub splice {
    my $object = shift;
    my ($offset, $length, @list) = @_;
    splice @$object, $offset, $length, @list;
}

sub Parse {
    my $object = shift;
    my $file = shift;

    open XM, $file || croak "unable to open file $file";
    while ($_ = <XM>) { last if (/\s*(\d+)(\s+\d+\.\d+){3}\s+\w+/); }
#
#   determine whether the file is "old" or "new" format
#
    my @fieldList = split;
    my $numFields = scalar @fieldList;
    my $lastField = $fieldList[13];
#    print STDERR "$numFields $lastField\n";
#    if (($numFields == 14) && ($lastField ne '*')) {
#        $format = 'new';
#    }
#    else {
#        $format = 'old';
#    }

    my ($querySeq, $targetSeq, $q);
    seek XM, 0, 0;              # return to top of file
    my $hitNum = -1;
    while (<XM>) {
        if (/\(\d+\)/) {
            $hitNum++;
	    if (/\)\s+\+\s/) {
                ($score, $mismatches, $ins, $del, $query, $qStart, $qEnd,
                 $qRem, $orient, $name, $type, $sRem, $tEnd, $tStart)
                    = split;
                $dir = 'F';
            }
            elsif (/\)\s+C\s+[^\(]/) {
		($score, $mismatches, $ins, $del, $query, $qStart, $qEnd,
		 $qRem, $orient, $name, $sRem, $tEnd, $tStart)
		    = split;
		$dir = 'R';
	    }
	    else {
		($score, $mismatches, $ins, $del, $query, $qStart, $qEnd,
		 $qRem, $name, $tStart, $tEnd, $sRem)
		    = split;
		$dir = 'F';
	    }

#
#   remove parenthesis from $sRem
#
            $qRem = $1 if ($qRem =~ /\(? (\d+) \)?/x);
            $sRem = $1 if ($sRem =~ /\(? (\d+) \)?/x);
            $queryLength = $qRem + $qEnd;
            $targetLength = $sRem + $tEnd;
            if ($name =~ /gi\|(\d+)/) {
                $targetId = $1;
            }
            else {
                $targetId = $name;
            }
            $object->score($hitNum, $score);
            $object->mismatches($hitNum, $mismatches);
            $object->ins($hitNum, $ins);
            $object->del($hitNum, $del);
            $object->queryName($hitNum, $query);
            $object->queryStart($hitNum, $qStart);
            $object->queryEnd($hitNum, $qEnd);
            $object->queryLength($hitNum, $queryLength);
            $object->targetId($hitNum, $targetId);
            $object->targetName($hitNum, $name);
            $object->targetType($hitNum, $type);
            $object->targetStart($hitNum, $tStart);
            $object->targetEnd($hitNum, $tEnd);
            $object->targetLength($hitNum, $targetLength);
            $object->direction($hitNum, $dir);
            $querySeq = '';
            $targetSeq = '';
        }
        elsif (/^C?\s*\S+\s+\d+\s+([A-Z\-]+)\s+\d+/) {
            if ($querySeq eq '') { # start of alignment
                $q = 1;         # indicates querySeq was just read
                $querySeq = $1;
                <XM>;           # read comparison line
            }
            elsif ($q) {
                $targetSeq .= $1;
                $q = 0;         # toggle off
            }
            elsif ($q == 0) {
                $querySeq .= $1;
                $q = 1;         # toggle back on
            }
        }
        elsif (/^Transitions/) {
            $object->querySeq($hitNum, $querySeq);
            $object->targetSeq($hitNum, $targetSeq);
            $querySeq = '';
            $targetSeq = '';
        }
    }
    close XM;
}
1;
