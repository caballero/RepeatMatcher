#!/usr/bin/perl

sub ReadMatrix {
    my $file = shift;
    my @alphabet;
    my %matrix;
    open (FILE, $file) || die "unable to open matrix file";
  LINE: while (<FILE>) {
      next if (/^ \s* \# /x);
      if (/[A-Za-z]/) {                       # line containing alphabet
          @alphabet = split;
          $i = 0;
          next LINE;
      }
      my @line = split;
      my $a = $alphabet[$i];
      for ($j = 0; $j <= $#alphabet; $j++) {
          my $b = $alphabet[$j];
          $matrix{$a, $b} = $line[$j];
      }
      $i++;
  }
    close FILE;
    return (\@alphabet, \%matrix);
}
1;
