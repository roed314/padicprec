#!/usr/bin/perl -w

use strict;

my @columns;
my $nlines = 0;

while(1) {
  my @lines;
  my $last = 1;
  while(<>) {
    chomp; s/\s*$//;
    last if /^\s*$/;
    push @lines, $_;
    $last = 0;
  }
  last if $last == 1;
  push @columns, \@lines;
  my $n = 1 + $#lines;
  $nlines = $n if $n > $nlines;
}

my $ncol = 1 + $#columns;
my $coldef = "|" . ("l|" x $ncol);

print "\\begin{tabular}{$coldef}\n";
for (my $i = 0; $i < $nlines; $i++) {
  for (my $j = 0; $j < $ncol; $j++) {
    my $t = $columns[$j][$i];
    if ($t =~ s/^\|\s*//) {
      print "\\hfill\\verb?" . $t . "?";
    } else {
      print "\\verb?" . $t . "?";
    }
    if ($j < $ncol-1) { print " & "; }
    else { print " \\\\\n"; }
  }
}

print "\\end{tabular}\n";
