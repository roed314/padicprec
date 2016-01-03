#!/usr/bin/perl -w
use strict;

my %results;

my $n; my $modulus; my $size;
my $valuation; my $jagged; my $newton;
my $lattice; my $notdiffused; my $diffused;
while(<>) {
  chomp; s/\s*$//;
  if (/\* n = (.*)/) {
    $n = $1;
  } elsif (/\* modulus = (.*)/) {
    $modulus = $1;
  } elsif (/\* Size: (.*)/) {
    $size = $1;
  } elsif (/\* Valuation: (.*)/) {
    $valuation = $1;
  } elsif (/\* Jagged: (.*)/) {
    $jagged = $1;
  } elsif (/\* Newton: (.*)/) {
    $newton = $1;
  } elsif (/\* Lattice: (.*) \((.*) \+ (.*)\)/) {
    $lattice = $1;
    $notdiffused = $2;
    $diffused = $3;
  } elsif (/^$/) {
    my $entry = [ $size, $valuation, $jagged, $newton, $lattice, $notdiffused, $diffused ];
    if (defined $results{$modulus}) {
      if (defined $results{$modulus}{$n}) {
        push @{$results{$modulus}{$n}}, $entry;
      } else {
        $results{$modulus}{$n} = [$entry];
      }
    } else {
      $results{$modulus} = { $n => [$entry] };
    }
  }
}

for my $modulus (keys(%results)) {
  my %resmod = %{$results{$modulus}};
  for my $n (keys(%resmod)) {
    my @resmodn = @{$resmod{$n}};
    $size = $valuation = $jagged = $newton = 0;
    $lattice = $notdiffused = $diffused = 0;
    for my $entry (@resmodn) {
      my ($s, $v, $j, $nw, $l, $nd, $d) = @$entry;
      $size += $s;
      $valuation += $v*$s;
      $jagged += $j*$s;
      $newton += $nw*$s;
      $lattice += $l*$s;
      $notdiffused += $nd*$s;
      $diffused += $d*$s;
    }
    print "* modulus = $modulus\n";
    print "* n = $n\n";
    print "* Size: $size\n";
    print "* Valuation: ", ($valuation/$size), "\n";
    print "* Jagged: ", ($jagged/$size), "\n";
    print "* Newton: ", ($newton/$size), "\n";
    print "* Lattice: ", ($lattice/$size), " (", ($notdiffused/$size), " + ", ($diffused/$size), ")\n";
    print "\n";
  }
}
