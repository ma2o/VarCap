#!/usr/bin/perl
#$-q all.q@cube[ab]*

use strict;
use warnings;

my $reads = $ARGV[0];
my $new_suffix = $ARGV[1];
open(my $handle, "<", $reads) || die "cannot open $reads: $!\n";
my $count = "1";

while( my $line = <$handle>){
  if ($count == "1" && $line =~ /^@/){
    if($line =~ /$new_suffix$/){
      # do not modify, as readnames suffix already exists
      print $line;
    } else {
      my @rname = split(' ',$line);
      print $rname[0].$new_suffix."\n";
    }
  } else {
    if ($count == "3"){
      $line = "+\n";
    }
    print $line;
  }
  if ($count == "4"){
    $count = "0";
  }
  $count++;
}
