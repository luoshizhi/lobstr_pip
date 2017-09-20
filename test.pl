#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin);

my $trfdir="111";
print "cat $trfdir/*dat  |perl -lane \'BEGIN{\$i=\"\"};(\$i)=\$_=~/Sequence:\\s+?(\\S+)/  if /Sequence: /;next unless /\^\\d/;print \"\$i\\t\$_\"\'|perl -lane \'s/\\s+/\\t/g;print \' \>$trfdir/GRCh38.p10.str.bed\n";
