#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
if (@ARGV<1) {
	print "perl $0 txt\n";
	exit;
}

open IN,"$ARGV[0]" or die;
my @lines=<IN>;
my ($achr,$as,$ae,$ag,$key);
my @filter;my $newflag=0;
for (my $i=1;$i<@lines;$i ++) {
	if (@filter<1) {
		($achr,$as,$ae,$ag)=(split/\s+/,$lines[$i-1])[0,1,2,8];
		$key=$lines[$i-1];
	}else{
		($achr,$as,$ae,$ag)=(split/\s+/,$key)[0,1,2,8];
	}
	my ($bchr,$bs,$be,$bg)=(split/\s+/,$lines[$i])[0,1,2,8];
	if ($achr ne $bchr) {
		push @filter,$key;
		$key=$lines[$i];
		$newflag=1;
		next;
	}

	if (($as<$be && $ae>$be) ||($as<$bs && $ae>$bs) || ($as<$bs&&$ae>$be) || ($as>$bs&&$ae<$be)) {   ## overlap  
		
		if ($ag>$bg) {
			$key=$key;
		}else{
			$key=$lines[$i];
		}
	}else{                                                            ## not overlap
		if(@filter<1||$newflag==1) {
			push @filter,$lines[$i-1];
			$newflag=0;
		}
		$key=$lines[$i];
	}
	push @filter,$key;
}

#
print  @filter;