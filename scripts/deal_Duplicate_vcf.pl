#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin $Script);
use Data::Dumper;

my ($vcf,$out,);

GetOptions(
	"vcf:s" => \$vcf,
	"out:s"=>\$out,

);

my $usage = <<USE;
Usage:
description:build trf index for Calling STR
author: Luoshizhi, luoshizhi\@genomics.cn
version :beta.1
date: 2017/9/14
usage: perl $Script [options]
	Common options:
	-vcf*		<str>	referent fasta for trf
	-out		<str>	out [./]
	-help|?			print help information
	
	Software options:
	
e.g.:
	perl $0 -i bam.list  -outdir ./outdir
USE
die $usage unless ($vcf &&$out);


open (IN,$vcf) or die;
open (OUT,">$out") or die;
while( <IN>){
	if (/^#/) {
		print OUT $_;
		next;
	}
	my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@sam)=split;
	if ($REF eq $ALT) {
		$ALT=".";
		for (my $i=0;$i<@sam ;$i++) {
			$sam[$i]=~s/\d\/\d:/0\/0/;
		}
	}
	print OUT join"\t",$CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@sam;
	print OUT "\n";
}
