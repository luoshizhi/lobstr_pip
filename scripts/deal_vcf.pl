#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin $Script);
use Data::Dumper;

my ($vcf,$outdir,$key);

GetOptions(
	"vcf:s" => \$vcf,
	"outdir:s"=>\$outdir,
	"key:s"=> \$key,
);

my $usage = <<USE;
Usage:
description:build trf index for Calling STR
author: Luoshizhi, luoshizhi\@genomics.cn
version :beta.1
date: 2017/9/14
usage: perl $Script [options]
	Common options:
	-vcf*		<str> vcf
	-outdir		<str>	outdir [./]
	-key	<str>	keyword
	-help|?			print help information
	
	Software options:
	
e.g.:
	perl $0 -i bam.list  -outdir ./outdir
USE
die $usage unless ($vcf);

$outdir ||="./";
&mkpath($outdir);
$outdir = File::Spec->rel2abs($outdir);

my $bchr;
my $bpos;
my $bt;
open(IN ,$vcf) or die;
while (<IN>) {
	if (/^#/) {
		print $_;
		next;
	}

	my ($chr,$pos)=split;
	unless  (defined $bchr && defined$bpos) {
		$bchr=$chr ;
		$bpos =$pos;
		$bt=$_;
		next;
	}
	
	if ($bchr eq $chr && $bpos > $pos) {
		$bchr=$bchr ;
		$bpos =$bpos;
		next;
	}elsif($bchr eq $chr && $bpos < $pos){
		$bchr=$chr ;
		$bpos =$pos;
		#unless (exists $cunzai{$bchr}{$bpos}){
			print $bt;
			#$cunzai{$bchr}{$bpos}=1;
		#}
		$bt=$_;
		next;
	}elsif($bchr ne $chr){
		$bchr=$chr ;
		$bpos =$pos;
		$bt=$_;
		print $bt;
	}
}
