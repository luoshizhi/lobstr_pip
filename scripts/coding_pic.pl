#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin $Script);
use Data::Dumper;


my ($vcf,$outdir,$key,$anno,$ref);

GetOptions(
	"vcf:s" => \$vcf,
	"outdir:s"=>\$outdir,
	"key:s"=> \$key,
	"anno"=>\$anno,
	"ref:s"=> \$ref,
);

my $usage = <<"	USE";
Usage:
description:build trf index for Calling STR
author: Luoshizhi
version :beta.1
date: 2017/9/27
usage: perl $Script [options]
	Common options:
	-vcf*		<str>	input vcf file
	-outdir		<str>	outdir [./]
	-anno		<str> if input vcf has not anno ,force to use -anno
	-key	<str>	keyword
	-help|?			print help information
	-ref default [GRCh38]  hg38/hg19/GRCh38
	Software options:
	
e.g.:
	perl $0 -vcf vcf  -outdir ./outdir -key key
	USE
die $usage unless ($vcf && $key &&$outdir);

$outdir ||="./";
$ref ||="GRCh38";
&mkpath($outdir);
$outdir = File::Spec->rel2abs($outdir);
$vcf = File::Spec->rel2abs($vcf);

if ($anno) {
	system("java -jar $Bin/../snpEff/snpEff.jar $ref $vcf > $outdir/$key.anno.vcf");
}else{
	system("ln -fs $vcf $outdir/$key.anno.vcf");
}

#stat noncoding and coding 
open IN, "$outdir/$key.anno.vcf" or die;
my %chr;
while (<IN>) {
	chomp;
	next if/^#/;
	my $chr=(split)[0];
	if (/intron|intergenic/i) {
		$chr{$chr}{'noncoding'} ||=0;
		$chr{$chr}{'noncoding'} ++;
	}else{
		$chr{$chr}{'coding'} ||=0;
		$chr{$chr}{'coding'} ++;
	}
}
close IN;

open OUT, ">$outdir/$key.stat" or die;
print OUT "#type\tcut\tvalue\n";
foreach my $chr (sort keys %chr) {
	my $total=$chr{$chr}{'coding'}+$chr{$chr}{'noncoding'};
	my $codingpre=$chr{$chr}{'coding'}/$total;
	my $noncodingpre=$chr{$chr}{'noncoding'}/$total;
	print OUT "$chr\tcoding\t$codingpre\n";
	print OUT "$chr\tnoncoding\t$noncodingpre\n";
}
close OUT;

###plot pic

my $cmd="export R_LIBS=$Bin/../lib/R_LIBS:\$R_LIBS  && Rscript  $Bin/fillBar.r --infile $outdir/$key.stat --outfile $outdir/$key.str.dis.png --x.col 1 --group.col 2 --y.col 3 --x.lab \" \" --group.lab \"$key STR\" --y.lab \"SRT percent\" --title.lab \"$key STR distribution\" --x.angle 45";
print "$cmd\n";
system($cmd);
