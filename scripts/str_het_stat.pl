#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin $Script);
use Data::Dumper;


my ($vcf,$outdir);

GetOptions(
	"vcf:s" => \$vcf,
	"outdir:s"=>\$outdir,
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
	-help|?			print help information
	Software options:
	
e.g.:
	perl $0 -vcf vcf  -outdir ./outdir 
	USE
die $usage unless ($vcf  &&$outdir);

$outdir ||="./";
&mkpath($outdir);
$outdir = File::Spec->rel2abs($outdir);
$vcf = File::Spec->rel2abs($vcf);


#stat noncoding and coding 
open IN, "$vcf" or die;
my@sam;
my %str;my %strlen;
while (<IN>) {
	chomp;
	next if/^##/;
	my@tmp=split/\t+/;
	if (/#CHROM/){
	@sam=@tmp[9..$#tmp];
		next;
	}
	my $motif_len=length(($tmp[7]=~/MOTIF=(\w+);/)[0]);
	my @allelotype_len=map length($_),($tmp[3],(split/,/,$tmp[4]));
	my@sam_info=@tmp[9..$#tmp];
	for (my $i=0;$i<@sam;$i++) {
		my $GT=(split/:/,$sam_info[$i])[0];
		my @GT=sort {$a <=> $b} (split/\//,$GT);
		if ($GT eq "." || $GT eq "0/0") {                                   #0/0:Hom.ref
			$str{$sam[$i]}{$motif_len}{'Hom.ref'} ||=0;
			$str{$sam[$i]}{$motif_len}{'Hom.ref'} ++;
			push @{$strlen{$sam[$i]}{'Hom.ref'}},($allelotype_len[$GT[0]],$allelotype_len[$GT[1]]);                             #stat str len for boxplot
		}elsif($GT[0] eq $GT[1] ){                                        #1/1,2/2:Hom.noref
			$str{$sam[$i]}{$motif_len}{'Hom.noref'} ||=0;
			$str{$sam[$i]}{$motif_len}{'Hom.noref'} ++;
			push @{$strlen{$sam[$i]}{'Hom.noref'}},($allelotype_len[$GT[0]],$allelotype_len[$GT[1]]);
		}elsif($GT[0] ne $GT[1]){                               #Het
			if	($GT[0] ==0){                                             #0/1:Het.ref/noref
				$str{$sam[$i]}{$motif_len}{'Het.ref/noref'} ||=0;             
				$str{$sam[$i]}{$motif_len}{'Het.ref/noref'} ++;
				push @{$strlen{$sam[$i]}{'Het.ref/noref'}},($allelotype_len[$GT[0]],$allelotype_len[$GT[1]]);
			}else{
				$str{$sam[$i]}{$motif_len}{'Het.noref/noref'} ||=0;             #1/2:Het.noref/noref
				$str{$sam[$i]}{$motif_len}{'Het.noref/noref'} ++;
				push @{$strlen{$sam[$i]}{'Het.noref/noref'}},($allelotype_len[$GT[0]],$allelotype_len[$GT[1]]);         
			}
		}
	}
}
close IN;


for (my $i=0;$i<@sam;$i++) {
	open OUT, ">$outdir/$sam[$i].str.unit.stat" or die;
		print OUT "#type\tcut\tvalue\n";
		foreach my $motif_len (sort keys %{$str{$sam[$i]}}) {
			my $totalcategory=$str{$sam[$i]}{$motif_len}{'Hom.ref'}+$str{$sam[$i]}{$motif_len}{'Hom.noref'}+$str{$sam[$i]}{$motif_len}{'Het.ref/noref'}+$str{$sam[$i]}{$motif_len}{'Het.noref/noref'};
			print OUT "$motif_len\tHom.ref\t",$str{$sam[$i]}{$motif_len}{'Hom.ref'}/$totalcategory,"\n";
			print OUT "$motif_len\tHom.noref\t",$str{$sam[$i]}{$motif_len}{'Hom.noref'}/$totalcategory,"\n";
			print OUT "$motif_len\tHet.ref/noref\t",$str{$sam[$i]}{$motif_len}{'Het.ref/noref'}/$totalcategory,"\n";
			print OUT "$motif_len\tHet.noref/noref\t",$str{$sam[$i]}{$motif_len}{'Het.noref/noref'}/$totalcategory,"\n";
	}
	close OUT;
	open CAT, ">$outdir/$sam[$i].str.category.stat" or die;
	print CAT (map "Hom.ref\t".$_."\n",@{$strlen{$sam[$i]}{'Hom.ref'}});
	print CAT (map "Hom.noref\t".$_."\n", @{$strlen{$sam[$i]}{'Hom.noref'}});
	print CAT (map "Het.ref/noref\t".$_."\n",@{$strlen{$sam[$i]}{'Het.ref/noref'}});
	print CAT (map "Het.noref/noref\t".$_."\n",@{$strlen{$sam[$i]}{'Het.noref/noref'}});
	close CAT;

}


###plot pic

for (my $i=0;$i<@sam;$i++) {
	my $cmd="export R_LIBS=$Bin/../lib/R_LIBS:\$R_LIBS  && Rscript  $Bin/fillBar.r --infile $outdir/$sam[$i].str.unit.stat --outfile $outdir/$sam[$i].str.unit.fraction.png --x.col 1 --group.col 2 --y.col 3 --x.lab \" STR repreat unit length(bp)\" --group.lab \"$sam[$i] STR\" --y.lab \"fraction of STRs\" --title.lab \"$sam[$i] STR fraction distribution\" ";
	print "$cmd\n";
	system($cmd);
	$cmd="export R_LIBS=$Bin/../lib/R_LIBS:\$R_LIBS  && Rscript $Bin/oneFactorBox.r --infile $outdir/$sam[$i].str.category.stat --outfile $outdir/$sam[$i].str.category.png --value.col 2 --x.col 1 --x.lab \"STR category\" --y.lab \"Total STR Length(bp)\" --title.lab \" \" --filter 100";
	print "$cmd\n";
	system($cmd);
}