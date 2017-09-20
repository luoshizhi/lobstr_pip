#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin $Script);

my ($fa,$outdir,$key);

GetOptions(
	"fa:s" => \$fa,
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
	-fa*		<str>	referent fasta for trf
	-outdir		<str>	outdir [./]
	-key	<str>	keyword
	-help|?			print help information
	
	Software options:
	
e.g.:
	perl $0 -i bam.list  -outdir ./outdir
USE
die $usage unless ($fa && $key &&$outdir);

$outdir ||="./";
&mkpath($outdir);
$outdir = File::Spec->rel2abs($outdir);

### Tandem Repeats Finder 
my $trf="$Bin/bin/trf409.legacylinux64";
my $trfdir="$outdir/trf";
my $cmd="cd $trfdir & $trf 2 7 7 80 10 16 6 -f -d -m -l 6";
&CMD($cmd);
$cmd="cat $trfdir/*dat  |perl -lane \'BEGIN{\$i=\"\"};(\$i)=\$_=~/Sequence:\\s+?(\\S+)/  if /Sequence: /;next unless /\^\\d/;print \"\$i\\t\$_\"\'|perl -lane \'s/\\s+/\\t/g;print \' \>$trfdir/$key.str.bed";
&CMD($cmd);
$cmd="perl -lane \'print if \$F[2]-\$F[1]>=10\' $trfdir/$key.str.bed > $trfdir/$key.str.bed2";
&CMD($cmd);
$cmd="perl $Bin/script/filter_score_redu_overlap.pl $trfdir/$key.str.bed2 $Bin/scripts/score_threshold  \|uniq >$trfdir/$key.str.bed3 ";
&CMD($cmd);
$cmd="perl $Bin/script/rmnear.pl $trfdir/$key.str.bed3 40 \|uniq >$trfdir/$key.str.bed4";
&CMD($cmd);
$cmd="mv $trfdir/$key.str.bed4 $outdir/$key.str.bed && rm  $trfdir/$key.str.bed2  $trfdir/$key.str.bed3";
&CMD($cmd);
$cmd="python $Bin/share/lobSTR/scripts/lobstr_index.py --str $outdir/$key.str.bed  --ref $fa --out $outdir/$key";
&CMD($cmd);
$cmd="python $Bin/share/lobSTR/scripts/GetSTRInfo.py $outdir/$key.str.bed  $fa >$outdir/$key.tab";
&CMD($cmd);
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#######################################################################################

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
###########################################################
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}
#############################################################



sub CMD()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}



sub mkpath
{
	my $dir=@_;
	system "mkdir  -p $dir" if  !-d $dir;
}
