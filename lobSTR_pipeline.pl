#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin/lib";
use AnaMethod;

my ($input, $build_index,$method, $outdir, $move, $monitorOption, $help, $qsubMemory, $pymonitor);

GetOptions(
	"input:s" => \$i,
	"build_index"=>\$build_index,
	"outdir:s" => \$outdir,
	"move:s" => \$move,
	"method" => \$method,
	"m:s" => \$monitorOption,
	"help|?" => \$help,
	"qsubMemory:s" => \$qsubMemory,
	"pymonitor:s" => \$pymonitor,
);

my $usage = <<USE;
Usage:
description:Calling STR by lobSTR
author: Luoshizhi, luoshizhi\@genomics.cn
version :beta.1
date: 2017-08-15
usage: perl $0 [options]
	Common options:
	-input*		<str>	allsample bam file list. format:sample_name bam
	-outdir		<str>	outdir.[./]
	-build_index	<str>	ref.fa
	-move		<str>	if this parameter is set,final result will be moved to it from output dir.
	-m		<str>	monitor options. will create monitor shell while defined this option
	-qsubMemory	<str>	"25G,10G,5G,1G"
	-help|?			print help information

	Software options:
	-pymonitor	<str>	monitor path [\$Bin/bin/monitor]
	
e.g.:
	perl $0 -i bam.list  -outdir ./outdir
USE
die $usage unless ($i && $sample_pair && $outdir);
mkpath($outdir);
$outdir = File::Spec->rel2abs($outdir);

$qsubMemory ||= "25G,5G,5G,1G";
my @qsubMemory = split /,/,$qsubMemory;
$qsubMemory[0] ||= "25G";
$qsubMemory[1] ||= "10G";
$qsubMemory[2] ||= "5G";
$qsubMemory[3] ||= "1G";
$pymonitor ||= "$Bin/bin/monitor";
my ($shell, $process, $list)=("$outdir/shell/", "$outdir/process/", "$outdir/list");
mkpath($shell);mkpath($process);mkpath($list);



my ($sample,$bam,$dependent) = &ReadInfo2($input);
my %sample = %$sample;
my %bam = %$bam;
my %dependent = %$dependent;

my $dependence = "$list/dependence.txt";
open TXT, ">$dependence" or die $!;

###step1 build_index
foreach my $sample (keys %sample) {
	my $content="$env";
	if ($build_index) {
		$content .="";
		if (defined $dependent{$sample}) {
		
		}
	}
}



my $pre_process = "$shell/pre_process_$sn.sh";
my $meerkat_c="$shell/meerkat_$sn.sh";
my $content="$env";
$content .="ln -s $normal $process/$sn/ && \\\n";
$content .="ln -s $normal.bai $process/$sn/ && \\\n";
$content .="ln -s $tumor $process/$sn/ && \\\n";
$content .="ln -s $tumor.bai $process/$sn/ &&\\\n";
$normal="$process/$sn/$normal_basename";
$tumor="$process/$sn/$tumor_basename";
$normal =~ s/.bam//g;
$tumor =~ s/.bam//g;
$content .="perl $meerkat/pre_process.pl -t 4 -s 20 -k 1500 -q 15 -b $normal.bam -I $hg19_bioDB/hg19.fasta -A $hg19_bioDB/hg19.fasta.fai -W $bwa -S $samtools && \\\n";
$content .="perl $meerkat/pre_process.pl -t 4 -s 20 -k 1500 -q 15 -b $tumor.bam -I $hg19_bioDB/hg19.fasta -A $hg19_bioDB/hg19.fasta.fai -W $bwa -S $samtools";
AnaMethod::generateShell($pre_process,$content);
print TXT "$pre_process:$qsubMemory[2]\t$meerkat_c:$qsubMemory[2]\n";

###step2 meerkat
my $mechanism="$shell/mechanism_$sn.sh";

$content="$env ";
$content .= "mv $tumor.blacklist.gz $tumor.blacklist.real.gz &&\\\n";
$content .= "ln $normal.blacklist.gz $tumor.blacklist.gz &&\\\n";

$content .= "perl $meerkat/meerkat.pl -s 20 -p 3 -o 1 -Q 10 -d 5 -t 8 -b $normal.bam -F $hg19_bioDB -W $bwa -S $samtools -B $blast &&\\\n";
$content .= "perl $meerkat/meerkat.pl -s 20 -p 3 -o 1 -Q 10 -d 5 -t 8 -b $tumor.bam -F $hg19_bioDB -W $bwa -S $samtools -B $blast ";
AnaMethod::generateShell($meerkat_c,$content);
print TXT "$meerkat_c:$qsubMemory[2]\t$mechanism:$qsubMemory[2]\n";

###step 3 mechanism
my $somatic_calling="$shell/somatic_calling_$sn.sh";
$content  ="$env ";
$content .= "perl $meerkat/mechanism.pl -b $normal.bam -R $meerkat_DB/hg19_rmsk.txt && \\\n";
$content .= "perl $meerkat/mechanism.pl -b $tumor.bam -R $meerkat_DB/hg19_rmsk.txt ";
AnaMethod::generateShell($mechanism,$content);
print TXT "$mechanism:$qsubMemory[2]\t$somatic_calling:$qsubMemory[2]\n";

###step4 somatic_calling
my $somatica="$tumor.somatica.variants";
my $somaticb="$tumor.somaticb.variants";
my $somaticc="$tumor.somaticc.variants";
my $somaticd="$tumor.somaticd.variants";
my $somatice="$tumor.somatice.variants";
my $somaticf="$tumor.somaticf.variants";
my $somaticg="$tumor.somaticg.variants";

my $germline_calling="$shell/germline_calling_$sn.sh";
$content ="$env ";
$content .= "perl $meerkat/somatic_sv.pl -i $tumor.variants -o $somatica -R $meerkat_DB/hg19_rmsk.txt -F re_try/ -l 1000 && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somatica -o $somaticb -R $meerkat_DB/hg19_rmsk.txt -n 1 -B $normal.bam -I $normal.isinfo  -D 5 -Q 10  -y 6 && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somaticb -o $somaticc -R $meerkat_DB/hg19_rmsk.txt -u 1 -B $normal.bam -Q 10 -S $samtools && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somaticc -o $somaticd -R $meerkat_DB/hg19_rmsk.txt -f 1 -B $normal.bam -Q 10 -S $samtools && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somaticd -o $somatice -R $meerkat_DB/hg19_rmsk.txt -e 1 -B $tumor.bam -I $tumor.isinfo -D 5 -Q 10 -S $samtools&& \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somatice -o $somaticf -R $meerkat_DB/hg19_rmsk.txt -z 1 -S $samtools && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somaticf -o $somaticg -R $meerkat_DB/hg19_rmsk.txt -d 40 -t 20 -S $samtools";
AnaMethod::generateShell($somatic_calling,$content);
print TXT "$somatic_calling:$qsubMemory[2]\t$germline_calling:$qsubMemory[2]\n";


###step 5 germline_calling

my $germa="$normal.germa.variants";
my $germb="$normal.germb.variants";
my $germc="$normal.germc.variants";
my $germd="$normal.germd.variants";
my $germe="$normal.germe.variants";

my $annotation="$shell/annotation_$sn.sh";
$content ="$env ";
$content .= "perl $meerkat/somatic_sv.pl -i $normal.variants -o $germa -R $meerkat_DB/hg19_rmsk.txt -E 0 && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $germa -o $germb -R $meerkat_DB/hg19_rmsk.txt -E 0 -e 1 -B $normal.bam -I $normal.isinfo  && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $germb -o $germc -R $meerkat_DB/hg19_rmsk.txt -E 0 -u 1 -B $normal.bam && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $germc -o $germd -R $meerkat_DB/hg19_rmsk.txt -E 0 -d 40 -t 40 && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $germd -o $germe -R $meerkat_DB/hg19_rmsk.txt -E 0 -z 1 -p 5 -P 10 ";
AnaMethod::generateShell($germline_calling,$content);
print TXT "$germline_calling:$qsubMemory[2]\t$annotation:$qsubMemory[3]\n";

###step 6 annotation
$content ="$env ";
$content .= "perl $meerkat/fusions.pl -i $somaticg -G  $meerkat_DB/hg19_refGene.sorted.txt";
AnaMethod::generateShell($annotation,$content);

	
}
close TXT;
close I;
if(defined $pymonitor && defined $monitorOption){
	`echo "$pymonitor $monitorOption -i $dependence " >$list/${method}_qsub.sh`;
}

sub mkpath
{
	my $dir=@_;
	system "mkdir  -p $dir" if  !-d $dir;
}

sub ABSOLUTE_DIR
#$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd `;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd `;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd `;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub ReadSampleInfo {
        my ($file) = @_;
        my (%hashSample,%hashbam,%hashDepend);
        open IN, "$file" or die $!;
        while (<IN>) {
                chomp;
                next if(/^\s*$/);
                s/\s*$//;
                s/^\s*//;
                my @tmp = split /\t+/;
                $hashSample{$tmp[0]}=$tmp[1];
				$hashbam{$tmp[0]}=$tmp[2];
                $hashDepend{$tmp[0]}=$tmp[3] if(@tmp >= 4);
        }
        close IN;
        return (\%hashSample,\%hashbam,\%hashDepend);
}

sub Readpair2 {
        my ($file) = @_;
        my ($control,$treatment,%pair,%C,%T);
        open IN, "$file" or die $!;
        while (<IN>) {
                next if(/^\s*$/);
                chomp;
                s/\s*$//;
                s/^\s*//;
                next if /^\s+#/;
                if(/(\S+)\t(\S+)/){
                        $control = $1;
                        $treatment = $2;
                }
                $T_N{$treatment}=control;
        }
        close IN;
        return (%T_N);
}