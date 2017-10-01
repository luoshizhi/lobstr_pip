#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin/lib";
use AnaMethod;
use Data::Dumper;
my ($input, $build_index,$allelotype_arg,$noise_model,$filter,$outdir,$env, $move, $monitorOption, $help, $qsubMemory, $pymonitor);

GetOptions(
	"input:s" => \$input,
	"build_index:s"=>\$build_index,
	"outdir:s" => \$outdir,
	"move:s" => \$move,
	"allelotype:s" => \$allelotype_arg,
	"noise_model:s" =>\$noise_model,
	"env:s"=>\$env,
	"filter:s"=>\$filter,
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
	-input*	<str>	allsample bam file list. format:sample_name bam
	-outdir	<str>	outdir.[./]
	-allelotype	<str>	allelotype_arg default["--command classify --filter-mapq0 --filter-clipped --max-repeats-in-ends 3 --min-read-end-match 10"]
	-noise_model	<str>	default[\$Bin/share/lobSTR/models/v3.pcrfree]
	-filter	<str>	default["--loc-cov 5 --loc-log-score 0.8 --loc-call-rate 0.8 --call-cov 5 --call-log-score 0.8"]
	-env	<str> 
	-build_index	<str>	ref.fa
	-move	<str>	if this parameter is set,final result will be moved to it from output dir.
	-m	<str>	monitor options. will create monitor shell while defined this option["taskmonitor -P common -p meekat_test  -q bc.q"]
	-qsubMemory <str>	"25G,10G,5G,1G"
	-help|?	print help information

	Software options:
	-pymonitor	<str>	monitor path [\$Bin/bin/monitor]

e.g.:
	perl $0 -input bam.list  -outdir ./outdir
USE
die $usage unless ($input && $outdir);
mkpath($outdir);
$outdir = File::Spec->rel2abs($outdir);

$noise_model ||="$Bin/share/lobSTR/models/illumina_v3.pcrfree";
$allelotype_arg ||="--command classify --min-het-freq 0.2 --filter-mapq0 --filter-clipped --max-repeats-in-ends 3 --min-read-end-match 10";
$filter ||="--loc-cov 5 --loc-log-score 0.8 --loc-call-rate 0.8 --call-cov 5 --call-log-score 0.8";

$monitorOption ||="taskmonitor -P common -p lobSTR  -q bc.q";
$qsubMemory ||= "25G,5G,5G,1G";
my @qsubMemory = split /,/,$qsubMemory;
$qsubMemory[0] ||= "25G";
$qsubMemory[1] ||= "10G";
$qsubMemory[2] ||= "5G";
$qsubMemory[3] ||= "1G";
$pymonitor ||= "$Bin/bin/monitor";
my ($shell, $process, $list)=("$outdir/shell/", "$outdir/process/", "$outdir/list");
mkpath($shell);mkpath($process);mkpath($list);
$env ||=" ";

my ($bam,$dependent) = &ReadInfo2($input);
my %bam = %$bam;
my %dependent = %$dependent;
#print Dumper %bam;
#die;
my $dependence = "$list/dependence.txt";
open TXT, ">$dependence" or die $!;

###step1 build_index
foreach my $sample (keys %bam) {
	
	my $shell_t="$shell/$sample";
	my $process_t="$process/$sample";
	mkpath($shell_t);mkpath($process_t);mkpath("$process_t/tmp");
	my $content="$env";
	
	my $build_index_sh="$shell/build_index.sh";
	my $lobSTR_sh="$shell_t/lobSTR.sh";
	if (defined $dependent{$sample}) {
		if (defined $build_index) {
			
			$content .="";
			AnaMethod::generateShell($build_index_sh,$content);
			print TXT "$dependent{$sample}\t$build_index_sh:$qsubMemory[2]\n";
			print TXT "$build_index_sh:$qsubMemory[2]\t$lobSTR_sh:$qsubMemory[2]\n";
		}else{
			print TXT "$dependent{$sample}\t$lobSTR_sh:$qsubMemory[2]\n";
		}
		 
	}else{
		print TXT "$lobSTR_sh:$qsubMemory[2]\n";
	}
	
	#source env
	$content="$env\n";
	
	$content .="export PERL5LIB=$Bin/lib/PERLLIB:\$PERL5LIB\n";
	$content .="export R_LIBS=$Bin/lib/R_LIBS:\$R_LIBS\n";
	$content .="export PATH=$Bin/bin:\$PATH\n";

	$content .="###allelotype\n";
	$content .="$Bin/bin/allelotype $allelotype_arg \\\n";   #allelotype arg
	$content .="--bam $bam{$sample} \\\n";
	$content .="--noise_model $noise_model \\\n";
	$content .="--out $process_t/$sample.STR \\\n";
	$content .="--strinfo $Bin/DB/database/GRCh38.p10.info.tab \\\n";
	$content .="--index-prefix  $Bin/DB/database/GRCh38.p10.ref/lobSTR_\n";
	$content .="if [ \$? -ne 0 ];then  echo ERRO:callstr ;exit 127 ;fi\n";


	$content .="vcf-sort -c $process_t/$sample.STR.vcf -t $process_t/tmp > $process_t/$sample.STR.sort.vcf \n";
	$content .="if [ \$? -ne 0 ];then  echo ERRO:vcf-sort ;exit 127 ;fi\n";
	$content .="mv $process_t/$sample.STR.sort.vcf  $process_t/$sample.STR.vcf  \n";


	$content .="###filter_vcf\n";
	$content .="python $Bin/share/lobSTR/scripts/lobSTR_filter_vcf.py --vcf $process_t/$sample.STR.vcf $filter > $process_t/$sample.STR.mark.vcf\n";
	$content .="perl -lane \'print if /^#/; print if /PASS.+PASS/\'   $process_t/$sample.STR.mark.vcf > $process_t/$sample.STR.filter.vcf \n";
	$content .="perl -lane \'if (/^#/){print; next};next if \$F[4] eq \".\";print if \$F[3] ne \$F[4] \'$process_t/$sample.STR.filter.vcf  > $process_t/$sample.STR.filter.final.vcf\n";
	$content .="###zbgip and tabxi\n";
	$content .="bgzip $process_t/$sample.STR.vcf && tabix $process_t/$sample.STR.vcf.gz \n";
	$content .="bgzip $process_t/$sample.STR.mark.vcf && tabix $process_t/$sample.STR.mark.vcf.gz";


	
	AnaMethod::generateShell($lobSTR_sh,$content);
	if (defined $dependent{$sample}) {
		print TXT "$dependent{$sample}\t$build_index_sh:$qsubMemory[2]\n";
	}
}

close TXT;
if(defined $pymonitor && defined $monitorOption){
	`echo "$pymonitor $monitorOption -i $dependence " >$list/lobSTR_qsub.sh`;
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

sub ReadInfo2 {
        my ($file) = @_;
        my (%hashbam,%hashDepend);
        open IN, "$file" or die $!;
        while (<IN>) {
                chomp;
                next if(/^\s*$/);
                s/\s*$//;
                s/^\s*//;
                my @tmp = split;
				$hashbam{$tmp[0]}=$tmp[1];
                $hashDepend{$tmp[0]}=$tmp[2] if(@tmp >= 3);
        }
        close IN;
        return (\%hashbam,\%hashDepend);
}

sub Readpair2 {
        my ($file) = @_;
        my ($control,$treatment,%pair,%C,%T,%T_N);
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
                $T_N{$treatment}=$control;
        }
        close IN;
        return (%T_N);
}
