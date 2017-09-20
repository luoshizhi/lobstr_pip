#!/bin/bash
echo hostname: `hostname`
echo ==========start at : `date` ==========
 
###allelotype
/home/luoshizhi/project/pipeline/lobstr_pip/bin/allelotype --command classify --filter-mapq0 --filter-clipped --max-repeats-in-ends 3 --min-read-end-match 10 \
--bam /zfssz3/BC_COO_P9/F17FTSECWLJ0893/HUMnccR/User/wanshq/refseq/gatk4.0/result/sample60003/sample60003.gatk4.recal.bam \
--noise_model /home/luoshizhi/project/pipeline/lobstr_pip/share/lobSTR/models/illumina_v3.pcrfree \
--out /home/luoshizhi/project/pipeline/lobstr_pip/example/process//sample60003/sample60003.STR \
--strinfo /home/luoshizhi/project/pipeline/lobstr_pip/DB/database/GRCh38.p10.info.tab \
--index-prefix  /home/luoshizhi/project/pipeline/lobstr_pip/DB/database/GRCh38.p10.ref/lobSTR_
###filter_vcf
python /home/luoshizhi/project/pipeline/lobstr_pip/share/lobSTR/scripts/lobSTR_filter_vcf.py --vcf /home/luoshizhi/project/pipeline/lobstr_pip/example/process//sample60003/sample60003.STR.vcf --loc-cov 5 --loc-log-score 0.8 --loc-call-rate 0.8 --call-cov 5 --call-log-score 0.8 > /home/luoshizhi/project/pipeline/lobstr_pip/example/process//sample60003/sample60003.STR.mark.vcf
perl -lane 'print if /^#/; print if /PASS.+PASS/'   /home/luoshizhi/project/pipeline/lobstr_pip/example/process//sample60003/sample60003.STR.mark.vcf > /home/luoshizhi/project/pipeline/lobstr_pip/example/process//sample60003/sample60003.STR.filter.vcf
 && \
echo ==========end at : `date` ========== && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /home/luoshizhi/project/pipeline/lobstr_pip/example/shell//sample60003/lobSTR.sh.sign
