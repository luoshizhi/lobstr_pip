#!/usr/bin/env python
"""
Filter lobSTR based on locus and call level filters.
Output VCF has filters listed for filtered loci, and a filter
  flag set for calls that are to be removed.

Locus level filters are applied before call level filters
"""

import argparse
import gzip
import math
import numpy as np
import os
import sys
import vcf
import tempfile

SMALLNUM = 10e-200
KEY_LOCLOGSCORE = "LocLogScore"
KEY_LOCCOV = "LocCov"
KEY_LOCMAXREFLEN = "LocMaxRefLen"
KEY_LOCCALLRATE = "LocCallRate"
KEY_CALLCOV = "SampleCov"
KEY_CALLLOGSCORE = "SampleLogScore"
KEY_CALLDISTEND = "SampleDistEnd"

def MSG(string):
    sys.stderr.write(string.strip() + "\n")

def LogScore(score):
    return -1*math.log10(1-score+SMALLNUM)

def GetWriter(reader, filters):
    """
    Get VCF Writer with the appropriate metadata
    """
    tmpdir = tempfile.mkdtemp(prefix="lobstr.")
    tmpfile = os.path.join(tmpdir, "header.vcf")
    f = open(tmpfile, "w")
    for line in reader._header_lines: f.write(line.strip() + "\n")
    for ft in filters.keys():
        name = ft + str(filters[ft]["Value"])
        desc = filters[ft]["Description"]
        f.write("##FILTER=<ID=%s,Description=\"%s\">\n"%(name, desc))
    f.write("##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Call-level filter.\">\n")
    f.write("#" + "\t".join(reader._column_headers + reader.samples) + "\n")
    f.close()
    writer = vcf.Writer(sys.stdout, vcf.Reader(open(tmpfile, "rb")))
    return writer

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="Input unfiltered VCF", type=str, required=True)
    parser.add_argument("--loc-log-score", help="Min mean log score cutoff to include a locus", type=float, default=0.0)
    parser.add_argument("--loc-cov", help="Min mean coverage to include a locus", type=float, default=0.0)
    parser.add_argument("--loc-max-ref-length", help="Max reference length of a locus to include", type=int, default=10000)
    parser.add_argument("--loc-call-rate", help="Min call rate to include a locus", type=float, default=0.0)
    parser.add_argument("--call-dist-end", help="Max mean absolute difference in distance from read ends to include a call", type=float, default=100)
    parser.add_argument("--call-log-score", help="Min log score cutoff to include a call", type=float, default=0.0)
    parser.add_argument("--call-cov", help="Min coverage to include a call", type=int, default=0)
    parser.add_argument("--ignore-samples", help="Ignore these samples when appying filters. File with one sample/line", type=str)
    args = parser.parse_args()

    # Set up filters:
    filters = {} # Name->{description, value)
    if args.loc_log_score > 0:
        filters[KEY_LOCLOGSCORE] = {
            "Description": "Average -log10(1-Q) less than %s"%args.loc_log_score,
            "Value": args.loc_log_score
        }
    if args.loc_cov > 0:
        filters[KEY_LOCCOV] = {
            "Description": "Average DP less than %s"%args.loc_cov,
            "Value": args.loc_cov
            }
    if args.loc_max_ref_length < 10000:
        filters[KEY_LOCMAXREFLEN] = {
            "Description": "Reference length greater than %s"%args.loc_max_ref_length,
            "Value": args.loc_max_ref_length
            }
    if args.loc_call_rate > 0:
        filters[KEY_LOCCALLRATE] = {
            "Description": "Call rate less than %s"%args.loc_call_rate,
            "Value": args.loc_call_rate
            }

    # Open VCF reader
    if args.vcf == "-":
        try:
            reader = vcf.Reader(sys.stdin)
        except (IOError, ValueError, StopIteration) as e:
            MSG("Problem reading VCF file. Is this a valid VCF?")
            sys.exit(1)
    else:
        try:
            reader = vcf.Reader(open(args.vcf, "rb"))
        except (IOError, ValueError) as e:
            MSG("Problem reading VCF file. Is this a valid VCF?")
            sys.exit(1)

    # Check that lobSTR fields are present
    fmtkeys = ["DP", "Q", "GT", "DISTENDS", "GB"]
    for fk in fmtkeys:
        if fk not in reader.formats:
            MSG("Required field %s not in VCF FORMAT. Was this file generated by allelotype?"%fk)
            sys.exit(1)
    infokeys = ["END"]
    for ik in infokeys:
        if ik not in reader.infos:
            MSG("Required field %s not in VCF INFO. Was this file generated by allelotype?"%ik)
            sys.exit(1)

    # Get any samples to ignore
    ignore_samples = []
    if args.ignore_samples:
        ignore_samples = map(lambda x: x.strip(), open(args.ignore_samples, "r").readlines())
        ignore_samples = [item for item in ignore_samples if item in reader.samples]
        MSG("Ignoring %s samples"%len(ignore_samples))

    # Set up VCF writer with new filter FORMAT field
    writer = GetWriter(reader, filters)

    # Keep track of data to print at end
    allfilters = []
    allcoverages = []
    passcoverages = []
    numallcalls = 0
    numpasscalls = 0

    # Test if there are any samples left
    if len(reader.samples)-len(ignore_samples) == 0:
        MSG("Exiting, no samples remaining")
        sys.exit(0)

    # Go through each record
    for record in reader:
        # Set up new sample info. only need to do this once per record
        if "FT" not in reader.formats:
            record.add_format("FT")
        samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':'))
        for fmt in samp_fmt._fields:
            if fmt == "FT":
                entry_type = "String"
                entry_num = 1
            else:
                entry_type = reader.formats[fmt].type
                entry_num = reader.formats[fmt].num
            samp_fmt._types.append(entry_type)
            samp_fmt._nums.append(entry_num)
        new_samples = []
        # Set call level filters and get info for locus filters
        coverages = []
        scores = []
        for sample in record:
            sample_filters = []
            if sample["GT"]:
                if sample.sample in ignore_samples:
                    sample_filters.append("IGNORE")
                else:
                    coverages.append(sample["DP"])
                    scores.append(sample["Q"])
                    if sample["DP"] < args.call_cov:
                        sample_filters.append(KEY_CALLCOV + str(args.call_cov))
                    if LogScore(sample["Q"]) < args.call_log_score:
                        sample_filters.append(KEY_CALLLOGSCORE + str(args.call_log_score))
                    if abs(sample["DISTENDS"]) > args.call_dist_end:
                        sample_filters.append(KEY_CALLDISTEND + str(args.call_dist_end))
                    allcoverages.append(sample["DP"])
                    numallcalls = numallcalls + 1
                    if len(sample_filters)==0:
                        passcoverages.append(sample["DP"])
                        numpasscalls = numpasscalls + 1
            if len(sample_filters) == 0:
                if sample["GT"]:
                    sample_filters.append("PASS")
                else:
                    sample_filters.append("NOCALL")
            # Make new sample with filter field
            sampdat = []
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key != "FT":
                    sampdat.append(sample[key])
                else: sampdat.append(",".join(sample_filters))
            call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
            new_samples.append(call)
        record.samples = new_samples
        if KEY_LOCMAXREFLEN in filters:
            str_length = record.INFO["END"]-record.POS+1
            if str_length > filters[KEY_LOCMAXREFLEN]["Value"]:
                record.add_filter(KEY_LOCMAXREFLEN+str(filters[KEY_LOCMAXREFLEN]["Value"]))
        if KEY_LOCCALLRATE in filters:
            if len(ignore_samples) > 0:
                call_rate = len(coverages)*1.0/(len(reader.samples)-len(ignore_samples))
            else: call_rate = record.num_called*1.0/len(reader.samples)
            if call_rate < filters[KEY_LOCCALLRATE]["Value"]:
                record.add_filter(KEY_LOCCALLRATE+str(filters[KEY_LOCCALLRATE]["Value"]))
        if KEY_LOCCOV in filters:
            mean_cov = np.mean(coverages)
            if mean_cov < filters[KEY_LOCCOV]["Value"]:
                record.add_filter(KEY_LOCCOV+str(filters[KEY_LOCCOV]["Value"]))
        if KEY_LOCLOGSCORE in filters:
            mean_score = np.mean(map(LogScore, scores))
            if mean_score < filters[KEY_LOCLOGSCORE]["Value"]:
                record.add_filter(KEY_LOCLOGSCORE+str(filters[KEY_LOCLOGSCORE]["Value"]))
        if record.FILTER is None or len(record.FILTER) == 0:
            record.add_filter("PASS")
        allfilters.append(":".join(record.FILTER))
        writer.write_record(record)
        
    if len(allfilters) == 0:
        MSG("No loci detected")
        sys.exit(0)
    if numallcalls == 0:
        MSG("No calls detected")
        sys.exit(0)
        
    MSG("### Locus Filter counts ###")
    for f in set(allfilters):
        MSG("%s: %s (%s)"%(f, allfilters.count(f), allfilters.count(f)*1.0/len(allfilters)))
    MSG("\n### Sample counts ###")
    MSG("Number of calls: %s"%numallcalls)
    MSG("  Mean coverage: %s"%np.mean(allcoverages))
    MSG("Number of calls passing filtering: %s (%s)"%(numpasscalls, numpasscalls*1.0/numallcalls))
    MSG("  Mean coverage: %s"%np.mean(passcoverages))