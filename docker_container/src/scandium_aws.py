#!/usr/bin/env python3
"Completes processes to run Scandium in AWS Batch."
import sys
import argparse
import subprocess
from pathlib import Path
from time import strftime

import boto3

from aws_tools import s3utils

TMP_OUTDIR = "Scandium"
REDIRECT_OUTPUT_CMD_SUFFIX = '> %s 2>&1'

def main():
    print("\nBegun execution %s\n" % (strftime("%Y-%m-%d %H:%M:%S")), flush=True)
    args = parse_args()
    handle_output_and_run(args)
    print("\nCompleted processes at %s\n" % (strftime("%Y-%m-%d %H:%M:%S")), flush=True)


def parse_args():
    """ CLI interface to identify Scandium parameters
    RETURNS: args
    """
    parser = argparse.ArgumentParser(
        description="""Identifies Scandium input as local path or s3uri
        for download.""",
        usage="""%(prog)s [options]. Use '%(prog)s --help' or refer to the 
        Scandium documentation for more""",
    )
    parser.add_argument(
        "-i", type=s3utils.conditionally_download_alignment, required=True, metavar=""
    )
    parser.add_argument("-o", type=conditional_tmp_output, required=True, metavar="")

    other_opt = parser.add_argument_group("Other Optionals")
    other_opt.add_argument("-R", type=s3utils.conditionally_download_alignment, metavar="")
    other_opt.add_argument("-b", metavar="")
    other_opt.add_argument("-f", type=s3utils.conditionally_download_s3uri, metavar="")
    other_opt.add_argument("-g", metavar="")
    other_opt.add_argument("-k", metavar="")
    other_opt.add_argument("-m", metavar="")
    other_opt.add_argument("-n", type=s3utils.conditionally_download_s3uri, metavar="")
    other_opt.add_argument("-p", metavar="")
    other_opt.add_argument("-r", type=s3utils.conditionally_download_s3uri, metavar="")
    other_opt.add_argument("-t", type=s3utils.conditionally_download_s3uri, metavar="")
    other_opt.add_argument("-B", metavar="")
    other_opt.add_argument("-D", metavar="")
    other_opt.add_argument("-H", metavar="")
    other_opt.add_argument("-L", metavar="")
    other_opt.add_argument("-P", metavar="")
    other_opt.add_argument("-T", metavar="")
    other_opt.add_argument("-U", metavar="")

    block = parser.add_argument_group(
        """The following generate block regions used for data smoothing in
        Coverage Uniformity Analysis"""
    )
    block.add_argument("-l", metavar="")
    block.add_argument("-u", metavar="")

    flags = parser.add_argument_group("Flags")
    flags.add_argument("-a", action="store_true")
    flags.add_argument("-C", action="store_true")
    flags.add_argument("-d", action="store_true")
    flags.add_argument("-s", action="store_true")
    flags.add_argument("-w", action="store_true")
    flags.add_argument("-G", action="store_true")
    flags.add_argument("-M", action="store_true")
    flags.add_argument("-O", action="store_true")
    flags.add_argument("-V", action="store_true")
    flags.add_argument("-W", action="store_true")

    try:
        args = parser.parse_args()
    except:
        sys.exit("""\nParse args failed. Verify input, ensure valid aws config,
            and provide valid scandium parameters.""")
    return args


def conditional_tmp_output(output_dir):
    """ Runs with temporary local output directory 
    if valid s3uri. Note: '-o'parameter produces many
    outfiles under the supplied output_dir 
        -- OR -- 
    Runs with un-altered local output_dir.
    RETURNS: TUPLE with tmp output_dir and s3uri
        -- OR --
    Un-altered output_dir
    """
    if output_dir.startswith("s3://"):
        Path(TMP_OUTDIR).mkdir(parents=True, exist_ok=True)
        print("Temporary local directory '%s' created" % (TMP_OUTDIR), flush=True)            
        # make sure s3utils handles this as a directory
        if not output_dir.endswith("/"):
            output_dir = output_dir + "/"
        return TMP_OUTDIR, output_dir
    else:
        return output_dir


def handle_output_and_run(args):
    """ If output is an s3uri, runs commands with tmp output
    then transfers to S3.
    """
    if type(args.o) is tuple:
        upload_location = args.o[1]
        # set tmp output location
        args.o = args.o[0]
        print("Running with temporary output to: %s" % (args.o), flush=True)
        run_tools(args)
        # upload output to s3
        out_bucket, out_key = s3utils.get_bucket_keyname(upload_location)
        s3utils.s3_upload_dir_recursively(out_bucket, out_key, TMP_OUTDIR)
    else:
        run_tools(args)


def run_tools(args):
    """ Submits subprocess commands for Scandium. 
    RETURNS: None
    """
    # parse argparse Namespace for command
    d_params = {
        "-{}".format(k): str(v)  # string type to handle bool
        for (k, v) in vars(args).items()
        if v is not None and v is not False
    }
    cmd = ["scandium"]
    for param in d_params.items():
        cmd.extend(param)
    # stdout to logfile
    logfile = TMP_OUTDIR + "/scandium_log_" + strftime("%Y%m%d_%H%M%S")
    cmd.append(REDIRECT_OUTPUT_CMD_SUFFIX % (logfile))
    # run Scandium
    print("Executing command: ")
    s3utils.run_cmd(" ".join(cmd))


if __name__ == "__main__":
    main()
