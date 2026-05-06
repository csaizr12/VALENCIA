import os
from pathlib import Path
import shutil

from subprocess import run

# this function runs gffcompare using protein or transcript evidence, 
# manages the output files and returns the result for each evidence
def run_gffcompare(outbase, protein_path, transcripts_path, 
                   anotation_path, results, kinds=[]):
    # define the command template and create a dedicated 'gffcompare_results' directory
    cmd = "gffcompare -r {} -o {} {}"
    outpath = outbase / "gffcompare_results"
    if not outpath.exists():
            outpath.mkdir(parents=True, exist_ok=True)
    # processing evidence
    for kind in kinds:
        # select the correct evidence file path 
        if kind == "proteins_evidence":
            evidence_path = protein_path
        elif kind == "transcripts_evidence":
            evidence_path = transcripts_path
        elif kind == "CDS_evidence":
            evidence_path = protein_path
        else:
            continue            
        # construct output file names and a list of expected suffixes
        outfile = "{}".format(kind)
        cmd_run = cmd.format(evidence_path.absolute(), outfile, Path(anotation_path.absolute()))
        # if output already exists, skip execution
        if (outpath / outfile).is_file():
             log_msg = "Gffcompare already done, skipping it"
             results["compare_"+kind] = {"outfile": outpath / outfile, "log_msg": log_msg,
                               "returncode": 0, "cmd": cmd_run}
        # otherwise, run the command via a systeam shell
        else:
            cmd_results = run(cmd_run, shell=True, capture_output=True, cwd=outpath)
            # on success
            if cmd_results.returncode == 0:
                # log the completion
                log_msg = "Gffcompare successfully done"
                results["compare_"+kind]= {"outfile": outpath / outfile, "log_msg": log_msg,
                                "returncode": cmd_results.returncode,
                                "cmd": cmd_run}
            # on failure, capture and log the specific error details
            else:
                log_msg = "Gffcompare error: {}".format(cmd_results.stderr.decode())
                results["compare_"+kind]= {"outfile": outpath / outfile, "log_msg": log_msg, 
                                "returncode": cmd_results.returncode,
                                "cmd": cmd_run}
 

    

