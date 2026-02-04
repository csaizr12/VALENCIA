from subprocess import run

def run_gffcompare(outbase, protein_evidence, transcripts_evidence, 
                   anotation_target, results, kinds=[]):
    cmd = "gffcompare -r {} -G -o {} {} {}"

    outpath = outbase / "gffcompare_results"
    if not outpath.exists():
            outpath.mkdir(parents=True, exist_ok=True)
    
    for kind in kinds:
         if kind == "protein_evidence":
              evidence_file = protein_evidence
         elif kind == "transcript_evidence":
              evidence_file = transcripts_evidence
         else:
              continue
    
    for kind in kinds:
        outfile = outpath/"{}.gffcompare".format(kind)
        cmd_run = cmd.format(anotation_target, evidence_file, outfile)
        if outfile.is_file():
             log_msg = "Gffcompare already done, skipping it"
             results[kinds] = {"outfile": outfile, "log_msg": log_msg, "returncode": 0,
                                 "cmd": cmd_run}
        else:
            cmd_results = run(cmd_run, shell=True, capture_output=True)
            if cmd_results.returncode == 0:
                log_msg = "Gffcompare successfully done"
                results[kinds]= {"outfile": outfile, "log_msg": log_msg, "returncode": cmd_results.returncode,
                                    "cmd": cmd_run}
            else:
                log_msg = "Gffcompare error: {}".format(cmd_results.stderr.decode())
                results[kinds]= {"outfile": outfile, "log_msg": log_msg, "returncode": cmd_results.returncode,
                                "cmd": cmd_run}
    return results    

    

