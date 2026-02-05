from subprocess import run

def run_gffcompare(outbase, protein_path, transcripts_path, 
                   anotation_path, results, kinds=[]):
    cmd = "gffcompare -r {} -o {} {}"

    outpath = outbase / "gffcompare_results"
    if not outpath.exists():
            outpath.mkdir(parents=True, exist_ok=True)

    for kind in kinds:
        if kind == "proteins_evidence":
            evidence_file = protein_path
        elif kind == "transcripts_evidence":
            evidence_file = transcripts_path
        else:
            continue
        outfile = (outpath/"{}".format(kind)).resolve()
        cmd_run = cmd.format(anotation_path, outfile, evidence_file)
        print(cmd_run)
        if outfile.is_file():
             log_msg = "Gffcompare already done, skipping it"
             results[kind] = {"outfile": outfile, "log_msg": log_msg, "returncode": 0,
                                 "cmd": cmd_run}
        else:
            cmd_results = run(cmd_run, shell=True, capture_output=True)
            if cmd_results.returncode == 0:
                log_msg = "Gffcompare successfully done"
                results[kind]= {"outfile": outfile, "log_msg": log_msg,
                                 "returncode": cmd_results.returncode,
                                    "cmd": cmd_run}
            else:
                log_msg = "Gffcompare error: {}".format(cmd_results.stderr.decode())
                results[kind]= {"outfile": outfile, "log_msg": log_msg, 
                                "returncode": cmd_results.returncode,
                                "cmd": cmd_run}
    return results    

    

