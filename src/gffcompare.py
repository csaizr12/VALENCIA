from subprocess import run
def run_gffcompare(outbase, protein_path, anotation_path, results):
    results = {}
    cmd = "gffcompare -r {} {} -o {}"
    
    outpath = outbase / "gffcompare_results"
    
    if not outpath.exists():
            outpath.mkdir(parents=True)
    outfile = outpath / "gffcompare"
    cmd_run = cmd.format(anotation_path,protein_path, outfile)
    if outfile.is_file():
        log_msg = "Gffcompare already done, skipping it"
        results["gffcompare"] = {"outfile": outfile, "log_msg": log_msg, "returncode": 0,
                                 "cmd": cmd_run}
    else:
        cmd_results = run(cmd_run, shell=True, capture_output=True)
        if cmd_results.returncode == 0:
            log_msg = "Gffcompare successfully done"
            results["gffcompare"]= {"outfile": outfile, "log_msg": log_msg, "returncode": cmd_results.returncode,
                                    "cmd": cmd_run}
        else:
            log_msg = "Gffcompare error: {}".format(cmd_results.stderr.decode())
            results["gffcompare"]= {"outfile": outfile, "log_msg": log_msg, "returncode": cmd_results.returncode,
                                "cmd": cmd_run}
    return results    

    

