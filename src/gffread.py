from subprocess import run


def run_gffread(outbase, genome_assembly, annotation, kinds=[]):
    results = {"transcript": {}, "protein": {}}
    gffread_modes = {"transcript": "w", "protein": "y"}
    cmd = "gffread -{} {} -g {} {}"

    for kind in kinds:
         if "target" not in kind:
            outpath = outbase / "evidence_sequences"
         else:
            outpath = outbase / "target_sequences"
            if not outpath.exists():
                outpath.mkdir(parents=True)
                outfile = outpath / "{}.fasta".format(kind)
                cmd_run = cmd.format(gffread_modes[kind], outfile, genome_assembly, annotation)
                if outfile.isfile():
                    log_msg = "Gffread in {} mode already, done, skipping it".format(kind)
                    return {"outfile": outfile, "log_msg": log_msg, "returncode": 0,
                            "cmd": cmd}
                else:
                    results = run(cmd_run, shell=True, capture_output=True)
                    if results.returncode == 0:
                        log_msg = "Gffread in {} mode successfully done".format(kind)
                    else:
                        log_msg = "Gffread in {} mode error: {}".format(kind, results.stderr.decode())
                    return {"outfile": outfile, "log_msg": log_msg, "returncode": results.returncode,
                            "cmd": cmd}
                    



        

    