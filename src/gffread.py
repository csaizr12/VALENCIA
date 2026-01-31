from subprocess import run


def run_gffread(outbase, genome_assembly, annotation_target, kinds=[]):
    results = {"transcripts": {}, "proteins": {}}
    gffread_modes = {"transcripts": "w", "proteins": "y"}
    cmd = "gffread -{} {} -g {} {}"

    if "target" not in kinds:
        outpath = outbase / "evidence_sequences"
    else:
        outpath = outbase / "target_sequences"
        kinds = ["transcripts", "proteins"]
    if not outpath.exists():
            outpath.mkdir(parents=True)
    for kind in kinds:
            kind = kind.split("_")[0]
            print(kind)
            outfile = outpath / "{}.fasta".format(kind)
            cmd_run = cmd.format(gffread_modes[kind], outfile, genome_assembly, annotation_target)
            print(cmd_run)
            if outfile.is_file():
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
                    



        

    