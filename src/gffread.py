from subprocess import run


def run_gffread(outbase, genome_assembly, annotation_path, results, kinds=[]):
    results = {}
    gffread_modes = {"transcripts": "w", "proteins": "y",
                    "transcripts_target": "w", "proteins_target": "y"}
    cmd = "gffread -{} {} -g {} {}"

    if "annotation_target" not in kinds:
        outpath = outbase / "evidence_annotation_sequences"
    else:
        outpath = outbase / "target_annotation_sequences"
    if "annotation_target" in kinds:
        kinds = ["transcripts_target", "proteins_target"]
    else:
         kinds = [kind.split("_")[0] for kind in kinds]
    if not outpath.exists():
            outpath.mkdir(parents=True)
    for kind in kinds:
            outfile = outpath / "{}.fasta".format(kind)
            cmd_run = cmd.format(gffread_modes[kind], outfile, genome_assembly, annotation_path)
            if outfile.is_file():
                log_msg = "Gffread in {} mode already, done, skipping it".format(kind)
                results[kind] = {"outfile": outfile, "log_msg": log_msg, "returncode": 0,
                                 "cmd": cmd_run}
            else:
                cmd_results = run(cmd_run, shell=True, capture_output=True)
                if cmd_results.returncode == 0:
                    log_msg = "Gffread in {} mode successfully done".format(kind)
                else:
                    log_msg = "Gffread in {} mode error: {}".format(kind, cmd_results.stderr.decode())
                results[kind]= {"outfile": outfile, "log_msg": log_msg, "returncode": cmd_results.returncode,
                                "cmd": cmd_run}
    return results

        

