from subprocess import run

# the function runs gffread to create FASTA files from a genome annotation,
# organizing outputs, skipping existing files, and returning the status
def run_gffread(outbase, genome_assembly, annotation_path, 
                results, kinds=[]):
    results = {}
    # define gffread_modes as:
    gffread_modes = {"transcripts": "w", "proteins": "y",
                    "transcripts_target": "w", 
                    "proteins_target": "y"}
    # define base command template:
    cmd = "gffread -{} {} -g {} {}"

    # determine output directory and adjust kings
    if "annotation_target" not in kinds:
        outpath = outbase / "evidence_annotation_sequences"
    else:
        outpath = outbase / "target_annotation_sequences"
    if "annotation_target" in kinds:
        kinds = ["transcripts_target", "proteins_target"]
    else:
         kinds = [kind.split("_")[0] for kind in kinds]

    # ensure output directory exists
    if not outpath.exists():
            outpath.mkdir(parents=True)
            
    # process each requested kind
    for kind in kinds:
            # build output file path
            outfile = outpath / "{}.fasta".format(kind)
            # build the command
            cmd_run = cmd.format(gffread_modes[kind], 
                                 outfile, genome_assembly,
                                  annotation_path)
            # if output already exists, skip execution
            if outfile.is_file():
                log_msg = "Gffread in {} mode already, done, skipping it".format(kind)
                results[kind] = {"outfile": outfile, "log_msg": log_msg, 
                                 "returncode": 0,"cmd": cmd_run}
            # otherwise, run the command via a systeam shell
            else:
                cmd_results = run(cmd_run, shell=True, capture_output=True)
                 # on success
                if cmd_results.returncode == 0:
                    # log the completion
                    log_msg = "Gffread in {} mode successfully done".format(kind)
                # on failure, capture and log the specific error details
                else:
                    log_msg = "Gffread in {} mode error: {}".format(kind,
                               cmd_results.stderr.decode())
                results[kind]= {"outfile": outfile, "log_msg": log_msg,
                                "returncode": cmd_results.returncode,
                                "cmd": cmd_run}
    return results

        

