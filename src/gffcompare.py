import os
import shutil

from subprocess import run

def run_gffcompare(outbase, protein_path, transcripts_path, 
                   anotation_path, results, kinds=[]):

    cmd = "gffcompare -r {} -o {} {}"

    outpath = outbase / "gffcompare_results"
    if not outpath.exists():
            outpath.mkdir(parents=True, exist_ok=True)

    for kind in kinds:
        if kind == "proteins_evidence":
            evidence_path = protein_path
        elif kind == "transcripts_evidence":
            evidence_path = transcripts_path
        else:
            continue
        outfile = outpath/"{}".format(kind)
        out_prefix = "{}.{}.{}"
        suffixes = ["tmap", "refmap"]
        cmd_run = cmd.format(evidence_path, outfile, anotation_path)
        if outfile.is_file():
             log_msg = "Gffcompare already done, skipping it"
             results[kind] = {"outfile": outfile, "log_msg": log_msg,
                               "returncode": 0, "cmd": cmd_run}
        else:
            cmd_results = run(cmd_run, shell=True, capture_output=True)
            if cmd_results.returncode == 0:
                log_msg = "Gffcompare successfully done"
                results[kind]= {"outfile": outfile, "log_msg": log_msg,
                                "returncode": cmd_results.returncode,
                                "cmd": cmd_run}
                for suffix in suffixes:
                     ref_fpath = evidence_path.parent / out_prefix.format(kind, anotation_path.name, suffix)
                     new_fpath = outpath/ "{}.{}.{}".format(kind, anotation_path.name, suffix)
                     shutil.move(ref_fpath, new_fpath)

            else:
                log_msg = "Gffcompare error: {}".format(cmd_results.stderr.decode())
                results[kind]= {"outfile": outfile, "log_msg": log_msg, 
                                "returncode": cmd_results.returncode,
                                "cmd": cmd_run}

    return results    

    

