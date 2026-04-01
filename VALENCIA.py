# import libraries 
import argparse
import os
import sys

from pathlib import Path

from src.gffcompare import run_gffcompare
from src.gffread import run_gffread
from src.parsers import add_tmap_info, get_gene_isoform_dict_from_target_annotation, add_tmap_info
from src.distance import edit_distance
from src.add_features_to_gff import add_features_to_gff
from src.VALENCIA_plotter import generate_quality_panel





#function generate parse_arguments:
def parse_arguments():
    desc = 'script to evaluated annotation with biological evidence'
    parser = argparse.ArgumentParser(description=desc)
    help_transcriptome = '(Required) Annotation with transcriptome evidence'
    parser.add_argument('--transcriptome_evidence','-t', 
                        type=str, help=help_transcriptome, required=True)
    help_protein = '(Required) Annotation with protein evidence'
    parser.add_argument('--protein_evidence', '-p',
                         type=str, help=help_protein, required=True)
    help_genome = '(Required) Genome assembly in fasta format'
    parser.add_argument('--genome_assembly', '-g',
                         type=str, help=help_genome, required=True)
    help_target = '(Required) Target annotation to be evaluated'
    parser.add_argument('--annotation_target', '-x',
                         type=str, help=help_target, required=True)
    help_outbase = '(Required) Outbase'
    parser.add_argument('--outbase','-o', 
                        type=str, help=help_outbase, required=True)
    return parser.parse_args()

# function to get arguments and return dictionary
def get_arguments():
    parser = parse_arguments()
    return {'transcripts_evidence': Path(parser.transcriptome_evidence).absolute(), 
            'proteins_evidence': Path(parser.protein_evidence).absolute(), 
            'genome_assembly': Path(parser.genome_assembly).absolute(), 
            'annotation_target': Path(parser.annotation_target).absolute(), 
            'outbase': Path(parser.outbase).absolute()}

# main function
def main():
    print(">>> Working directory:", os.getcwd())
    print("Command: {}".format(" ".join(sys.argv)))
    args = get_arguments()
    outbase = args["outbase"]
    results = {}
    if not outbase.exists():
        outbase.mkdir(parents=True, exist_ok=True)
    #generate_sequences for evidence
    for option, path in args.items():
       if "evidence" in option or "target" in option:
            # determine categories to run based on the option
            if option == "transcripts_evidence":
                kinds_to_run = ["transcripts_evidence"]
            elif option == "proteins_evidence":
                kinds_to_run = ["proteins_evidence", "CDS_evidence"]
            elif option == "annotation_target":
                kinds_to_run = ["annotation_target"]
            run_gffread(outbase, args["genome_assembly"],
                        path, results, kinds=kinds_to_run)
    for kind, result in results.items():
        if result["returncode"] != 0:
            print("Error in {}: {}".format(kind, result["log_msg"]))
    # run gffcompare for evidence vs target annotation and manage outputs 
    for option, path in args.items():
        if "evidence" in option:
            # determine categories to run based on the option
            if option == 'transcripts_evidence':
                 categories = ["transcripts_evidence"]
            else:
                 categories = ["proteins_evidence","CDS_evidence"]
            for category in categories:                
        #run compare(evidence, target)
                 run_gffcompare(outbase,args["proteins_evidence"], 
                                args["transcripts_evidence"],
                                args["annotation_target"], 
                                results, kinds=[category]) 
    for kind, result in results.items():
        if result["returncode"] != 0:
            print("Error in {}: {}".format(kind, result["log_msg"]))
    # dictionary with gene and isoforms from target annotation 
    with open(args["annotation_target"], "r") as target_annotation:
        gene_dict = get_gene_isoform_dict_from_target_annotation(target_annotation)
        
    # search tmap files obteined to gffcompare and edit_distance
    results_dir = Path(outbase) / 'gffcompare_results'

    for tmap_file in results_dir.glob('*.tmap'):
        # add info to gene_dict
        gene_dict = add_tmap_info(gene_dict, str(tmap_file))
        # add results gffread
        gene_dict = edit_distance(gene_dict, results["transcripts_target"]["outfile"],
                                  results["transcripts_evidence"]["outfile"],
                                  results["proteins_target"]["outfile"], 
                                  results['proteins_evidence']["outfile"], 
                                  results["CDS_target"]["outfile"],
                                  results["CDS_evidence"]["outfile"])
        print(gene_dict)
    # add features to gff
    add_features_to_gff(outbase, args["annotation_target"], gene_dict) 
    # generate quality panel
    if (outbase / "Athaliana_447_Araport11.gene_exons_with_evidence_features.gff3").exists():
        generate_quality_panel(outbase / "Athaliana_447_Araport11.gene_exons_with_evidence_features.gff3",
        output_png = outbase / "VALENCIA_Quality_Report.png")       
# run main function 
if __name__ == '__main__':
    main()