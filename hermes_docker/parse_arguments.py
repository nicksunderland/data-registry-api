import argparse
from pathlib import Path


def parse_colmap(colmap_str):
    colmap_dict = {}
    key_value_pairs = colmap_str.split(',')
    for pair in key_value_pairs:
        key, value = pair.split('=')
        colmap_dict[key.strip()] = value.strip()
    return colmap_dict


def parse_args():
    opts = argparse.ArgumentParser()
    opts.add_argument('--gwaspath', type=str, required=True)
    opts.add_argument('--refdir', type=str, required=True)
    opts.add_argument('--ancestry', type=str, required=True, choices=['AFR', 'AMR', 'EAS', 'EUR', 'SAS'])
    opts.add_argument('--phenotype', type=str, required=True)
    opts.add_argument('--colmap', type=str, required=True)
    opts.add_argument('--dichotomous', action='store_true', required=True)
    opts.add_argument('--cases', type=int, required=False)
    opts.add_argument('--controls', type=int, required=False)
    opts.add_argument('--subjects', type=int, required=False)
    opts.add_argument('--outputdir', type=str, required=True)

    # parse arguments
    args = opts.parse_args()
    setattr(args, 'base_filename', str(Path(args.gwaspath).with_suffix("").stem))
    setattr(args, 'column_map', parse_colmap(args.colmap))
    output_dir = Path(args.outputdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    setattr(args, 'outputdir', str(output_dir))
    setattr(args, 'processed_file', f'{args.outputdir}/processed_{args.base_filename}.json')
    setattr(args, 'output_file', f'{args.outputdir}/output_{args.base_filename}.tsv')
    setattr(args, 'log_file', f'{args.outputdir}/log_{args.base_filename}.log')

    # checks
    if args.dichotomous and (args.cases is None or args.controls is None) and ('n' not in args.column_map or 'ncases' not in args.column_map):
        raise ValueError("For dichotomous traits, provide either 'cases' & 'controls' as arguments OR ensure columns 'n' and 'ncases' exist in the file.")
    if not args.dichotomous and 'n' not in args.column_map and args.subjects is None:
        raise ValueError("For continuous traits, provide either 'subjects' as an argument OR ensure column 'n' exists in the file.")

    return args
