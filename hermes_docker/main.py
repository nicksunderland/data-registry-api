import os
from parse_arguments import parse_args
from data_processing import IntakeUtilities, IntakeDebug, DataIntake
from data_filtering import FilterIntake
from qc_plots import manhattan_plot, qq_plot, pz_plot, eaf_plot


def main():
    """
    usage: main.py[-h]
     --gwaspath GWASPATH
     --refdir REFDIR
     --ancestry {AFR, AMR, EAS, EUR, SAS}
     --phenotype PHENOTYPE
     --colmap COLMAP
     --dichotomous
    [--cases CASES]
    [--controls CONTROLS]
    [--subjects SUBJECTS]
     --outputdir OUTPUTDIR

     example: see example_run.sh
    """

    # parse the input args and make additions (parse_arguments.py)
    args = parse_args()

    # get reference data and store meta-data
    utils = IntakeUtilities(refdir=args.refdir,
                            phenotype=args.phenotype,
                            ancestry=args.ancestry,
                            column_map=args.column_map,
                            dichotomous=args.dichotomous,
                            dataset=args.base_filename,
                            cases=args.cases,
                            controls=args.controls,
                            subjects=args.subjects)

    # logs
    intake_debug = IntakeDebug(args.log_file)

    # read in the data and process to a temporary .json file in output dir
    data_intake = DataIntake(args.gwaspath, utils)
    data_intake.process_file(args.processed_file, intake_debug)
    intake_debug.write_log()

    # qc filtering and writing final file to output dir
    filter_intake = FilterIntake(inputfile=args.processed_file, outputfile=args.output_file)
    filter_intake.filter()

    # clean up temporary processed .json file
    os.remove(args.processed_file)

    # create plots
    manhattan_plot(args.output_file, args.outputdir, f'manhattan_{args.base_filename}')
    qq_plot(args.output_file, args.outputdir, f'qqplot_{args.base_filename}')
    pz_plot(args.output_file, args.outputdir, f'pzplot_{args.base_filename}')
    eaf_plot(args.output_file, args.outputdir, f'eafplot_{args.base_filename}')


if __name__ == '__main__':
    main()
