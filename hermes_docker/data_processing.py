from genome_references import FaFinder
from genome_references import G1000Reference
import math
from scipy.stats import norm
import numpy as np
import json
import gzip


class IntakeUtilities:
    def __init__(self, refdir, phenotype, ancestry, column_map, dichotomous, dataset, cases, controls, subjects):
        self.metadata = {'dichotomous': dichotomous,
                         'ancestry': ancestry,
                         'phenotype': phenotype,
                         'column_map': column_map,
                         'dataset': dataset,
                         'cases': cases,
                         'controls': controls,
                         'subjects': subjects}
        self.fa_finder = FaFinder(refdir)
        self.g1000_reference = G1000Reference(refdir, self.metadata['ancestry'])


class LineFlipper:
    AF_CHECK_THRESHOLD = 0.3
    NEVER_COMPLIMENT_THRESHOLD = 0.1
    ALWAYS_COMPLIMENT_THRESHOLD = 1 - NEVER_COMPLIMENT_THRESHOLD

    def __init__(self, line, utils, dataset_compliment_fraction, intake_debug):
        self.line = line
        self.utils = utils
        self.dataset_compliment_fraction = dataset_compliment_fraction
        self.intake_debug = intake_debug

    def format_line(self):
        if self.line.var_id.is_valid():
            if self.line.var_id.is_unambiguous():
                flip_output = self.build_unambiguous_line()
            else:
                flip_output = self.build_ambiguous_line()
            if flip_output is not None:
                self.intake_debug.count(flip_output)
                return self.line.line_string(flip_output)
        else:
            raw_line = self.line.line_string(FlipOutput.raw())
            self.intake_debug.num_skipped += 1
            if len(self.intake_debug.skipped_lines) < 100000:
                self.intake_debug.skipped_lines.append(raw_line)
                print(f"Invalid line: {raw_line}")

    def build_unambiguous_line(self):
        # Note that with ATGC required these four cases cover all ref/alt/actual_ref combinations
        if self.line.var_id.ref_match():
            return FlipOutput(flip=False, compliment=False, null_beta=False, missing_af=False, is_ambiguous=False)
        elif self.line.var_id.ref_compliment():
            return FlipOutput(flip=False, compliment=True, null_beta=False, missing_af=False, is_ambiguous=False)
        elif self.line.var_id.alt_match():
            return FlipOutput(flip=True, compliment=False, null_beta=False, missing_af=False, is_ambiguous=False)
        elif self.line.var_id.alt_compliment():
            return FlipOutput(flip=True, compliment=True, null_beta=False, missing_af=False, is_ambiguous=False)
        else:
            raise Exception(f'Invalid ref/alt/actual_ref '
                            f'({self.line.var_id.ref}/{self.line.var_id.alt}/{self.line.var_id.actual_ref})')

    def build_ambiguous_line(self):
        if self.line.var_id.ref_match():
            return self.build_ambiguous_line_with_correct_ref()
        elif self.line.var_id.alt_match():
            return self.build_ambiguous_line_with_incorrect_ref()
        else:
            raise Exception(f'Invalid ref/alt/actual_ref '
                            f'({self.line.var_id.ref}/{self.line.var_id.alt}/{self.line.var_id.actual_ref})')

    def af_low_match(self, af):
        return af is not None and af < self.AF_CHECK_THRESHOLD and self.line.eaf < self.AF_CHECK_THRESHOLD

    def af_high_match(self, af):
        return af is not None and af > 1 - self.AF_CHECK_THRESHOLD and self.line.eaf > 1 - self.AF_CHECK_THRESHOLD

    def af_match(self, af):
        return self.af_low_match(af) or self.af_high_match(af)

    def af_low_not_match(self, af):
        return af is not None and af < self.AF_CHECK_THRESHOLD and self.line.eaf > 1 - self.AF_CHECK_THRESHOLD

    def af_high_not_match(self, af):
        return af is not None and af > 1 - self.AF_CHECK_THRESHOLD and self.line.eaf < self.AF_CHECK_THRESHOLD

    def af_not_match(self, af):
        return self.af_low_not_match(af) or self.af_high_not_match(af)

    def build_ambiguous_line_with_correct_ref(self):
        if self.dataset_compliment_fraction < self.NEVER_COMPLIMENT_THRESHOLD:
            return FlipOutput(flip=False, compliment=False, null_beta=False, missing_af=False, is_ambiguous=True)
        elif self.dataset_compliment_fraction > self.ALWAYS_COMPLIMENT_THRESHOLD:
            return FlipOutput(flip=True, compliment=True, null_beta=False, missing_af=False, is_ambiguous=True)
        else:
            af = self.utils.g1000_reference.get(self.line.var_id.format_var_id(self.line.var_id.ref, self.line.var_id.alt))
            return self.build_ambiguous_line_with_correct_ref_with_af(af)

    def build_ambiguous_line_with_correct_ref_with_af(self, af):
        if af is None or self.line.eaf is None:
            return FlipOutput(flip=False, compliment=False, null_beta=True, missing_af=True, is_ambiguous=True)
        elif self.af_match(af):
            return FlipOutput(flip=False, compliment=False, null_beta=False, missing_af=False, is_ambiguous=True)
        elif self.af_not_match(af):
            return FlipOutput(flip=True, compliment=True, null_beta=False, missing_af=False, is_ambiguous=True)
        else:
            return FlipOutput(flip=False, compliment=False, null_beta=True, missing_af=False, is_ambiguous=True)

    def build_ambiguous_line_with_incorrect_ref(self):
        if self.dataset_compliment_fraction < self.NEVER_COMPLIMENT_THRESHOLD:
            return FlipOutput(flip=True, compliment=False, null_beta=False, missing_af=False, is_ambiguous=True)
        elif self.dataset_compliment_fraction > self.ALWAYS_COMPLIMENT_THRESHOLD:
            return FlipOutput(flip=False, compliment=True, null_beta=False, missing_af=False, is_ambiguous=True)
        else:
            af = self.utils.g1000_reference.get(self.line.var_id.format_var_id(self.line.var_id.alt, self.line.var_id.ref))
            return self.build_ambiguous_line_with_incorrect_ref_with_af(af)

    def build_ambiguous_line_with_incorrect_ref_with_af(self, af):
        if af is None or self.line.eaf is None:
            return FlipOutput(flip=False, compliment=False, null_beta=True, missing_af=True, is_ambiguous=True)
        elif self.af_match(af):
            return FlipOutput(flip=False, compliment=True, null_beta=False, missing_af=False, is_ambiguous=True)
        elif self.af_not_match(af):
            return FlipOutput(flip=True, compliment=False, null_beta=False, missing_af=False, is_ambiguous=True)
        else:
            return FlipOutput(flip=False, compliment=False, null_beta=True, missing_af=False, is_ambiguous=True)


class VarId:
    compliment = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    @staticmethod
    def normalize_chromosomes(value):
        if value == "23":
            return "X"
        elif value == "24":
            return "Y"
        elif value == "25":
            return "XY"
        elif value == "26" or value == "M" or value == "m":
            return "MT"
        else:
            return value

    def __init__(self, ref, alt, formatted_line, multiallelic, col_map, utils):
        self.ref = ref
        self.alt = alt
        self.compliment_ref = self.get_compliment(self.ref)
        self.compliment_alt = self.get_compliment(self.alt)
        self.formatted_line = formatted_line
        self.multiallelic = multiallelic
        self.col_map = col_map
        self.norm_chr = VarId.normalize_chromosomes(formatted_line[self.col_map['chromosome']])
        self.position = int(formatted_line[self.col_map['position']])
        self.utils = utils
        self.actual_ref = self.get_actual_ref(self.norm_chr, self.position)

    def get_compliment(self, bases):
        if all([base in self.compliment for base in bases]):
            return ''.join([self.compliment[base] for base in bases])

    def get_actual_ref(self, chromosome, position):
        max_ref_length = max([len(self.ref), len(self.alt)])  # will be equal in the case of ambiguous line
        return self.utils.fa_finder.get_actual_ref(chromosome, position, max_ref_length)

    def is_unambiguous(self):
        return self.compliment_ref is not None and \
               self.compliment_alt is not None and \
               not self.alt == self.compliment_ref

    def is_complimentary(self):
        return self.actual_ref is not None and \
               (self.compliment_ref == self.actual_ref[:len(self.ref)] or
                self.compliment_alt == self.actual_ref[:len(self.alt)])

    def is_valid(self):
        return self.ref != self.alt and \
               self.compliment_ref is not None and \
               self.compliment_alt is not None and \
               self.actual_ref is not None

    def ref_match(self):
        return self.ref == self.actual_ref[:len(self.ref)]

    def ref_compliment(self):
        return self.compliment_ref == self.actual_ref[:len(self.ref)]

    def alt_match(self):
        return self.alt == self.actual_ref[:len(self.alt)]

    def alt_compliment(self):
        return self.compliment_alt == self.actual_ref[:len(self.alt)]

    def flip_and_or_compliment(self, flip, compliment):
        if flip:
            return (self.compliment_alt, self.compliment_ref) if compliment else (self.alt, self.ref)
        else:
            return (self.compliment_ref, self.compliment_alt) if compliment else (self.ref, self.alt)

    def format_var_id(self, ref, alt):
        return '{}:{}:{}:{}'.format(self.norm_chr, self.position, ref, alt)


class Line:
    null_link = ['', 'null', 'na', 'NaN', 'NA']

    def __init__(self, var_id, formatted_line, col_map, utils):
        self.var_id = var_id
        self.formatted_line = formatted_line
        self.col_map = col_map
        self.utils = utils

        self.pValue = self.optional_float_column('pValue')
        self.oddsRatio = self.optional_float_column('oddsRatio')
        self.beta = self.optional_float_column('beta')
        self.maf = None  # derived from eaf
        self.eaf = self.optional_float_column('eaf')
        self.stdErr = self.optional_float_column('stdErr')
        self.zScore = None  # derived from beta and standard error
        self.n = self.optional_float_column('n')
        self.ncases = self.optional_float_column('ncases')
        self.derive_optional_columns()

    def is_valid(self):
        return self.var_id.is_valid()

    def optional_float_column(self, column):
        if column in self.col_map:
            col_str = self.formatted_line[self.col_map[column]]
            if col_str not in self.null_link:
                return float(col_str)

    def infer_eff_n(self):
        if self.utils.metadata['dichotomous']:
            return 4.0 / ((1.0 / self.utils.metadata['cases']) + (1.0 / self.utils.metadata['controls']))
        else:
            return self.utils.metadata['subjects']

    def derive_optional_columns(self):
        if self.beta is None and self.oddsRatio is not None:
            self.beta = math.log(self.oddsRatio)
        if self.stdErr is None and self.beta is not None and self.pValue is not None:
            self.stdErr = -abs(self.beta) / norm.ppf(self.pValue / 2) if self.beta != 0.0 else 1.0
        if self.zScore is None and self.beta is not None and self.stdErr is not None:
            self.zScore = self.beta / self.stdErr
        if self.maf is None and self.eaf is not None:
            self.maf = self.eaf if self.eaf <= 0.5 else 1 - self.eaf
        if self.n is None:
            self.n = self.infer_eff_n()

    def get_beta(self, flip_output):
        if not flip_output.null_beta and self.beta is not None:
            return self.beta if not flip_output.flip else -self.beta

    def get_odds_ratio(self, flip_output):
        if self.oddsRatio is not None:
            return self.oddsRatio if not flip_output.flip else 1 / self.oddsRatio

    def get_eaf(self, flip_output):
        if self.eaf is not None:
            return self.eaf if not flip_output.flip else 1 - self.eaf

    def get_zScore(self, flip_output):
        if self.zScore is not None:
            return self.zScore if not flip_output.flip else -self.zScore

    def line_string(self, flip_output):
        ref, alt = self.var_id.flip_and_or_compliment(flip_output.flip, flip_output.compliment)
        ancestry_eaf = self.utils.g1000_reference.get(self.var_id.format_var_id(ref, alt))
        return json.dumps({
            'varId': self.var_id.format_var_id(ref, alt),
            'chromosome': self.var_id.norm_chr,
            'position': self.var_id.position,
            'reference': ref,
            'alt': alt,
            'multiAllelic': self.var_id.multiallelic,
            'dataset': self.utils.metadata['dataset'],
            'phenotype': self.utils.metadata['phenotype'],
            'ancestry': self.utils.metadata['ancestry'],
            'g1000_eaf': ancestry_eaf,
            'pValue': self.pValue if self.pValue > 0.0 else np.nextafter(0, 1),
            'beta': self.get_beta(flip_output),
            'oddsRatio': self.get_odds_ratio(flip_output),
            'eaf': self.get_eaf(flip_output),
            'maf': self.maf,
            'stdErr': self.stdErr,
            'zScore': self.get_zScore(flip_output),
            'n': self.n,
            'ncases': self.ncases
        })


class FlipOutput:
    def __init__(self, flip, compliment, null_beta, missing_af, is_ambiguous):
        self.flip = flip
        self.compliment = compliment
        self.null_beta = null_beta
        self.missing_af = missing_af
        self.is_ambiguous = is_ambiguous

    @staticmethod
    def raw():
        return FlipOutput(False, False, False, False, False)


class IntakeDebug:
    def __init__(self, log_name):
        self.log_name = log_name
        self.num_total = 0
        self.num_skipped = 0
        self.skipped_lines = []
        self.num_unambiguous = 0
        self.num_unambiguous_flipped = 0
        self.num_unambiguous_compliment = 0
        self.num_ambiguous = 0
        self.num_ambiguous_flipped = 0
        self.num_ambiguous_compliment = 0
        self.num_null_beta = 0
        self.num_missing_af = 0

    def __str__(self):
        skipped_line_list = '\n\t'.join(self.skipped_lines)
        return f'Total: {self.num_total}\n' \
               f'Unambiguous: {self.num_unambiguous}\n' \
               f'\tFlipped: {self.num_unambiguous_flipped}\n' \
               f'\tCompliment: {self.num_unambiguous_compliment}\n' \
               f'Ambiguous: {self.num_ambiguous}\n' \
               f'\tFlipped: {self.num_ambiguous_flipped}\n' \
               f'\tCompliment: {self.num_ambiguous_compliment}\n' \
               f'Null Beta: {self.num_null_beta}\n' \
               f'Missing AF: {self.num_missing_af}\n' \
               f'Skipped: {self.num_skipped}\n' \
               f'\t{skipped_line_list}\n'

    def count(self, flip_output):
        self.num_total += 1
        if flip_output.is_ambiguous:
            self.num_ambiguous += 1
            if flip_output.flip:
                self.num_ambiguous_flipped += 1
            if flip_output.compliment:
                self.num_ambiguous_compliment += 1
        if not flip_output.is_ambiguous:
            self.num_unambiguous += 1
            if flip_output.flip:
                self.num_unambiguous_flipped += 1
            if flip_output.compliment:
                self.num_unambiguous_compliment += 1
        if flip_output.null_beta:
            self.num_null_beta += 1
        if flip_output.missing_af:
            self.num_missing_af += 1

    def write_log(self):
        with open(self.log_name, 'w') as f:
            f.write(str(self))


class LineSplitter:
    def __init__(self, header, utils):
        self.utils = utils
        self.col_map = self.get_col_map(header)

    @staticmethod
    def normalize_spaces(value):
        return value.strip().replace(' ', '_')

    @staticmethod
    def split_line(line):
        return [LineSplitter.normalize_spaces(column) for column in line.split('\t')]

    def get_col_map(self, header):
        split_header = LineSplitter.split_line(header)
        header_map = {v: k for k, v in self.utils.metadata['column_map'].items() if v is not None}
        return {header_map[column]: idx for idx, column in enumerate(split_header) if column in header_map}

    def var_id_iterator(self, line):
        formatted_line = LineSplitter.split_line(line)
        ref = formatted_line[self.col_map['reference']].upper()
        alt = formatted_line[self.col_map['alt']].upper()
        multiallelic = ',' in ref or ',' in alt
        for split_ref in ref.split(','):
            for split_alt in alt.split(','):
                yield split_ref, split_alt, formatted_line, multiallelic

    def number_of_compliments(self, raw_line):
        num_unambiguous = 0
        num_unambiguous_and_complimentary = 0
        for split_ref, split_alt, formatted_line, multiallelic in self.var_id_iterator(raw_line):
            try:
                var_id = VarId(split_ref, split_alt, formatted_line, multiallelic, self.col_map, self.utils)
                if var_id.is_unambiguous():
                    num_unambiguous += 1
                    if var_id.is_complimentary():
                        num_unambiguous_and_complimentary += 1
            except Exception as err:
                print(f"Error: {err}")
                print(f"Unable to parse line (func. number_of_compliments): {raw_line.strip()}")
        return num_unambiguous, num_unambiguous_and_complimentary

    def write_lines(self, raw_line, dataset_compliment_fraction, f_out, intake_debug):
        for split_ref, split_alt, formatted_line, multiallelic in self.var_id_iterator(raw_line):
            try:
                var_id = VarId(split_ref, split_alt, formatted_line, multiallelic, self.col_map, self.utils)
                line = Line(var_id, formatted_line, self.col_map, self.utils)
                line_flipper = LineFlipper(line, self.utils, dataset_compliment_fraction, intake_debug)
                line_str = line_flipper.format_line()
                if line_str is not None:
                    f_out.write(f'{line_str}\n')
            except Exception as err:
                intake_debug.num_skipped += 1
                if len(intake_debug.skipped_lines) < 100000:
                    print(err)
                    intake_debug.skipped_lines.append(raw_line.strip())
                    print(f"Unable to parse line (func. write_lines): {raw_line.strip()}")


class DataIntake:
    def __init__(self, file, utils):
        self.file = file
        self.utils = utils
        self.dataset_compliment_fraction = None

    def generate_dataset_compliment_fraction(self):
        num_unambiguous = 0
        num_unambiguous_and_complimentary = 0

        if self.file.endswith('.gz'):
            with gzip.open(self.file, 'r') as f_in:
                line_splitter = LineSplitter(f_in.readline().decode(), self.utils)
                new_line = f_in.readline().decode()
                while len(new_line) > 0:
                    new_unambiguous, new_unambiguous_and_complimentary = line_splitter.number_of_compliments(new_line)
                    num_unambiguous += new_unambiguous
                    num_unambiguous_and_complimentary += new_unambiguous_and_complimentary
                    new_line = f_in.readline().decode()
        else:
            with open(self.file, 'r') as f_in:
                line_splitter = LineSplitter(f_in.readline(), self.utils)
                new_line = f_in.readline()
                while len(new_line) > 0:
                    new_unambiguous, new_unambiguous_and_complimentary = line_splitter.number_of_compliments(new_line)
                    num_unambiguous += new_unambiguous
                    num_unambiguous_and_complimentary += new_unambiguous_and_complimentary
                    new_line = f_in.readline()

        if num_unambiguous > 0:
            self.dataset_compliment_fraction = num_unambiguous_and_complimentary / num_unambiguous

    def process_file(self, output_name, intake_debug):
        if self.dataset_compliment_fraction is None:
            self.generate_dataset_compliment_fraction()

        if self.file.endswith('.gz'):
            with gzip.open(self.file, 'r') as f_in:
                with open(output_name, 'w') as f_out:
                    line_splitter = LineSplitter(f_in.readline().decode(), self.utils)
                    new_line = f_in.readline().decode()
                    while len(new_line) > 0:
                        line_splitter.write_lines(new_line, self.dataset_compliment_fraction, f_out, intake_debug)
                        new_line = f_in.readline().decode()
        else:
            with open(self.file, 'r') as f_in:
                with open(output_name, 'w') as f_out:
                    line_splitter = LineSplitter(f_in.readline(), self.utils)
                    new_line = f_in.readline()
                    while len(new_line) > 0:
                        line_splitter.write_lines(new_line, self.dataset_compliment_fraction, f_out, intake_debug)
                        new_line = f_in.readline()