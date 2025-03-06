import gzip
import json
from Bio import SeqIO


class FaFinder:
    def __init__(self, data_path):
        self.fasta_path = f'{data_path}/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz'
        self.fa_dict = None
        with gzip.open(self.fasta_path, 'rt') as fasta_file:
            fasta_sequences = SeqIO.parse(fasta_file, 'fasta')
            self.fa_dict = {fasta.id: str(fasta.seq) for fasta in fasta_sequences}

    def get_actual_ref(self, chromosome, position, ref_length):
        if chromosome == 'XY':  # For reference treat 'XY' as 'X'
            return self.fa_dict['X'][int(position) - 1:int(position) - 1 + ref_length]
        elif chromosome in self.fa_dict:
            return self.fa_dict[chromosome][int(position) - 1:int(position) - 1 + ref_length]


class G1000Reference:
    def __init__(self, data_path, ancestry):
        self.data_path = data_path
        self.var_to_af = self.get_var_to_af(ancestry)

    def get_var_to_af(self, ancestry):
        out = {}
        with open(f'{self.data_path}/var_to_af_{ancestry}.json', 'r') as f:
            data = json.load(f)
            for variant, af in data.items():
                chromosome, position, ref, alt = variant.split(':')[:4]
                if chromosome not in out:
                    out[chromosome] = {}
                if position not in out[chromosome]:
                    out[chromosome][position] = []
                out[chromosome][position].append((ref, alt, af))
        return out

    @staticmethod
    def ref_alt_match(ref_alt_af, ref, alt):
        return ref_alt_af[0] == ref and ref_alt_af[1] == alt

    @staticmethod
    def ref_alt_prefix_match(ref_alt_af, ref, alt):
        return ref_alt_af[0][:len(ref)] == ref and \
               ref_alt_af[1][:len(alt)] == alt and \
               ref_alt_af[0][len(ref):] == ref_alt_af[1][len(alt):]

    @staticmethod
    def either_match(ref_alt_af, ref, alt):
        return G1000Reference.ref_alt_match(ref_alt_af, ref, alt) or \
               G1000Reference.ref_alt_prefix_match(ref_alt_af, ref, alt)

    def get(self, variant):
        chromosome, position, ref, alt = variant.split(':')
        ref_alt_afs = self.var_to_af.get(chromosome, {}).get(position, [])
        for ref_alt_af in ref_alt_afs:
            if G1000Reference.either_match(ref_alt_af, ref, alt):
                return ref_alt_af[2]

