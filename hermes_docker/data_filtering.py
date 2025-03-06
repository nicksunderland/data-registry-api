import shutil
from pathlib import Path
from pyspark.sql import SparkSession
from pyspark.sql.functions import lit
from pyspark.sql.types import StructType, StructField, StringType, IntegerType, BooleanType, FloatType, DoubleType


class VariantColumnFilter:
    def __init__(self, column_name, regex_pattern, nullable=True, value_range=None):
        self.column_name = column_name
        self.regex_pattern = regex_pattern
        self.nullable = nullable
        self.value_range = value_range

    def apply_value_range(self, spark_df):
        if self.value_range is None:
            return spark_df[self.column_name].rlike(self.regex_pattern)
        else:
            return spark_df[self.column_name].rlike(self.regex_pattern) & \
                   (spark_df[self.column_name].cast('float') >= self.value_range[0]) & \
                   (spark_df[self.column_name].cast('float') <= self.value_range[1])

    def apply_nullable(self, spark_df):
        if self.nullable:
            return (self.apply_value_range(spark_df)) | \
                   spark_df[self.column_name].isNull()
        else:
            return (self.apply_value_range(spark_df)) & \
                   ~spark_df[self.column_name].isNull()

    def split(self, spark_df):
        condition = self.apply_nullable(spark_df)
        return spark_df.filter(~condition), spark_df.filter(condition)


class FilterIntake:
    variants_schema = StructType([
        StructField('varId', StringType(), nullable=False),
        StructField('chromosome', StringType(), nullable=False),
        StructField('position', IntegerType(), nullable=False),
        StructField('reference', StringType(), nullable=False),
        StructField('alt', StringType(), nullable=False),
        StructField('multiAllelic', BooleanType(), nullable=False),
        StructField('dataset', StringType(), nullable=False),
        StructField('phenotype', StringType(), nullable=False),
        StructField('ancestry', StringType(), nullable=True),
        StructField('g1000_eaf', FloatType(), nullable=True),
        StructField('pValue', DoubleType(), nullable=False),
        StructField('beta', FloatType(), nullable=True),
        StructField('oddsRatio', DoubleType(), nullable=True),
        StructField('eaf', FloatType(), nullable=True),
        StructField('maf', FloatType(), nullable=True),
        StructField('stdErr', FloatType(), nullable=True),
        StructField('zScore', FloatType(), nullable=True),
        StructField('n', FloatType(), nullable=True),
        StructField('ncases', FloatType(), nullable=True)
    ])

    good_chromosome = "^([1-9]{1}|1[0-9]{1}|2[0-4]{1}|X|Y|XY|MT)$"  # 1-24 + X + Y + XY + MT are valid
    good_positive_integer = "^([1-9]{1}[0-9]*|0)$"  # positive or zero only
    good_float = "^-?[0-9]+.?[0-9]*[eE]?-?[0-9]*$"  # includes scientific notation and signed
    good_positive_float = "^[0-9]+.?[0-9]*[eE]?-?[0-9]*$"  # includes scientific notation
    good_base = "^[atcgATCG]+$"  # case insensitive, only ATCG, no multialleles (commas)

    filters_to_run = [
        VariantColumnFilter("chromosome", good_chromosome, nullable=False),
        VariantColumnFilter("position", good_positive_integer, nullable=False),
        VariantColumnFilter("reference", good_base, nullable=False),
        VariantColumnFilter("alt", good_base, nullable=False),
        VariantColumnFilter("pValue", good_positive_float, nullable=False, value_range=[0, 1]),
        VariantColumnFilter("oddsRatio", good_positive_float),
        VariantColumnFilter("beta", good_float),
        VariantColumnFilter("stdErr", good_positive_float),
        VariantColumnFilter("eaf", good_positive_float, value_range=[0, 1]),
        VariantColumnFilter("n", good_positive_float)
    ]

    def __init__(self, inputfile, outputfile):
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.spark = None

    def read_variants_json(self, spark):
        return spark.read.json(self.inputfile, schema=self.variants_schema).repartition(100)

    @staticmethod
    def write_variant_json(df, file):
        p = Path(file)
        extensions = "".join(p.suffixes)
        spark_out_dir = Path(str(p).removesuffix(extensions))
        if spark_out_dir.is_dir():
            shutil.rmtree(spark_out_dir)
        Path(spark_out_dir).mkdir(parents=True, exist_ok=True)

        df.write \
            .mode('overwrite') \
            .option("delimiter", "\t") \
            .option("header", "true") \
            .csv(str(spark_out_dir))

        part_files = list(Path(spark_out_dir).glob("part-*"))

        with open(file, 'w') as outfile:
            first_file = True
            for part_file in part_files:
                with open(part_file, 'r') as infile:
                    lines = infile.readlines()
                    if first_file:
                        outfile.writelines(lines)
                        first_file = False
                    else:
                        outfile.writelines(lines[1:])

        shutil.rmtree(spark_out_dir)

    def filter(self):
        spark = SparkSession.builder.appName('qc').getOrCreate()
        df = self.read_variants_json(spark)
        bad_df = None

        for filter_to_run in self.filters_to_run:
            new_bad_df, df = filter_to_run.split(df)
            new_bad_df = new_bad_df.withColumn('filter_reason', lit(filter_to_run.column_name))
            bad_df = new_bad_df if bad_df is None else bad_df.union(new_bad_df)

        base_path = str(Path(self.outputfile).parent / Path(self.outputfile).with_suffix("").stem)
        good_filepath = base_path + ".tsv"
        bad_filepath = base_path + "_bad.tsv"
        self.write_variant_json(bad_df, bad_filepath)
        self.write_variant_json(df, good_filepath)

        spark.stop()
