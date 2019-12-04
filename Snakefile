
# import python libraries
import csv
import os

INPUT_DIR = "input"
DATA_DIR = "data"
OUTPUT_DIR = "output"
REPORTS_DIR = "output/reports"

INTERSECT_DIR = "/mnt/output/intersect"
INDIVIDUAL_DIR = "/mnt/output/individual"


INPUT_FILES = glob_wildcards("{INPUT_DIR}/{files}.bam").files
DATA_FILES = glob_wildcards("{DATA_DIR}/{files}.fa").files

# variant callers methods and pairs to analize
VC_METHODS = ["freebayes", "haplotype", "samtools"]
VC_PAIRS = [["freebayes", "haplotype"], ["freebayes", "samtools"], ["haplotype", "samtools"]]

rule all:
    input:
        expand("output/reports/{file}", REPORTS_DIR=REPORTS_DIR, file=INPUT_FILES)
    output:
        touch(".status")

rule report:
    input:
        count_individual = [INDIVIDUAL_DIR+"/"+method+"_"+"{inputfile}.FILTER.summary" for method in VC_METHODS],
        count_intersected = [INTERSECT_DIR+"/"+pair[0]+"_"+pair[1]+"_{inputfile}.FILTER.summary" for pair in VC_PAIRS],

    output:
         "output/reports/{inputfile}"
    run:
        def get_N_VARIANTS(file):
            with open(str(file), "r", newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for r in reader:
                    result =r["N_VARIANTS"]
            return result

        result = [[None, None, None],
                  [None, None, None],
                  [None, None, None]]

        for i, individual_file in enumerate(input.count_individual):
            result[i][0] = get_N_VARIANTS(individual_file)

        for i, intersected_file in enumerate(input.count_intersected):
            result[int(i==2)][int(i>0)+1] = get_N_VARIANTS(intersected_file)

        with open(str(output), "w+") as fout:
            for s in result:
                fout.write("\t".join(str(s) for s in s) + "\n")

rule count_intersected_variants:
    input:
        "{INTERSECT_DIR}/{file1}_{file2}_{inputfile}.vcf.gz",
    output:
        "{INTERSECT_DIR}/{file1}_{file2}_{inputfile}.FILTER.summary",
    shell: """
        vcftools --gzvcf {input} --FILTER-summary --stdout > {output}
        """



rule intersect_vcf:
    input:
        expand("{OUTPUT_DIR}/{dir}/{file}.vcf.gz.tbi", OUTPUT_DIR=OUTPUT_DIR, dir=VC_METHODS, file=INPUT_FILES),
    output:
        "{INTERSECT_DIR}/{file1}_{file2}_{inputfile}.vcf.gz",
    shell:  """
        vcf-isec -f -n +2 "{OUTPUT_DIR}/{wildcards.file1}/{wildcards.inputfile}.vcf.gz" \
                          "{OUTPUT_DIR}/{wildcards.file2}/{wildcards.inputfile}.vcf.gz" \
                          | bgzip -c > {output}
        """


rule count_individual_variants:
    input:
        "output/{method}/{file}.vcf.gz",
    output:
        "{INDIVIDUAL_DIR}/{method}_{file}.FILTER.summary",
    shell:   """
        vcftools --gzvcf {input} --FILTER-summary --stdout > {output}
        """

rule index_vcf:
    input:
        expand("{INPUT_DIR}/{file}.bam.bai", INPUT_DIR=INPUT_DIR, file=INPUT_FILES),
        expand("{DATA_DIR}/{DATA_FILES}.{ext}", DATA_DIR=DATA_DIR, DATA_FILES=DATA_FILES, ext=["fa.fai", "dict"]),
        source = "{OUTPUT_DIR}/{method}/{file}.vcf",
    output:
        "{OUTPUT_DIR}/{method}/{file}.vcf.gz",
        "{OUTPUT_DIR}/{method}/{file}.vcf.gz.tbi",
    shell:   """
        bgzip {input.source}
        tabix  -p vcf "{input.source}.gz"
        """




rule vc_haplotypecaller:
    output:
        "{OUTPUT_DIR}/haplotype/{inputfile}.vcf"

    shell: """
        java -jar $GATK -R "{DATA_DIR}/{DATA_FILES}.fa" \
                        -T HaplotypeCaller \
                        -I "{INPUT_DIR}/{wildcards.inputfile}.bam" -o {output}
        """

rule vc_freebayes:
    output:
        "{OUTPUT_DIR}/freebayes/{inputfile}.vcf"

    shell:  """
        freebayes -f "{DATA_DIR}/{DATA_FILES}.fa" \
                     "{INPUT_DIR}/{wildcards.inputfile}.bam" > {output}
        """


rule vc_samtools:
    output:
        "{OUTPUT_DIR}/samtools/{inputfile}.vcf"
    shell:  """
        samtools mpileup -uf "{DATA_DIR}/{DATA_FILES}.fa" \
                             "{INPUT_DIR}/{wildcards.inputfile}.bam" | bcftools view -vcg - > {output}
        """


rule index_bam:
    input: "{input_dir}/{inputfile}.bam"
    output: "{input_dir}/{inputfile}.bam.bai"
    shell:   """
        samtools index {input} > {output}
        """

rule reference_dict:
   input: "{data_dir}/{inputfile}.fa"
   output: "{data_dir}/{inputfile}.dict"
   shell:   """
        java -jar $PICARD CreateSequenceDictionary R={input} O={output}
        """

rule reference_index:
   input: "{data_dir}/{inputfile}.fa"
   output: "{data_dir}/{inputfile}.fa.fai"
   shell:   """
        samtools faidx {input} > {output}
        """
