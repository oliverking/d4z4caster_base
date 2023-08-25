## conda activate snakemake
## https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html
## snakemake -p -j2  -s ~/Dropbox/bss_2020/ok_snakemake_v7.py -C proj='.'
## snakemake  -dag -s ~/Dropbox/bss_2020/ok_snakemake_v5.py
## add log files!
## https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html

import os
import glob
import re
import pandas as pd

## in docker these should all be in PATH so can likely be ommited

samtools_path = "/usr/bin"
hisat2_path = "/usr/bin"
bowtie2_path = "/usr/bin"
bismark_path = "/usr/bin"
rscript_path = "/usr/bin"
fastqc_path = "/usr/bin"
### log all version numbers?

hisat2_mode = "local" ## hardcoded below
## local or global?
## bowtie2 or hisat2?

## default outpath is Rmd dir rather than calling dir!
## put in setup chunk?: knitr::opts_knit$set(root.dir = getwd())
rmd_path = "/scripts"
rmd_file = "bss_summary_0.1.1.Rmd"
rmd_local = "script_" + rmd_file


## new indices won't be persistent --- keep in host rather than Docker?
ref_folder = "/refs/ref_v2"
ref_file = ref_folder + "/ref_v2_pad100.fa" ## used by R

bismark_temp = "bismark_temp"


## add function to clean up fasta?
## all upper-case, unix-style line breaks



# Note: glob_wildcards returns a tuple so comma is needed
# BAMFILES, = glob_wildcards("2020-09-22_BC-fusion_library/BAM/{ids}.bam")
# ls 2020-09-22_BC-fusion_library/BAM/*.bam
# ...MethylationAnalysis.bam ...-Fusion_Prim.bam  ...-Fusion_Primers_230.bam

proj_name = ""
# proj_dir = config["proj"]
bam_dir = "UBAM"

# os.path.expanduser(path)
# os.path.expandvars(path)

BAMFILES = glob.glob(bam_dir + "/*.bam") # recursive=T?
print(BAMFILES)

## allow override?: info_file = config["info_file"]
# skip / if proj_dir = ""
info_file = glob.glob("*sample_info*.csv")

print(info_file)

if len(info_file) > 0 : ## should warn if len > 1 !

    # proj_name = re.sub("sample_info.*$", "", os.path.basename(info_file[0]))
    print("sample info file" , os.path.basename(info_file[0]))
    print("project name:" , proj_name)
    # need to remove spaces?

    samp_info = pd.read_csv(info_file[0], skiprows=1)
    if "Barcode" in samp_info.columns:
        samp_info = samp_info.set_index("Barcode", drop=False)
    else :
        samp_info = samp_info.set_index(samp_info.columns[0], drop=False)
    print(samp_info)
    if "Sample Name (required)" in samp_info.columns :
        name_col = "Sample Name (required)"
    elif "Sample Name" in samp_info.columns :
        name_col = "Sample Name"
    else:
        name_col = samp_info.columns[2]
        print("'Sample Name (required)' and 'Sample Name' not found in samp_info, using", names_col)
    # need to remove spaces?

    BAMDICT = {}
    for i in BAMFILES:
        ## TODO: match arbitrary strings from samp_info?
        BAMSTEM = re.search("(IonXpress_.*)_rawlib",i).group(1)
        if BAMSTEM in samp_info.index:
            BAMDICT[i] = samp_info[name_col][BAMSTEM]
        else:
            print("no match for ", BAMSTEM, " in sample_info")
            BAMDICT[i] = BAMSTEM
        print(i, ' -> ', BAMSTEM, '->', BAMDICT[i])

else :
    print("no sample_info.csv file found in", bam_dir)
    ## continue but don't substitute names?
    quit()

# print(BAMDICT.values())

INVDICT = {value: key for key, value in BAMDICT.items()}
for id in INVDICT:
    print(id, '->', INVDICT[id])


rule all:
    input:
        expand("bismark_bam/local.{id}_bismark_hisat2.sorted.bam.bai", id=INVDICT.keys()),
        ## files not made if there are no reads in that orientation?
        # expand("CpG_OT_local.{id}_bismark_hisat2.txt.gz", id=INVDICT.keys()),
        expand(f"r_out/{proj_name}report_{{id}}.html", id=INVDICT.keys()),
        expand("fastqc_out/{id}_fastqc.html",  id=INVDICT.keys())
        # expand("{bam}.bai", bam=BAMFILES) # <-- not needed for UBAM

rule index_bam:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    shell:
        "{samtools_path}/samtools index {input}"

rule sort_bam:
    input:
        "{sample}.bam"
    output:
        "{sample}.sorted.bam"
    shell:
        "{samtools_path}/samtools sort {input} -o {output}"

rule fastqc:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        "fastqc_out/{sample}_fastqc.html"
    shell:  ## -o for other dir, but needs to be created in advance? Or will snakemake do this?
        "{fastqc_path}/fastqc {input} -o fastqc_out"

## BAM path no longer relevant after here
rule bam2fq:
    input:  ## https://snakemake.readthedocs.io/en/stable/project_info/faq.html
        lambda wildcards: INVDICT[wildcards.sample]
    output:
        "fastq/{sample}.fastq.gz"
    shell:
        "{samtools_path}/samtools fastq -F 0x900 {input} | gzip > {output}"

## output not explictly used by name in shell command so used at dummy ix2 in bismark_bt
rule index_bt:
    input:
        ref_folder
    output:
        f"{ref_folder}/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"
    shell:
        "{bismark_path}/Bismark/bismark_genome_preparation --bowtie2 --verbose {input}"

## output not explictly used by name in shell command so used at dummy ix2 in bismark_hs
rule index_hs:
    input:
        ref_folder
    output:
        f"{ref_folder}/Bisulfite_Genome/GA_conversion/BS_GA.1.ht2"
        #expand("{rf}/Bisulfite_Genome/GA_conversion/BS_GA.1.ht2", rf=ref_folder)
    shell:
        "{bismark_path}/bismark_genome_preparation --path_to_aligner {hisat2_path} --hisat2 --verbose {input}"

### pre-create since there can be an error if two calls to bismark_hs try to autocontruct it at the same time
rule create_temp:
    output:
        directory(bismark_temp)
    shell:
        "mkdir -p {output}"

rule bismark_hs:
    input:
        fq="fastq/{sample}.fastq.gz",
        ix=ref_folder,
        ix2=f"{ref_folder}/Bisulfite_Genome/GA_conversion/BS_GA.1.ht2",
        temp=bismark_temp
    output:
        "bismark_bam/local.{sample}_bismark_hisat2.bam"
    shell:
        "{bismark_path}/bismark --output_dir bismark_bam --temp_dir {input.temp} --non_directional --hisat2 --path_to_hisat2 {hisat2_path} --samtools_path {samtools_path} --local --prefix local --genome {input.ix} -q {input.fq}"

rule bismark_me:
    input:
        "bismark_bam/local.{sample}_bismark_hisat2.bam"
    output:
        "bismark_txt/any_C_context_local.{sample}_bismark_hisat2.txt.gz"
    shell:
        "{bismark_path}/bismark_methylation_extractor --output bismark_txt --gz --samtools_path {samtools_path} --yacht {input}"

rule bismark_me2:
    input:
        "local.{sample}_bismark_hisat2.bam"
    output:
        "bismark_txt/CpG_OT_local.{sample}_bismark_hisat2.txt.gz" ## etc
    shell:
        "{bismark_path}/bismark_methylation_extractor --output bismark_txt -gz --samtools_path {samtools_path}  {input}"

# Rmarkdown sets directory of script as working directory, which is inconvenient
# when using relative paths. So here we'll copy a version locally (also for archiving)

rule copy_rmd:
    input:
        f"{rmd_path}/{rmd_file}"
    output:
        rmd_local # scripts/rmd_local ## what would need to change in Rmd file for this?
    shell:
        "cp {input} {output}"

# will script command work instead?
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html

## "gzcat {input.txt} | head > {output}" ### gives error since gzcat doesn't finish
## (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

rule cpg_matrix:
    input:
        script=rmd_local,
        txt="bismark_txt/any_C_context_local.{sample}_bismark_hisat2.txt.gz"
    output:
        html=f"r_out/{proj_name}report_{{sample}}.html",
        # pdf=f"r_out/{proj_name}report_{{sample}}.pdf",
        txt=f"r_out/{proj_name}report_{{sample}}.txt"
    params:
        stub=f"r_out/{proj_name}report_{{sample}}"
    shell:
        "{rscript_path}/Rscript -e " # Note lowercase s in Rscript on Unix!
        #"'Sys.setenv(RSTUDIO_PANDOC=\"/Applications/RStudio.app/Contents/MacOS/pandoc\"); "
        "'rmarkdown::render(\"{input.script}\", "
        "params=list(ref_file=\"{ref_file}\", in_file=\"{input.txt}\", out_txt=\"{output.txt}\"), "
        "output_format=\"html_document\", "
        "output_file=\"{params.stub}\")' "
        ## can one specify mutiple output names? not making stub.pdf and stub.html
