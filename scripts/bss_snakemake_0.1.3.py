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

# allow overide with config["odk"] and config["mapper"]?
bismark_name = "bismark_odk"
# bismark_name = "bismark"
mapper = "bt2"
# mapper = "hisat2"

### log all version numbers?

# hisat2_mode = "local" ## hardcoded below
## local or global?
## bowtie2 or hisat2?

## default outpath is Rmd dir rather than calling dir!
## put in setup chunk?: knitr::opts_knit$set(root.dir = getwd())
rmd_path = "/scripts"
rmd_file = "bss_summary_0.1.3.Rmd"
rmd_local = "script_" + rmd_file

## new indices won't be persistent --- keep in host rather than Docker?
ref_folder = "/refs/ref_v3"
ref_file = ref_folder + "/ref_v3_pad100.fa" ## used by R

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
        expand(f"bismark_bam/local.{{id}}_bismark_{mapper}.sorted.ds.bam.bai", id=INVDICT.keys()),
        expand(f"bismark_bam/local.{{id}}_bismark_{mapper}.sorted.ds.100.bam.bai", id=INVDICT.keys()),
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
    log:
        err = "logs/index_bam/{sample}.err",
        out = "logs/index_bam/{sample}.out"
    shell:
        "{samtools_path}/samtools index {input} 2> {log.err} 1> {log.out}"

rule sort_bam:
    input:
        "{sample}.bam"
    output:
        "{sample}.sorted.bam"
    log:
        err="logs/sort_bam/{sample}.err",
        out="logs/sort_bam/{sample}.out"
    shell:
        "{samtools_path}/samtools sort {input} -o {output} 2> {log.err} 1> {log.out}"


## Flag EXACT duplicates? (Not just same start/end coords)
## What about same cigar string (PICARD)
#### not enough, as ID of mismatches may differ
## May want to do this later, since it may flatten out distintion between
## true sequence and sequencing errors.
rule dedup_bam:
    input:
        "{sample}.bam"
    output:
        "{sample}.dedup.bam"
    log:
        err="logs/dedup_bam/{sample}.err",
        out="logs/dedup_bam/{sample}.out"
    shell:
        "{samtools_path}/samtools sort {input} -o {output} 2> {log.err} 1> {log.out}"

max_per_chrom = 10000
## target max number of reads to retain per chromsome, before dedup

rule downsample_bam:
    input:
        bam="bismark_bam/{sample}.bam",
        bai="bismark_bam/{sample}.bam.bai"
    output:
        "bismark_bam/{sample}.ds.bam"
    params:
        "bismark_bam/{sample}"
    log:
        err="logs/downsample_bam/{sample}.err",
        out="logs/downsample_bam/{sample}.out"
    shell:
        "python /scripts/bss_downsample_bam.py {params} {max_per_chrom} 2> {log.err} 1> {log.out}"


## downsample again, max_per_chrom 100

rule downsample_bam_100:
    input: 
        bam="bismark_bam/{sample}.bam",
        bai="bismark_bam/{sample}.bam.bai"
    output:
        "bismark_bam/{sample}.100.bam"
    params:
        "bismark_bam/{sample}"
    log:
        err="logs/downsample_bam_again/{sample}.err",
        out="logs/downsample_bam_again/{sample}.out"
    shell:
        "python /scripts/bss_downsample_bam.py {params} 100 100 2> {log.err} 1> {log.out}"

rule fastqc:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        "fastqc_out/{sample}_fastqc.html"
    log:
        err="logs/fastqc/{sample}.err",
        out="logs/fastqc/{sample}.out"
    shell:  ## -o for other dir, but needs to be created in advance? Or will snakemake do this?
        "{fastqc_path}/fastqc {input} -o fastqc_out 2> {log.err} 1> {log.out}"

## BAM path no longer relevant after here
rule bam2fq:
    input:  ## https://snakemake.readthedocs.io/en/stable/project_info/faq.html
        lambda wildcards: INVDICT[wildcards.sample]
    output:
        "fastq/{sample}.fastq.gz"
    log:
        err="logs/bam2fq/{sample}.err",
        out="logs/bam2fq/{sample}.out"
    shell:
        "{samtools_path}/samtools fastq -F 0x900 {input} | gzip > {output} 2> {log.err}"

## output not explictly used by name in shell command so used at dummy ix2 in bismark_bt
rule index_bt:
    input:
        ref_folder
    output:
        f"{ref_folder}/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"
    log:
        err = "logs/index_bt2.err",
        out = "logs/index_bt2.out"
    shell:
        "{bismark_path}/bismark_genome_preparation --bowtie2 --verbose {input} 2> {log.err} 1> {log.out}"

## output not explictly used by name in shell command so used at dummy ix2 in bismark_hs
rule index_hs:
    input:
        ref_folder
    output:
        f"{ref_folder}/Bisulfite_Genome/GA_conversion/BS_GA.1.ht2"
        #expand("{rf}/Bisulfite_Genome/GA_conversion/BS_GA.1.ht2", rf=ref_folder)
    log:
        err = "logs/index_hisat2.err",
        out = "logs/index_hisat2.out"
    shell:
        "{bismark_path}/bismark_genome_preparation --path_to_aligner {hisat2_path} --hisat2 "
        "--verbose {input} 2> {log.err} 1> {log.out}"

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
        ix2=ancient(f"{ref_folder}/Bisulfite_Genome/GA_conversion/BS_GA.1.ht2"),
        temp=bismark_temp
    output:
        "bismark_bam/local.{sample}_bismark_hisat2.bam"
    log:
        err="logs/bismark_hisat2/{sample}.err",
        out="logs/bismark_hisat2/{sample}.out"
    shell:
        "{bismark_path}/{bismark_name} --output_dir bismark_bam --temp_dir {input.temp} --non_directional "
        "--hisat2 --path_to_hisat2 {hisat2_path} --samtools_path {samtools_path} --local --prefix local "
        "--genome {input.ix} -q {input.fq} 2> {log.err} 1> {log.out}"


rule bismark_bt:
    input:
        fq="fastq/{sample}.fastq.gz",
        ix=ref_folder,
        ix2=ancient(f"{ref_folder}/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"),
        temp=bismark_temp
    output:
        "bismark_bam/local.{sample}_bismark_bt2.bam"
    log:
        err = "logs/bismark_bt2/{sample}.err",
        out = "logs/bismark_bt2/{sample}.out"
    shell:
        "{bismark_path}/{bismark_name} --output_dir bismark_bam --temp_dir {input.temp} --non_directional "
        "--bowtie2 --path_to_bowtie2 {bowtie2_path} --samtools_path {samtools_path} --local --prefix local "
        "--genome {input.ix} -q {input.fq} 2> {log.err} 1> {log.out}"


rule bismark_me:
    input:
        "bismark_bam/{stem}.bam" ## run on just sorted.ds.bams based on rule all
    output:
        "bismark_txt/any_C_context_{stem}.txt.gz" ## etc
    log:
        err = "logs/bismark_me2/{stem}.err",
        out = "logs/bismark_me2/{stem}.out"
    shell:
        "{bismark_path}/bismark_methylation_extractor --output bismark_txt --gz "
        "--samtools_path {samtools_path} --yacht {input} 2> {log.err} 1> {log.out}"


rule bismark_me2:
    input:
        "bismark_bam/{stem}.bam"
    output:
        "bismark_txt/CpG_OT_{stem}.txt.gz" ## etc
    log:
        err = "logs/bismark_me2/{stem}.err",
        out = "logs/bismark_me2/{stem}.out"
    shell:
        "{bismark_path}/bismark_methylation_extractor --output bismark_txt -gz "
        "--samtools_path {samtools_path} {input} 2> {log.err} 1> {log.out}"

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
        txt=f"bismark_txt/any_C_context_local.{{id}}_bismark_{mapper}.sorted.ds.txt.gz"
    output:
        html=f"r_out/{proj_name}report_{{id}}.html",
        # pdf=f"r_out/{proj_name}report_{{sample}}.pdf",
        txt=f"r_out/{proj_name}report_{{id}}.txt"
    params:
        stub=f"r_out/{proj_name}report_{{id}}"
    log:
        err=f"logs/r/{proj_name}{{id}}.err",
        out=f"logs/r/{proj_name}{{id}}.out"
    shell:
        "{rscript_path}/Rscript -e " # Note lowercase s in Rscript on Unix!
        #"'Sys.setenv(RSTUDIO_PANDOC=\"/Applications/RStudio.app/Contents/MacOS/pandoc\"); "
        "'rmarkdown::render(\"{input.script}\", "
        "params=list(ref_file=\"{ref_file}\", in_file=\"{input.txt}\", out_txt=\"{output.txt}\"), "
        "output_format=\"html_document\", "
        "output_file=\"{params.stub}\")' 2> {log.err} 1> {log.out}"
        ## can one specify mutiple output names? not making stub.pdf and stub.html
