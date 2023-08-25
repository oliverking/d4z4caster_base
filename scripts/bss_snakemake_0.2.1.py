## https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html
## https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html

# can override some or all with snakemake --configfile or --config
configfile: "/scripts/config_0.2.1.yaml"

import os
import glob
import re
import pandas as pd

################################################################################
#################### Moved these to config.yaml file ###########################
################################################################################

verbose = config["verbose"]
use_csv = config["use_csv"]
no_match = config["no_match"]
input_type = config["input_type"]
paired_end = config["paired_end"]
bam_id_regex = config["bam_id_regex"]
bam_id_regex_fallback = config["bam_id_regex_fallback"]
fastq_id_regex = config["fastq_id_regex"]
fastq_id_regex_fallback = config["fastq_id_regex_fallback"]
fasta_id_regex = config["fasta_id_regex"]
fasta_id_regex_fallback = config["fasta_id_regex_fallback"]
info = config["info"]
info_glob = config["info_glob"]
info_glob_fallback = config["info_glob_fallback"]
info_skip_rows = config["info_skip_rows"]
id_col = config["id_col"]
id_idx = config["id_idx"]
new_col = config["new_col"]
new_col_fallback = config["new_col_fallback"]
new_idx = config["new_idx"]
report_format = config["report_format"]
ref_path = config["ref_path"]
ref_subdir = config["ref_subdir"]
ref_fasta = config["ref_fasta"]
ref_meta = config["ref_meta"]
max_per_ref = config["max_per_ref"]

################################################################################
################################################################################

if verbose:
    print("Options from yaml(s) and/or command line:")
    print("# verbose:", verbose)
    print("# use_csv:", use_csv)
    print("# no_match:", no_match)
    print("# input_type:",  input_type)
    print("# paired_end:", paired_end)
    print("# bam_id_regex:", bam_id_regex)
    print("# bam_id_regex_fallback:", bam_id_regex_fallback)
    print("# fastq_id_regex:", fastq_id_regex)
    print("# fastq_id_regex_fallback:", fastq_id_regex_fallback)
    print("# fasta_id_regex:", fasta_id_regex)
    print("# fasta_id_regex_fallback:", fasta_id_regex_fallback)
    print("# info:", info)
    print("# info_glob:", info_glob)
    print("# info_glob_fallback:", info_glob_fallback)
    print("# info_skip_rows:", info_skip_rows)
    print("# id_col:", id_col)
    print("# id_idx:", id_idx)
    print("# new_col:", new_col)
    print("# new_col_fallback:", new_col_fallback)
    print("# new_idx:", new_idx)
    print("# report_format:", report_format)
    print("# ref_path:", ref_path)
    print("# ref_subdir:", ref_subdir)
    print("# ref_fasta:", ref_fasta)
    print("# ref_meta:", ref_meta)
    print("# max_per_ref:", max_per_ref)

################################################################################
################################################################################

# in Docker container these should all be in PATH so can likely be omitted
samtools_path = "/usr/bin"
hisat2_path = "/usr/bin"
bowtie2_path = "/usr/bin"
bismark_path = "/usr/bin"
rscript_path = "/usr/bin"
fastqc_path = "/usr/bin"
bbmap_path = "/usr/bin"

bismark_name = "bismark_odk" # has some added functionality vs standard bismark
mapper = "bt2"               # "hisat2" option had issues and is not supported

# log all version numbers somewhere?

# Note: Rscript default outpath is Rmd dir rather than calling dir!
# put in setup chunk?: knitr::opts_knit$set(root.dir = getwd())
rmd_path = "/scripts"
rmd_file = "bss_make_report_0.2.1.Rmd"
rmd_local = "script_" + rmd_file
r_file = "bss_summarize_run_0.2.1.R"
downsample_script = "/scripts/bss_downsample_bam_0.1.7.py"
summary_pdf = "run_summary.pdf"               # add project name or version number?
summary_forced_pdf = "run_summary_forced.pdf" # add project name or version number?

ref_folder = ref_path + "/" + ref_subdir
ref_fasta_full = ref_folder + "/" + ref_fasta  # name gets passed to Rscript
ref_meta_full = ref_meta   #  an already include same or different path

bismark_temp = "bismark_temp"

# Note: glob_wildcards returns a tuple so comma is needed
# INFILES, = glob_wildcards("2020-09-22_BC-fusion_library/BAM/{ids}.bam")
# ls 2020-09-22_BC-fusion_library/BAM/*.bam
# ...MethylationAnalysis.bam ...-Fusion_Prim.bam  ...-Fusion_Primers_230.bam

proj_name = ""
bam_dir = "UBAM"
bam_dir_fallback= "BAM" # if bam_dir not found
fastq_dir = "fastq"     # not to be changed
fasta_dir = "fasta"     # not to be changed

if report_format == "html_document":
    report_ext = "html"
elif report_format == "pdf_document":
    report_ext = "pdf"
elif report_format == "ioslides_presentation" :
    report_ext = "slides.html"
elif report_format == "word_document" :
    report_ext = "docx"
else:
    print("unknown report_format; quitting")
    quit()

# os.path.expanduser(path); os.path.expandvars(path); use recursive=T?
# note: dir and file names should be case-sensitive if underlying OS filesystem is
if input_type == "bam" :
    if not os.path.isdir(bam_dir) :
        bam_dir = bam_dir_fallback
        if not os.path.isdir(bam_dir) :
            print("directories", bam_dir, "and", bam_dir_fallback, "not found")
            quit()
    INFILES = glob.glob(bam_dir + "/*.bam")
    id_regex = bam_id_regex
    id_regex_fallback = bam_id_regex_fallback
elif input_type == "fastq" :
    if not os.path.isdir(fastq_dir) :
        print("directory", fastq_dir, "not found")
        quit()
    INFILES = glob.glob(fastq_dir + "/*.fastq.gz")
    id_regex = fastq_id_regex
    id_regex_fallback = fastq_id_regex_fallback
    if verbose and use_csv: print("use_csv set to False for fastq input")
    use_csv = False

elif input_type == "fasta" :
    if not os.path.isdir(fasta_dir) :
        print("directory", fasta_dir, "not found")
        quit()
    INFILES1 = glob.glob(fasta_dir + "/*.fasta")
    INFILES2 = glob.glob(fasta_dir + "/*.fasta.gz")
    INFILES3 = glob.glob(fasta_dir + "/*.fa")
    INFILES4 = glob.glob(fasta_dir + "/*.fa.gz")
    # disjoint so no dups when cconcatenated
    INFILES = INFILES1 + INFILES2 + INFILES3 + INFILES4
    # or: MAYBES = glob.glob(fasta_dir + "/*.fa*")
    # INFLILES = [x for x in MAYBES if re.match("\\.fa(sta)?(.gz)?$)",x)]
    id_regex = fasta_id_regex
    id_regex_fallback = fasta_id_regex_fallback
else:
    print("unknown input file type", input_type, "; quitting")
    quit()


MYDICT = {} # to associate new names with input filenames

if verbose: print(INFILES)

if use_csv:
    if os.path.isfile(info):
        if verbose: print("located", info)
        info_file = [info] #list of length 1
    else:
        if verbose: print("did not locate", info, "; trying info_glob")
        info_file = glob.glob(info_glob)
        if (len(info_file) == 0):
            if verbose: print("no matches to", info_glob, "; trying info_glob_fallback")
            info_file = glob.glob(info_glob_fallback)

    if verbose: print("info_file:", info_file)

    if len(info_file) > 0 :
        if len(info_file) > 1 :
            if verbose: print(len(info_file), "files found matching pattern, using first")
        info_file = info_file[0] # even if length is 1, want string rather than list
        if verbose:
            print("sample info file" , os.path.basename(info_file))
            print("project name:" , proj_name)
        samp_info = pd.read_csv(info_file, skiprows=info_skip_rows)
        if id_col in samp_info.columns :
            samp_info = samp_info.set_index(id_col, drop=False)
        else :
            samp_info = samp_info.set_index(samp_info.columns[id_idx], drop=False)
            if verbose: print(id_col, "column not found in samp_info, using", samp_info.columns[0])
        if verbose: print(samp_info)
        if new_col in samp_info.columns :
            name_col = new_col
        elif new_col_fallback in samp_info.columns :
            name_col = new_col_fallback
        else : # TODO: test for errors here?
            name_col = samp_info.columns[new_idx]
            if verbose: print(new_col, "and", new_col_fallback, "not found in samp_info, using", name_col)
        use_samp_info = True
    else:
        print("sample info file not found in project directory; quitting")
        quit()
else :
    use_samp_info = False

if verbose: print ("Mappings from filenames -> pattterns from regex -> name stems from csv:")
for i in INFILES :
    MYMATCH = re.search(id_regex, i)
    if MYMATCH is None :
        MYMATCH = re.search(id_regex_fallback, i)
    if MYMATCH is None:
        if no_match == "quit" :
            print ("## sample", i, "does not have expected file name pattern; quitting")
            quit()
        elif no_match == "skip" :
            print ("## sample", i, "does not have expected file name pattern; skipping")
        else :
            print ("## unknown no_match; quitting")
            quit()
    else:
        MYSTEM = MYMATCH.group(1)
        if use_samp_info and (MYSTEM in samp_info.index) :
            MYDICT[i] = samp_info[name_col][MYSTEM]
        else :
            if verbose: print("## no match for", MYSTEM, "in sample_info")
            MYDICT[i] = MYSTEM
        if verbose: print("#", i, '->', MYSTEM, '->', MYDICT[i])


# print(MYDICT.values())
INVDICT = {value: key for key, value in MYDICT.items()}

if verbose:
    print ("Mappings from name stems to input filenames:")
    for id in INVDICT:
        print("#", id, '->', INVDICT[id])


# can make some output files conditional on e.g. input_type, if desired
# https://stackoverflow.com/questions/64949149/is-it-possible-to-add-a-conditional-statement-in-snakemakes-rule-all
rule all:
    input:
        summary_pdf,
        expand(f"bismark_bam/local.{{id}}_bismark_{mapper}.sorted.ds.bam.bai", id=INVDICT.keys()),
        expand(f"bismark_bam/local.{{id}}_bismark_{mapper}.sorted.ds.100.bam.bai", id=INVDICT.keys()),
        expand(f"r_out/{proj_name}report_{{id}}.{report_ext}", id=INVDICT.keys()),
        expand("fastqc_out/{id}_fastqc.html", id=INVDICT.keys())
        ## files not made if there are no reads in that orientation?
        # expand("CpG_OT_local.{id}_bismark_hisat2.txt.gz", id=INVDICT.keys()),


rule index_bam:
    input:
        "bismark_bam/{sample}.bam"
    output:
        "bismark_bam/{sample}.bam.bai",
    log:
        err = "logs/index_bam/{sample}.err",
        out = "logs/index_bam/{sample}.out"
    shell:
        "{samtools_path}/samtools index {input} 2> {log.err} 1> {log.out}"


rule index_stats:
    input:
        bam="bismark_bam/{sample}.bam",
        bai="bismark_bam/{sample}.bam.bai"
    output:
        "bismark_bam/{sample}.bam.idxstats"
    log:
        err = "logs/index_stats/{sample}.err"
        # out = "logs/index_stats/{sample}.out"
    shell:
        "{samtools_path}/samtools idxstats {input.bam} > {output} 2> {log.err} "


rule sort_bam:
    input:
        "bismark_bam/{sample}.bam"
    output:
        "bismark_bam/{sample}.sorted.bam"
    log:
        err="logs/sort_bam/{sample}.err",
        out="logs/sort_bam/{sample}.out"
    shell:
        "{samtools_path}/samtools sort {input} -o {output} 2> {log.err} 1> {log.out}"


## Flag EXACT duplicates? (Not just same start/end coords)
## What about same cigar string (PICARD)
## -- not enough, as ID of mismatches may differ
## May want to do this later, since it may flatten out distinction between
## true sequence and sequencing errors.
rule dedup_bam:
    input:
        "bismark_bam/{sample}.bam"
    output:
        "bismark_bam/{sample}.dedup.bam"
    log:
        err="logs/dedup_bam/{sample}.err",
        out="logs/dedup_bam/{sample}.out"
    shell:
        "{samtools_path}/samtools sort {input} -o {output} 2> {log.err} 1> {log.out}"


# target max number of reads to retain per chromsome, before deduplication
# max_per_ref = 10000  ## moved to config.yaml

rule downsample_bam:
    input:
        bam="bismark_bam/{sample}.bam",
        bai="bismark_bam/{sample}.bam.bai",
        idx="bismark_bam/{sample}.bam.idxstats"
    output:
        bam="bismark_bam/{sample}.ds.bam",
        bai="bismark_bam/{sample}.ds.bam.bai",
        idx="bismark_bam/{sample}.ds.bam.idxstats"
    params:
        "bismark_bam/{sample}"
    log:
        err="logs/downsample_bam/{sample}.err",
        out="logs/downsample_bam/{sample}.out"
    shell:
        "python {downsample_script} {params} {max_per_ref} 2> {log.err} 1> {log.out}"


## downsample again, with max_per_ref 100, for genome browser plots
## renamed placeholder {sample} to {sampleds} just for clarity
rule downsample_bam_100:
    input:
        bam="bismark_bam/{sampleds}.bam",
        bai="bismark_bam/{sampleds}.bam.bai",
        idx="bismark_bam/{sampleds}.bam.idxstats"
    output:
        bam="bismark_bam/{sampleds}.100.bam",
        bai="bismark_bam/{sampleds}.100.bam.bai",
        idx="bismark_bam/{sampleds}.100.bam.idxstats"
    params:
        "bismark_bam/{sampleds}"
    log:
        err="logs/downsample_bam_again/{sampleds}.err",
        out="logs/downsample_bam_again/{sampleds}.out"
    shell:
        "python {downsample_script} {params} 100 100 2> {log.err} 1> {log.out}"

## downsample_bam also indexes the bam, so give priority
ruleorder: downsample_bam > downsample_bam_100 > index_bam
ruleorder: downsample_bam > downsample_bam_100 > index_stats


rule fastqc:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        "fastqc_out/{sample}_fastqc.html"
    log:
        err="logs/fastqc/{sample}.err",
        out="logs/fastqc/{sample}.out"
    shell:  ## -o for other dir. Needs to be created in advance? Or will snakemake do this?
        "{fastqc_path}/fastqc {input} -o fastqc_out 2> {log.err} 1> {log.out}"


# Convert (U)BAM to fastq, for input to BISMARK. Optionally change sample names, too,
#  with wildcards: ttps://snakemake.readthedocs.io/en/stable/project_info/faq.html
if input_type == "bam" :
    rule bam2fq:
        input:
            lambda wildcards: INVDICT[wildcards.sample]
        output:
            "fastq/{sample}.fastq.gz"
        log:
            err="logs/bam2fq/{sample}.err",
            out="logs/bam2fq/{sample}.out"
        shell:
            "{samtools_path}/samtools fastq -F 0x900 {input} | gzip > {output} 2> {log.err}"

# BISMARK should work with fasta input, but we'd have to skip the fastqc step.
# Instead, convert fasta files to fastq.gz files (with fake QUALS) and procede from there
# FIXME: allow different extensions (e.g. fa, fa.gz, fasta, fasta.gz) and new names using lambda wildcards?
if input_type == "fasta" :
    rule fa2fq:
        input:
            lambda wildcards: INVDICT[wildcards.sample]
            # "fasta/{sample}.fasta.gz"
        output:
            "fastq/{sample}.fastq.gz"
        log:
            err="logs/fa2fq/{sample}.err",
            out="logs/fa2fq/{sample}.out"
        shell:
            "{bbmap_path}/bbmap/reformat.sh qfake=30 ow=t in={input} out={output} 2> {log.err} 1> {log.out}"


# TODO: add abi2fa or abi2fq?
# https://www.biostars.org/p/309727/
# https://rdrr.io/bioc/CrispRVariants/man/abifToFastq.html

## output not explictly used by name in shell command so used as dummy ix2 in bismark_bt
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


### pre-create to avoid error if two calls to bismark_hs try to autoconstruct it at the same time?
rule create_temp:
    output:
        directory(bismark_temp)
    shell:
        "mkdir -p {output}"


rule bismark_bt:
    input:
        fq = "fastq/{sample}.fastq.gz",
        ix = ref_folder,
        ix2 = ancient(f"{ref_folder}/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"),
        temp = bismark_temp
    output:
        bam = "bismark_bam/local.{sample}_bismark_bt2.bam",
        ser = "bismark_bam/local.{sample}_bismark_bt2_SE_report.txt"
    log:
        err = "logs/bismark_bt2/{sample}.err",
        out = "logs/bismark_bt2/{sample}.out"
    shell:
        "{bismark_path}/{bismark_name} --output_dir bismark_bam --temp_dir {input.temp} --non_directional "
        "--bowtie2 --path_to_bowtie2 {bowtie2_path} --samtools_path {samtools_path} --local --prefix local "
        "--genome {input.ix} -q {input.fq} 2> {log.err} 1> {log.out}"

# Note: it is possible that the output file any_C_context file will legitimately be
# empty,in which case it is automatically deleted by bismark_methylation_extractor.
# Adding "&& touch {output}" is so that an empty file is re-created in this case.
rule bismark_me:
    input:
        "bismark_bam/{stem}.bam" ## run on just sorted.ds.bams based on rule all
    output:
        "bismark_txt/any_C_context_{stem}.txt.gz" ## etc. May be others
    log:
        err = "logs/bismark_me2/{stem}.err",
        out = "logs/bismark_me2/{stem}.out"
    shell:
        "( {bismark_path}/bismark_methylation_extractor --output bismark_txt --gz "
        "--samtools_path {samtools_path} --yacht {input} && touch {output} ) 2> {log.err} 1> {log.out}"


rule bismark_me2:
    input:
        "bismark_bam/{stem}.bam"
    output:
        "bismark_txt/CpG_OT_{stem}.txt.gz" ## etc -- list others?
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
# "gzcat {input.txt} | head > {output}"  # gives error since gzcat doesn't finish
# (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
rule cpg_matrix:
    input:
        script=rmd_local,
        txt=f"bismark_txt/any_C_context_local.{{id}}_bismark_{mapper}.sorted.ds.txt.gz",
        bam=f"bismark_bam/local.{{id}}_bismark_{mapper}.sorted.ds.100.bam",
        bai=f"bismark_bam/local.{{id}}_bismark_{mapper}.sorted.ds.100.bam.bai", # used implicitly?
        ser=f"bismark_bam/local.{{id}}_bismark_{mapper}_SE_report.txt",
        idx=f"bismark_bam/local.{{id}}_bismark_{mapper}.sorted.bam.idxstats"
    output:
        ext=f"r_out/{proj_name}report_{{id}}.{report_ext}",
        txt=f"r_out/{proj_name}report_{{id}}.txt"
    params:
        stub=f"r_out/{proj_name}report_{{id}}"
    log:
        err=f"logs/r/{proj_name}{{id}}.err",
        out=f"logs/r/{proj_name}{{id}}.out"
    shell:
        "{rscript_path}/Rscript -e " # Note lowercase s in Rscript on UNIX!
        #"'Sys.setenv(RSTUDIO_PANDOC=\"/Applications/RStudio.app/Contents/MacOS/pandoc\"); "
        "'rmarkdown::render(\"{input.script}\", "
        "params=list(ref_fasta=\"{ref_fasta_full}\", ref_meta=\"{ref_meta_full}\", "
        "max_reads=\"{max_per_ref}\", in_file=\"{input.txt}\", ser_file=\"{input.ser}\", "
        "idx_file=\"{input.idx}\", bam_file=\"{input.bam}\", out_txt=\"{output.txt}\", "
        "format=\"{report_format}\"), "
        "output_format=\"{report_format}\", "
        "output_file=\"{output.ext}\")' 2> {log.err} 1> {log.out}"
        # can one specify mutiple output names? not making stub.pdf and stub.html
        # had "output_file=\"{params.stub}\")', but it looks like .html does not get automatically added
        #  as file extension if params stub already has a dot in it


rule summary_plot:
    input:
        expand(f"r_out/{proj_name}report_{{id}}.txt", id=INVDICT.keys()),
        expand(f"bismark_bam/local.{{id}}_bismark_{mapper}_SE_report.txt", id=INVDICT.keys()),
        expand(f"bismark_bam/local.{{id}}_bismark_{mapper}.sorted.bam.idxstats", id=INVDICT.keys())
    output:
        summary_pdf
    log:
        err = "logs/summary_plot.err",
        out = "logs/summary_plot.out"
    shell:
        "{rscript_path}/Rscript --vanilla {rmd_path}/{r_file} . {output} 2> {log.err} 1> {log.out}"


# make summary plot from existing output files, even if some failed
# snakemake --forcerun summary_plot_force
rule summary_plot_force:
    input:
    output:
        summary_forced_pdf
    log:
        err = "logs/summary_plot_forced.err",
        out = "logs/summary_plot_forced.out"
    shell:
        "{rscript_path}/Rscript --vanilla {rmd_path}/{r_file} . {output} 2> {log.err} 1> {log.out}"

## END
