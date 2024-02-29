#  Instructions for D4Z4caster program

---

## Installation

The simplest way to obtain the program, including all its dependencies, is with Docker.

It is first necessary to install and start Docker; see instructions at https://docs.docker.com/get-docker/

Then open a terminal window on your computer. Your computer will be referred to as the "host computer", and the Docker container its "guest". In this document we use `host$` as prefix for commands that are run from the terminal of the host computer, as opposed to commands run from within the Docker container; `host$` is *not* part of the command that is typed. Also, long commands may be wrapped across multiple lines in this document but should be entered as single lines at the command-line (or with a backslash `\` at the end of each line before the last).

**Option 1**: pull image from DockerHub. From the command-line of the host computer:

*(currently a private repository)*

`host$ docker pull kingod/d4z4caster` 

Rename local copy so that subsequent commands use d4z4caster rather than kingod/d4z4caster

`host$ docker tag kingod/d4z4caster  d4z4caster`   
  

**Option 2**: get docker file and supporting files from GitHub and build docker image locally. (Building image may take 30 minutes or longer.) From the command-line of the host computer:


*(currently a private repository)*

`host$ git clone oliverking/d4z4caster_base`

`host$ docker build -t d4z4caster d4z4caster_base`

---

## Basic usage

### Starting Docker container in interactive mode

Basic command to launch docker container. From the command-line of the host computer:

`host$ docker run --rm -it --entrypoint /bin/bash d4z4caster`

Increasing memory to 10g and mounting project directory `/path/of/host/project/` on host computer as `/bss`

`host$ docker run --memory=10g --rm -it --entrypoint /bin/bash -v /path/of/host/project/:/bss d4z4caster`

Here `/path/of/host/project/` should be changed to the **full** (not relative) path the the project directory on the host computer.  

Note: on Windows host computers, Windows-style paths may be needed, e.g. `C:\path\of\host\project\`

The starting working directory in the Docker container is `/bss`, and the user should now see a command-line prompt *within* the Docker container that should end in `:/bss#` . By mounting your host directory to `/bss` the files you mounted should be visible with the command `ls -la`.


### Running program from prompt within Docker container

These commands are run from the command-line **within** the docker container.

Use `-n` (or `--dry-run`) to see jobs that *would* be run:

`snakemake -n -p -c1 -s /scripts/bss_snakefile.py`

Run workflow using `n` cores (`-cn` or `--cores n`), here with `n = 1`:

`snakemake -p -c1 -s /scripts/bss_snakefile.py`

### Exiting Docker container

`exit`


### Alternative: starting Docker container and running program in one command

The default entrypoint for the Docker container is the command

`snakemake -p -c1 -s scripts/bss_snakefile.py`

but in the examples above this was overridden with `/bin/bash` to give an interactive prompt.

One can use the default entrypoint to directly run the snakemake command non-interactively as follows:

`host$ docker run --memory=10g --rm -t -v /path/of/host/project/:/bss d4z4caster`

Here for a dry-run, the `-n` option can be added to the end of the command, as follows:

`host$ docker run --memory=10g --rm -t -v /path/of/host/project/:/bss d4z4caster -n`

Similarly, the other configuration options for snakemake discussed below can be added to the end of this command.


### Demo data

A small demo dataset is included in the `/demo` directory of the Docker image.  This can be used to test the installation:

`host$ docker run --memory=10g --rm -it --entrypoint /bin/bash d4z4caster`

`cd /demo`

`snakemake -p -c1 -s /scripts/bss_snakefile.py`

Results will be in the `/demo` directory of the Docker container, and not by default visible in the host filesystem.  There are several ways around this. One is mounting an *empty* host directory as `/bss` when starting the Docker container and copying the files from `/demo` into `/bss` before running the snakemake workflow

`host$ docker run --memory=10g --rm -it --entrypoint /bin/bash -v /path/of/empty/host/folder/:/bss d4z4caster`

`cp -i -r /demo/* .`

`snakemake -p -c1 -s /scripts/bss_snakefile.py`

*CAUTION*: For the above, please make sure the host directory you mount as `/bss` is empty, so that the demo files do not overwrite or get intermingled with your real files! The `-i` flag in `cp -i` should warn you if you are about to overwrite any files, but checking that the working directory `/bss` is empty with `ls -lah` before the copy command may be prudent.

Note that the two input BAMS in the `/demo/UBAM` are identical, so results for both samples should be the same since the same seeds are used for randomized steps (e.g. downsampling).  

`diff r_out/report_demoA.txt r_out/report_demoB.txt`

Also, in these demo bam files all reads uniquely map to one of the reference sequences. This will not typically be the case for real data. Here we have filtered out all reads that do not map to one of the reference sequences to avoid any potential privacy concerns (e.g., if enough SNPs could be called from reads from non-targeted regions scattered throughout the genome there might be a risk of subject re-identification).

---

## Input Files

Default options are based on file types and naming conventions used in our sequencing workflow; see section "Options" for how to change these options.  

### BAM files

Reads for each sample are in unaligned BAM (uBAM) files, one file per sample, with single-end reads. In the default workflow each filename is of the form `IonXpress_XYZ_rawlib.basecaller.bam` for a three digit number XYZ.  These files must be in a top-level subdirectory of the project directory, named `UBAM` (first checked) or `BAM` (fallback).  

Note: these directory names are case-sensitive, if the operating system's filesystem is case sensitive. The same is true of other file names discussed below.

Note: if the reads in the input BAMs have *already* been aligned to some reference, the alignments will be ignored: the first step of the pipeline is converting the BAMS to fastq files, which will be mapped by Bismark/Bowtie2. It is important that unmapped reads from any earlier alignment are included in the BAM (unless the same reference sequences and same bisulfite-aware aligner had been used) and that multi-mapping reads are included only once. Users may consider using `java -jar picard.jar RevertSam -I in.bam -O out.bam` for preprocessing if reads in input BAMs have already been aligned.

### Sample info file

Then snakemake workflow generates many output files in many subfolders, and for the most part the files derived from a given input file have the same basename as the input file but with different prefixes and suffixes/extensions.  There is however the opportunity at the first step – the conversion of BAM files to fastq file - to change the basename of all subsequent files to either a specific substring of the original filename or any other name of the user's choosing.

So for example `UBAM/IonXpress_001_rawlib.basecaller.bam` could give rise to `fastq/IonXpress_001.fastq` instead of `fastq/SampleA.fastq.gz` instead of `fastq/IonXpress_001_rawlib.basecaller.fastq.gz`, and this change carries through the remainder of the pipeline, e.g. generating `r_out/report_IonXpress_001.html` or `r_out/report_SampleA.html` instead of `r_out/report_IonXpress_001_rawlib.basecaller.html`

The mapping between original filenames and new names relies on two things:

First, a regular expression (`bam_id_regex`) defining what portion of the initial file names to extract.  The default is `"(IonXpress_[0-9]+)_"`, which will extract just the matching pattern in parentheses from a filename, so will for example extract `IonXpress_001` from `UBAM/IonXpress_001_rawlib.basecaller.bam`. (Note that the folder name `UBAM/` is considered part of the filename here).

Second, a file (with default name `sample_info.csv`) that includes a column (default name `Barcode`) for extracted IDs and a column for new names (default name `Sample Name (requires)`). If the extracted ID is found in the former of these columns it is replaced by the corresponding value in the latter of these columns, so e.g. can be used to replace `IonXpress_001` by `SampleA`; if not, the extracted ID is used for output files.  If an ID cannot be extracted from a filename because it does not match the `bam_id_regex` pattern, a more permissive pattern `bam_id_regex_fallback` is checked. If that also fails then the sample is skipped, though there is also the option to quit in such cases (`no_match="quit"`).  

The default options for describing the layout of the `sample_info.csv` file are designed for the csv files exported by Ion Torrent software, which have the following structure:

```
CSV Version (required),2.2,,,,,,,
Barcode,Control Type,Sample Name (required),Sample ID,Sample Description,DNA/RNA/Fusions,Reference library,Target regions BED file,Hotspot regions BED file
IonXpress_001,,SampleA,,Saliva,DNA,,,
IonXpress_002,,SampleB,,Saliva,DNA,,,
...
```

These options can be changed to describe other file layouts, as described in the next section. The options include fallback values that can also be changed; e.g., if a file `sample_info.csv` is not found in the top-level project directory then the default is to look for a filename containing the string `sample` and ending in `.csv`.

The sample info file is not strictly required and can be ignored by setting `use_csv` to `False`; in this case one can still use the `bam_id_regex` to extract portions of filenames to be used for output files.


### Reference files

The references sequences and indices to which reads are aligned to are in the `/refs` directory; this contains:

- a subdirectory (default name `ref_v3/`) containing:
  - a single fasta file (default name `pad100.fa`) with reference sequences that have been padded with 100 Ns on both sides.  (The padding is to avoid problems with Bismark extracting CpG context for alignments that include the first or last base of a reference sequence.)  
  - an optional yaml file, by default having the same name as the fasta file but with `.yaml` added. This can be used to specify what range of CpG sites in each reference sequence to use for scoring, which can be useful for example if coverage tends to be poor toward the 5' end.
- a subdirectory `Bisulfite_Genome/` with Bowtie2 index files generated by Bismark.  

Notes:

- The format of the yaml is as below, consisting of lines with (unindented) seqnames that exactly match names in the fasta file followed by lines with (indented) indices of first and last CpGs in the range to be used for scoring. The first CpG in a sequence has index `1` and the index `.inf` can be used for the last CpG if one doesn't want to bother counting the total number. If a seqname from the fasta does not appear in the yaml at all, all CpGs for it will be used for scoring. If a seqname in the yaml does not appear in the fasta file, that entry in the yaml is ignored.

```
seqname_1
  - 10
  - 40
seqname_2
  - 4
  - 50
...  
```

- If the user changes the reference fasta,  new index files should be automatically generated by the snakemake workflow, though these will exist only transiently within the current container. They could be made permanent by e.g. snapshotting/checkpointing the container or by mounting a directory on the host directory as `/refs`.  

- This workflow is designed for mapping of reads to a small number of *targeted* sequences (targeted by PCR amplification), rather than for genome-scale mapping. Some steps, for example, make plots of every CpG in every reference sequence, or expect a high proportion of reads to particular reference sequences to cover a particular CpGs, and will not work well if there are very many or very long reference sequences.

---

## Options

Various options for the snakemake workflow are specified in a file scripts/config.yaml that is read in by the snakemake.py file. This config file can be edited by the user, but as it resides within the docker image it may be simpler for the user to copy the config file into the project space, edit that version, and pass that filename when calling snakemake

`snakemake -p -c1 -s /scripts/bss_snakefile.py --configfile=myconfig.yaml`

The parameters from scripts/config.yaml are still read in, but their values will be overriden by any new values defined in `myconfig.yaml`. One can also override one or more parameters directly at the command line with the `--config` flag

`snakemake -p -c1 -s /scripts/bss_snakefile.py --config report_format=pdf_document verbose=True`

See the `config.yaml` file for a list of configurable options.

### Optional other input file formats:

**FASTQ**:

- Input files in fastq format can be used if option `input_type=fastq` is specified.
- This essentially just skips the first step of the default Snakemake workflow, so for compatibility with this, files should be in a directory named `/fastq` and should be gzipped with extensions `fastq.gz`,
and further renaming of output based on pattern matching and `sample_info.csv` is not done.

**FASTA**:

- Input files in fasta format can be used if option `input_type=fasta` is specified.  
- This is an experimental feature, to allow processing of low throughput Sanger-sequencing data, but the program is primarily designed with high-throughput data with hundreds or thousands of reads per reference sequence in mind.
- files should be in a directory named `/fasta` and can optionally be gzipped, with extensions `.fa`, `.fa.gz`, `.fasta`, or `.fasta.gz`
- renaming of output based on `sample_info.csv` proceeds just as for BAM input, though with pattern matching defined by `fasta_id_regex` rather than `bam_id_regex`.
- files are converted to fastq files (with fake Phred quality scores of 30) for compatibility with default Snakemake workflow; thus the quality scores in the `/fastqc` reports should be ignored, but other FastQC plots, e.g. of read-length distribution, may be useful.

---

## Output

The following directories are created, each containing one or more files corresponding to each of the input files in the `UBAM/` or `BAM/` directory

- `fastq/`
 - UBAM converted to fastq format (as Bismark does not directly accept BAM input).
- `fastqc_out/`
 - html report from running FastQC on fastq file
 - zip file containing this html file as well as some text summary files
- `bismark_bam/`
 - bam files with alignments of reads from fastq to reference sequences, using modified version of Bismark, both unsorted (`_bt2.bam`) and sorted by ref sequence and position (`_bt2.sorted.bam`), along with:
   - Bismark summary report for mapping (`_report.txt`)
   - downsampled versions of the sorted bams (`.ds.bam`) capped at 10000 reads per reference sequence (selected randomly).
   - further downsampled versions of the sorted bams (`.ds.100.bam`), capped at 100 reads per reference sequences (selected randomly).
   - index files (`.bai`) for sorted bams (including downsampled versions).
   - index stat files (`.idxstat`) for sorted bams (including downsampled versions).
   - fastq files (`_unmapped_reads.fq.gz`) with unmapped reads.
   - bam files (`.ambig.bam`) with randomly chosen mappings for ambiguously mapping reads .

- `bismark_txt/`
 - Output from bismark_methylation_extractor program on `.ds.bam` files:
   - table of methylation calls at every CpG and CH context for every read (`any_C_context_.*.txt.gz`)
   - methylation bias table (`.M-bias.txt`)
   - summary statistics (`_splitting_report.txt`)
- `r_out/`
 - report with plots and tables in format specified by `report_format` parameter:
   - `ioslides_presentation` for html slideshow (`.slide.html`),
   - `html_document` for single-page html (`.html`),
   - `pdf_document` for pdf (`.pdf`)
   - `word_document` for Word document (`.docx`) (table layout is clunky, though)
 - tab-delimited text file (`.txt`) with summary statistics
 - an R data file (`.Rdata`) with more extensive output that can be loaded into R (see `.Rmd` file for details)

There is also a `logs/` directory, with subdirectories for various steps of the snakemake workflow, and a `bismark_temp/` for temporary files created during mapping that should be automatically deleted once mapping is complete.  

The top-level project directory also includes a copy of the Rmarkdown script (`.Rmd`) that was used for generating the reports in the `r_out/` directory, and a file named `run_summary.pdf` that plots  per-sample summary statistics after all samples have been processed.

---

## Troubleshooting

- If the snakemake workflow fails with an error, look in the log files for clues as to what went wrong. The output in the terminal window, which is also also written to a file in the hidden folder `.snakmake/log/` often points to a log file for a specific sample at a specific step that causes the error.
- Rerunning the same snakemake command should pick up where the initial run left off. Since jobs that don't depend on one another may be run in a different order, the failure point may come latter or perhaps not at all. If rerunning, one might try the following options:
  - `-c1` (or `--cores 1`) to use a single core, if failure was with a higher number of cores (which may e.g. use more RAM).
  - `-k` (or `--keepon`) to run jobs for other samples even if jobs for some samples fail
  - `-w 20` (or `--latency-wait 20`), in case errors were caused by filesystem latency
  - `--config verbose=True`, if this option been set to False in the `config.yaml` file
  - `--ri` (or `--rerun-incomplete`) if snakemake cautions that some output is incomplete

Usually after all sample have been successfully processed a report named `run_summary.pdf` is generated. User can force generation of a summary plot named `run_summary_forced.pdf` from whichever files have already been successfully processed with this command:

`snakemake -p -c1 --force summary_plot_force -s /scripts/bss_snakefile.py`

---

### License

The custom code (R, python and perl) is licensed under a GNU General Public License v3.0 license (see file `LICENSE.txt` or `COPYING.txt`). This includes a modified version of Bismark perl code, which is itself licensed under GNU General Public License v3.0.

The workflow also uses unmodified versions of Bowtie2, samtools, BBmap, FastQC, and R with various libraries (Biostrings, bbmle, ggplot2, GViz, etc), some of which are subject to other licenses. See documentation for those packages for details.

### Copyright

Custom code and documentation are copyright (C) 2023 University of Massachusetts

### Citing

Until a manuscript detailing this NGS-based extension of our earlier work is published, please cite one or both of the following, as appropriate.

For assay and quartile-based scoring:

- Jones TI, Yan C, Sapp PC, McKenna-Yasek D, Kang PB, Quinn C, Salameh JS, King OD, Jones PL. Identifying diagnostic DNA methylation profiles for facioscapulohumeral muscular dystrophy in blood and saliva using bisulfite sequencing. Clin Epigenetics. 2014; 6(1):23. PMID: 25400706.

For mixture-of-beta-binomial scoring:

- Jones TI, King OD, Himeda CL, Homma S, Chen JC, Beermann ML, Yan C, Emerson CP, Miller JB, Wagner KR, Jones PL. Individual epigenetic status of the pathogenic D4Z4 macrosatellite correlates with disease in facioscapulohumeral muscular dystrophy. Clin Epigenetics. 2015; 7:37. PMID: 25904990.


###  Patent

Although this software is being made available under a GNU General Public License, note that there is a patent related to the DZ4Z methylation assay and analysis methodology:

- Jones PL, Jones T, Salameh J, Quinn C, King OD, inventors;  University of Massachusetts Medical Center, assignee. Molecular Diagnosis of FSHD By Epigenetic Signature.  Patent no. US 10,870,886 B2, issued 12/22/2020


### Contacts

- For software:
  - Oliver King (oliver.king@umassmed.edu)
- For assay:
  - Takako Jones (takakojones@med.unr.edu)
  - Peter Jones (peterjones@med.unr.edu)
