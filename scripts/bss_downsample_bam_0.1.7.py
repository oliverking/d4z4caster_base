import subprocess
import sys
import shutil
import os
import pandas

## could do this all as shell script: while read COL_1 COL_2 ...

# bam_stub = "local.PLJ10378_bismark_bt2.sorted"
if len(sys.argv) > 1:
    bam_stub = sys.argv[1]

max_per_chrom = 10000
if len(sys.argv) > 2:
    max_per_chrom = int(sys.argv[2])

new_ext = ".ds.bam"
if len(sys.argv) > 3:
    new_ext = "." + sys.argv[3] + ".bam"

## Refinement done in memory, not good for very large max
exact_count = max_per_chrom <= 100000

bam_name = bam_stub + ".bam"
bam_index = bam_stub + ".bam.bai"
bam_stats = bam_name + ".idxstats"
sam_temp = bam_stub + ".temp.sam"
sub_temp = bam_stub + ".temp.sub.sam"
bam_new = bam_stub + new_ext
bam_new_stats = bam_stub + new_ext + ".idxstats"

print("## input bam:",  bam_name)
print("## target max reads per ref:", max_per_chrom)
print("## NOTE: actual max reads may vary by approx sqrt of target number")
## Can also overshoot by a bit, then do exact downsampling of those,
## e.g. with shuf or BBmap shuffle.sh or UNIX sort -R | head -n X | sort -k3

if not os.path.exists(bam_name):
    print("could not locate", bam_name, file=sys.stderr)
    sys.exit(1)

if not os.path.exists(bam_index):
    print("bam index not found, so indexing", bam_name)
    subprocess.run(["samtools", "index", bam_name])

if not os.path.exists(bam_stats):
    print("idxstats not found, so creating for", bam_name)
    stats_file = open(bam_stats, 'w')
    subprocess.run(["samtools", "idxstats", bam_name], stdout=stats_file)
    stats_file.close()

idxstats = pandas.read_table(bam_stats, names=["seqname","length","count","other"], index_col=0)

## assumes samtools_path is in path. samtools_path = "/usr/bin"

#  exclude "*" here? has zero counts for these files
if max(idxstats["count"]) <= max_per_chrom :
    print("## no downsampling needed, so just copying bam, bai, idxstats")
    shutil.copyfile(bam_name, bam_new)
    shutil.copyfile(bam_index, bam_new + ".bai")
    shutil.copyfile(bam_stats, bam_new + ".idxstats")
    sys.exit(0)

temp_file = open(sam_temp, 'w')
# could write everything to bam, use samtools merge...
print("## writing bam header to temp file", sam_temp )
subprocess.run(["samtools", "view", "-H", bam_name], stdout=temp_file)
##subprocess.run(["samtools", "view", "-H", bam_name, ">>", sam_temp], shell=True)

print("## iterating over", len(idxstats.index), "ref sequences")
for sn in idxstats.index :
    # includes sn = "*" but it appears that gets escaped or quoted properly below
    num = idxstats["count"][sn]
    if num > max_per_chrom :
        s = round(max_per_chrom/num, 5)
        spp = round((max_per_chrom + 5*(max_per_chrom**0.5))/num, 5)
        ## don't want spp>1 since init digits are seed; 1.00 is okay, but not 1?
        print("###", sn, ", n =", num, ", s =", s, ", spp = ", spp)
        ## type(index_stats["count"][4]/max_per_chrom) is <class 'numpy.float64'>
        if exact_count:
            if spp >= 1 :
                sam = subprocess.Popen(["samtools", "view", bam_name, sn], stdout=subprocess.PIPE)
            else :
                sam = subprocess.Popen(["samtools", "view", "-s", str(spp), bam_name, sn], stdout=subprocess.PIPE)
            sort1 = subprocess.Popen(["sort", "-R"], stdin=sam.stdout, stdout=subprocess.PIPE)
            head = subprocess.Popen(["head", "-n",  str(max_per_chrom)], stdin=sort1.stdout, stdout=subprocess.PIPE)
            sort2 = subprocess.Popen(["sort", "-k3"], stdin=head.stdout, stdout=temp_file)
            sort2.communicate()
            ## should double-check that this works properly!
            ### kill subprocess now??
        else:
            subprocess.run(["samtools", "view", "-s", str(s), bam_name, sn], stdout=temp_file)
    else:
        print("###", sn, ", n = ", num, ", s = NA", )
        subprocess.run(["samtools", "view", bam_name, sn], stdout=temp_file)

temp_file.close()

out_file = open(bam_new, 'w')
print("## converting ", sam_temp, "to", bam_new)
subprocess.run(["samtools", "view", "-b", sam_temp], stdout=out_file)
out_file.close()

print("## removing ", sam_temp)
os.remove(sam_temp)

print("## indexing", bam_new, flush=True)
subprocess.run(["samtools", "index", bam_new])

print("## writing idxstats", bam_new, flush=True)
stats_file = open(bam_new_stats, 'w')
subprocess.run(["samtools", "idxstats", bam_new], stdout=stats_file)
stats_file.close()

## could just cat from bam_stats!
print("## idxstats before downsampling:", flush=True)
subprocess.run(["samtools", "idxstats", bam_name])

## could just cat from bam_new_stats!
print("## idxstats after downsampling:", flush=True)
subprocess.run(["samtools", "idxstats", bam_new])
