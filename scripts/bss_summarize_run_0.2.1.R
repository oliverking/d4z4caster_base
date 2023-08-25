## Collects up text files with summary statistics from r_out and bismark_bam folders 
## and plots these for all samples found. In snakemake workflow this gets called when all 
## processes are complete, but in case some jobs fail this can also be run manually.
## Note that in this case it make be that some bismark_bam summary files do not 
## have corresponding r_out summary files; if require_all == T then only samples with 
## all three summary files found are included.

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

require_all = T

args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) proj_dir = args[1] else proj_dir = "."
if (!dir.exists(proj_dir)) stop(proj_dir, "not found")

## v0.1.7: changed x-axis font size from 6 to 4 or 5, to accommodate more samples
## v0.2.0: changed stack2m from left_join to inner_join, and added require_all logic
##         (both for forcing summary plot to be made even if some samples failed)

# RScript --vanilla ~/bss_docker_v6/scripts/bss_summarize_run_0.2.0.R  ~/bisulfite_2020/JFBSS_0019_2022_01_26
# proj_dir = "~/bisulfite_2020/JFBSS_0019_2022_01_26"

## in current working dir
output_file = "run_summary_plot.pdf"
if (length(args) > 1) output_file = args[2] 

min_reads = -1  ## for plots requiring min read count per reference sequence; if -1 is adaptive
if (length(args) > 2) output_file = as.integer(args[3]) 

files1 = dir(path=file.path(proj_dir, "r_out"), pattern="^report.*.txt$", full=T)
dat1 = lapply(files1, read.table, header=T)
# sapply(dat1, dim)

## TODO: error-checking in case some are incomplete?

names(dat1) = sub(".txt","", basename(files1))
stack1 = do.call(rbind, dat1)
stack1$name = sub("report_", "", sub("\\.\\d+$", "", rownames(stack1))) 
stack1 = stack1 %>% mutate(filt = factor(ifelse(filtered=="YES", "filtered", "unfiltered"),
                                       levels=c("unfiltered", "filtered"))) 

files2 = dir(path=file.path(proj_dir, "bismark_bam"), pattern="sorted.bam.idxstats$", full=T)
dat2 = lapply(files2,  read.table, header=F, sep="\t")
names(dat2) = sub("local.", "", sub("_bismark_(bt2|hisat2).sorted.bam.idxstats","", basename(files2)))
stack2 = do.call(rbind, dat2) %>% filter(V1 != "*")
stack2$name =sub("report_", "", sub("\\.\\d+$", "", rownames(stack2))) 
colnames(stack2)[1:3] = c("ref_name", "length", "reads")

stack1$ref_name = factor(stack1$ref_name)
stack2$ref_name = factor(stack2$ref_name)
stack1$reads[which(is.na(stack1$reads))] = 0

max_reads = max(stack1$reads[stack1$filtered=="YES"])

if (min_reads == -1) {
    if (max_reads > 100) {
      min_reads = 100
    } else if (max_reads > 10) {
      min_reads = 10;
    } else {
      min_reads = 2
    }
}

stack1w = stack1 %>% dplyr::select(name, ref_name, filtered, reads) %>% 
  pivot_wider(names_from=filtered, values_from=reads) %>% mutate(frac=YES/pmax(NO,1))
  
# stack2m = left_join(stack2, stack1w %>% select(name, ref_name, frac))
# even without require_all, require both being joined here
stack2m = inner_join(stack2, stack1w %>% select(name, ref_name, frac))

files3 = dir(path=file.path(proj_dir, "bismark_bam"), pattern="SE_report.txt$", full=T)
dat3 = sapply(files3, FUN=function(x) grep("Mapping efficiency:", readLines(x), value=T)[1])
dat3b = lapply(files3, FUN=function(x) { y=readLines(x);
              y[grep("Mapping efficiency:", y)[1] + c(-2:3)]})
dat3c = lapply(dat3b, FUN=function(x) fread(text=sub("%", "", x)))
names(dat3c) = sub("local.", "", sub("_bismark_(bt2|hisat2)_SE_report.txt","", basename(files3)))
stack3c = dplyr::bind_rows(dat3c, .id="name") %>% 
  reshape2::melt(id.vars=c("name", "V1"), measure.vars="V2")

stack3c$mapping = ifelse(stack3c$V1=="Sequences did not map uniquely:", "non-uniquely mapped",
                         ifelse(stack3c$V1=="Number of alignments with a unique best hit from the different alignments:", 
                                "uniquely mapped",
                                ifelse(stack3c$V1=="Sequences with no alignments under any condition:", "not mapped", NA)))
stack3c = stack3c[!is.na(stack3c$mapping),]
stack3c$mapping = factor(stack3c$mapping, levels=c("not mapped", "non-uniquely mapped", "uniquely mapped")) 

## TODO: also include non-unique?
names(dat3) = sub("local.", "", sub("_bismark_(bt2|hisat2)_SE_report.txt","", basename(files3)))
stack3 = data.frame(do.call(rbind, strsplit(dat3, split="\t")))
stack3$name = rownames(stack3)
stack3$Mapping.efficiency = as.numeric(sub("%", "", stack3$X2, fixed=T))
options(scipen=10)

ref_palette = palette()[-1]; 
# darken grey, or change color of grid lines in plots.
# seems to be gray on some devices, gray62 on others?
ref_palette = sub("gray.*$","#9999AA",ref_palette)
cat("palette:", ref_palette,"\n")

if (require_all) {
  keep_samp = intersect(intersect(stack1$name, stack2$name), stack3$name)
  stack1 = stack1 %>% filter(name %in% keep_samp)
  stack2 = stack2 %>% filter(name %in% keep_samp)
  stack2m = stack2m %>% filter(name %in% keep_samp)
  stack3 = stack3 %>% filter(name %in% keep_samp)
  stack3c = stack3c %>% filter(name %in% keep_samp)
  
  if (F) {
    names_split = strsplit(keep_samp, "_")
    if (all(sapply(names_split, length)==2)) {
      n1 = sapply(strsplit(keep_samp, "_"),"[[",1)
      n2 = sapply(strsplit(keep_samp, "_"),"[[",2)
      name_perm = order(n2, n1)
      names_sorted = keep_samp[name_perm]
      stack1$name = factor(stack1$name, levels=names_sorted)
      stack2$name = factor(stack2$name, levels=names_sorted)
      stack3$name = factor(stack3$name, levels=names_sorted)
      stack2m$name = factor(stack2m$name, levels=names_sorted)
      stack3c$name = factor(stack3c$name, levels=names_sorted)
    }
  }
}


num_samp = length(unique(stack3c$name))
fig_width = ifelse(num_samp <= 80, 9, 1 + round(num_samp/10))
fig_height = 4 
asp_ratio = 1.5*fig_height/fig_width # for half-width plots 

pdf(output_file, width=fig_width, height=fig_height)

ggplot(stack3c, aes(x=name, y=value, fill=mapping)) + 
  geom_bar(stat="identity", position="fill" ) + theme_light() + theme(aspect.ratio=asp_ratio) + 
  scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=4)) + 
  labs(y="Percent of reads") + 
  xlab(NULL) + coord_cartesian(ylim=c(0,1), expand=F) + 
  ggtitle("Mapping statistics")

ggplot(stack2, aes(x=name, fill=ref_name, y=reads)) + 
  geom_bar(stat="identity") + theme_light() + theme(aspect.ratio=asp_ratio) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=4)) + 
  labs(y="number of reads") + 
  xlab(NULL) + scale_fill_manual(values=ref_palette, drop=F) +
  ggtitle("All mapped reads with no downsampling or filtering")

ggplot(stack2m, aes(x=name, fill=ref_name, y=reads*frac)) + 
  geom_bar(stat="identity") + theme_light() + theme(aspect.ratio=asp_ratio) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=4)) + 
  labs(y="number of reads") + 
  xlab(NULL) + scale_fill_manual(values=ref_palette, drop=F) +
  ggtitle("Approximate number of reads with filtering but without downsampling")

ggplot(stack1, aes(x=name, fill=ref_name, y=reads)) + 
  geom_bar(stat="identity") + theme_light() + theme(aspect.ratio=asp_ratio) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=4)) + 
  facet_wrap(~filt) + xlab(NULL) + scale_fill_manual(values=ref_palette, drop=F) +
  labs(y="number of reads") + 
  ggtitle("Downsampled to max of 10000 reads per ref_name, with or without additional filtering")


stack1x = stack1 %>% filter(reads >= min_reads)

gg1 = ggplot(stack1x, aes(x=name, col=ref_name, y=Q1)) + 
  geom_point(shape=2) + theme_light() + theme(aspect.ratio=asp_ratio) + 
  ylim(c(0,100)) + labs(y="Q1 of percent methylation") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=4)) + 
  xlab(NULL) + scale_color_manual(values=ref_palette, drop=F) + 
  ggtitle(paste0("Q1 for ref_name with at least ", min_reads, " reads in the sample"))

if (nrow(stack1x) > 0) gg1 = gg1 + facet_wrap(~filt)

gg1

stack1x = stack1 %>% filter(reads >= min_reads & filtered=="YES")

gg2 = ggplot(stack1x, aes(x=name, col=ref_name, y=Q2)) + 
  geom_linerange(aes(ymin=pct10, ymax=pct90), linewidth=0.4) + 
  geom_linerange(aes(ymin=Q1, ymax=Q3), linewidth=1) + 
  geom_point(shape=1, size=1.8) + theme_light() +  ylim(c(0,100)) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=4)) + 
  xlab(NULL) + scale_color_manual(values=ref_palette, drop=F) + 
  ggtitle(paste0("Quantiles for ref_seqs with at least ", min_reads, " filtered reads:")) + 
  labs(subtitle="circle is Q2; bar is Q1-Q3; whiskers are PCT10-PCT90", y="percent methylation")

if (nrow(stack1x) > 0) gg2 = gg2 + facet_wrap(~ref_name, ncol=2)

gg2

dev.off()

options(scipen=0)

