library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) proj_dir = args[1] else proj_dir = "."
if (!dir.exists(proj_dir)) stop(proj_dir, "not found")

# RScript --vanilla ~/Dropbox/bss_summarize_run.R  ~/bisulfite_2020/JFBSS_0004_2021_03_24
# RScript --vanilla ~/Dropbox/bss_summarize_run.R .
# proj_dir = "~/bisulfite_2020/JFBSS_0005_2021_03_26"

## in current working dir
output_file = "run_summary_plot.pdf"
if (length(args) > 1) output_file = args[2] 

files1 = dir(path=file.path(proj_dir, "r_out"), pattern="*.txt", full=T)
dat1 = lapply(files1, read.table, header=T)
names(dat1) = sub(".txt","", basename(files1))
stack1 = do.call(rbind, dat1)
stack1$name = sub("report_", "", sub("\\.\\d+$", "", rownames(stack1))) 
stack1 = stack1 %>% mutate(filt=factor(ifelse(filtered=="YES", "filtered", "unfiltered"),
                                       levels=c("unfiltered", "filtered"))) 

files2 = dir(path=file.path(proj_dir, "bismark_bam"), pattern="sorted.bam.idxstats", full=T)
dat2 = lapply(files2,  read.table, header=F, sep="\t")

names(dat2) = sub("local.", "", sub("_bismark_(bt2|hisat2).sorted.bam.idxstats","", basename(files2)))
stack2 = do.call(rbind, dat2) %>% filter(V1 != "*")
stack2$name =sub("report_", "", sub("\\.\\d+$", "", rownames(stack2))) 
colnames(stack2)[1:3]=c("ref_name", "length", "reads")

stack1$ref_name = factor(stack1$ref_name)
stack2$ref_name = factor(stack2$ref_name)
stack1$reads[which(is.na(stack1$reads))] = 0


stack1w = stack1 %>% select(name, ref_name, filtered, reads) %>% 
  pivot_wider(names_from=filtered, values_from=reads) %>% mutate(frac=YES/pmax(NO,1))
  
stack2m = left_join(stack2, stack1w %>% select(name, ref_name, frac))

options(scipen=10)

# ref_palette = RColorBrewer::brewer.pal(n=length(unique(stack2$ref_name)), name="Set2")
ref_palette = palette()[-1]; 
# darken grey, or change color of grid lines in plots
# grey on some devices, gray 62 on others!
ref_palette = sub("gray.*$","#9999AA",ref_palette)
#plot(1:7, col=ref_palette, pch=16, cex=2)
cat("palette:", ref_palette,"\n")

## basename removes trailing slash
pdf(output_file, width=9, height=4)

ggplot(stack2, aes(x=name, fill=ref_name, y=reads)) + 
  geom_bar(stat="identity") + theme_light() + theme(aspect.ratio=0.618) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=6)) + 
  labs(y="number of reads") + 
  xlab(NULL) + scale_fill_manual(values=ref_palette, drop=F) +
  ggtitle("All mapped reads with no downsampling or filtering")

ggplot(stack2m, aes(x=name, fill=ref_name, y=reads*frac)) + 
  geom_bar(stat="identity") + theme_light() + theme(aspect.ratio=0.618) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=6)) + 
  labs(y="number of reads") + 
  xlab(NULL) + scale_fill_manual(values=ref_palette, drop=F) +
  ggtitle("Approximate number of reads with filtering but without downsampling")


ggplot(stack1, aes(x=name, fill=ref_name, y=reads)) + 
  geom_bar(stat="identity") + theme_light() + theme(aspect.ratio=0.618) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=6)) + 
  facet_wrap(~filt) + xlab(NULL) + scale_fill_manual(values=ref_palette, drop=F) +
  labs(y="number of reads") + 
  ggtitle("Downsampled to max of 10000 reads per ref_name, with or without additional filtering")

ggplot(stack1 %>% filter(reads >= 100), aes(x=name, col=ref_name, y=Q1)) + 
  geom_point(shape=2) + theme_light() + theme(aspect.ratio=0.618) + 
  ylim(c(0,100)) + labs(y="Q1 of percent methylation") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=6)) + 
  facet_wrap(~filt) + xlab(NULL) + scale_color_manual(values=ref_palette, drop=F) + 
  ggtitle("Q1 for ref_name with at least 100 reads in the sample")

ggplot(stack1 %>% filter(reads >= 100 & filtered=="YES"), aes(x=name, col=ref_name, y=Q2)) + 
  geom_linerange(aes(ymin=pct10, ymax=pct90), size=0.4) + 
  geom_linerange(aes(ymin=Q1, ymax=Q3), size=1) + 
  geom_point(shape=1, size=1.8) + theme_light() +  ylim(c(0,100)) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size=6)) + 
  facet_wrap(~ref_name, ncol=2) + xlab(NULL) + scale_color_manual(values=ref_palette, drop=F) + 
  ggtitle("Quantiles for ref_seqs with at least 100 filtered reads:") + 
  labs(subtitle="circle is Q2; bar is Q1-Q3; whiskers are PCT10-PCT90", y="percent methylation")

dev.off()

options(scipen=0)

