import glob
import os
import pandas as pd
import itertools
basedir = ""

genome="genomes/hg19.fa"
def get_group_files(sample,dm):
	g = dm[dm["sample"] == sample]["group"].values[0]
#	print(g)
	return dm[dm["group"] == g]["sample"].values


# IDR permutations
def idr_perm(sample, suffix=".narrowPeak",sep="___",prefix="."):
	
	perm_file = []

	for i, s in enumerate(sample[:-1]):
		for s2 in sample[i+1:]:
			perm_file.append("%s/%s___%s"%(prefix,s,s2))

	return perm_file
		
	
# General config

design_matrix = pd.read_csv("atacseq_design_ndhlovu_mar_17_202.csv")


design_matrix["group"] = design_matrix.apply(lambda x: ("%s.%s.%s"%(x.cell_type, x.tissue, x.subsetname)).lower() ,axis=1)
design_matrix["sample"] = design_matrix.bam.apply(lambda x:x.replace(".bam",""))


group_map = pd.Series(design_matrix["group"].values,index=design_matrix["bam"].values)
sample_map = pd.Series(map(lambda x:"%s"%x,design_matrix["bam"].values), index=design_matrix["group"].values)

groups = design_matrix["group"].drop_duplicates().values

print(sample_map)
print(groups)
print(design_matrix.columns)

rule all:
	input:
		genrich_peaks=expand("peaks/genrich.peaks/{group}.narrowPeak.gz", group=groups),
		macs2_peaks=design_matrix.apply(lambda x:"peaks/macs2/{group}/{sample}_peaks.narrowPeak".format(**x),1),
		macs2_bedpe_peaks=design_matrix.apply(lambda x:"peaks/macs2.bedpe/{group}/{sample}_peaks.narrowPeak".format(**x),1),
		macs2_bedpe_idr = expand("peaks/{macs2_method}_idr/{group}_idr_{optimal}.bed",macs2_method=["macs2","macs2.bedpe"],group=groups,
					 optimal=["relaxed","optimal"]),
		macs2_counts = expand("peak_counts/macs2/{macs2_method}/{optimal}/{sample}_counts.bed", macs2_method=["macs2","macs2.bedpe"], sample=design_matrix["sample"].values, optimal=["relaxed","optimal"]),
	        macs2_idr_merge = expand("peaks/{macs2_method}_idr/{group}_idr_{optimal}.bed",group=groups, macs2_method=["macs2","macs2.bedpe"], optimal=["relaxed","optimal"]),
		
		macs2_counts_dedup = expand("peak_counts.dedup/macs2/{macs2_method}/{optimal}/{sample}_counts.bed", macs2_method=["macs2","macs2.bedpe"], sample=design_matrix["sample"].values, optimal=["relaxed","optimal"]),
		macs2_genrich_overlap = expand("peaks/macs2.genrichoverlap/{group}.bed",group=groups),
		#rose=expand("rose/{group}/{group}_Plot_points.png",group=groups),
		#bigwig=expand("bigwig/{group}.bw",group=groups),
		nucleoatac=expand("nucleoatac/{group}/{group}.nfrpos.bed.gz",group=groups),
		nucleoatac_combinded=expand("nucleoatac/{group}/{group}.nucmap_combined.bed.gz",group=groups),
		nucleoatac_smooth_bw=expand("nucleoatac/{group}/{group}.nucleoatac_signal.smooth.bw", group=groups),
		tobias=expand("tobias.bindetect/genrich.peaks/norm/{group}/bindetect_figures.pdf",group=groups)

rule idr:
	input:
		lambda wildcards: expand("peaks/{macs2_method}/{group}/{sample}_peaks.narrowPeak", sample=design_matrix[design_matrix['group'] == wildcards.group]["sample"].values,
		group=wildcards.group, macs2_method=wildcards.macs2_method)
	output:
		relaxed="peaks/{macs2_method}_idr/{group}_idr_relaxed.bed",
		optimal="peaks/{macs2_method}_idr/{group}_idr_optimal.bed"


	
	params:
		idr=0.1
	threads: 7

	shell:
		"""
		FOLDER=peaks/{wildcards.macs2_method}/{wildcards.group}
		FILES="{input}"
		mkdir ${{FOLDER}}/idr

		echo {input} | tr " " "\n" | head -n-1 |
		
		while read file

			do 
				BFILES=$(echo ${{BFILES}} | tr " " "\n" | tail -n + 2)
				BASENAME1=$(basename ${{file}})
				echo ${{BFILES}} | tr " " "\n" |
					while read file2
					do
						BASENAME2=$(basename ${{file2}})
						echo "idr --idr {params.idr} --samples <( sort -k8,8nr ${{file}}) <( sort -k8,8nr ${{file2}}) --input-file-type narrowPeak --rank p.value 
							--output-file ${{FOLDER}}/idr/${{BASENAME1}}___${{BASENAME2}}.bed" 
					done
			done | parallel -j {threads}

			optimal_file=$(wc -l ${{FOLDER}}/idr/*.bed | head -n-1 | sort -n | rev | cut -f 1 -d\ | rev | tail -n 1)
			
			bedtools merge -i <(bedtools sort -i <(cat ${{FOLDER}}/idr/*.bed)) > {output.relaxed}
			cp ${{optimal_file}} {output.optimal}


		"""
				
		


#Genrich requires name sorted .....
rule namesort:
	input:
		IN="input.bam/{sample}"
	output:
		OUT="align/namesorted/{group}/{sample}.namesorted.bam"
		#OUT="align/namesorted/{group}/{sample}.namesorted.bam"
	threads: 8
	shell:
		"samtools sort -T /dev/shm -n {input.IN} -@ {threads} | samtools view -b -o {output.OUT}"
rule groupbam:
	input:
		IN="input.bam/{sample}"
	output:
		OUT="align/inital/{group}/{sample}"
	shell:
		"ln -s -f {input.IN} {output.OUT}"


rule genrichpeaks:
	input:
		#expand("align/namesorted/{group}/{sample}.namesorted.bam",sample=design_matrix["sample"],group=group_map[design_matrix["sample"].values])
		#rules.namesort.output.OUT
		#lambda wildcards: sample_map[wildcards.group].values
		#"align/namesorted/{group}/{sample}.namesorted.bam"
		lambda wildcards: expand("align/namesorted/{{group}}/{sample}.namesorted.bam",sample=sample_map[wildcards.group].values)
	params: 
		blacklist="assets/hg19-blacklist.bed",
		files=lambda wildcards: ",".join(expand("align/namesorted/{group}/{sample}.namesorted.bam",sample=sample_map[wildcards.group],group=wildcards.group))
	threads: 1
	resources:
	    mem=40 		
	output:
		narrowpeak="peaks/genrich.peaks/{group}.narrowPeak.gz",
		log="peaks/genrich.peaks/{group}.log"
	shell:
		"Genrich -t {params.files} -o {output.narrowpeak} -j -r -E {params.blacklist}  -q 0.05 -e chrM -f {output.log} -z"

rule merge_genrich:
	input:
		lambda wildcards: expand("peaks/genrich.peaks/{group}.narrowPeak.gz",group=groups)
	output:
		"peaks/merged.peaks/genrich.all.narrowPeak"
	shell:
		"bedtools merge -i <(bedtools sort -i <(zcat {input})) > {output}"


rule bam2bedpe:
	input:
		"input.bam/{sample}.bam"

	params:
		blacklist="assets/hg19-blacklist.bed"
	
	output:
		"align/bedpe/{group}/{sample}/{sample}.bedpe.gz"

	threads: 1
	shell:
		"bedtools intersect -v -a {input} -b {params.blacklist} | samtools sort -n | bedtools bamtobed -bedpe -i - 2> /dev/zero | gzip -c -1 > {output}"
 
rule macs2peaks_bedpe:
	input:
		"align/bedpe/{group}/{sample}/{sample}.bedpe.gz"
	params:
		sample=lambda wildcards:wildcards.sample,
		group=lambda wildcards:wildcards.group

	output:
		"peaks/macs2.bedpe/{group}/{sample}_peaks.narrowPeak"
	shell:
		"macs2 callpeak -t {input} -f BEDPE --keep-dup auto -p 0.01  -g hs -n {params.sample} --nomodel --outdir peaks/macs2.bedpe/{params.group}"

rule macs2peaks:
	input:
		#expand("align/initial/{group}/{sample}.bam",sample=design_matrix["sample"],group=group_map[design_matrix["sample"].values])
		#rules.groupbam.output
		"input.bam/{sample}.bam"
	params:
		sample=lambda wildcards: wildcards.sample,
		group=lambda wildcards: wildcards.group
	log: "log/macs2/{group}/{sample}.log"
	threads: 1
	output:
		#macs_narrowpeak="peaks/macs2/{group}/{sample}.narrowPeak",
		macs_narrowpeak="peaks/macs2/{group}/{sample}.narrowPeak"
	shell:
		"macs2 callpeak -t {input} -f BAMPE --keep-dup all -g hs -n {params.sample} --nomodel --outdir peaks/macs2/{wildcards.group}/  2> {log}"

rule macs2merged:
	input:
                "align/merged.dedup/{group}.bam"

	params:
		group=lambda wildcards: wildcards.group
	output:
		"peaks/macs2.merged_q0.1/{group}/{group}_peaks.narrowPeak"
	threads: 1
	shell:
		"macs2 callpeak -t {input} -f BAMPE --keep-dup all -g hs -n {params.group} --nomodel --outdir peaks/macs2.merged/{params.group}/"


rule genrichmacsoverlap:
	input:
		genrich="peaks/genrich.peaks/{group}.narrowPeak.gz",
		macs2="peaks/macs2.merged/{group}/{group}_peaks.narrowPeak"
	output:
		"peaks/macs2.genrichoverlap/{group}.bed"
	shell:
		"bedtools merge -i <(bedtools sort -i <(cat <(bedtools intersect -v -a {input.macs2} -b {input.genrich}) <(bedtools intersect -a {input.macs2} -b {input.genrich}))) > {output}"


rule tempbamdedup:
	input:
		"input.bam/{sample}"
	output:
		BAM="align/dedup/{sample}",
		METRICS="align/dedup/{sample}.txt"
		
	threads: 5
	shell:
		"java -XX:-UseTransparentHugePages -XX:+AlwaysPreTouch  -XX:ParallelGCThreads=2  -jar /apps/picard-2.17.11/picard.jar MarkDuplicates MAX_RECORDS_IN_RAM=350000 VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true ASSUME_SORTED=true M={output.METRICS} I={input} O={output.BAM} && samtools index {output.BAM}"

		
rule tempbammerge:
	input:
		lambda wildcards: expand("align/dedup/{sample}",  sample=sample_map[wildcards.group])
	output:
		"align/merged.dedup/{group}.bam"
	threads: 7
	shell:
		"samtools merge -1 -@ {threads} {output} {input}"

rule bigwigsignalmerge:
	input:
		"align/merged.dedup/{group}.bam"
	output:
		"bigwig/{group}.bw"
	threads: 4
	shell:
		"bamCoverage --bam {input} --effectiveGenomeSize 2864785220 -o {output} -p {threads} -bs 10 --blackListFileName assets/hg19-blacklist-merged.bed --normalizeUsing BPM"
	
rule rosemacs2genrich:
	input:
		BAM="align/merged.dedup/{group}.bam",
                PEAK="peaks/macs2.genrichoverlap/{group}.bed"
	output:
		BEDTOGFF="rose/{group}/{group}.gff",
		ROSE="rose/{group}/{group}_Plot_points.png"
	threads: 1
	shell:
                """
                cat {input.PEAK} | awk '{{print $1"\t"$4"\t"$4"\t"$2"\t"$3"\t"$7"\t"$6"\t"$4"\t"$4}}' | sort -u -k2,2 > {output.BEDTOGFF}
                python2 ROSE_main.py -g hg19 -i {output.BEDTOGFF} -r {input.BAM} -o rose/{wildcards.group} -n {wildcards.group}
                """

rule nucleoatac_slop:
	input:
		BED="peaks/macs2.genrichoverlap/{group}.bed"
	output:
		BED="peaks/macs2.genrichoverlap.slop/{group}.bed"
	shell:
		"bedtools slop -i {input.BED} -g assets/chrom.sizes.hg19 -b 500 > {output.BED}"


rule nucleoatac_occ:
	input:
		BED="peaks/macs2.genrichoverlap.slop/{group}.bed",
		BAM="align/merged.dedup/{group}.bam"
	params: 
		BASENAME="{wildcards.group}"
	output:
		"nucleoatac/{group}/{group}.occ.bedgraph.gz", "nucleoatac/{group}/{group}.occ.lower_bound.bedgraph.gz", 
		"nucleoatac/{group}/{group}.occ.upper_bound.bedgraph.gz","nucleoatac/{group}/{group}.fragmentsizes.txt",
		"nucleoatac/{group}/{group}.occpeaks.bed.gz", "nucleoatac/{group}/{group}.nuc_dist.txt"
	threads: 14
	shell:
		"nucleoatac occ --fasta assets/genome_hg19.fa --bed {input.BED} --bam {input.BAM} --out nucleoatac/{wildcards.group}/{wildcards.group} --cores {threads}"

rule nucleoatac_vprocess:
	input:
		"nucleoatac/{group}/{group}.nuc_dist.txt"
	output:
		"nucleoatac/{group}/{group}.VMat"
	shell:
		"nucleoatac vprocess --sizes {input} --out nucleoatac/{wildcards.group}/{wildcards.group}"

rule nucleoatac_nuc:
	input:
		#BED="peaks/macs2.genrichoverlap.slop/{group}.bed",
		BED="peaks/genrich.peaks/allpeaks.nofilter.merged.slop.bed",
		BAM="align/merged.dedup/{group}.bam",
		VMAT="nucleoatac/{group}/{group}.VMat",
		FRAGMENTSIZES="nucleoatac/{group}/{group}.fragmentsizes.txt"
	output:
		"nucleoatac/{group}/{group}.nucpos.bed.gz",
		"nucleoatac/{group}/{group}.nucpos.redundant.bed.gz",
		"nucleoatac/{group}/{group}.nucleoatac_signal.bedgraph.gz",
		"nucleoatac/{group}/{group}.nucleoatac_signal.smooth.bedgraph.gz",
	log:    "log/{group}.nuc.log"
	threads: 28
	shell:
		"OMP_NUM_THREADS=1 nucleoatac nuc --fasta assets/genome_hg19.fa --bed {input.BED} --vmat {input.VMAT} --bam {input.BAM} --out nucleoatac/{wildcards.group}/{wildcards.group} --sizes {input.FRAGMENTSIZES}  --cores {threads} > {log}"


rule nucleoatac_merge:
	input:
		OCC="nucleoatac/{group}/{group}.occpeaks.bed.gz",
		NUCPOS="nucleoatac/{group}/{group}.nucpos.bed.gz"
	output:
		"nucleoatac/{group}/{group}.nucmap_combined.bed.gz"
	shell:
		"nucleoatac merge --occpeaks {input.OCC} --nucpos {input.NUCPOS} --out nucleoatac/{wildcards.group}/{wildcards.group}"

rule nucleoatac_nfr:
	input:
		BED="peaks/macs2.genrichoverlap.slop/{group}.bed",
		OCC="nucleoatac/{group}/{group}.occ.bedgraph.gz",
		NUCPOS="nucleoatac/{group}/{group}.nucpos.bed.gz",
		BAM="align/merged.dedup/{group}.bam"

	output:
		"nucleoatac/{group}/{group}.nfrpos.bed.gz"
	shell:
		"nucleoatac nfr --bam {input.BAM} --fasta assets/genome_hg19.fa --bed {input.BED} --occ_track {input.OCC} --calls {input.NUCPOS} --out nucleoatac/{wildcards.group}/{wildcards.group}"

rule nucleoatac_bedgraph_intermediate:
	input:
		SMOOTH_BEDGRAPH="nucleoatac/{group}/{group}.nucleoatac_signal.smooth.bedgraph.gz"
	output:
		temp("nucleoatac/{group}/{group}.nucleoatac_signal.smooth.bedgraph")
	threads: 4
	shell:
		"gunzip -c {input.SMOOTH_BEDGRAPH} > {output}"

rule nucleoatac_sig_bw:
	input:
		SMOOTH_BEDGRAPH="nucleoatac/{group}/{group}.nucleoatac_signal.smooth.bedgraph"
	output:
		"nucleoatac/{group}/{group}.nucleoatac_signal.smooth.bw"
	threads: 4
	shell:
		"bedGraphToBigWig {input.SMOOTH_BEDGRAPH} assets/chrom.sizes.hg19 {output}"


rule tobias_atacorrect_single:
	input:
		BAM="align/dedup/{sample_nobam}.bam",
		PEAKS="peaks/genrich.peaks/{group}.narrowPeak.gz"
	output:
		SIGNAL="tobias.single/atacorrect/{sample_nobam}_corrected.bw"
	
	params: 
		blacklist="assets/hg19-blacklist.bed",
		name="{sample_nobam}"
	threads: 14
	shell:
		"TOBIAS ATACorrect --split 1000 --blacklist {params.blacklist} --cores {threads} --prefix {params.name} --outdir tobias.single/atacorrect --genome {genome} --peaks {input.PEAKS} --bam {input.BAM}"

rule tobias_atacorrect:
	input:
		BAM="align/merged.dedup/{group}.bam",
		PEAKS="peaks/genrich.peaks/{group}.narrowPeak.gz"
	output:
		ATACCORRECT_NORM_REP="tobias/{group}/atacorrect.norm/{group}_corrected.bw",
	params:
		ATACORRECT_NORM="tobias/{group}/atacorrect.norm",
		blacklist="assets/hg19-blacklist.bed",
		name="{group}"
	
	threads: 14
	shell:
		"TOBIAS ATACorrect --blacklist {params.blacklist} --cores {threads} --bam {input.BAM} "
		" --prefix {params.name} --genome {genome} --peaks <( zcat {input.PEAKS} ) --outdir {params.ATACORRECT_NORM} "

	
rule tobias_fp_score:
	input:
		SIGNAL_NORM="tobias/{group}/atacorrect.norm/{group}_corrected.bw",
		PEAKS="peaks/genrich.peaks/{group}.narrowPeak.gz"
	output:
		FP_NORM="tobias/{group}/footprints.norm/{group}.bw",

	threads: 14
	shell:
		"TOBIAS FootprintScores --cores {threads} -r <( zcat {input.PEAKS} ) -s {input.SIGNAL_NORM} -o {output.FP_NORM}"




rule tobias_bindetect:
	input:
		FP_NORM=expand("tobias/{group}/footprints.norm/{group}.bw",group=groups),
		#PEAKS=expand("peaks/genrich.peaks/{group}.narrowPeak.gz",group=groups)
		PEAKS="peaks/merged.peaks/genrich.all.narrowPeak"
	output:
		FP_COMP_NORM="tobias.bindetect/genrich.peaks/norm/{group}/bindetect_figures.pdf",
	params:
		motifs="tfmotifs/JASPAR2020.txt",
	threads: 20
	shell:
		"TOBIAS BINDetect  --cores {threads} --skip-excel "
		"--outdir tobias.bindetect/genrich.peaks/norm --signals {input.FP_NORM} "
		"--motifs {params.motifs} --genome {genome} --peaks {input.PEAKS} "

		


rule merge_idr_peaks:
	input:
		expand("peaks/{{macs2_method}}_idr/{group}_idr_{{optimal}}.bed",group=groups)
	
	output:
		"peaks/merged_peaks/macs2/{macs2_method}_idr_{optimal}.bed"
	

	shell:
		"bedtools merge -i <(bedtools sort -g assets/chrom.sizes.hg19 -i <(cat {input})) > {output}"

rule merge_macs2_genrich:
	input:
		expand("peaks/macs2.genrichoverlap/{group}.bed",group=groups)

	output:
		expand("peaks/merged/macs2.genrichoverlap.bed")

	shell:
		"bedtools merge -i <(bedtools sort -i <(cat {input})) > {output}"
	
rule peakcutcounts_macs:
	input:
		PEAKFILE="peaks/merged_peaks/macs2/{macs2_method}_idr_{optimal}.bed",
		BAM="input.bam/{sample}.bam"

	
	output:
		"peak_counts/macs2/{macs2_method}/{optimal}/{sample}_counts.bed"

	threads: 2	
	shell: 
		"bedtools multicov -bams {input.BAM} -bed {input.PEAKFILE} > {output}"

		

rule peakcutcounts_macs_dedup:
	input:
		PEAKFILE="peaks/merged_peaks/macs2/{macs2_method}_idr_{optimal}.bed",
		BAM="align/dedup/{sample}.bam"

	
	output:
		"peak_counts.dedup/macs2/{macs2_method}/{optimal}/{sample}_counts.bed"

	threads: 2	
	shell: 
		"bedtools multicov -bams {input.BAM} -bed {input.PEAKFILE} > {output}"
