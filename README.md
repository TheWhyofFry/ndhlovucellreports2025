# ndhlovucellreports2025
Description of a subset analysis performed for the manuscript `Spatial regulation of CD8+ T cells at the HLA-E-NKG2A axis drives HIV persistence in Lymph Node B Cell Follicles (Cell Reports, 2025)`


# ATAC-Seq


Raw reads were passed through the pipeline in `pipelines/atacseq`. Simply, reads were trimmed and mapped to the `hg19` reference genome, afterwhich the peaks were called using Genrich and MACS2 and finally an intersected (overlapping peaks) and union (inclusive of all peaks in all samples) were generated. These counds can be found in the R data structure (RDS) file, `assets/dds.atac.counts.trim.analyzed_cell_reports.rds`. A corresponding table with annotations can be found in `assets/atac.counts.peak.anno.df.rds`.  While the analysis was performed using a full design matrix that includes additional samples, these are not part of the publication and have been omitted. There was very little patient overlap between the lymph node and PBMC samples, and thus the quantitative peak analysis was only performed using the tissues (LN and Blood) as covariates. This, together with cell subtype variation, has an inflationary impact on variance, and it is for this reason we did not correct for multiple testing. ATAC-Seq plots were generated with custom functions (`plot.genetrack`) found in `scripts/plots.r`. BigWig files are available upon request.




# RNA-Seq


RNA-Seq data was passed through the RNA-Seq Snakemake pipeline (in `pipelines/rnaseq`). From the pipeline, the bootstrapped and imported with the R package _sleuth_. The design, unlike the ATAC-Seq analysis, were paired with the PID files. The _sleuth_ model is saved under `assets/sleuth_model_ln_vs_pb.rds`. A simple contrast between LN and PBMC samples, while accounting for inter-participant variation can be found under `assets/sleuth_test_ln_vs_pb.rds`. Since the amount of samples are small, for visualization purposes we transformed the data with VST and included 5 surrogate variables. This was merely to aid in visualization and was not included in the official analysis, where statistics were pulled from the orignal _sleuth_ model. The _VST-adjusted_ samples is saved under `assets/dds.lymphnode_study.vst.rb_cell_reports_2025.rds`. The corresponding function in the `scripts/atacplot.R` file is `boxplot.gene`.




