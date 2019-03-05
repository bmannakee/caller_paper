---
output:
    pdf_document:
        fig_caption: yes

bibliography: caller_paper.bib
---
# Targets
- Bioinformatics
- Plus Comp. Bio
- BMC bioinformatics
- Nature Scientific Reports.
- Genome Biology (published MuSE [@Fan2016])
# Introduction

Cancer is an evolutionary process, and understanding initiation, progression, and metastasis will require applications of evolutionary theory.
One of the major tools in the evolutionary theory toolbox is the allele frequency spectrum.
This allele frequency spectrum is constructed from 

Tumor heterogeneity has been associated with prognosis (1-4 in chuang paper) and the evolutionary trajectory helps identify the number of tumor subclones and their selective advantage.

The variant allele frequency spectrum that is currently used most often in cancer is truncated at a level above 5-10% because of difficulties in identifying low frequency variants.

- There are two main tracks in variant calling.
    - Heuristic filters
    - Statistical models of sequencing error
- We focus here on a model of mutation probability, including but not limited to sequencing error.
- Many types of callers, all assume there is no biological preference for mutation at a given site. Any site specific estimates are site specific sequencing/alignment error models[@Xu2018]. 
- Mutect2, FreeBayes and others are haplotype based callers
- Callers with site specific variant probabilities generate them either from other samples or through deep sequencing (deepSNV,EBCall,LoLoPicker). They are essentially generating a site specific sequencing error model, not a site specific probability of mutation
- Need to think about how the method applies to UMI (barcode) based sequencing, which are mostly deep targeted
MuSE is continuous time markov evolutionary model, still assuming no biological difference in site specific mutation probability[@Fan2016]
- very little attention to the statistical model, either in competition or development
- there is useful biology.....
    - [@Temko2018] links between mutational processes and driver mutations
    - [@VandenEynden2017] mutational signature critical for estimating selection
    - [@Kandoth2013,@Alexandrov2013a] Underlying mutational processes generate tumor and tumor type specific mutation signatures
- Rather than using a constant probability for mutation, as other variant callers do, we convert that to an average or expected mutation probability, and compute the probability conditional on context and genome composition
- Poisson models make similar assumptions about the probability of an allele at a site. (Illumina technical note https://www.illumina.com/Documents/products/technotes/technote_somatic_variant_caller.pdf).
- we simulate neutral tumor evolution, and assign vafs using a Beta(1,6) distribution
    - if M(f) is proportional to 1/f, then an exponential distribution is implied [@tarabichi2017,@Williams2017](and the answering note by De, which also has a strong argument about why we need lower frequencies to do evolutionary inference). We choose a beta distribution to draw vafs and tuned to achieve a slightly fatter distribution in the 2-5% range in which we are most interested.
- Need a list of why evolutionary inference on tumors is important. Resistance, virulence(heterogeneity), biology (mutation rate/signature/micro-environment).

# Results

## Sensitivity and specificity in simulated data
We can look at interactions now that we have settled on an experimental design.

1. Brief description of simulations, see methods
2. Words about the figure
    - What is the linear regime in the Mutect ROC curves about?
    - Is it related to the uniform prior, and does it give a good explanation of the performance difference?

Experiment 2 is a 100X whole genome with ~29000 spiked variants, most of which are under 2% because of the way the simulation works.


Experiment 10 is a 100X whole genome with the same variants as Experiment 2, but with a uniform vaf distribution. Prior method is still better, but in the uniform scenario there are only a small fraction of the total variants that are challenging to call.

Experiment 9 is a whole exome that was supposed to have 1,7,11, but instead has a random set of signatures due to a bug

## Sensitivity in real data
We examined two validation datasets from real tumors. An acute myeloid leukemia whole genome was sequenced to average coverage of 365X, and over 200,000 mutations validated by deep sequencing, generating a set of "platinum" consensus calls for the tumor. In addition to the full dataset we also called mutations on two downsample datasets, one retaining 50% of the original reads and one retaining 25%. ROC curves were generated using the "platinum" calls as cases, and sites where validation sequencing depth was greater than 100X and no variant reads were found as controls. Both algorithms perform similarly and nowhere along the curve is the {what is the name of this thing} method below raw mutect calls. The {method} calls a higher fraction of platinum calls at every odds threshold, and is especially effective at the common threshold of 2:1 odds in favor of the mutation.



***Going to need a table of AUROCs in the supplement for this***

## Effect of odds threshold
This has very little effect, even in an exome, as the figure inserted shows.
1. As threshold goes to infinity you get mutect.
2. As threshold goes to zero you should also get mutect.
3. Observe very little difference in the middle

## Effect of number of mutations
We will have this from the difference between exome and wgs on the same vaf distribution and signature.
This is likely to have some signature dependence.
1. How to approach this?
   - At what point does the empirical make more sense than the dirichlet.
   - I think never, they will converge
   - What is the stopping point with a low number of high confidence mutations
   - Implementation of the dirichlet should let us create an estimation of total error between the final empirical at a given threshold and the dirichlet at every point in the process. Maybe a plot of this?

## Effect of variant allele frequency distribution
1. TCGA data for different distributions.
   - Different cancer types?
   - Hypermutators vs. not?
   - This should only be related to the number of mutations that are confident and contribute to the prior
   - If that is the case, is there an analytical way to better describe this?
   THE ONLY EFFECT IS ON THE ROC. EASIER DISTRIBUTIONS SHRINK THE EFFECT BECAUSE SO FEW ARE NEAR THE CRITICAL POINT



# Methods

## 
100X whole genome and 500X whole exome for each of three signatures

1,7,11 UV (Very concentrated at C>T)
1,4,5 Tobacco (Slight concentation at C>A and C>T)
1,3,5 Breast (diffuse)

All vafs will be from the beta(1,6) which is a fat exponential


# Figures

![roc curve figure experiment 9](figures/roc_and_called_curves.png)

![ figure experiment 9](figures/WES_thresholds_exp9.png)

![roc curve figure experiment 9](figures/beta_1_6.png)
Figure 1 - aml31 no downsample roc

![Figure 1 - aml31 no downsample roc](figures/aml31_no_downsample_roc.png)

Figure 2 - aml31 no downsample fraction called

![Figure 2 - aml31 no downsample fraction called](figures/aml31_no_downsample_fraction_called.png)

Figure 2a - aml31 no downsample vaf

![Figure 2 - aml31 no downsample vaf](figures/aml31_no_downsample_vaf.png)

Figure 3 - aml31 50 percent downsample roc

![Figure 3 - aml31 50 percent downsample roc](figures/aml31_downsampled_50_percent_roc.png)

Figure 4 - aml31 50 percent downsample fraction called

![Figure 4 - aml31 50 percent downsample fraction called](figures/aml31_downsampled_50_percent_fraction_called.png)

Figure 4a - aml31 50 percent downsample vaf

![Figure 4 - aml31 50 percent downsample vaf](figures/aml31_downsampled_50_percent_vaf.png)

Figure 5 - aml31 25 percent downsample roc

![Figure 5 - aml31 25 percent downsample roc](figures/aml31_downsampled_25_percent_roc.png)

Figure 6 - aml31 25 percent downsample fraction called

![Figure 6 - aml31 25 percent downsample fraction called](figures/aml31_downsampled_25_percent_fraction_called.png)

Figure 6a - aml31 25 percent downsample vaf

![Figure 6 - aml31 25 percent downsample vaf](figures/aml31_downsampled_25_percent_vaf.png)

Figure 7 - cell paper roc

![Figure 7 - cell paper roc](figures/cell_paper_roc.png)

Figure 8 - cell paper fraction called

![Figure 8 - cell paper fraction called](figures/cell_paper_fraction_called.png)

Figure 8a - cell paper vaf

![Figure 8 - cell paper vaf](figures/cell_paper_vaf.png)


<!-- Figure 9 - experiment 2 roc

![Figure 7 - experiment 2 roc](figures/experiment2_roc.png)

Figure 10 - experiment 2 fraction called

![Figure 8 - aml31 experiment 2 fraction called](figures/experiment2_fraction_called.png)

Figure 10a - experiment 2 vaf

![Figure 8 - aml31 experiment 2 vaf](figures/experiment2_vaf.png)


Figure 11 - experiment 10 (uniform vaf) roc

![Figure 11 - experiment 10 (uniform vaf) roc](figures/experiment10_roc.png)

Figure 12 - experiment 10 (uniform vaf) fraction called

![Figure 12 - experiment 10 (uniform vaf) fraction called](figures/experiment10_fraction_called.png)

Figure 12a - experiment 10 (uniform vaf) vaf

![Figure 12a - experiment 10 (uniform vaf) vaf](figures/experiment10_vaf.png) -->

# References