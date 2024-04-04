# rtPRS-CS
**rtPRS-CS** is a python-based command line tool that performs real-time online updating of polygenic risk score (PRS) weights in a target dataset one sample at-a-time. For each target sample, rtPRS-CS first calculates a PRS using the most recent set of SNP weights, then use stochastic gradient descent to update the SNP weights given the genotype, phenotype, and covariates of the target sample.

## Getting Started

- Clone this repository using the following git command:

    `git clone https://github.com/getian107/rtPRS.git`

    Alternatively, download the source files from the github website (`https://github.com/getian107/rtPRS`)
    
- Download the LD reference panels and extract files:

    LD reference panels constructed using the 1000 Genomes Project phase 3 samples:
    
     [AFR reference](https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=0 "AFR reference") (~4.44G);
     `tar -zxvf ldblk_1kg_afr.tar.gz`
     
     [AMR reference](https://www.dropbox.com/s/uv5ydr4uv528lca/ldblk_1kg_amr.tar.gz?dl=0 "AMR reference") (~3.84G);
     `tar -zxvf ldblk_1kg_amr.tar.gz`
        
     [EAS reference](https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=0 "EAS reference") (~4.33G);
     `tar -zxvf ldblk_1kg_eas.tar.gz`
        
     [EUR reference](https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0 "EUR reference") (~4.56G);
     `tar -zxvf ldblk_1kg_eur.tar.gz`
     
     [SAS reference](https://www.dropbox.com/s/hsm0qwgyixswdcv/ldblk_1kg_sas.tar.gz?dl=0 "SAS reference") (~5.60G);
     `tar -zxvf ldblk_1kg_sas.tar.gz`
    
    LD reference panels constructed using the UK Biobank data ([Notes](https://www.dropbox.com/s/y3hsc15kwjxwjtd/UKBB_ref.txt?dl=0 "Notes")):
    
     [AFR reference](https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz?dl=0 "AFR reference") (~4.93G);
     `tar -zxvf ldblk_ukbb_afr.tar.gz`
     
     [AMR reference](https://www.dropbox.com/s/y7ruj364buprkl6/ldblk_ukbb_amr.tar.gz?dl=0 "AMR reference") (~4.10G);
     `tar -zxvf ldblk_ukbb_amr.tar.gz`
    
     [EAS reference](https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=0 "EAS reference") (~5.80G);
     `tar -zxvf ldblk_ukbb_eas.tar.gz`
    
     [EUR reference](https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=0 "EUR reference") (~6.25G);
     `tar -zxvf ldblk_ukbb_eur.tar.gz`
    
     [SAS reference](https://www.dropbox.com/s/nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz?dl=0 "SAS reference") (~7.37G);
     `tar -zxvf ldblk_ukbb_sas.tar.gz`
    
    Note that these files are identical to the reference panels used in **PRS-CS**.  
    Therefore, there is no need to download again if you are already using **PRS-CS**.
    
    For regions that don't have access to Dropbox, reference panels can be downloaded from the
    [alternative download site](https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference).

- Download the SNP information file and put it in the same folder containing the reference panels:

    1000 Genomes reference: [SNP info](https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3?dl=0 "SNP info") (~106M)
    
    UK Biobank reference: [SNP info](https://www.dropbox.com/s/oyn5trwtuei27qj/snpinfo_mult_ukbb_hm3?dl=0 "SNP info") (~108M)
    
- PRScsx requires Python packages **scipy** (https://www.scipy.org/) and **h5py** (https://www.h5py.org/) installed.
 
- Once Python and its dependencies have been installed, running

    `./rtPRScs.py --help` or `./rtPRScs.py -h`

    will print a list of command-line options.
    

## Using rtPRS-CS

`
python rtPRScs.py --ref_dir=PATH_TO_REFERENCE --n_gwas=GWAS_SAMPLE_SIZE --pst_eff=POSTERIOR_EFFECTS --psi_est=PSI_ESTIMATES
                  --vld_prefix=VALIDATION_DATASET_PREFIX --vld_phn_cov=VALIDATION_DATASET_PHENOTYPE_COVARIATES --vld_frq_prefix=VALIDATION_FRQ_FILE_PREFIX
                  --tst_prefix=TESTING_DATASET_PREFIX --tst_phn_cov=TESTING_DATASET_PHENOTYPE_COVARIATES --out_file=OUTPUT_FILENAME
                  [--rate=LEARNING_RATE --imp=IMPLICIT_SGD --order=ORDER_OF_ALGORITHM --chrom=CHROM] 
`
 - PATH_TO_REFERENCE (required): Full path to the directory that contains the SNP information file and LD reference panels. If the 1000 Genomes reference is used, the folder would contain the SNP information file `snpinfo_mult_1kg_hm3` and one or more of the LD reference files: `ldblk_1kg_afr`, `ldblk_1kg_amr`, `ldblk_1kg_eas`, `ldblk_1kg_eur`, `ldblk_1kg_sas`; if the UK Biobank reference is used, the folder would contain the SNP information file `snpinfo_mult_ukbb_hm3` and one or more of the LD reference files: `ldblk_ukbb_afr`, `ldblk_ukbb_amr`, `ldblk_ukbb_eas`, `ldblk_ukbb_eur`, `ldblk_ukbb_sas`.

 - GWAS_SAMPLE_SIZE (required): Sample sizes of the GWAS, in the same order of the GWAS summary statistics files, separated by comma.
   
 - POSTERIOR_EFFECTS (required): Starting values for SNP effect sizes, output from running **PRS-CS** on a set of baseline GWAS summary statistics.
   
 - PSI_ESTIMATES (required): SNP local shrinkage parameters output from running **PRS-CS** (with the option --write_psi) on a set of baseline GWAS summary statistics.
   
 - VALIDATION_DATASET_PREFIX (required): Full path and the prefix of the bim file for the validation dataset. This file is used to provide a list of SNPs that are available in the target dataset.

 - VALIDATION_DATASET_PHENOTYPE_COVARIATES (required): Full path to the space- or tab-delimited text file containing the validation sample phenotype and covariates values. The file shoud following the PLINK "pheno" file format. i.e. column 1 contains the family ID, column 2 contains the individual ID, column 3 contains the phenotype, and any remaining columns contain covariates.

 - VALIDATION_FRQ_FILE_PREFIX (required): Full path and the prefix of the .frq file output from running PLINK's "--freq" command on the validation dataset samples.

 - TESTING_DATASET_PREFIX (required): Full path and the prefix of the bim file for the target (test) validation dataset.

 - TESTING_DATASET_PHENOTYPE_COVARIATES (required): Full path to the space- or tab-delimited text file containing the test (target) sample phenotype and covariates values. The file shoud following the PLINK "pheno" file format. i.e. column 1 contains the family ID, column 2 contains the individual ID, column 3 contains the phenotype, and any remaining columns contain covariates.

 - OUTPUT_FILENAME (required): Output filename prefix of the calculated PRS and final posterior effect size estimates.

 - LEARNING_RATE (optional): The learning rate parameter, controlling how each incoming sample is weighted when updating SNP effect sizes. Default is 1.

 - IMPLICIT_SGD (optional): A boolean value indicating whether implicit (True) or explicit (False) stochastic gradient descent should be employed when updating SNP effect sizes. Default is True.

 - ORDER_OF_ALGORITHM (optional): A string indicating whether the 1st or 2nd order derivative be used in SGD optimization. Default is 2nd.

 - CHROM (optional): A numeric value or range of values indicating which chromosome(s) should be run. Deafult is 1-23.

## Output

rtPRS-CS writes two output files to the user-specified directory. The first is a .txt file containing the standardized PRS for each individual. The second is a .pst_eff file containing the vector of SNP weights after incorporating all target (training) samples.

## Support

Please direct any problems or questions to Tian Ge (tge1@mgh.harvard.edu).

