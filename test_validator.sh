#!/bin/sh
python scoring_harness/TESLA_validation.py --patientId 1 examples/1/TESLA_OUT_2.csv examples/1/TESLA_OUT_3.csv examples/1/TESLA_OUT_4.csv examples/1/TESLA_VCF.vcf
python scoring_harness/TESLA_validation.py --patientId 1 examples/1/TESLA_OUT_2.csv examples/1/TESLA_OUT_3.csv examples/1/TESLA_OUT_4.csv examples/1/TESLA_VCF.vcf examples/1/1_EXOME_N.bam examples/1/1_EXOME_T.bam examples/1/1_RNA_T.bam --validatingBAM --patientId 1  
python scoring_harness/TESLA_validation.py --patientId 1 examples/1/TESLA_OUT_2.csv examples/1/TESLA_OUT_3.csv examples/1/TESLA_OUT_4.csv examples/1/TESLA_MAF.maf
