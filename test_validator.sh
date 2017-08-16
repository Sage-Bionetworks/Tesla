#!/bin/sh
echo "Test VCF validation"
python scoring_harness/TESLA_validation.py --patientId 10 examples/10/TESLA_OUT_1.csv examples/10/TESLA_OUT_2.csv examples/10/TESLA_OUT_3.csv examples/10/TESLA_OUT_4.csv examples/10/TESLA_OUT_5.csv examples/10/TESLA_VCF.vcf
echo "Test TESLA_ranking_method validation"
python scoring_harness/TESLA_validation.py --patientId 10 examples/10/TESLA_OUT_1.csv examples/10/TESLA_OUT_2.csv examples/10/TESLA_OUT_3.csv examples/10/TESLA_OUT_4.csv examples/10/TESLA_OUT_5.csv examples/10/TESLA_VCF.vcf examples/10/TESLA_ranking_method.txt 
echo "Test bam name validation"
python scoring_harness/TESLA_validation.py --patientId 10 examples/10/TESLA_OUT_1.csv examples/10/TESLA_OUT_2.csv examples/10/TESLA_OUT_3.csv examples/10/TESLA_OUT_4.csv examples/10/TESLA_OUT_5.csv examples/10/TESLA_VCF.vcf examples/10/10_EXOME_N.bam examples/10/10_EXOME_T.bam examples/10/10_RNA_T.bam --validatingBAM 
echo "Test MAF instead of VCF"
python scoring_harness/TESLA_validation.py --patientId 10 examples/10/TESLA_OUT_1.csv examples/10/TESLA_OUT_2.csv examples/10/TESLA_OUT_3.csv examples/10/TESLA_OUT_4.csv examples/10/TESLA_OUT_5.csv examples/10/TESLA_MAF.maf
echo "Test missing files 2 and 4"
python scoring_harness/TESLA_validation.py --patientId 10 examples/10/TESLA_OUT_1.csv examples/10/TESLA_OUT_3.csv examples/10/TESLA_OUT_5.csv examples/10/TESLA_VCF.vcf
