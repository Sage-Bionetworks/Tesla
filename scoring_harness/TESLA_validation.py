import os
import re
import argparse
import sys
import math
import operator
try:
	import pandas as pd
except ImportError:
	raise ImportError("Please make sure you have pandas installed: pip install pandas")

def checkType(submission, cols, colType, optional=False):
	for col in cols:
		if optional:
			assert all(submission[col].apply(lambda x: isinstance(x, colType) or math.isnan(x))), "All values in %s column must be type: %s [%s]" % (col, re.sub(".+['](.+)['].+","\\1",str(colType)), submission[col])
		else:
			assert all(submission[col].apply(lambda x: isinstance(x, colType))), "All values in %s column must be type: %s [%s]" % (col, re.sub(".+['](.+)['].+","\\1",str(colType)), submission[col])

def validate_1(submission_filepath):
	"""
	Validates first TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_1.csv
	"""
	print("VALIDATING %s" % submission_filepath)
	required_cols = pd.Series(["VAR_ID","CHROM","POS","OA_CALLER"])
	integer_cols = ['VAR_ID','POS']
	#NO duplicated VAR_ID
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from TESLA_OUT_1.csv: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
	#CHECK: CHROM must be 1-22 or X
	chroms = range(1,23)
	chroms = [str(i) for i in chroms]
	chroms.extend(["X","Y","MT"])
	submission.CHROM = submission.CHROM.astype(str)
	assert all(submission.CHROM.isin(chroms)), "CHROM values must be 1-22, X, Y or MT. You have: %s" % ", ".join(set(submission.CHROM[~submission.CHROM.isin(chroms)]))
	#CHECK: integer, string and float columns are correct types
	checkType(submission, integer_cols, int)

	return(True,"Passed Validation!")

def validate_2(submission_filepath):
	"""
	Validates second TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_2.csv
	"""
	#VAR_ID have to check out with first file
	print("VALIDATING %s" % submission_filepath)
	required_cols = pd.Series(["RANK","VAR_ID","PROT_POS","HLA_ALLELE","HLA_ALLELE_MUT","HLA_ALT_BINDING","HLA_REF_BINDING","PEP_LEN","ALT_EPI_SEQ","REF_EPI_SEQ","RANK_METRICS","RANK_DESC","ADDN_INFO"])

	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from TESLA_OUT_2.csv: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])

	integer_cols = ['VAR_ID','PROT_POS','PEP_LEN',"RANK"]
	string_cols = ['HLA_ALLELE','ALT_EPI_SEQ','REF_EPI_SEQ']
	#CHECK: RANK must be ordered from 1 to nrows
	assert all(submission.RANK == range(1, len(submission)+1)), "RANK column must be sequencial and must start from 1 to the length of the data"
	#CHECK: integer, string and float columns are correct types
	checkType(submission, integer_cols, int)
	checkType(submission, string_cols, str)
	checkType(submission, ['HLA_ALLELE_MUT',"RANK_DESC","ADDN_INFO"], str, optional=True)
	checkType(submission, ['HLA_ALT_BINDING','HLA_REF_BINDING'], float, optional=True)

	return(True,"Passed Validation!")

def validate_3(submission_filepath):
	"""
	Validates third TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_3.csv
	"""
	print("VALIDATING %s" % submission_filepath)
	required_cols = pd.Series(["VAR_ID","PROT_POS","HLA_ALLELE","HLA_ALLELE_MUT","HLA_ALT_BINDING","HLA_REF_BINDING","PEP_LEN","ALT_EPI_SEQ","REF_EPI_SEQ","STEP_ID"])
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from TESLA_OUT_3.csv: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
	integer_cols = ['VAR_ID','PROT_POS','PEP_LEN',"STEP_ID"]
	string_cols = ['HLA_ALLELE',"ALT_EPI_SEQ","REF_EPI_SEQ"]

	#CHECK: integer, string and float columns are correct types
	checkType(submission, integer_cols, int)
	checkType(submission, string_cols, str)
	checkType(submission, ['HLA_ALLELE_MUT'], str, optional=True)
	checkType(submission, ['HLA_ALT_BINDING','HLA_REF_BINDING'], float, optional=True)

	return(True,"Passed Validation!")

#Validate workflow
def validate_4(submission_filepath):
	"""
	Validates fourth TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_4.csv
	"""
	print("VALIDATING %s" % submission_filepath)
	required_cols = pd.Series(["STEP_ID","PREV_STEP_ID","DESC"])
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from TESLA_OUT_4.csv: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])

	checkType(submission, ["STEP_ID"], int)
	prevStepIds = [i.split(";") for i in submission['PREV_STEP_ID']]
	prevStepIds = reduce(operator.add, prevStepIds)
	stepIds = submission['STEP_ID'].tolist()
	stepIds.append(-1)
	assert all([int(i) in stepIds for i in prevStepIds]), "PREV_STEP_IDs must be -1 or existing STEP_IDs"
	checkType(submission, ["DESC"], str)

	return(True,"Passed Validation!")

### VALIDATING VCF
def contains_whitespace(x):
	"""
	Helper function for validateVCF.  No whitespace is allowed in VCF files

	:returns:     Sum of the the amount of whitespace in a string
	"""
	return(sum([" " in i for i in x if isinstance(i, str)]))

# Resolve missing read counts
def validateVCF(filePath):
	"""
	This function validates the VCF file to make sure it adhere to the genomic SOP.

	:params filePath:     Path to VCF file

	:returns:             Text with all the errors in the VCF file
	"""
	print("VALIDATING %s" % filePath)
	required_cols = pd.Series(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","TUMOR","NORMAL"])
	#FORMAT is optional
	with open(filePath,"r") as foo:
		for i in foo:
			if i.startswith("#CHROM"):
				headers = i.replace("\n","").split("\t")
	if headers is not None:
		submission = pd.read_csv(filePath, sep="\t",comment="#",header=None,names=headers)
	else:
		raise ValueError("Your vcf must start with the header #CHROM")

	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from TESLA_VCF.vcf: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])

	#Require that they report variants mapped to either GRCh37 or hg19 without
	#the chr-prefix. variants on chrM are not supported
	chroms = range(1,23)
	chroms = [str(i) for i in chroms]
	chroms.append("X")
	submission['#CHROM'] = submission['#CHROM'].astype(str)
	assert all(submission['#CHROM'].isin(chroms)), "CHROM values must be 1-22, or X. You have: %s" % ", ".join(set(submission.CHROM[~submission.CHROM.isin(chroms)]))
	#No white spaces
	temp = submission.apply(lambda x: contains_whitespace(x), axis=1)
	assert sum(temp) == 0, "Your vcf file should not have any white spaces in any of the columns"
	#I can also recommend a `bcftools query` command that will parse a VCF in a detailed way,
	#and output with warnings or errors if the format is not adhered too
	return(True,"Passed Validation!")

def validate_VAR_ID(submission1_filepath, submission2_filepath, submission3_filepath):
	submission1 = pd.read_csv(submission1_filepath)
	submission2 = pd.read_csv(submission2_filepath)
	submission3 = pd.read_csv(submission3_filepath)

	assert all(submission2['VAR_ID'].isin(submission1['VAR_ID'])), "TESLA_OUT_2.csv VAR_ID's must be part of TESLA_OUT_1.csv's VAR_IDs"
	assert all(submission3['VAR_ID'].isin(submission1['VAR_ID'])), "TESLA_OUT_3.csv VAR_ID's must be part of TESLA_OUT_1.csv's VAR_IDs"

	return(True, "Passed Validation!")

def validate_STEP_ID(submission3_filepath, submission4_filepath):
	submission3 = pd.read_csv(submission3_filepath)
	submission4 = pd.read_csv(submission4_filepath)

	assert all(submission3['STEP_ID'].isin(submission4['STEP_ID'].append(pd.Series([-1])))), "TESLA_OUT_3.csv STEP_ID's must be part of TESLA_OUT_4.csv's STEP_IDs"

	return(True, "Passed Validation!")

validation_func = {"TESLA_OUT_1.csv":validate_1,
				   "TESLA_OUT_2.csv":validate_2,
				   "TESLA_OUT_3.csv":validate_3,
				   "TESLA_OUT_4.csv":validate_4,
				   "TESLA_VCF.vcf":validateVCF}

def validate_files(filelist, patientId, validatingBAM=False):
	required=["TESLA_OUT_1.csv","TESLA_OUT_2.csv","TESLA_OUT_3.csv","TESLA_OUT_4.csv","TESLA_VCF.vcf"]
	if validatingBAM:
		print("VALIDATING BAMS")
		required.extend(["%s_EXOME_N.bam" % patientId ,"%s_EXOME_T.bam"% patientId,"%s_RNA_T.bam"% patientId])
	requiredFiles = pd.Series(required)
	basenames = [os.path.basename(name) for name in filelist]
	assert all(requiredFiles.isin(basenames)), "All %d submission files must be present and submission files must be named %s" % (len(required), ", ".join(required))
	for filepath in filelist:
		if not os.path.basename(filepath).endswith(".bam"):
			validation_func[os.path.basename(filepath)](filepath)
	onlyTesla = [i for i in filelist if "TESLA_OUT_" in i]
	order = pd.np.argsort(onlyTesla)
	print("VALIDATING THAT VARID EXISTS IN TESLA_OUT_{1,2,3}.csv")
	validate_VAR_ID(onlyTesla[order[0]],onlyTesla[order[1]],onlyTesla[order[2]])
	print("VALIDATING THAT STEPID EXISTS IN TESLA_OUT_{3,4}.csv")
	validate_STEP_ID(onlyTesla[order[2]],onlyTesla[order[3]])
	return(True, "Passed Validation!")

def perform_validate(args):
	validate_files(args.file, args.patientId, validatingBAM=args.validatingBAM)
	print("Passed Validation")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Validate TESLA files per sample')
	parser.add_argument("file", type=str, nargs="+",
						help='path to TESLA files (Must have TESLA_OUT_{1..4}.csv and TESLA_VCF.vcf), bam files are optional (include --validatingBAM and --patientId parameters if you include the BAM files)')
	parser.add_argument("--validatingBAM",action="store_true")
	parser.add_argument("--patientId",type=str,default=None,
						help='Patient Id')
	if ("--validatingBAM" in sys.argv) and ("--patientId" not in sys.argv):
	    parser.error("--patientId must be given if --validatingBAM is used")
	args = parser.parse_args()
	perform_validate(args)
