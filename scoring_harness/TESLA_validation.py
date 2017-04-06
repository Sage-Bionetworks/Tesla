import os
import re
import argparse
import sys
import math
import string
import getpass
try:
	import synapseclient
except ImportError:
	print("Please Install Synapse Client Library")
	print(">>> pip install synapseclient")
	sys.exit(1)
try:
	import pandas as pd
except ImportError:
	print("Please Pandas")
	print(">>> pip install pandas")
	sys.exit(1)

def synapse_login():
	try:
		syn = synapseclient.login()
	except Exception as e:
		print("Please provide your synapse username/email and password (You will only be prompted once)")
		try: 
			input = raw_input
		except NameError: 
			pass
		Username = input("Username: ")
		Password = getpass.getpass()
		syn = synapseclient.login(email=Username, password=Password,rememberMe=True)
	return syn

def configureHLA(i):
	return(str(i).replace("*","").split("(")[0])

def checkType(submission, cols, colType, fileName, optional=False,vcf=False):
	for col in cols:
		if optional:
			assert all(submission[col].apply(lambda x: isinstance(x, colType) or math.isnan(x))), "%s: All values in %s column must be type or blank/na: %s" % (fileName, col, re.sub(".+['](.+)['].+","\\1",str(colType)))
		elif vcf:
			assert all(submission[col].apply(lambda x: isinstance(x, colType) or x == ".")), "%s: All values in %s column must be type or .: %s" % (fileName, col, re.sub(".+['](.+)['].+","\\1",str(colType)))
		else:
			if col == "REF_EPI_SEQ":
				message = "%s: Please fill blank values with - that are equivalent to the length of the associated PEP_LEN.  (ie. if PEP_LEN is 5, then the REF_EPI_SEQ should be -----). All values in %s column must be type: %s"
			else:
				message = "%s: No blank values allowed and all values in %s column must be type: %s"
			assert all(submission[col].apply(lambda x: isinstance(x, colType))), message % (fileName, col, re.sub(".+['](.+)['].+","\\1",str(colType)))

def checkDelimiter(submission, cols, fileName, allowed=[';']):
	for col in cols:
		assert all(submission[col].apply(lambda x: all([i not in x if i not in allowed else True for i in string.punctuation]))),  "%s: All values in %s column must only have punctuation: [%s].  No other punctuation allowed in the string." % (fileName, col, " or ".join(allowed))

def validate_1(submission_filepath, validHLA):
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
	assert all(required_cols.isin(submission.columns)), "TESLA_OUT_1.csv: These column headers are missing: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
	#CHECK: CHROM must be 1-22 or X
	chroms = range(1,23)
	chroms = [str(i) for i in chroms]
	chroms.extend(["X","Y","MT"])
	strchroms = ["chr"+i for i in chroms]
	strchroms.pop()
	strchroms.append("chrM")
	submission.CHROM = submission.CHROM.astype(str)
	if submission.CHROM[0].startswith("chr"):
		assert all(submission.CHROM.isin(strchroms)), "TESLA_OUT_1.csv: CHROM values must be chr1-22, chrX, chrY, or chrM. You have: %s" % "or ".join(set(submission.CHROM[~submission.CHROM.isin(strchroms)]))
	else:
		assert all(submission.CHROM.isin(chroms)), "TESLA_OUT_1.csv: CHROM values must be 1-22, X, Y, or MT. You have: %s" % ", ".join(set(submission.CHROM[~submission.CHROM.isin(chroms)]))
	#CHECK: integer, string and float columns are correct types
	checkType(submission, integer_cols, int, "TESLA_OUT_1.csv")
	assert all(submission.VAR_ID == range(1, len(submission)+1)), "TESLA_OUT_1.csv: VAR_ID column must be sequencial and must start from 1 to the length of the data"
	checkType(submission, ['OA_CALLER'], str, "TESLA_OUT_1.csv",optional=True)
	submission['OA_CALLER'] = submission['OA_CALLER'].fillna('')
	checkDelimiter(submission, ['OA_CALLER'], "TESLA_OUT_1.csv")
	return(True,"Passed Validation!")

def validate_2(submission_filepath, validHLA):
	"""
	Validates second TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_2.csv
	"""
	#VAR_ID have to check out with first file
	print("VALIDATING %s" % submission_filepath)
	required_cols = pd.Series(["RANK","VAR_ID","PROT_POS","HLA_ALLELE","HLA_ALLELE_MUT","HLA_ALT_BINDING","HLA_REF_BINDING","PEP_LEN","ALT_EPI_SEQ","REF_EPI_SEQ","RANK_METRICS","RANK_DESC","ADDN_INFO"])

	submission = pd.read_csv(submission_filepath,na_values="n/a")
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "TESLA_OUT_2.csv: These column headers are missing: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])

	integer_cols = ['VAR_ID','PEP_LEN',"RANK"]
	string_cols = ['HLA_ALLELE','ALT_EPI_SEQ','REF_EPI_SEQ','RANK_METRICS']
	checkType(submission, integer_cols, int, 'TESLA_OUT_2.csv')
	#CHECK: RANK must be ordered from 1 to nrows
	assert all(submission.RANK == range(1, len(submission)+1)), "TESLA_OUT_2.csv: RANK column must be sequencial and must start from 1 to the length of the data"
	#CHECK: integer, string and float columns are correct types
	checkType(submission, string_cols, str, 'TESLA_OUT_2.csv')
	submission['RANK_DESC'] = submission['RANK_DESC'].fillna('').apply(str)
	checkType(submission, ['HLA_ALLELE_MUT',"RANK_DESC","ADDN_INFO"], str, 'TESLA_OUT_2.csv', optional=True)
	checkType(submission, ['HLA_ALT_BINDING','HLA_REF_BINDING'], float, 'TESLA_OUT_2.csv', optional=True)
	checkDelimiter(submission, ['RANK_METRICS'], "TESLA_OUT_2.csv",allowed=[';',':',".","_","-"])

	final_PROT_POS = []
	PROT_POS = [str(i).split(";") for i in submission['PROT_POS']]
	for i in PROT_POS:
		final_PROT_POS.extend(i)
	#PROT_POS = reduce(operator.add, PROT_POS)
	try:
		[int(i) for i in final_PROT_POS]
	except ValueError as e:
		raise AssertionError("TESLA_OUT_4.csv: PROT_POS must be semi-colon separated and must all be integers.")

	assert all(submission[['PEP_LEN','REF_EPI_SEQ']].apply(lambda x: len(x['REF_EPI_SEQ']) == x['PEP_LEN'], axis=1)), "TESLA_OUT_2.csv: Length of REF_EPI_SEQ values must be equal to the PEP_LEN"
	assert all(submission[['PEP_LEN','ALT_EPI_SEQ']].apply(lambda x: len(x['ALT_EPI_SEQ']) == x['PEP_LEN'], axis=1)), "TESLA_OUT_2.csv: Length of ALT_EPI_SEQ values must be equal to the PEP_LEN"
	assert all(submission['HLA_ALLELE'].apply(lambda x: configureHLA(x) in validHLA)), "TESLA_OUT_2.csv: HLA_ALLELE must be part of this list for this patient: %s" % ", ".join(validHLA)
	return(True,"Passed Validation!")

def validate_3(submission_filepath, validHLA):
	"""
	Validates third TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_3.csv
	"""
	print("VALIDATING %s" % submission_filepath)
	required_cols = pd.Series(["VAR_ID","PROT_POS","HLA_ALLELE","HLA_ALLELE_MUT","HLA_ALT_BINDING","HLA_REF_BINDING","PEP_LEN","ALT_EPI_SEQ","REF_EPI_SEQ","STEP_ID"])
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "TESLA_OUT_3.csv: These column headers are missing: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
	integer_cols = ['VAR_ID','PEP_LEN']
	string_cols = ['HLA_ALLELE',"ALT_EPI_SEQ","REF_EPI_SEQ"]
	PROT_POS = [str(i).split(";") for i in submission['PROT_POS']]
	final_PROT_POS = []
	PROT_POS = [str(i).split(";") for i in submission['PROT_POS']]
	for i in PROT_POS:
		final_PROT_POS.extend(i)
	try:
		[int(i) for i in final_PROT_POS]
	except ValueError as e:
		raise AssertionError("TESLA_OUT_4.csv: PROT_POS must be semi-colon separated and must all be integers.")

	#CHECK: integer, string and float columns are correct types
	checkType(submission, integer_cols, int, 'TESLA_OUT_3.csv')
	checkType(submission, string_cols, str, 'TESLA_OUT_3.csv')
	checkType(submission, ['STEP_ID'], int, 'TESLA_OUT_3.csv', optional=True)
	checkType(submission, ['HLA_ALLELE_MUT'], str, 'TESLA_OUT_3.csv', optional=True)
	checkType(submission, ['HLA_ALT_BINDING','HLA_REF_BINDING'], float, 'TESLA_OUT_3.csv', optional=True)

	assert all(submission[['PEP_LEN','REF_EPI_SEQ']].apply(lambda x: len(x['REF_EPI_SEQ']) == x['PEP_LEN'], axis=1)), "TESLA_OUT_3.csv: Length of REF_EPI_SEQ values must be equal to the PEP_LEN"
	assert all(submission[['PEP_LEN','ALT_EPI_SEQ']].apply(lambda x: len(x['ALT_EPI_SEQ']) == x['PEP_LEN'], axis=1)), "TESLA_OUT_3.csv: Length of ALT_EPI_SEQ values must be equal to the PEP_LEN"
	assert all(submission['HLA_ALLELE'].apply(lambda x: configureHLA(x) in validHLA)), "TESLA_OUT_3.csv: HLA_ALLELE must be part of this list for this patient: %s" % ", ".join(validHLA)

	return(True,"Passed Validation!")

def turnInt(i):
	try:
		i = int(i)
	except ValueError:
		i = -2
	return(i)

#Validate workflow
def validate_4(submission_filepath, validHLA):
	"""
	Validates fourth TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_4.csv
	"""
	print("VALIDATING %s" % submission_filepath)
	required_cols = pd.Series(["STEP_ID","PREV_STEP_ID","DESC"])
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "TESLA_OUT_4.csv: These column headers are missing: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
	
	checkType(submission, ["STEP_ID"], int, 'TESLA_OUT_4.csv')
	assert all(~submission['PREV_STEP_ID'].isnull()), "TESLA_OUT_4.csv: There must not be any NULL values in PREV_STEP_ID.  NULL values must be -1."
	prevStepIds = [str(i).split(";") for i in submission['PREV_STEP_ID']]
	final_prevStepIds = []
	for i in prevStepIds:
		final_prevStepIds.extend(i)
	stepIds = submission['STEP_ID'].tolist()
	stepIds.append(-1)
	assert all([turnInt(i) in stepIds for i in final_prevStepIds]), "TESLA_OUT_4.csv: PREV_STEP_IDs must be -1 or existing STEP_IDs"
	checkType(submission, ["DESC"], str, 'TESLA_OUT_4.csv')

	return(True,"Passed Validation!")

### VALIDATING VCF
def contains_whitespace(x):
	"""
	Helper function for validateVCF.  No whitespace is allowed in VCF files

	:returns:     Sum of the the amount of whitespace in a string
	"""
	return(sum([" " in i for i in x if isinstance(i, str)]))

# Resolve missing read counts
def validateVCF(filePath, validHLA):
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
		raise ValueError("TESLA_VCF.vcf: This file must start with the header #CHROM")

	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "TESLA_VCF.vcf: These column headers are missing: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])

	#Require that they report variants mapped to either GRCh37 or hg19 without
	#the chr-prefix. variants on chrM are not supported
	chroms = range(1,23)
	chroms = [str(i) for i in chroms]
	chroms.extend(["X","Y","MT"])
	strchroms = ["chr"+i for i in chroms]
	strchroms.pop()
	strchroms.append("chrM")
	submission['#CHROM'] = submission['#CHROM'].astype(str)
	if submission['#CHROM'][0].startswith("chr"):
		assert all(submission['#CHROM'].isin(strchroms)), "TESLA_VCF.vcf: CHROM values must be chr1-22, chrX, chrY, or chrM. You have: %s" % ", ".join(set(submission['#CHROM'][~submission['#CHROM'].isin(strchroms)]))
	else:
		assert all(submission['#CHROM'].isin(chroms)), "TESLA_VCF.vcf: CHROM values must be 1-22, X, Y, or MT. You have: %s" % ", ".join(set(submission['#CHROM'][~submission['#CHROM'].isin(chroms)]))
	
	#CHECK: integer, string columns are correct types and the missing value used must be .
	checkType(submission, ['POS'], int, 'TESLA_VCF.vcf')
	try:
		submission.QUAL = [float(i) if i != "." else i for i in submission.QUAL]
	except ValueError as e:
		raise AssertionError("TESLA_VCF.vcf: QUAL values must be numeric or .")
	checkType(submission, ['QUAL'], float, 'TESLA_VCF.vcf',vcf=True)

	required_string_cols = ['REF',"ALT",'FORMAT','TUMOR','NORMAL']
	opt_string_cols = ["ID",'FILTER','INFO']
	checkType(submission, required_string_cols, str, 'TESLA_VCF.vcf')
	checkType(submission, opt_string_cols, str, 'TESLA_VCF.vcf',vcf=True)

	#No white spaces
	temp = submission.apply(lambda x: contains_whitespace(x), axis=1)
	assert sum(temp) == 0, "TESLA_VCF.vcf: This file should not have any white spaces in any of the columns"
	#I can also recommend a `bcftools query` command that will parse a VCF in a detailed way,
	#and output with warnings or errors if the format is not adhered too
	return(True,"Passed Validation!")

def validateMAF(filePath, validHLA):
	#print("VALIDATING %s" % filePath)
	return(True, "Passed Validation!")

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

	submission3['STEP_ID'] = submission3['STEP_ID'].fillna(-1)
	assert all(submission3['STEP_ID'].isin(submission4['STEP_ID'].append(pd.Series([-1])))), "TESLA_OUT_3.csv STEP_ID's must be part of TESLA_OUT_4.csv's STEP_IDs"

	return(True, "Passed Validation!")

validation_func = {"TESLA_OUT_1.csv":validate_1,
				   "TESLA_OUT_2.csv":validate_2,
				   "TESLA_OUT_3.csv":validate_3,
				   "TESLA_OUT_4.csv":validate_4,
				   "TESLA_VCF.vcf":validateVCF,
				   "TESLA_MAF.maf":validateMAF}

def validate_files(filelist, patientId, validHLA, validatingBAM=False):
	required=["TESLA_OUT_1.csv","TESLA_OUT_2.csv","TESLA_OUT_3.csv","TESLA_OUT_4.csv"]
	vcfmaf = ["TESLA_VCF.vcf","TESLA_MAF.maf"]
	if validatingBAM:
		print("VALIDATING BAMS")
		required.extend(["%s_EXOME_N.bam" % patientId ,"%s_EXOME_T.bam"% patientId,"%s_RNA_T.bam"% patientId])
	requiredFiles = pd.Series(required)
	vcfmafFiles = pd.Series(vcfmaf)
	basenames = [os.path.basename(name) for name in filelist]
	assert all(requiredFiles.isin(basenames)), "All %d submission files must be present and submission files must be named %s" % (len(required), ", ".join(required))
	assert sum(vcfmafFiles.isin(basenames)) == 1, "Must have TESLA_VCF.vcf or TESLA_MAF.maf file"
	for filepath in filelist:
		if not os.path.basename(filepath).endswith(".bam"):
			validation_func[os.path.basename(filepath)](filepath, validHLA)
	onlyTesla = [i for i in filelist if "TESLA_OUT_" in i]
	order = pd.np.argsort(onlyTesla)
	print("VALIDATING THAT VARID EXISTS IN TESLA_OUT_{1,2,3}.csv")
	validate_VAR_ID(onlyTesla[order[0]],onlyTesla[order[1]],onlyTesla[order[2]])
	print("VALIDATING THAT STEPID EXISTS IN TESLA_OUT_{3,4}.csv")
	validate_STEP_ID(onlyTesla[order[2]],onlyTesla[order[3]])
	return(True, "Passed Validation!")


def perform_validate(args):
	syn = synapse_login()
	metadataPath = syn.get("syn8371011").path
	metadata = pd.read_csv(metadataPath)
	HLA = metadata[['patientId','classIHLAalleles']][~metadata['classIHLAalleles'].isnull()]
	assert args.patientId in metadata['patientId'], "Patient Id must be in the metadata"
	listHLA = HLA['classIHLAalleles'][HLA['patientId'] == args.patientId]
	validHLA = [i.replace("*","").split(";") for i in listHLA]
	final_validHLA = []
	for i in validHLA:
		final_validHLA.extend(i)
	final_validHLA = set([i.split("(")[0] for i in final_validHLA])
	validate_files(args.file, args.patientId, final_validHLA, validatingBAM=args.validatingBAM)
	print("Passed Validation")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Validate TESLA files per sample')
	parser.add_argument("file", type=str, nargs="+",
						help='path to TESLA files (Must have TESLA_OUT_{1..4}.csv and TESLA_VCF.vcf), bam files are optional (include --validatingBAM and --patientId parameters if you include the BAM files)')
	parser.add_argument("--patientId",type=int, required=True,
						help='Patient Id')
	parser.add_argument("--validatingBAM",action="store_true")
	args = parser.parse_args()
	perform_validate(args)
