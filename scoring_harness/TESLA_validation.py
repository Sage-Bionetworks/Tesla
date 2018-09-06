import os
import re
import argparse
import subprocess
from subprocess import Popen, PIPE, STDOUT
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
try:
	import yaml
except ImportError:
	print("Please Install pyyaml")
	print(">>> pip install pyyaml")
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
	return(str(i).replace("*","").split("(")[0].replace("HLA-",''))

def checkType(submission, cols, colType, fileName, optional=False,vcf=False):
	for col in cols:
		if optional:
			assert all(submission[col].apply(lambda x: isinstance(x, colType) or math.isnan(x))), "%s: All values in %s column must be type: %s or blank/NA" % (fileName, col, re.sub(".+['](.+)['].+","\\1",str(colType)))
		elif vcf:
			assert all(submission[col].apply(lambda x: isinstance(x, colType) or x == ".")), "%s: All values in %s column must be type: %s or ." % (fileName, col, re.sub(".+['](.+)['].+","\\1",str(colType)))
		else:
			if col == "REF_EPI_SEQ":
				message = "%s: Please fill blank/NA values with - that are equivalent to the length of the associated PEP_LEN.  (ie. if PEP_LEN is 5, then the REF_EPI_SEQ should be -----). All values in %s column must be type: %s"
			else:
				message = "%s: No blank/NA values allowed and all values in %s column must be type: %s"
			if col == "VAR_ID":
				assert all([i != "nan" for i in submission[col]]),  message % (fileName, col, re.sub(".+['](.+)['].+","\\1",str(colType)))
			assert all(submission[col].apply(lambda x: isinstance(x, colType))), message % (fileName, col, re.sub(".+['](.+)['].+","\\1",str(colType)))

def checkDelimiter(submission, cols, fileName, allowed=[';']):
	for col in cols:
		assert all(submission[col].apply(lambda x: all([i not in x if i not in allowed else True for i in string.punctuation]))),  "%s: All values in %s column must only have punctuation: [%s].  No other punctuation allowed in the string." % (fileName, col, " or ".join(allowed))

def intSemiColonListCheck(submission, fileName, col):
	allResults = []
	results = [str(i).split(";") for i in submission[col]]
	for i in results:
			allResults.extend(i)
	try:
		[int(i) for i in allResults]
	except ValueError as e:
		raise AssertionError("%s: %s can be semi-colon separated but all values must be integers." %(fileName, col))
	return(pd.Series(allResults).astype(int))

def semiColonIdListCheck(submission, fileName, col, vcfIDs):
	allResults = []
	results = [str(i).replace(" ","").split(";") for i in submission[col]]
	for i in results:
		allResults.extend(i)
	assert all(pd.Series(allResults).isin(vcfIDs)), "%s: %s can be semi-colon separated list but all values must be in your VCF ID column." %(fileName, col)

def validate_1_2(submission_filepath, validHLA):
	"""
	Validates first and second TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_{1..2}.csv
	"""
	#VAR_ID have to check out with first file
	print("VALIDATING %s" % submission_filepath)
	basename = os.path.basename(submission_filepath)
	required_cols = pd.Series(["RANK","VAR_ID","PROT_POS","HLA_ALLELE","HLA_ALLELE_MUT","HLA_ALT_BINDING","HLA_REF_BINDING","PEP_LEN","ALT_EPI_SEQ","REF_EPI_SEQ","RANK_METRICS","RANK_DESC","ADDN_INFO","SCORE",'REF_ALLELE_EXP','ALT_ALLELE_EXP'])

	submission = pd.read_csv(submission_filepath,na_values="n/a")
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "%s: These column headers are missing: %s" % (basename,", ".join(required_cols[~required_cols.isin(submission.columns)]))
	submission['VAR_ID'] = submission['VAR_ID'].astype(str)
	integer_cols = ['PEP_LEN',"RANK"]
	string_cols = ['HLA_ALLELE','ALT_EPI_SEQ','REF_EPI_SEQ','RANK_METRICS','VAR_ID']
	checkType(submission, integer_cols, int, basename)
	#CHECK: RANK must be ordered from 1 to nrows
	assert all(submission.RANK == range(1, len(submission)+1)), "%s: RANK column must be sequencial and must start from 1 to the length of the data" % basename
	#CHECK: integer, string and float columns are correct types
	checkType(submission, string_cols, str, basename)
	submission['RANK_DESC'] = submission['RANK_DESC'].fillna('').apply(str)
	checkType(submission, ['HLA_ALLELE_MUT',"RANK_DESC","ADDN_INFO"], str, basename, optional=True)
	checkType(submission, ['HLA_ALT_BINDING','HLA_REF_BINDING','SCORE','REF_ALLELE_EXP','ALT_ALLELE_EXP'], float, basename, optional=True)
	checkDelimiter(submission, ['RANK_METRICS'], basename,allowed=[';',':',".","_","-"])
	intSemiColonListCheck(submission, basename, 'PROT_POS')

	assert all(submission[['PEP_LEN','REF_EPI_SEQ']].apply(lambda x: len(x['REF_EPI_SEQ']) == x['PEP_LEN'], axis=1)), "%s: Length of REF_EPI_SEQ values must be equal to the PEP_LEN" % basename
	assert all(submission[['PEP_LEN','ALT_EPI_SEQ']].apply(lambda x: len(x['ALT_EPI_SEQ']) == x['PEP_LEN'], axis=1)), "%s: Length of ALT_EPI_SEQ values must be equal to the PEP_LEN" % basename
	assert all(submission['HLA_ALLELE'].apply(lambda x: configureHLA(x) in validHLA)), "%s: HLA_ALLELE must be part of this list for this patient: %s" % (basename,", ".join(validHLA))
	return(True,"Passed Validation!")

def validate_3_4(submission_filepath, validHLA):
	"""
	Validates third and fourth TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_{3..4}.csv
	"""
	print("VALIDATING %s" % submission_filepath)
	basename = os.path.basename(submission_filepath)
	required_cols = pd.Series(["VAR_ID","PROT_POS","HLA_ALLELE","HLA_ALLELE_MUT","HLA_ALT_BINDING","HLA_REF_BINDING","PEP_LEN","ALT_EPI_SEQ","REF_EPI_SEQ","STEP_ID",'SCORE','REF_ALLELE_EXP','ALT_ALLELE_EXP'])
	submission = pd.read_csv(submission_filepath,na_values="n/a")
	assert all(required_cols.isin(submission.columns)), "%s: These column headers are missing: %s" % (basename,", ".join(required_cols[~required_cols.isin(submission.columns)]))
	if not submission.empty:
		#CHECK: Required headers must exist in submission
		integer_cols = ['PEP_LEN']
		string_cols = ['HLA_ALLELE',"ALT_EPI_SEQ","REF_EPI_SEQ",'VAR_ID']
		submission['VAR_ID'] = submission['VAR_ID'].astype(str)
		intSemiColonListCheck(submission, basename, 'PROT_POS')
		#intSemiColonListCheck(submission, "TESLA_OUT_2.csv", 'VAR_ID')

		#CHECK: integer, string and float columns are correct types
		checkType(submission, integer_cols, int, basename)
		checkType(submission, string_cols, str, basename)
		#Fill STEP_ID na's with an integer and change the entire column to int
		submission['STEP_ID'] = submission['STEP_ID'].fillna(-1)
		submission['STEP_ID'] = submission['STEP_ID'].apply(int)
		checkType(submission, ['STEP_ID'], int, basename, optional=True)
		checkType(submission, ['HLA_ALLELE_MUT'], str, basename, optional=True)
		checkType(submission, ['HLA_ALT_BINDING','HLA_REF_BINDING','REF_ALLELE_EXP','ALT_ALLELE_EXP'], float, basename, optional=True)

		assert all(submission[['PEP_LEN','REF_EPI_SEQ']].apply(lambda x: len(x['REF_EPI_SEQ']) == x['PEP_LEN'], axis=1)), "%s: Length of REF_EPI_SEQ values must be equal to the PEP_LEN" % basename
		assert all(submission[['PEP_LEN','ALT_EPI_SEQ']].apply(lambda x: len(x['ALT_EPI_SEQ']) == x['PEP_LEN'], axis=1)), "%s: Length of ALT_EPI_SEQ values must be equal to the PEP_LEN" % basename
		assert all(submission['HLA_ALLELE'].apply(lambda x: configureHLA(x) in validHLA)), "%s: HLA_ALLELE must be part of this list for this patient: %s" % (basename,", ".join(validHLA))
	return(True,"Passed Validation!")

def turnInt(i):
	try:
		i = int(i)
	except ValueError:
		i = -2
	return(i)



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

	required_string_cols = ['REF',"ALT",'FORMAT','TUMOR','NORMAL','ID']
	opt_string_cols = ['FILTER','INFO']
	submission['ID'] = submission['ID'].apply(str)
	checkType(submission, required_string_cols, str, 'TESLA_VCF.vcf')
	checkType(submission, opt_string_cols, str, 'TESLA_VCF.vcf',vcf=True)
	assert sum(submission.ID.duplicated()) == 0, "TESLA_VCF.vcf: The ID column must not have any duplicates"
	#No white spaces
	temp = submission.apply(lambda x: contains_whitespace(x), axis=1)
	assert sum(temp) == 0, "TESLA_VCF.vcf: This file should not have any white spaces in any of the columns"
	#I can also recommend a `bcftools query` command that will parse a VCF in a detailed way,
	#and output with warnings or errors if the format is not adhered too
	try:
		cmd = ['docker','run','-v','%s:/TESLA_VCF.vcf' % os.path.abspath(filePath), "thomasvyu/vcf-validator:0.6"]
		#cmd = ['./%s' % vcfValidator, "-i", filePath]
		p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
	except OSError as e:
		raise ValueError('Please make sure docker is installed.')
	output = p.stdout.read()
	result = '\nAccording to the VCF specification, the input file is valid\n' in output
	assert result, output
	return(True,"Passed Validation!")


def validate_yaml(submission_filepath, template_filepath="../TESLA_YAML.yaml"):

    """
    Validate the TESLA submission yaml file.
    :param submission_filepath: Path of submission file TESLA_YAML.yaml
    :param template_filepath: Path of TESLA_YAML.yaml template
    """
    print("VALIDATING %s" % submission_filepath)

    template = yaml.load(open(template_filepath))
    submission = yaml.load(open(submission_filepath))
    required_steps = template.keys()
    submited_steps = submission.keys()
    parameter_steps = [key for key in required_steps
                       if "key_parameters" in template[key].keys()]

    # check steps
    missing_steps = [key for key in required_steps
                     if key not in submited_steps]
    assert len(missing_steps) == 0, (
        "Step(s) missing from TESLA_YAML.yaml: [" +
        ", ".join(missing_steps) + "]")

    extra_steps = [key for key in submited_steps if key not in required_steps]
    assert len(extra_steps) == 0, (
        "Extra step(s) in TESLA_YAML.yaml: [" +
        ", ".join(extra_steps) + "]")

    # check used fields
    missing_used_fields = [step for step in submited_steps
                           if 'used' not in submission[step].keys()]
    assert len(missing_used_fields) == 0, (
        "Step(s) in TESLA_YAML.yaml are missing used field: "
        "[" + ", ".join(missing_used_fields) + "]")

    bad_used_fields = [step for step in submited_steps
                       if type(submission[step]['used']) is not bool]
    assert len(bad_used_fields) == 0, (
        "Step(s) in TESLA_YAML.yaml have used fields that are not boolean: "
        "[" + ", ".join(bad_used_fields) + "]")

    steps_used = [step for step in submited_steps if submission[step]['used']]
    assert len(steps_used) > 0, "No steps indicated in TESLA_OUT.yaml"

    # check changed fields
    missing_changed_fields = [step for step in submited_steps
                              if 'changed' not in submission[step].keys()]
    assert len(missing_changed_fields) == 0, (
        "Step(s) in TESLA_YAML.yaml are missing changed field: "
        "[" + ", ".join(missing_changed_fields) + "]")

    bad_changed_fields = [step for step in submited_steps
                          if type(submission[step]['changed']) is not bool]
    assert len(bad_changed_fields) == 0, (
        "Step(s) in TESLA_YAML.yaml have changed fields that are not boolean: "
        "[" + ", ".join(bad_changed_fields) + "]")

    # check comment fields
    missing_comment_fields = [step for step in submited_steps
                              if 'comment' not in submission[step].keys()]
    assert len(missing_comment_fields) == 0, (
        "Step(s) in TESLA_YAML.yaml are missing comment field: "
        "[" + ", ".join(missing_comment_fields) + "]")

    # check key_parameters fields
    missing_key_parameters = [step for step in parameter_steps
                              if 'key_parameters' not in
                              submission[step].keys()]

    assert len(missing_key_parameters) == 0, (
        "Step(s) in TESLA_YAML.yaml are missing key_parameters field: "
        "[" + ", ".join(missing_key_parameters) + "]")

    used_key_parameters = [step for step in parameter_steps
                           if submission[step]['key_parameters'] is not None]

    VALUE_RELATIONSHIPS = [">", ">=", "=", "<=", "<"]
    CATEGORY_RELATIONSHIPS = ["in", "not"]
    ALL_RELATIONSHIPS = VALUE_RELATIONSHIPS + CATEGORY_RELATIONSHIPS

    for step in used_key_parameters:
        for param in submission[step]['key_parameters']:

            assert 'name' in param.keys(), (
                "Step: " + step + " key parmater missing name field.")

            assert 'relationship' in param.keys(), (
                "Step: " + step + " key parameter: " + param['name'] +
                " missing relationship field.")

            assert param['relationship'] in ALL_RELATIONSHIPS, (
                "Step: " + step + " key parameter: " + param['name'] +
                " not one of [" + " ".join(ALL_RELATIONSHIPS) + "].")

            if param['relationship'] in VALUE_RELATIONSHIPS:

                assert 'unit' in param.keys(), (
                    "Step: " + step + " key parameter: " + param['name'] +
                    " missing unit field.")

                assert 'value' in param.keys(), (
                    "Step: " + step + " key parameter: " + param['name'] +
                    " missing value field.")

            if param['relationship'] in CATEGORY_RELATIONSHIPS:

                assert 'values' in param.keys(), (
                    "Step: " + step + " key parameter: " + param['name'] +
                    " missing values field.")

    return(True, "Passed Validation!")


def validate_VAR_ID(submission1_filepath, submission3_filepath, submissionvcf_filepath, submission2_filepath=None, submission4_filepath=None,patientVCFDf=None):
	with open(submissionvcf_filepath,"r") as foo:
		for i in foo:
			if i.startswith("#CHROM"):
				headers = i.replace("\n","").split("\t")
	submissionvcf = pd.read_csv(submissionvcf_filepath, sep="\t",comment="#",header=None,names=headers)
	submission1 = pd.read_csv(submission1_filepath)
	submission3 = pd.read_csv(submission3_filepath)
	submissionvcf['ID'] = submissionvcf['ID'].apply(str)
	semiColonIdListCheck(submission3, "TESLA_OUT_3.csv", 'VAR_ID', submissionvcf['ID'])
	semiColonIdListCheck(submission1, "TESLA_OUT_1.csv", 'VAR_ID', submissionvcf['ID'])

	#assert all(submission3['VAR_ID'].apply(str).isin(submissionvcf['ID'])), "TESLA_OUT_3.csv VAR_ID's must be part of TESLA_VCF.vcf's IDs"
	#assert all(submission1['VAR_ID'].apply(str).isin(submissionvcf['ID'])), "TESLA_OUT_1.csv VAR_ID's must be part of TESLA_VCF.vcf's IDs"
	if submission2_filepath is not None and submission4_filepath is not None:
		submission2 = pd.read_csv(submission2_filepath)
		submission4 = pd.read_csv(submission4_filepath)
		semiColonIdListCheck(submission2, "TESLA_OUT_2.csv", 'VAR_ID', submissionvcf['ID'])
		#assert all(submission2['VAR_ID'].apply(str).isin(patientVCFDf[2])), "TESLA_OUT_2.csv VAR_ID's must be part of the patient VCF's ID's"
		semiColonIdListCheck(submission4, "TESLA_OUT_2.csv", 'VAR_ID', submissionvcf['ID'])
		#assert all(submission4['VAR_ID'].apply(str).isin(patientVCFDf[2])), "TESLA_OUT_4.csv VAR_ID's must be part of the patient VCF's ID's"
	return(True, "Passed Validation!")

def validate_STEP_ID(submission3_filepath, submission5_filepath, submission4_filepath=None):
	submission3 = pd.read_csv(submission3_filepath)
	submission5 = pd.read_csv(submission5_filepath)
	submission3['STEP_ID'] = submission3['STEP_ID'].fillna(-1)
	assert all(submission3['STEP_ID'].isin(submission5['STEP_ID'].append(pd.Series([-1])))), "TESLA_OUT_3.csv STEP_ID's must be part of TESLA_OUT_5.csv's STEP_IDs"

	if submission4_filepath is not None:
		submission4 = pd.read_csv(submission4_filepath)
		submission4['STEP_ID'] = submission4['STEP_ID'].fillna(-1)
		assert all(submission4['STEP_ID'].isin(submission5['STEP_ID'].append(pd.Series([-1])))), "TESLA_OUT_4.csv STEP_ID's must be part of TESLA_OUT_5.csv's STEP_IDs"
	return(True, "Passed Validation!")

validation_func = {"TESLA_OUT_1.csv":validate_1_2,
				   "TESLA_OUT_2.csv":validate_1_2,
				   "TESLA_OUT_3.csv":validate_3_4,
				   "TESLA_OUT_4.csv":validate_3_4,
				   "TESLA_OUT_YAML.yaml":validate_yaml,
				   "TESLA_VCF.vcf":validateVCF}

def validate_files(syn, filelist, patientId, validHLA, validatingBAM=False):
	required=["TESLA_OUT_1.csv","TESLA_OUT_3.csv","TESLA_OUT_5.csv","TESLA_VCF.vcf"]
	optional = ["TESLA_OUT_2.csv", "TESLA_OUT_4.csv"]
	if validatingBAM:
		print("VALIDATING BAMS")
		required.extend(["%s_EXOME_N.bam" % patientId ,"%s_EXOME_T.bam"% patientId,"%s_RNA_T.bam"% patientId])
	requiredFiles = pd.Series(required)
	optionalFiles = pd.Series(optional)
	basenames = [os.path.basename(name) for name in filelist]

	useOptional = all(optionalFiles.isin(basenames))
	assert all(requiredFiles.isin(basenames)), "All %d submission files must be present and submission files must be named %s" % (len(required), ", ".join(required))
	assert useOptional or sum(optionalFiles.isin(basenames)) == 0, "TESLA_OUT_2.csv, TESLA_OUT_4.csv.  Both files MUST either be present or missing.  If missing, you are missing predictions from VCF. If this is not as intended, please submit again."
	for filepath in filelist:
		if os.path.basename(filepath) in ['TESLA_OUT_1.csv','TESLA_OUT_2.csv','TESLA_OUT_3.csv','TESLA_OUT_4.csv']:
			validation_func[os.path.basename(filepath)](filepath, validHLA)
		elif not os.path.basename(filepath).endswith(".bam") and os.path.basename(filepath) != "TESLA_ranking_method.txt":
			validation_func[os.path.basename(filepath)](filepath)
	onlyTesla = [i for i in filelist if "TESLA_OUT_" in i or "TESLA_VCF" in i]
	order = pd.np.argsort(onlyTesla)

	if useOptional:
		patientFiles = syn.tableQuery('SELECT * FROM syn8292741 where patientId = "%s" and fileFormat = "vcf"' % patientId)
		patientFilesDf = patientFiles.asDataFrame()
		patientVCFEnt = syn.get(patientFilesDf['id'][0])
		patientVCFDf = pd.read_csv(patientVCFEnt.path,sep="\t",comment="#",header=None)
		print("VALIDATING THAT VARID EXISTS IN TESLA_OUT_{1,3}.csv and maps to ID in TESLA_VCF.vcf and TESLA_OUT_{2,4}.csv maps to ID in %s" % patientVCFEnt.name)
		validate_VAR_ID(onlyTesla[order[0]],onlyTesla[order[2]],onlyTesla[order[5]],submission2_filepath=onlyTesla[order[1]],submission4_filepath=onlyTesla[order[3]],patientVCFDf=patientVCFDf)
	else:
		print("VALIDATING THAT VARID EXISTS IN TESLA_OUT_{1,3}.csv and maps to ID in TESLA_VCF.vcf")
		validate_VAR_ID(onlyTesla[order[0]],onlyTesla[order[1]],onlyTesla[order[3]])

	if useOptional:
		print("VALIDATING THAT STEPID EXISTS IN TESLA_OUT_{3,4,5}.csv")
		validate_STEP_ID(onlyTesla[order[2]],onlyTesla[order[4]],submission4_filepath=onlyTesla[order[3]])
	else:
		print("VALIDATING THAT STEPID EXISTS IN TESLA_OUT_{3,5}.csv")
		validate_STEP_ID(onlyTesla[order[1]],onlyTesla[order[2]])

	return(True, useOptional, "Passed Validation!")


def perform_validate(args):
	syn = synapse_login()
	#metadataPath = syn.get("syn8371011").path
	#metadata = pd.read_csv(metadataPath)
	metadataTable = syn.tableQuery('SELECT * FROM syn8292741')
	metadata = metadataTable.asDataFrame()
	HLA = metadata[['patientId','classIHLAalleles']][~metadata['classIHLAalleles'].isnull()]
	HLA.drop_duplicates("classIHLAalleles",inplace=True)
	patientIds = set(metadata['patientId'][~metadata['patientId'].isnull()].apply(int))
	assert args.patientId in patientIds, "Patient Id must be in the metadata"
	listHLA = HLA['classIHLAalleles'][HLA['patientId'] == args.patientId]
	validHLA = [i.replace("*","").split(";") for i in listHLA]
	final_validHLA = []
	for i in validHLA:
		final_validHLA.extend(i)
	final_validHLA = set([i.split("(")[0] for i in final_validHLA])
	validate_files(syn, args.file, args.patientId, final_validHLA, validatingBAM=args.validatingBAM)
	print("Passed Validation")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Validate TESLA files per sample')
	parser.add_argument("file", type=str, nargs="+",
						help='path to TESLA files (Must have TESLA_OUT_{1..4}.csv, TESLA_YAML.yaml and TESLA_VCF.vcf), bam files are optional (include --validatingBAM and --patientId parameters if you include the BAM files)')
	parser.add_argument("--patientId",type=int, required=True,
						help='Patient Id')
	parser.add_argument("--validatingBAM",action="store_true")
	args = parser.parse_args()
	perform_validate(args)
