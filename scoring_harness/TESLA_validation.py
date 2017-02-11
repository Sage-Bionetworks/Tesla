import os
import re
import json
import argparse
from multiprocessing.dummy import Pool as ThreadPool 
from functools import partial
try:
	import httplib2 as http
except ImportError:
	raise ImportError("Please make sure you have httplib2 installed: pip install httplib2")
try:
	import pandas as pd
except ImportError:
	raise ImportError("Please make sure you have pandas installed: pip install pandas")
try:
	from urlparse import urlparse
except ImportError:
	from urllib.parse import urlparse

pool = ThreadPool(4)

#Validate genes
def hgncRestCall(path):
	"""
	This function does the rest call to the genenames website

	:params path:     The gene symbol url path to add to the base uri

	:returns:         If the symbol exists, returns True and the corrected symbol, otherwise returns False and None.
	"""
	headers = {'Accept': 'application/json',}

	uri = 'http://rest.genenames.org'

	target = urlparse(uri+path)
	method = 'GET'
	body = ''
	h = http.Http()
	response, content = h.request(target.geturl(),
								  method,
								  body,
								  headers)
	if response['status'] == '200':
		data = json.loads(content)
		if len(data['response']['docs']) == 0:
			return(False, [None])
		else:
			#print(len(data['response']['docs']))
			mapped = [symbol['symbol'] for symbol in data['response']['docs']]
			return(True, mapped)
	else:
		return(False, [None])

# Validation of gene names
def validateSymbol(gene, returnMapping=False, entrez=False):
	"""
	This function does validation of symbols

	:params gene:               Gene symbol
	:params returnMapping:      Return mapping of old gene to new gene

	:returns:                   Check if the provided gene name is a correct symbol and print out genes 
								that need to be remapped or can't be mapped to anything
	"""
	path = '/fetch/symbol/%s' %  gene
	verified, symbol = hgncRestCall(path)
	if not verified and not entrez:
		path = '/fetch/prev_symbol/%s' %  gene
		verified, symbol = hgncRestCall(path)
	if not verified and not entrez:
		path = '/fetch/alias_symbol/%s' %  gene
		verified, symbol = hgncRestCall(path)
	if not verified and entrez:
		path = '/fetch/ensembl_gene_id/%s' %  gene
		verified, symbol = hgncRestCall(path)
		if verified:
			return(True)
		else:
			return(False)
	if gene in symbol:
		return(True)
	else:
		if symbol[0] is None:
			print("%s cannot be remapped. Please correct." % gene)
		else:
			#if "MLL4", then the HUGO symbol should be KMT2D and KMT2B
			if len(symbol) > 1:
				print("%s can be mapped to different symbols: %s. Please correct." % (gene, ", ".join(symbol)))
			else:
				print("%s should be remapped to %s" % (gene, symbol[0]))
				if returnMapping:
					return({gene, symbol[0]})
				else:
					return(True)
		return(False)

def checkType(submission, cols, colType):
	for col in cols:
		assert all(submission[col].apply(lambda x: isinstance(x, colType))), "All values in %s column must be type: %s" % (col, re.sub(".+['](.+)['].+","\\1",str(colType)))

def validate_1(submission_filepath):
	"""
	Validates first TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_1.csv
	"""
	required_cols = pd.Series(["VAR_ID","GENE","SYMBOL","CHROM","START","END","CLASS","REF","ALT","REFCOUNT_T","ALTCOUNT_T",
					 "REFCOUNT_N","ALTDEPTH_N","REF_FPKM_T","REF_HTSEQ_T","ALT_FPKM_T","ALT_HTSEQ_T","QUAL","VARID","INFO"])
	integer_cols = ['VAR_ID','START','END','REFCOUNT_T','ALTCOUNT_T','REFCOUNT_N','ALTDEPTH_N','REF_HTSEQ_T','ALT_HTSEQ_T','QUAL']
	string_cols = ['REF','ALT','VARID','INFO']
	float_cols = ['ALT_FPKM_T','REF_FPKM_T']
	#CLASS_categories = ['intron','missense','silent','splice_site','in_frame_deletion','in_frame_insertion',
	#					'frame_shift_deletion','frame_shift_insertion',
	#					'nonsense_mutation','structural_variant','splice_isoform','other']
	CLASS_categories = ['transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained','frameshift_variant',
						'stop_lost','start_lost','transcript_amplification','inframe_insertion','inframe_deletion','missense_variant',
						'protein_altering_variant','splice_region_variant','incomplete_terminal_codon_variant','stop_retained_variant',
						'synonymous_variant','coding_sequence_variant','mature_miRNA_variant','5_prime_UTR_variant','3_prime_UTR_variant',
						'non_coding_transcript_exon_variant','intron_variant','NMD_transcript_variant','non_coding_transcript_variant',
						'upstream_gene_variant','downstream_gene_variant','TFBS_ablation','TFBS_amplification','TF_binding_site_variant',
						'regulatory_region_ablation','regulatory_region_amplification','feature_elongation','regulatory_region_variant',
						'feature_truncation','intergenic_variant']
	#NO duplicated VAR_ID
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from your file: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
	#CHECK: CHROM must be 1-22 or X
	assert all(submission.CHROM.isin(range(1,23) + ["X"])), "CHROM values must be 1-22, or X. You have: %s" % ", ".join(set(submission.CHROM[~submission.CHROM.isin(range(1,23) + ["X"])])) 
	#CHECK: CLASS must be in CLASS_categories
	assert all(submission.CLASS.isin(CLASS_categories)), "CLASS values must be one these: %s.  You have: %s" % (", ".join(CLASS_categories),", ".join(set(submission.CLASS[~submission.CLASS.isin(CLASS_categories)])))
	#CHECK: STRAND must be +/-
	#assert all(submission.STRAND.isin(['+','-'])), "STRAND values must be + or -.  You have: %s" %(", ".join(set(submission.STRAND[~submission.STRAND.isin(['+','-'])])))
	#CHECK: HUGO SYMBOL validation
	assert all(pool.map(validateSymbol, submission.SYMBOL.drop_duplicates())), "All gene symbols in GENE column must be up to date (hgnc standards)"
	#CHECK: ENTREZ GENE ID VALIDATION
	validateGene = partial(validateSymbol, entrez=True)
	assert all(pool.map(validateGene, submission.GENE.drop_duplicates())), "All gene symbols in GENE column must be up to date (hgnc standards)"

	#CHECK: integer, string and float columns are correct types
	checkType(submission, integer_cols, int)
	checkType(submission, string_cols, str)
	checkType(submission, float_cols, float)

	return(True,"Passed Validation!")


def validate_2(submission_filepath):
	"""
	Validates second TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_2.csv
	"""
	#VAR_ID have to check out with first file
	required_cols = pd.Series(["RANK","VAR_ID","PROT_POS","HLA_ALLELE","PEP_LEN","ALT_EPI_SEQ","WT_EPI_SEQ"])
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from your file: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])

	integer_cols = ['VAR_ID','PROT_POS','PEP_LEN']
	string_cols = ['HLA_ALLELE','ALT_EPI_SEQ','WT_EPI_SEQ']
	#CHECK: RANK must be ordered from 1 to nrows
	assert all(submission.RANK == range(1, len(submission)+1)), "RANK column must be sequencial and must start from 1 to the length of the data"
	#CHECK: integer, string and float columns are correct types
	checkType(submission, integer_cols, int)
	checkType(submission, string_cols, str)

	return(True,"Passed Validation!")


def validate_3(submission_filepath):
	"""
	Validates second TESLA file

	:param submission_filepath: Path of submission file TESLA_OUT_2.csv
	"""
	required_cols = pd.Series(["VAR_ID","PROT_POS","HLA_ALLELE","PEP_LEN","ALT_EPI_SEQ","WT_EPI_SEQ","STEP_ID"])
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from your file: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
	integer_cols = ['VAR_ID','PROT_POS','PEP_LEN','STEP_ID']
	string_cols = ['HLA_ALLELE','ALT_EPI_SEQ','WT_EPI_SEQ']

	#CHECK: integer, string and float columns are correct types
	checkType(submission, integer_cols, int)
	checkType(submission, string_cols, str)

	return(True,"Passed Validation!")

#Validate workflow
def validate_4(submission_filepath):
	required_cols = pd.Series(["STEP_ID","PREV_STEP_ID","DESC"])
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from your file: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
	
	checkType(submission, ["STEP_ID","PREV_STEP_ID"], int)
	checkType(submission, ["DESC"], str)

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

	assert all(submission3['STEP_ID'].isin(submission4['STEP_ID'])), "TESLA_OUT_3.csv STEP_ID's must be part of TESLA_OUT_4.csv's STEP_IDs"

	return(True, "Passed Validation!")

validation_func = {"TESLA_OUT_1.csv":validate_1,
				   "TESLA_OUT_2.csv":validate_2,
				   "TESLA_OUT_3.csv":validate_3,
				   "TESLA_OUT_4.csv":validate_4}

def validate_files(filelist):
	requiredFiles = pd.Series(["TESLA_OUT_1.csv","TESLA_OUT_2.csv","TESLA_OUT_3.csv","TESLA_OUT_4.csv"])
	basenames = [os.path.basename(name) for name in filelist]
	assert all(requiredFiles.isin(basenames)), "All four submission file must be present and submission files must be named TESLA_OUT_1.csv, TESLA_OUT_2.csv, TESLA_OUT_3.csv, or TESLA_OUT_4.csv"
	for filepath in filelist: 
		validation_func[os.path.basename(filepath)](filepath)
	order = pd.np.argsort(basenames)
	validate_VAR_ID(filelist[order[0]],filelist[order[1]],filelist[order[2]])
	validate_STEP_ID(filelist[order[2]],filelist[order[3]])
	return(True, "Passed Validation!")

def perform_validate(args):
	validate_files(args.file)
	print("Passed Validation")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Validate GENIE files')
	parser.add_argument("file", type=str, nargs="+",
						help='path to TESLA files')
	args = parser.parse_args()

	perform_validate(args)