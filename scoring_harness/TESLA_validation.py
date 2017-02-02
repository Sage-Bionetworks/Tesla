import pandas as pd
import os
import re
import httplib2 as http
import json

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

from multiprocessing.dummy import Pool as ThreadPool 

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
            return(False, None)
        else:
            #print(len(data['response']['docs']))
            mapped = [symbol['symbol'] for symbol in data['response']['docs']]
            return(True, mapped)
    else:
        return(False, None)

# Validation of gene names
def validateSymbol(gene, returnMapping=False):
    """
    This function does validation of symbols

    :params gene:               Gene symbol
    :params returnMapping:      Return mapping of old gene to new gene

    :returns:                   Check if the provided gene name is a correct symbol and print out genes 
                                that need to be remapped or can't be mapped to anything
    """
    path = '/fetch/symbol/%s' %  gene
    verified, symbol = hgncRestCall(path)
    if not verified:
        path = '/fetch/prev_symbol/%s' %  gene
        verified, symbol = hgncRestCall(path)
    if not verified:
        path = '/fetch/alias_symbol/%s' %  gene
        verified, symbol = hgncRestCall(path)       
    if gene in symbol:
        return(True)
    else:
        if symbol is None:
            print("%s cannot be remapped. Please correct." % gene)
        else:
            #if "MLL4", then the HUGO symbol should be KMT2D and KMT2B
            if len(symbol) > 1:
                print("%s can be mapped to different symbols: %s. Please correct." % (gene, ", ".join(symbol)))
            else:
                print("%s will be remapped to %s" % (gene, symbol[0]))
                if returnMapping:
                    return({gene, symbol[0]})
                else:
                    return(True)
        return(False)

def checkType(submission, cols, colType):
	for col in cols:
		assert all(submission[col].apply(lambda x: isinstance(x, colType))), "All values in %s column must be type: %s" % (col, re.sub(".+['](.+)['].+","\\1",str(colType)))

def validate_1(submission_filepath):
	#CHECK: Correct filename
	assert os.path.basename(submission_filepath) == "TESLA_OUT_1.csv", "Submission file must be named TESLA_OUT_1.csv"
	required_cols = pd.Series(["VAR_ID","GENE","CHROM","START","END","STRAND","CLASS","REFALLELE","MUTALLELE","REFCOUNT_T","MUTCOUNT_T",
					 "REFCOUNT_N","MUTDEPTH_N","REF_FPKM_T","REF_HTSEQ_T","MUT_FPKM_T","MUT_HTSEQ_T","QUAL","VARID","INFO"])
	integer_cols = ['VAR_ID','START','END','REFCOUNT_T','MUTCOUNT_T','REFCOUNT_N','MUTDEPTH_N','REF_HTSEQ_T','MUT_HTSEQ_T','QUAL']
	string_cols = ['REFALLELE','MUTALLELE','VARID','INFO']
	float_cols = ['MUT_FPKM_T','REF_FPKM_T']
	CLASS_categories = ['intron','missense','silent','splice_site','in_frame_deletion','in_frame_insertion',
						'frame_shift_deletion','frame_shift_insertion',
						'nonsense_mutation','structural_variant','splice_isoform','other']
	submission = pd.read_csv(submission_filepath)
	#CHECK: Required headers must exist in submission
	assert all(required_cols.isin(submission.columns)), "These column headers are missing from your file: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
	#CHECK: CHROM must be 1-22 or X
	assert all(submission.CHROM.isin(range(1,23) + ["X"])), "CHROM values must be 1-22, or X. You have: %s" % ", ".join(set(submission.CHROM[~submission.CHROM.isin(range(1,23) + ["X"])])) 
	#CHECK: CLASS must be in CLASS_categories
	assert all(submission.CLASS.isin(CLASS_categories)), "CLASS values must be one these: %s.  You have: %s" % (", ".join(CLASS_categories),", ".join(set(submission.CLASS[~submission.CLASS.isin(CLASS_categories)])))
	#CHECK: STRAND must be +/-
	assert all(submission.STRAND.isin(['+','-'])), "STRAND values must be + or -.  You have: %s" %(", ".join(set(submission.STRAND[~submission.STRAND.isin(['+','-'])])))
	#CHECK: gene validation
    assert all(pool.map(validateSymbol, submission.GENE.drop_duplicates())), "All gene symbols in GENE column must be correct"

	#CHECK: integer, string and float columns are correct types
	checkType(submission, integer_cols, int)
	checkType(submission, string_cols, str)
	checkType(submission, float_cols, float)

    return(True,"Passed Validation")


def validate_2(submission_filepath, goldstandard_path):
	assert os.path.basename(submission_filepath) == "TESLA_OUT_2.csv", "Submission file must be named TESLA_OUT_1.csv"
	required_cols = pd.Series(["RANK","VAR_ID","PROT_POS","HLA_ALLELE","PEP_LEN","MUT_EPI_SEQ","WT_EPI_SEQ"])
	submission = pd.read_csv(submission_filepath)
	assert required_cols.isin(submission.columns), "These column headers are missing from your file: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])
    return(True,"Passed Validation")


def validate_3(submission_filepath, goldstandard_path):
	assert os.path.basename(submission_filepath) == "TESLA_OUT_3.csv", "Submission file must be named TESLA_OUT_1.csv"
	required_cols = pd.Series(["VAR_ID","PROT_POS","HLA_ALLELE","PEP_LEN","MUT_EPI_SEQ","WT_EPI_SEQ","STEP_ID"])
	submission = pd.read_csv(submission_filepath)
	assert required_cols.isin(submission.columns), "These column headers are missing from your file: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])

    return(True,"Passed Validation")

#Validate workflow
def validate_4(submission_filepath, goldstandard_path):
	assert os.path.basename(submission_filepath) == "TESLA_OUT_4.csv", "Submission file must be named TESLA_OUT_4.csv"
    return(True,"Passed Validation")