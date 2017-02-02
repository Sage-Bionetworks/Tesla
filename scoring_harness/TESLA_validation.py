import pandas as pd

def validate_1(submission_filepath, goldstandard_path):
	assert os.path.basename(submission_filepath) == "TESLA_OUT_1.csv", "Submission file must be named TESLA_OUT_1.csv"
	required_cols = pd.Series(["VAR_ID","GENE","CHROM","START","END","STRAND","CLASS","REFALLELE","MUTALLELE","REFCOUNT_T","MUTCOUNT_T",
					 "REFCOUNT_N","MUTDEPTH_N","REF_FPKM_T","REF_HTSEQ_T","MUT_FPKM_T","MUT_HTSEQ_T","QUAL","VARID","INFO"])
	submission = pd.read_csv(submission_filepath)
	assert required_cols.isin(submission.columns), "These column headers are missing from your file: %s" % ", ".join(required_cols[~required_cols.isin(submission.columns)])

	assert submission.CHROM.isin(range(1,23) + ["X"]), "CHROM values must be 1-22, or X. You have: %s" % ", ",join(submission.CHROM[~submission.CHROM.isin(range(1,23) + ["X"])]) 

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