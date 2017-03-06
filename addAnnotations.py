import synapseclient
import pandas as pd
import argparse
import math

def addAnnotationHelper(x, data_round):
	result = syn.query('select id from file where parentId =="%s" and name == "%s"' %(x['uploadAccount'],x['teslaName']))
	entityId = result['results'][0]['file.id']
	print(entityId)
	ent = syn.get(entityId,downloadFile=False)
	annotations = x.to_dict()
	annotations.pop('teslaName')
	annotations.pop('uploadAccount')
	annotations.pop('rawName')
	toRemove = []
	for key in annotations:
		if not isinstance(annotations[key],str):
			if math.isnan(annotations[key]) or annotations[key] == "Not Applicable ":
				toRemove.append(key)
	for key in toRemove:
		annotations.pop(key)
	if annotations.get('tumorPurity(percent) ') is not None:
		annotations['tumorPurityPercent'] = annotations.pop('tumorPurity(percent) ')
	ent.annotations = annotations
	ent.round = data_round
	syn.store(ent)
	return(ent.id)

def addAnnotation(syn, data_round):
	metadataPath = syn.get("syn8371011").path
	metadata = pd.read_csv(metadataPath)
	metadata['exomePulldownFile'][metadata['exomePulldownFile'] == "TESLA_EXOME_REGIONS.bed.gz"] = "syn8313637"
	metadata['exomePulldownFile'][metadata['exomePulldownFile'] == "TESLA_EXOME_REGIONS2.bed"] = "syn8348425"
	metadata.apply(lambda x: addAnnotationHelper(x, data_round), axis=1)

def perform_add(syn, args):
	addAnnotation(syn, args.round)
	print("Annotations updated")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Validate TESLA files per sample')
	parser.add_argument("round", type=int,
						help='Round of data')
	args = parser.parse_args()
	syn = synapseclient.login()
	perform_add(syn, args)




