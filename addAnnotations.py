import synapseclient
import pandas as pd
import argparse



def addAnnotationHelper(x, data_round):
	result = syn.query('select id from file where parentId =="%s" and name == "%s"' %(x['uploadAccount'],x['teslaName']))
	entityId = result['results'][0]['file.id']
	print(entityId)
	ent = syn.get(entityId,downloadFile=False)
	annotations = x.to_dict()
	annotations.pop('teslaName')
	annotations.pop('uploadAccount')
	annotations.pop('rawName')
	annotations['tumorPurityPercent'] = str(annotations.pop('tumorPurity(percent) '))
	ent.annotations = annotations
	ent.round = data_round
	syn.store(ent)
	return(ent.id)

def addAnnotation(syn, data_round):
	metadataPath = syn.get("syn8303410").path
	metadata = pd.read_csv(metadataPath)
	metadata['qcFileName'] = metadata['qcFileName'].fillna('')
	metadata['checkpointInhibitor'] = metadata['checkpointInhibitor'].fillna('')
	metadata['classIHLAalleles'] = metadata['classIHLAalleles'].fillna('')
	metadata['isTreated'] = metadata['isTreated'].fillna('')
	metadata['organ'] = metadata['organ'].fillna('')
	metadata['sex'] = metadata['sex'].fillna('')
	metadata['tumorPurity(percent) '] = metadata['tumorPurity(percent) '].fillna('')
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




