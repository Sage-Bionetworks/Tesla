import synapseclient
import pandas as pd
syn = synapseclient.login()

def addAnnotations(x, data_round):
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

metadataPath = syn.get("syn8303410").path
metadata = pd.read_csv(metadataPath)
metadata['qcFileName'] = metadata['qcFileName'].fillna('')
metadata['checkpointInhibitor'] = metadata['checkpointInhibitor'].fillna('')
metadata['classIHLAalleles'] = metadata['classIHLAalleles'].fillna('')
metadata['collectionDate'] = metadata['collectionDate'].fillna('')
metadata['isTreated'] = metadata['isTreated'].fillna('')
metadata['organ'] = metadata['organ'].fillna('')
metadata['sex'] = metadata['sex'].fillna('')
metadata['tumorPurity(percent) '] = metadata['tumorPurity(percent) '].fillna('')
data_round = 1

metadata.apply(lambda x: addAnnotations(x, data_round), axis=1)





