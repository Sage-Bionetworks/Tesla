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
	#annotations['pairedEndId'] = str(annotations['pairedEndId'])
	ent.annotations = annotations
	ent.round = data_round
	syn.store(ent)
	return(ent.id)

metadataPath = syn.get("syn8290709").path
metadata = pd.read_excel(metadataPath)
metadata['qcFileName'] = metadata['qcFileName'].fillna('')
metadata['checkpointInhibitor'] = metadata['checkpointInhibitor'].fillna('')

data_round = 1

metadata.apply(lambda x: addAnnotations(x, data_round), axis=1)





