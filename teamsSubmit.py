temp = syn.getSubmissionBundles(syn.getEvaluation(8116290))

submitted = set()
invalidteams = set()
validteams = set()
for i, stat in temp:
	getTeam = filter(lambda x: x.get("key")=="team", stat.annotations['stringAnnos'])
	if len(getTeam) > 0:
		submitted.add(getTeam[0]['value'])
		if stat.status == "INVALID":
			invalidteams.add(getTeam[0]['value'])
		else:
			validteams.add(getTeam[0]['value'])

teams = syn.tableQuery('SELECT * FROM syn8220615')
teamsDf = teams.asDataFrame()


temp = syn.getSubmissionBundles(syn.getEvaluation(8116290), status=  "VALIDATED")

zippedSubmissions = dict()
for i, stat in temp:
	getTeam = filter(lambda x: x.get("key")=="team", stat.annotations['stringAnnos'])[0]['value']
	entity = syn.getSubmission(i,downloadFile=False).entity
	#if entity.name.endswith(".zip"):
	if zippedSubmissions.get(getTeam) is None:
		zippedSubmissions[getTeam] = [entity.name]
	else:
		zippedSubmissions[getTeam].append(entity.name)

for key in zippedSubmissions:
	zippedSubmissions[key] = set(zippedSubmissions[key])

# patients = dict()

# for key in zippedSubmissions:
# 	for i in zippedSubmissions[key]:
# 		if patients.get(i) is None:
# 			patients[i] = 1
# 		else:
# 			patients[i] = patients[i]+1