import synapseclient
import pandas as pd


def checkIfTesla(newMember):
	teams = syn.tableQuery('SELECT * FROM syn8220615')
	teamNames = teams.asDataFrame()
	total_members = set()
	for i in teamNames['realTeam']:
		members = syn.getTeamMembers(syn.getTeam(i))
		total_members.update(set([member['member']['userName'] for member in members]))
	tesla_part = syn.getTeamMembers(syn.getTeam(3346701))
	tesla_mem = [i['member']['userName'] for i in tesla_part]
	print([i for i in total_members if i not in tesla_mem])

def setPermissions():
	teams = syn.tableQuery('SELECT * FROM syn8220615')
	teamNames = teams.asDataFrame()
	total_members = set()
	for i in teamNames['realTeam']:
		team = syn.getTeam(i)
		perms = syn.setPermissions("syn7362874", team['id'], accessType=[u'READ'], warn_if_inherits=False, overwrite=False)


