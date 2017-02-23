import synapseclient
import pandas as pd


def checkIfTesla(newMember):
	teams = syn.tableQuery('SELECT * FROM syn8220615')
	teamNames = teams.asDataFrame()
	total_members = set()
	for i in teamNames['realTeam']:
		members = syn.getTeamMembers(syn.getTeam(i))
		total_members.update(set([member['member']['userName'] for member in members]))
	print(newMember in  total_members)