##-----------------------------------------------------------------------------
##
## challenge specific code and configuration
##
##-----------------------------------------------------------------------------

import TESLA_validation as TESLA_val
import zipfile
import os
import re
## A Synapse project will hold the assetts for your challenge. Put its
## synapse ID here, for example
CHALLENGE_SYN_ID = "syn7362874"

## Name of your challenge, defaults to the name of the challenge's project
CHALLENGE_NAME = "TESLA_consortium"

## Synapse user IDs of the challenge admins who will be notified by email
## about errors in the scoring script
ADMIN_USER_IDS = [3324230]



evaluation_queues = [
    {
        'id':8116290,
        'validation_func':TESLA_val.validate_files
    }
]
evaluation_queue_by_id = {q['id']:q for q in evaluation_queues}


## define the default set of columns that will make up the leaderboard
LEADERBOARD_COLUMNS = [
    dict(name='objectId',      display_name='ID',      columnType='STRING', maximumSize=20),
    dict(name='userId',        display_name='User',    columnType='STRING', maximumSize=20, renderer='userid'),
    dict(name='entityId',      display_name='Entity',  columnType='STRING', maximumSize=20, renderer='synapseid'),
    dict(name='versionNumber', display_name='Version', columnType='INTEGER'),
    dict(name='name',          display_name='Name',    columnType='STRING', maximumSize=240),
    dict(name='team',          display_name='Team',    columnType='STRING', maximumSize=240)]

## Here we're adding columns for the output of our scoring functions, score,
## rmse and auc to the basic leaderboard information. In general, different
## questions would typically have different scoring metrics.
leaderboard_columns = {}
for q in evaluation_queues:
    leaderboard_columns[q['id']] = LEADERBOARD_COLUMNS + [
        dict(name='score',         display_name='Score',   columnType='DOUBLE'),
        dict(name='rmse',          display_name='RMSE',    columnType='DOUBLE'),
        dict(name='auc',           display_name='AUC',     columnType='DOUBLE')]

## map each evaluation queues to the synapse ID of a table object
## where the table holds a leaderboard for that question
leaderboard_tables = {}


def validate_submission(syn, evaluation, submission, team_mapping, patientIds):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    assert 'teamId' in submission, "Must submit as part of a team and not as an individual"
    team = syn.getTeam(submission.teamId)
    teamIndex = team_mapping['realTeam'] == team['name']
    assert sum(teamIndex) == 1, "Must submit as one of these teams: %s" % ", ".join(team_mapping['realTeam'])
    teamDict = {'team':team_mapping['alias'][teamIndex].values[0]}
    submissionName = submission.entity.name
    if not submission.entity.name.endswith("bam"):
        submission = syn.getSubmission(submission.id)

    patientId = re.sub("(\d+).+","\\1",submissionName)
    assert patientId.isdigit(), "Wrong filenaming convention"
    assert int(patientId) in patientIds, "Patient Id must be part of the Id list"
    if submissionName.endswith(".zip"):
        assert submissionName == "%s.zip" % patientId, "Zip file must be named patientId.zip"
        #Unzip files here
        dirname = submission.entity.cacheDir
        try:
            zfile = zipfile.ZipFile(submission.filePath)
        except zipfile.BadZipfile as e:
            raise AssertionError("Must submit a zipped file containing TESLA_OUT_1.csv, TESLA_OUT_2.csv, TESLA_OUT_3.csv, and TESLA_OUT_4.csv")

        for name in zfile.namelist():
          zfile.extract(name, dirname)

        tesla_out_1 = os.path.join(dirname,'TESLA_OUT_1.csv')
        tesla_out_2 = os.path.join(dirname,'TESLA_OUT_2.csv')
        tesla_out_3 = os.path.join(dirname,'TESLA_OUT_3.csv')
        tesla_out_4 = os.path.join(dirname,'TESLA_OUT_4.csv')
        tesla_vcf = os.path.join(dirname,'TESLA_VCF.vcf')

        filelist = [tesla_out_1,tesla_out_2,tesla_out_3,tesla_out_4,tesla_vcf]
        assert all([os.path.exists(i) for i in filelist]), "TESLA_OUT_1.csv, TESLA_OUT_2.csv, TESLA_OUT_3.csv, TESLA_OUT_4.csv, and TESLA_VCF.vcf must all be in the zipped file"
        TESLA_val.validate_files(filelist,patientId,validatingBAM=False)
    else:
        assert submissionName in ["%s_EXOME_N.bam" % patientId,"%s_EXOME_T.bam" % patientId,"%s_RNA_T.bam" % patientId], "Bam files must be named patientId_EXOME_N.bam, patientId_EXOME_T.bam or patientId_RNA_T.bam"
    teamDict['patientId'] = patientId
    return True, "Validation passed!", teamDict


def score_submission(evaluation, submission):
    """
    Find the right scoring function and score the submission

    :returns: (score, message) where score is a dict of stats and message
              is text for display to user
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    score = config['scoring_func'](submission, config['goldstandard_path'])
    #Make sure to round results to 3 or 4 digits
    return (dict(score=round(score[0],4), rmse=score[1], auc=score[2]), "You did fine!")


