##-----------------------------------------------------------------------------
##
## challenge specific code and configuration
##
##-----------------------------------------------------------------------------

import TESLA_validation as TESLA_val
import zipfile
import os
import sys
import re
import operator
import pandas as pd
try:
    reduce
except NameError:
    from functools import reduce
## A Synapse project will hold the assetts for your challenge. Put its
## synapse ID here, for example
CHALLENGE_SYN_ID = "syn7362874"

## Name of your challenge, defaults to the name of the challenge's project
CHALLENGE_NAME = "TESLA_consortium"

## Synapse user IDs of the challenge admins who will be notified by email
## about errors in the scoring script
ADMIN_USER_IDS = [3360851]


CONFIG_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))



evaluation_queues = [
    {
        'id': 9614266,
        'validation_func': TESLA_val.validate_files,
        'patients': [str(i) for i in [1,2,3,4,5,10,11,12,14,15,16]]
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
        dict(name='score',         display_name='Score',   columnType='DOUBLE')]

## map each evaluation queues to the synapse ID of a table object
## where the table holds a leaderboard for that question
leaderboard_tables = {}

def validate_teamname(syn, evaluation, submission, team_mapping):
    assert 'teamId' in submission, "Must submit as part of a team and not as an individual"
    team = syn.getTeam(submission.teamId)
    teamIndex = team_mapping['realTeam'] == team['name']
    assert sum(teamIndex) == 1, "Must submit as: %s" % ",".join(team_mapping['realTeam'])
    teamDict = {'team':team_mapping['alias'][teamIndex].values[0]}
    return True, "Validation passed!", teamDict

def validate_submission(syn, evaluation, submission, patientIds, HLA):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    submissionName = submission.entity.name
    if not submission.entity.name.endswith("bam"):
        submission = syn.getSubmission(submission.id)

    config = evaluation_queue_by_id[int(evaluation.id)]
    allowed_patients = config['patients']
    patientId = re.sub("(\d+).+", "\\1", submissionName)
    assert patientId.isdigit(), "Wrong filenaming convention"
    assert int(patientId) in patientIds, "Patient Id must be part of the Id list"
    assert submissionName.endswith(".zip"), "Must submit a zip file"
    assert patientId in allowed_patients, (
        "Patient: " +
        str(patientId) +
        " not allowed in submission queue: " +
        str(evaluation.id) +
        ". Only patients: " +
        ", ".join(allowed_patients) +
        " are allowed.")
    hasVCF = False
    #if submissionName.endswith(".zip"):
    assert submissionName == "%s.zip" % patientId, "Zip file must be named patientId.zip"
    #Unzip files here
    dirname = os.path.dirname(submission.entity.path)
    try:
        zfile = zipfile.ZipFile(submission.filePath)
    except zipfile.BadZipfile as e:
        raise AssertionError("Must submit a zipped file containing TESLA_OUT_1.csv,TESLA_OUT_3.csv, TESLA_YAML.yaml and TESLA_VCF.vcf")

    for name in zfile.namelist():
      zfile.extract(name, dirname)

    tesla_out_1 = os.path.join(dirname,'TESLA_OUT_1.csv')
    tesla_out_2 = os.path.join(dirname,'TESLA_OUT_2.csv')
    tesla_out_3 = os.path.join(dirname,'TESLA_OUT_3.csv')
    tesla_out_4 = os.path.join(dirname,'TESLA_OUT_4.csv')
    tesla_yaml = os.path.join(dirname,'TESLA_YAML.yaml')
    tesla_ranking = os.path.join(dirname, 'TESLA_ranking_method.txt')
    tesla_vcf = os.path.join(dirname,'TESLA_VCF.vcf')
    filelist = [tesla_out_1,tesla_out_3,tesla_yaml,tesla_vcf]
    optionalFiles = [tesla_out_2,tesla_out_4]
    optionalExists = all([os.path.exists(i) for i in optionalFiles])
    assert all([os.path.exists(i) for i in filelist]), "TESLA_OUT_1.csv, TESLA_OUT_3.csv, TESLA_VCF.vcf, and TESLA_YAML.yaml must all be in the zipped file.  Please do NOT put your files in a folder prior to zipping them."
    assert optionalExists or sum([os.path.exists(i) for i in optionalFiles]) == 0, "TESLA_OUT_2.csv, TESLA_OUT_4.csv.  Both files MUST either be present or missing.  If missing, you are missing predictions from VCF. If this is not as intended, please submit again."
    if optionalExists:
        filelist.extend(optionalFiles)
    if os.path.exists(tesla_ranking):
        filelist.append(tesla_ranking)
    listHLA = HLA['classIHLAalleles'][HLA['patientId'] == int(patientId)]
    validHLA = [i.split(";") for i in listHLA]
    validHLA = reduce(operator.add, validHLA)
    validHLA = set([i.split("(")[0] for i in validHLA])
    validated, hasVCF, message = TESLA_val.validate_files(syn, filelist,patientId,validHLA,validatingBAM=False)
    # else:
    #     assert submissionName in ["%s_EXOME_N.bam" % patientId,"%s_EXOME_T.bam" % patientId,"%s_RNA_T.bam" % patientId], "Bam files must be named patientId_EXOME_N.bam, patientId_EXOME_T.bam or patientId_RNA_T.bam"
    teamDict = {'patientId':patientId,
                'hasVCF':hasVCF}
    return True, "Validation passed!", teamDict

