##-----------------------------------------------------------------------------
##
## challenge specific code and configuration
##
##-----------------------------------------------------------------------------

import TESLA_validation as TESLA_val

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
        'goldstandard_path':'path/to/sc1gold.txt'
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


def validate_submission(syn, evaluation, submission):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    validated, validation_message = config['validation_func'](submission, config['goldstandard_path'])

    if 'teamId' in submission:
        team = syn.getTeam(submission.teamId)
        if 'name' in team:
            teamDict = {'team':team['name']}
        else:
            teamDict = {'team':submission.teamId}
    else:
        raise AssertionError("Must submit as part of a team and not as an individual")

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

    filelist = [tesla_out_1,tesla_out_2,tesla_out_3,tesla_out_4]
    assert all([os.path.exists(i) for i in filelist]), "TESLA_OUT_1.csv, TESLA_OUT_2.csv, TESLA_OUT_3.csv, and TESLA_OUT_4.csv must all be in the zipped file"
    TESLA_val.validate_files(filelist)
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


