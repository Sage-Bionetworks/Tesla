# Make sure you create a directory where you want to keep your log files in this case (~/log)
cd ~/log && mv score.log score`date +"%Y_%m_%d"`.log && touch score.log

python ~/Tesla/createSubmissionQueueStatDf.py 8116290 --databaseSynId syn104076944
python ~/Tesla-Analysis/Scripts/synapse-IO/tesla-submissions.py