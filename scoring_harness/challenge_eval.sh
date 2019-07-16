# Automation of validation and scoring
# Make sure you point to the directory where challenge.py belongs and a log directory must exist for the output
cd ~/Tesla/scoring_harness
#---------------------
#Validate submissions
#---------------------
#Remove --send-messages to do rescoring without sending emails to participants
python3 challenge.py --send-messages --notifications --acknowledge-receipt validate 9614265 >> ~/log/score.log 2>&1

#--------------------
#Score submissions
#--------------------
#python3 challenge.py --send-messages --notifications score 9614265 >> ~/log/score.log 2>&1
