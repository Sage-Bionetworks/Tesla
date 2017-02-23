# Automation of validation and scoring
# Make sure you point to the directory where challenge.py belongs and a log directory must exist for the output
cd ~/Tesla/scoring_harness
#---------------------
#Validate submissions
#---------------------
#Remove --send-messages to do rescoring without sending emails to participants
python challenge.py -u thomas.yu --send-messages --notifications --acknowledge-receipt validate --all >> ~/log/score.log 2>&1

#--------------------
#Score submissions
#--------------------
#python challenge.py -u thomas.yu --send-messages --notifications score --all >> log/score.log 2>&1