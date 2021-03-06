library(synapseClient)
synapseLogin()
#Change DATE_START
#Update the roundNum to just be space when doing new rounds
#Change the roundNum to round%d when end of a round to save the plots
DATE_START = as.Date("2017-7-07")
getSubmissionCount = function(evalId, status) {
  LIMIT = 10
  OFFSET = 0
  subs = synGetSubmissions(evalId,status = status)
  challenge_stats = list()
  while (OFFSET < subs@totalNumberOfResults) {
    subs = synGetSubmissions(evalId,status = status, limit=LIMIT, offset = OFFSET)
    sub_stats = lapply(subs@results, function(x) {
      stat = synGetSubmissionStatus(x$id)
      for (i in stat@annotations@stringAnnos@content) {
        if (i$key == "team") {
          team = i$value
        } else if (i$key == "submissionName") {
          name = i$value 
        } else if (i$key == "patientId") {
          patientId = i$value
        } else if (i$key == "round") {
          rounds = i$value
        }
      }
      if (!exists("team")) {
        team = "None"
      }
      if (!exists("name")) {
        name = "None"
      }
      if (!exists("patientId")) {
        patientId = "None"
      }
      if (!exists("rounds")) {
        rounds = "None"
      }
      c(x$id, x$createdOn, team, name, patientId, rounds)
    })
    OFFSET = OFFSET + LIMIT
    challenge_stats = append(challenge_stats, sub_stats)
  }
  challenge_stats_df = data.frame(do.call(rbind, challenge_stats))
  colnames(challenge_stats_df) <- c("ID","time","team","fileName", "patientId","round")
  challenge_stats_df$Date = as.Date(challenge_stats_df$time)
  challenge_stats_df
}

submissionsPerWeek = function(challenge_stats_df, patientId, challengeSynId, roundNum) {
  weeks = seq(DATE_START, Sys.Date(), "weeks")
  numSubs = sapply(weeks, function(x) x= 0)
  print(challenge_stats_df)
  names(numSubs) <- weeks
  for (i in as.character(weeks)) {
    start = as.Date(i)
    for (y in as.character(challenge_stats_df$Date)) {
      yDate <- as.Date(y)
      if (yDate >= start && yDate < start+7) {
        numSubs[i]= numSubs[i]+1
      }
    }
  }
  for (i in seq_along(numSubs)[-1]-1) {
    numSubs[i+1] = numSubs[i] + numSubs[i+1]
  }
  png(sprintf("%s%s_submissions.png", roundNum, patientId),width = 600, height = 600)
  par(mar=c(8,4,4,2)+0.1)
  barplot(numSubs, main=sprintf("Number of Complete Submissions for patient: %s",patientId), ylab="Number of Submissions", las=3, ylim=c(0,25))
  mtext("Date", side=1, line=6)
  dev.off()
  synStore(File(sprintf("%s%s_submissions.png", roundNum, patientId),parentId = challengeSynId))
}

numTeamsOverTime = function(challenge_stats_df, challengeSynId, roundNum) {
  #weeks <- seq(min(challenge_stats_df$Date), max(challenge_stats_df$Date)+6, "weeks")
  if (max(challenge_stats_df$Date) > DATE_START) {
    weeks <- seq(DATE_START, Sys.Date()+6, "weeks")
    weekSegment = sapply(challenge_stats_df$Date, function(x) {
      corDate = as.character(tail(weeks[x >= weeks],n=1))
      if (length(corDate) >0) {
        corDate
      } else {
        NA
      }
    })
    dontKeep = is.na(weekSegment)
    submissions <- table(as.Date(unlist(weekSegment))[!dontKeep],challenge_stats_df$team[!dontKeep])
    for (i in seq(2, nrow(submissions))) {
      submissions[i,] = submissions[i-1,] + submissions[i,]
    }
    numberOfTeams = rowSums(submissions > 0)
    dates = as.Date(names(numberOfTeams))
    png(sprintf("%stotalTeamsSubmitted.png",roundNum),width=600, height=400)
    plot(dates,numberOfTeams, xaxt="n",xlab = "Dates",ylab = "Number of Teams",main="Cumulative Number of Teams Submitted",ylim = c(0, max(numberOfTeams)),type = "l")
    axis.Date(1, at = seq(min(dates), Sys.Date()+6, "weeks"))
    dev.off()
    synStore(File(sprintf("%stotalTeamsSubmitted.png",roundNum),parentId = challengeSynId))
  }
}

plotStats <- function(patientId, challenge_stats_df, roundNum) {
  allFiles = c(sprintf("%s_RNA_T.bam",patientId),
               sprintf("%s_EXOME_N.bam",patientId),
               sprintf("%s_EXOME_T.bam",patientId),
               sprintf("%s.zip", patientId))
  
  temp = challenge_stats_df[challenge_stats_df$patientId == patientId,]
  teams = as.character(unique(temp$team))
  submitAll = sapply(teams, function(x) {
    all(allFiles %in% temp$fileName[temp$team == x])
  })
  if (length(submitAll) == 0) {
    noDups = data.frame()
  } else {
    submitted = challenge_stats_df[challenge_stats_df$team %in% names(submitAll[submitAll]),]
    submitted = submitted[order(submitted$Date,decreasing = F),]
    noDups = submitted[!duplicated(submitted$team, fromLast = T),]
  }
  submissionsPerWeek(noDups, patientId, "syn7801079", roundNum=roundNum)
}

#TESLA STATS
challenge_stats_df = getSubmissionCount(8116290, "VALIDATED")
metadata = synTableQuery('SELECT * FROM syn8292741 where round = "2"')
metadataDf = metadata@values

for (i in unique(metadataDf$patientId[!is.na(metadataDf$patientId)])) {
  #plotStats(i, challenge_stats_df,"round1_")
  plotStats(i, challenge_stats_df, "")
}
#Get cumulative teams submitted over time
#numTeamsOverTime(challenge_stats_df, "syn7801079",roundNum='round1_')
numTeamsOverTime(challenge_stats_df, "syn7801079",roundNum='round2_')


## STATS
#only look at invalid only
#Regenerate stats for teams that have submitted
#How many of the teams submitted
#which teams have valid submissions
#which teams have invalid submissions
#teams that submitted 2/4
challenge_stats_df = getSubmissionCount(8116290, "VALIDATED")
challenge_invalid_df = getSubmissionCount(8116290, "INVALID")

round2_valid = challenge_stats_df[challenge_stats_df$round == 2,]
zippedFiles = apply(round2_valid, 1, function(x) {
  if (endsWith(x['fileName'], ".zip")) {
    return(x[c('ID','team','fileName')])
  }
})
allZipped = do.call(rbind,zippedFiles)
validZippedFiles = table(allZipped[,"team"],allZipped[,'fileName'])
write.csv(validZippedFiles,"round2_valid_zipped_files.csv")

hasTESLA2 = c()
for (i in allZipped[,"ID"]) {
  sub = synGetSubmission(i)
  foo = unzip(sub@filePath,list=T)
  if (sum(foo$Name == "TESLA_OUT_2.csv") == 1) {
    hasTESLA2 = c(hasTESLA2, T)
  } else{
    hasTESLA2 = c(hasTESLA2, F)
  }
}
VCF_preds_zipped = allZipped[hasTESLA2,]
validVCFZippedFiles = table(VCF_preds_zipped[,"team"],VCF_preds_zipped[,'fileName'])

write.csv(validVCFZippedFiles,"round2_valid_VCF_zipped_files.csv")

bamFiles = apply(round2_valid, 1, function(x) {
  if (endsWith(x['fileName'], ".bam")) {
    return(x[c('ID','team','fileName')])
  }
})
allBams = do.call(rbind,bamFiles)
bamFiles = table(allBams[,"team"],allBams[,'fileName'])
write.csv(bamFiles,"round2_bam_files.csv")
challenge_stats_df

database = synTableQuery('SELECT * FROM syn10407694')
database@values$submissionId <- challenge_stats_df$ID
database@values$team <- challenge_stats_df$team
database@values$fileName <- challenge_stats_df$fileName
database@values$patientId <- challenge_stats_df$patientId
database@values$round <- challenge_stats_df$round
database@values$dateTime <- challenge_stats_df$time
