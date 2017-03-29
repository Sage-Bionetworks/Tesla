library(synapseClient)
synapseLogin()
DATE_START = as.Date("2017-02-20")
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
        }
      }
      c(x$createdOn, team, name, patientId)
    })
    OFFSET = OFFSET + LIMIT
    challenge_stats = append(challenge_stats, sub_stats)
  }
  challenge_stats_df = data.frame(do.call(rbind, challenge_stats))
  colnames(challenge_stats_df) <- c("time","team","fileName", "patientId")
  challenge_stats_df$Date = as.Date(challenge_stats_df$time)
  challenge_stats_df
}

submissionsPerWeek = function(challenge_stats_df, patientId, challengeSynId) {
  weeks = seq(DATE_START, Sys.Date()+6, "weeks")
  numSubs = sapply(weeks, function(x) x= 0)
  names(numSubs) <- weeks
  for (i in as.character(weeks)) {
    start = as.Date(i)
    for (y in challenge_stats_df$Date) {
      if (y > start && y < start+6) {
        numSubs[i]= numSubs[i]+1
      }
    }
  }
  png(sprintf("%s_submissions.png", patientId),width = 600, height = 600)
  par(mar=c(8,4,4,2)+0.1)
  barplot(numSubs, main=sprintf("Number of Complete Submissions Per Week for patient: %s",patientId), ylab="Number of Submissions", las=3, ylim=c(0,nrow(challenge_stats_df)+1))
  mtext("Date", side=1, line=6)
  dev.off()
  synStore(File(sprintf("%s_submissions.png", patientId),parentId = challengeSynId))
}

numTeamsOverTime = function(challenge_stats_df, challengeSynId) {
  weeks <- seq(min(challenge_stats_df$Date), max(challenge_stats_df$Date)+6, "weeks")
  weekSegment = sapply(challenge_stats_df$Date, function(x) {
    as.character(tail(weeks[x >= weeks],n=1))
  })
  submissions <- table(as.Date(weekSegment),challenge_stats_df$team)
  for (i in seq(2, nrow(submissions))) {
    submissions[i,] = submissions[i-1,] + submissions[i,]
  }
  numberOfTeams = rowSums(submissions > 0)
  dates = as.Date(names(numberOfTeams))
  png("totalTeamsSubmitted.png",width=600, height=400)
  plot(dates,numberOfTeams, xaxt="n",xlab = "Dates",ylab = "Number of Teams",main="Cumulative Number of Teams Submitted",ylim = c(0, max(numberOfTeams)),type = "l")
  axis.Date(1, at = seq(DATE_START, Sys.Date()+6, "weeks"))
  dev.off()
  synStore(File("totalTeamsSubmitted.png",parentId = challengeSynId))
}

plotStats <- function(patientId, challenge_stats_df) {
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
  submissionsPerWeek(noDups, patientId, "syn7801079")
}

#TESLA STATS
challenge_stats_df = getSubmissionCount(8116290, "VALIDATED")
metadata = synGet("syn8371011")
metadataDf = read.csv(getFileLocation(metadata))
for (i in unique(metadataDf$patientId)) {
  plotStats(i, challenge_stats_df)
}
#Get cumulative teams submitted over time
numTeamsOverTime(challenge_stats_df, "syn7801079")


