install.packages("RMySQL")
library(RMySQL)
library(synapseClient)
synapseLogin()
mydb = dbConnect(MySQL(), user='tyu', password="sAGE()dw!", host='warehouse.c95bbsvwbjlu.us-east-1.rds.amazonaws.com')

downloadReport1 = dbSendQuery(mydb, "SELECT PAR.ENTITY_ID, AR.USER_ID ,COUNT(*)
                                     FROM warehouse.PROCESSED_ACCESS_RECORD PAR,
                                     warehouse.ACCESS_RECORD AR,
                                     (SELECT ID, CREATED_BY FROM warehouse.NODE_SNAPSHOT WHERE (ID, TIMESTAMP) IN (SELECT ID, TIMESTAMP FROM tyu.syn7362874)) AS NODE
                                     WHERE PAR.SESSION_ID = AR.SESSION_ID
                                     AND PAR.TIMESTAMP = AR.TIMESTAMP
                                     AND PAR.ENTITY_ID = NODE.ID
                                     AND PAR.NORMALIZED_METHOD_SIGNATURE IN ('GET /entity/#/file','GET /entity/#/version/#/file')
                                     AND AR.TIMESTAMP BETWEEN unix_timestamp(curdate())*1000 - (30*24*60*60*1000) AND  unix_timestamp(curdate())*1000
                                     GROUP BY PAR.ENTITY_ID, AR.USER_ID ;")

downloadReport1Data = fetch(downloadReport1, n=-1)

downloadReport2 = dbSendQuery(mydb, 'SELECT NODE.ID, FDR.USER_ID ,COUNT(*)
                                    FROM warehouse.FILE_DOWNLOAD_RECORD FDR,
                                    (SELECT ID, CREATED_BY FROM warehouse.NODE_SNAPSHOT WHERE (ID, TIMESTAMP) IN (SELECT ID, TIMESTAMP FROM tyu.syn7362874)) AS NODE
                                    WHERE FDR.ASSOCIATION_OBJECT_ID = NODE.ID
                                    AND FDR.ASSOCIATION_OBJECT_TYPE = "FileEntity"
                                    AND FDR.TIMESTAMP BETWEEN unix_timestamp(curdate())*1000 - (7*24*60*60*1000) AND  unix_timestamp(curdate())*1000
                                    GROUP BY NODE.ID, FDR.USER_ID;')

downloadReport2Data = fetch(downloadReport2, n=-1)
downloadReport2Data$ENTITY_ID <- downloadReport2Data$ID
downloadReport2Data$ID <- NULL
total <- rbind(downloadReport1Data,
          downloadReport2Data)


write.csv(total,"Tesla_download_Stats.txt") #uploaded to syn8442870
downloadStatsEnt <- synGet("syn8442870")
downloadStats <- read.csv(getFileLocation(downloadStatsEnt))

teams = synTableQuery('SELECT * FROM syn8220615')

teamMembers <- sapply(teams@values$realTeam, function(teamName) {
  teamId <- synRestGET(paste0('/teams?fragment=',teamName))$results[[1]]$id
  results <- synRestGET(paste0('/teamMembers/',teamId))
  members <- unlist(lapply(results$results, function(x){
    x$member$ownerId
  }))
  members[members != "3324230"]
})
teamMembers$TESLA_Consortium_Admins <- NULL

usersDownloadFile <- aggregate(USER_ID ~ ENTITY_ID, downloadStats, c)
teamDownloadStats <- apply(usersDownloadFile[!usersDownloadFile$ENTITY_ID %in% c(8303327,8303410),], 1, function(info){
  users = unlist(info$USER_ID)
  entId = info$ENTITY_ID
  ent = synGet(paste0("syn",entId), downloadFile=F)
  exist = lapply(teamMembers, function(members) {
    members <- unlist(members)
    if (length(members) >0) {
      if (any(members %in% users)) {
        TRUE
      } else {
        FALSE
      }
    }
  })
  data <- data.frame(unlist(exist))
  colnames(data) <- ent@properties$name
  data
  #data.frame("fileName" = ent@properties$name, "notDownloaded"= paste(names(exist)[unlist(exist)==F],collapse = ","),"downloaded"= paste(names(exist)[unlist(exist)==T],collapse = ","))
})

total = do.call(cbind,teamDownloadStats)
teamsDownloaded <- t(total)
write.table(teamsDownloaded, "DownloadStats.csv",sep="\t",quote = F)
#Number of teams downloaded data
#Which teams have downloaded the data
#Which ones haven't

sumFastqDownloads <- apply(teamsDownloaded[grepl("*fastq.gz",row.names(teamsDownloaded)),],2, sum)
names(sumFastqDownloads[sumFastqDownloads == 0])

sumNotFastqDownloads <- apply(teamsDownloaded[!grepl("*fastq.gz",row.names(teamsDownloaded)),],1, sum)

NotFastqDownloadStats <- data.frame("downloadStats" = sumNotFastqDownloads)
write.csv(NotFastqDownloadStats, "vcfbam_downloadStats.csv",quote=F)


#Of the other fileTypes, how many teams have downloaded (fileType)
