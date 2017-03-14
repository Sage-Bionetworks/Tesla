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
head(downloadReport1Data)
head(downloadReport2Data)
downloadReport2Data$ENTITY_ID <- downloadReport2Data$ID
downloadReport2Data$ID <- NULL
total <- rbind(downloadReport1Data,
          downloadReport2Data)
write.csv(total,"Tesla_download_Stats.txt") #uploaded to syn8442870
table(total$ENTITY_ID)

teams = synTableQuery('SELECT * FROM syn8220615')
teams@values$realTeam
syn.getTeam()

teamMembers <- sapply(teams@values$realTeam, function(teamName) {
  teamId <- synRestGET(paste0('/teams?fragment=',teamName))$results[[1]]$id
  results <- synRestGET(paste0('/teamMembers/',teamId))
  members <- unlist(lapply(results$results, function(x){
    x$member$ownerId
  }))
  members[members != "3324230"]
})

usersDownloadFile <- aggregate(USER_ID ~ ENTITY_ID, total, c)
downloadStats <- lapply(usersDownloadFile$USER_ID, function(users){
  users = unlist(users)
  lapply(teamMembers, function(members) {
    members <- unlist(members)
    if (length(members) >0) {
      if (any(members %in% users)) {
        TRUE
      }
    }
  }) 
})

unlist(downloadStats)
#Number of teams downloaded data
#Which teams have downloaded the data
#Which ones haven't

#group by ID, print out teams that have not downloaded data
#Teams that haven't downloaded the fastqs

#Of the other fileTypes, how many teams have downloaded (fileType)
