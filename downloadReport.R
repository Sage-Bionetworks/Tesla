# library(RMySQL)
#
# START = as.Date("2019-07-01")
# DAYS_BEFORE = as.numeric(Sys.Date() - START)
# #set this environment variable in .Renviron
# datawareHousePW = Sys.getenv("SAGEDATAWAREHOUSEPW")
# # Create node snapshot first
# # http://build-system-dw.sagebase.org:8080/job/Launch%20warehouse%20snapshot/
# # Fill out user and password with user and password set.
# # host and db is found in console output
# mydb = dbConnect(MySQL(),
#                  user = 'thomas',
#                  password = datawareHousePW,
#                  host = 'thomas-tesla-40.c95bbsvwbjlu.us-east-1.rds.amazonaws.com',
#                  db = 'warehouse')
# # Tries to list the tables
# dbListTables(mydb)
# # dbSendQuery(mydb, "DROP TABLE syn7362874;")
# 
# 
# getQueryInfo <- function(mydb, DAYS_BEFORE) {
#   downloadReport1 = dbSendQuery(mydb, sprintf("SELECT PAR.ENTITY_ID, NODE.CREATED_BY, COUNT(*)
#                                              FROM PROCESSED_ACCESS_RECORD PAR,
#                                              (SELECT ID, CREATED_BY FROM NODE_SNAPSHOT WHERE (ID, TIMESTAMP) IN (SELECT ID, TIMESTAMP FROM syn7362874)) AS NODE
#                                              WHERE PAR.ENTITY_ID = NODE.ID
#                                              AND PAR.NORMALIZED_METHOD_SIGNATURE IN ('GET /entity/#/file','GET /entity/#/version/#/file')
#                                              AND TIMESTAMP BETWEEN unix_timestamp(curdate())*1000 - (%d*24*60*60*1000) AND  unix_timestamp(curdate())*1000 - (%d*24*60*60*1000)
#                                              GROUP BY PAR.ENTITY_ID, NODE.CREATED_BY;", DAYS_BEFORE, DAYS_BEFORE-7))
#   
#   downloadReport1Data = fetch(downloadReport1, n=-1)
#   
#   
#   downloadReport2 = dbSendQuery(mydb, sprintf('SELECT NODE.ID, FDR.USER_ID ,COUNT(*)
#                                             FROM FILE_DOWNLOAD_RECORD FDR,
#                                             (SELECT ID, CREATED_BY FROM NODE_SNAPSHOT WHERE (ID, TIMESTAMP) IN (SELECT ID, TIMESTAMP FROM syn7362874)) AS NODE
#                                             WHERE FDR.ASSOCIATION_OBJECT_ID = NODE.ID
#                                             AND FDR.ASSOCIATION_OBJECT_TYPE = "FileEntity"
#                                             AND FDR.TIMESTAMP BETWEEN unix_timestamp(curdate())*1000 - (%d*24*60*60*1000) AND  unix_timestamp(curdate())*1000 - (%d*24*60*60*1000)
#                                             GROUP BY NODE.ID, FDR.USER_ID;',DAYS_BEFORE, DAYS_BEFORE-7))
#   
#   downloadReport2Data = fetch(downloadReport2, n=-1)
#   downloadReport2Data$ENTITY_ID <- downloadReport2Data$ID
#   downloadReport2Data$ID <- NULL
#   downloadStats <- rbind(downloadReport1Data,
#                          downloadReport2Data)
#   return(downloadStats)
# }
# NOT_FOUND_FILES = c(8303327,8303410,8077800,8077817,8077818,10156068,10166143,10166144,10166145,10166146,10166147,10166148)
# downloadStats = downloadStats[!downloadStats$ENTITY_ID %in% NOT_FOUND_FILES,]
# downloaded = table(downloadStats$ENTITY_ID)
# names(downloaded) = paste0("syn",names(downloaded))
# 
# downloadedDf = data.frame(downloaded)
# synNames = sapply(downloadedDf$Var1, function(x) {
#   temp = synGet(as.character(x),downloadFile=F)
#   temp@properties$name
# })
# colnames(downloadedDf) <- c("synId","numberOfDownloads")
# downloadedDf$synName = synNames
#write.csv(total,"Tesla_download_Stats.txt") #uploaded to syn8442870
#downloadStatsEnt <- synGet("syn8442870")
#downloadStats <- read.csv(getFileLocation(downloadStatsEnt))


### Above workflow is obsolete

library(synapser)
synLogin()

# Run sql queries first in mysql workbench

download_report1 = read.csv("download_report1.csv")
download_report2 = read.csv("download_report2.csv")
colnames(download_report1) = c("ENTITY_ID", "USER_ID", "COUNT")
colnames(download_report2) = c("ENTITY_ID", "USER_ID", "COUNT")
downloadStats = rbind(download_report1, download_report2)
teams = synTableQuery("SELECT * FROM syn8220615 where alias <> 'butcherbirds'")
teamsdf = teams$asDataFrame()
teamMembers <- sapply(teamsdf$realTeam, function(teamName) {
  print(teamName)
  teamId <- synRestGET(paste0('/teams?fragment=',teamName))$results[[1]]$id
  results <- synRestGET(paste0('/teamMembers/',teamId))
  members <- unlist(lapply(results$results, function(x){
    x$member$ownerId
  }))
  members[members != "3324230"]
})
teamMembers$TESLA_Consortium_Admins <- NULL

usersDownloadFile <- aggregate(USER_ID ~ ENTITY_ID, downloadStats, c)
teamDownloadStats <- apply(usersDownloadFile, 1, function(info){
  users = unlist(info$USER_ID)
  entId = info$ENTITY_ID
  ent = synGet(paste0("syn",entId), downloadFile = F)
  exist = lapply(teamMembers, function(members) {
    members <- unlist(members)
    if (length(members) > 0) {
      if (any(members %in% users)) {
        TRUE
      } else {
        FALSE
      }
    }
  })
  data <- data.frame(unlist(exist))
  colnames(data) <- ent$properties$name
  data
})

total = do.call(cbind, teamDownloadStats)
teamsDownloaded <- t(total)
write.csv(teamsDownloaded, "round3_download_statistics.csv", quote = F)

#Number of teams downloaded data
#Which teams have downloaded the data
#Which ones haven't

# metadata = synTableQuery('SELECT * FROM syn8292741 where round = "2"')
# metadataDf = metadata@values
# #only download round 2 stuff
# teamsDownloaded = teamsDownloaded[metadataDf$name,]
# write.table(teamsDownloaded,
#             "round2_DownloadStats.csv",
#             sep = "\t",
#             quote = F)
# fastqDownloads <- apply(teamsDownloaded[grepl("*fastq.gz",row.names(teamsDownloaded)),],2, sum)
# notFastqDownloads <- apply(teamsDownloaded[!grepl("*fastq.gz",row.names(teamsDownloaded)),],2, sum)
# vcfDownloads <- apply(teamsDownloaded[!grepl("*vcf.gz",row.names(teamsDownloaded)),],2, sum)
# bamDownloads <- apply(teamsDownloaded[!grepl("*bam",row.names(teamsDownloaded)),],2, sum)
# metadataDownloads <- apply(teamsDownloaded[!grepl("*csv|*txt|*xlsx",row.names(teamsDownloaded)),],2, sum)
# png("round2_dataTypeDownloaded.png", width = 600, height=400)
# DLs = c("fastq" = sum(fastqDownloads > 0), "vcf" = sum(vcfDownloads>0), "bam" = sum(bamDownloads>0), "metadata" = sum(metadataDownloads>0))
# barplot(DLs,main = "Number of Teams Downloaded File Type", xlab = "File Type",ylab = "Number of Teams")
# dev.off()
# synStore(File("./round2_dataTypeDownloaded.png",parentId = "syn8082860"))
# 
# NotFastqDownloadStats <- data.frame("downloadStats" = notFastqDownloads)
# write.csv(NotFastqDownloadStats, "vcfbam_downloadStats.csv",quote=F)

#Of the other fileTypes, how many teams have downloaded (fileType)
