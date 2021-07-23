library(bipartite)

azov_whole <- as.matrix(read.table("azov_whole.txt", header = T, row.names = 1))
japan_whole<- as.matrix(read.table("japan_whole.txt", header = T, row.names = 1))

ap_whole<- specieslevel(azov_whole, level = "higher", index = c("betweenness"))

jp_whole<- specieslevel(japan_whole, level = "higher", index = c("betweenness"))


azov_active<- as.matrix(read.table("azov_active.txt", header = T, row.names = 1))
japan_active<- as.matrix(read.table("japan_active.txt", header = T, row.names = 1))

ap_active<- specieslevel(azov_active, level = "higher", index = c("betweenness"))

jp_active<- specieslevel(japan_active, level = "higher", index = c("betweenness"))


azov_ecto<- as.matrix(read.table("azov_ecto.txt", header = T, row.names = 1))
japan_ecto<- as.matrix(read.table("japan_ecto.txt", header = T, row.names = 1))

ap_ecto<- specieslevel(azov_ecto, level = "higher", index = c("betweenness"))

jp_ecto<- specieslevel(japan_ecto, level = "higher", index = c("betweenness"))

azov_pasive<- as.matrix(read.table("azov_pasive.txt", header = T, row.names = 1))
japan_pasive<- as.matrix(read.table("japan_pasive.txt", header = T, row.names = 1))

ap_pasive<- specieslevel(azov_pasive, level = "higher", index = c("betweenness"))

jp_pasive<- specieslevel(japan_pasive, level = "higher", index = c("betweenness"))


write.csv2(ap_whole, file = "ap_whole.csv")
write.csv2(ap_active, file = "ap_active.csv")
write.csv2(ap_pasive, file = "ap_pasive.csv")
write.csv2(ap_ecto, file = "ap_ecto.csv")

write.csv2(jp_whole, file = "jp_whole.csv")
write.csv2(jp_active, file = "jp_active.csv")
write.csv2(jp_pasive, file = "jp_pasive.csv")
write.csv2(jp_ecto, file = "jp_ecto.csv")

