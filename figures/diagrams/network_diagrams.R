
library(ggplot2)
library(ggbipart)
library(network)
library(GGally)

azov_whole <- read.table("azov_whole.txt", header = T, row.names = 1)
japan_whole <- read.table("japan_whole.txt", header = T, row.names = 1)

azov_active <- read.table("azov_active.txt", header = T, row.names = 1)
japan_active <- read.table("japan_active.txt", header = T, row.names = 1)

azov_pasive <- read.table("azov_pasive.txt", header = T, row.names = 1)
japan_pasive <- read.table("japan_pasive.txt", header = T, row.names = 1)

azov_ecto <- read.table("azov_ecto.txt", header = T, row.names = 1)
japan_ecto <- read.table("japan_ecto.txt", header = T, row.names = 1)


# WHOLE COMMUNITY ---------------------------------------------------------------------------------
## SEA OF AZOV

#tiff(filename = "whole_azov.tiff", width = 9000, height = 5000, units = "px", res = 800)
bip1 = azov_whole
bip.net1 <- bip_init_network(as.matrix(azov_whole))
col= c("P"= "#D7E4F2", "A"= "#FF66B2")
bip_ggnet(bip.net1, as.matrix(bip1),
          size=2, shape= 18, #label= T,
          color= "mode", palette= col, 
          layout.exp= 0) + 
  geom_point(shape=23, fill="#D7E4F2", size = 4, color="grey40") +
  geom_text(aes(label= network.vertex.names(bip.net1)),
            color= "grey70", size= 1) +
  theme(legend.position="none")             
#dev.off()

## SEA OF JAPAN

#tiff(filename = "whole_japan.tiff", width = 9000, height = 5000, units = "px", res = 600)
bip2 = japan_whole
bip.net2 <- bip_init_network(as.matrix(japan_whole))
col= c("P"= "#D7E4F2", "A"= "#FF66B2")
bip_ggnet(bip.net2, as.matrix(bip2),
          size=2, shape= 18, #label= T,
          color= "mode", palette= col, 
          layout.exp= 0) + 
  geom_point(shape=23, fill="#D7E4F2", size = 4, color="grey40") +
  geom_text(aes(label= network.vertex.names(bip.net2)),
            color= "grey70", size= 1) +
  theme(legend.position="none")             
#dev.off()

#ACTIVELY TRANSMITTED COMMUNITY---------------------------------------------------------------
## SEA OF AZOV

#tiff(filename = "active_azov.tiff", width = 9000, height = 5000, units = "px", res = 800)
bip3 = azov_active
bip.net3 <- bip_init_network(as.matrix(azov_active))
col= c("P"= "#D7E4F2", "A"= "#FF66B2")
bip_ggnet(bip.net3, as.matrix(bip3),
          size=2, shape= 18, #label= T,
          color= "mode", palette= col, 
          layout.exp= 0) + 
  geom_point(shape=23, fill="#D7E4F2", size = 4, color="grey40") +
  geom_text(aes(label= network.vertex.names(bip.net3)),
            color= "grey70", size= 1) +
  theme(legend.position="none") 
#dev.off()

## SEA OF JAPAN

#tiff(filename = "active_japan.tiff", width = 9000, height = 5000, units = "px", res = 800)
bip4 = japan_active
bip.net4 <- bip_init_network(as.matrix(japan_active))
col= c("P"= "#D7E4F2", "A"= "#FF66B2")
bip_ggnet(bip.net4, as.matrix(bip4),
          size=2, shape= 18, #label= T,
          color= "mode", palette= col, 
          layout.exp= 0) + 
  geom_point(shape=23, fill="#D7E4F2", size = 4, color="grey40") +
  geom_text(aes(label= network.vertex.names(bip.net4)),
            color= "grey70", size= 1) +
  theme(legend.position="none") 
#dev.off()


# TROPHICALLY TRANSMITTED COMMUNITY---------------------------------------------------------------
# SEA OF AZOV

#tiff(filename = "pasive_azov.tiff", width = 9000, height = 5000, units = "px", res = 700)
bip5 = azov_pasive
bip.net5 <- bip_init_network(as.matrix(azov_pasive))
col= c("P"= "#D7E4F2", "A"= "#FF66B2")
bip_ggnet(bip.net5, as.matrix(bip5),
          size=2, shape= 18, #label= T,
          color= "mode", palette= col, 
          layout.exp= 0) + 
  geom_point(shape=23, fill="#D7E4F2", size = 4, color="grey40") +
  geom_text(aes(label= network.vertex.names(bip.net5)),
            color= "grey70", size= 1) +
  theme(legend.position="none")          
#dev.off()

## SEA OF JAPAN 

#tiff(filename = "pasive_japan.tiff", width = 9000, height = 5000, units = "px", res = 800)
bip6 = japan_pasive
bip.net6 <- bip_init_network(as.matrix(japan_pasive))
col= c("P"= "#D7E4F2", "A"= "#FF66B2")
bip_ggnet(bip.net6, as.matrix(bip6),
          size=2, shape= 18, #label= T,
          color= "mode", palette= col, 
          layout.exp= 0) + 
  geom_point(shape=23, fill="#D7E4F2", size = 4, color="grey40") +
  geom_text(aes(label= network.vertex.names(bip.net6)),
            color= "grey70", size= 1) +
  theme(legend.position="none")         
#dev.off()


# ECTOPARASITE COMMUNITY------------------------------------------------------------------------------------
## SEA OF AZOV

#tiff(filename = "ecto_azov.tiff", width = 9000, height = 5000, units = "px", res = 800)
bip7 = azov_ecto
bip.net7 <- bip_init_network(as.matrix(azov_ecto))
col= c("P"= "#D7E4F2", "A"= "#FF66B2")

bip_ggnet(bip.net7, as.matrix(bip7),
          size=2, shape= 18, #label= T,
          color= "mode", palette= col, 
          layout.exp= 0) + 
  geom_point(shape=23, fill="#D7E4F2", size = 4, color="grey40") +
  geom_text(aes(label= network.vertex.names(bip.net7)),
            color= "grey70", size= 1) +
  theme(legend.position="none") 
#dev.off()


## SEA OF JAPAN

#tiff(filename = "ecto_japan.tiff", width = 9000, height = 5000, units = "px", res = 800)
bip8 = japan_ecto
bip.net8 <- bip_init_network(as.matrix(japan_ecto))
col= c("P"= "#D7E4F2", "A"= "#FF66B2")
bip_ggnet(bip.net8, as.matrix(bip8),
          size=2, shape= 18, #label= T,
          color= "mode", palette= col, 
          layout.exp= 0) + 
  geom_point(shape=23, fill="#D7E4F2", size = 4, color="grey40") +
  geom_text(aes(label= network.vertex.names(bip.net8)),
            color= "grey70", size= 1) +
  theme(legend.position="none")              
#dev.off()

