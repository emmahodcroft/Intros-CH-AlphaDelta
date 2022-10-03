# We want our focal set to be primarily sequences from Switzerland that are Alpha/Delta prior to 31 March
# or 31 July (respectively).

setwd("2021_Delta/clusters")

### Find out how many from Metadata, if curious. But use method below for final files.

#meta <- read.csv("C:/Users/Emma/wsl/corona/ncov/data/downloaded_gisaid.tsv", sep="\t", as.is=T)
#this was downloaded on 15 Feb
meta <- read.csv("C:/Users/Emma/wsl/corona/Projects/Swiss-Wengen/data/metadata.tsv", sep="\t", as.is=T)

#Alpha
swiss_alpha <- meta[which(meta$country == "Switzerland" & meta$Nextstrain_clade == "20I (Alpha, V1)"),]
swiss_alpha$date <- as.Date(swiss_alpha$date)

length(which(swiss_alpha$date <= as.Date("2021-03-31")))

#Delta
swiss_delta <- meta[which(meta$country == "Switzerland" & 
    (meta$Nextstrain_clade == "21A (Delta)" | meta$Nextstrain_clade == "21I (Delta)" | meta$Nextstrain_clade == "21J (Delta)")),]
swiss_delta$date <- as.Date(swiss_delta$date)

length(which(swiss_delta$date <= as.Date("2021-07-31")))








####### Run allClusterDynamics_faster.py from CoVariants.org with correct dated_cluster & dated_limit
# move the resulting dated cluster files into this folder
# Now we strip out just the Swiss ones.

setwd("2021_Delta/clusters")

alphaclus <- read.table("input/cluster_20I.Alpha.V1-2021-03-31.txt", header=F)

swissalpha <- grep("Switzerland", alphaclus$V1)
length(swissalpha)
# 8083 - 1 Apr 22

write.table(alphaclus$V1[swissalpha], "cluster_20I.Alpha.V1-swiss-2021-03-31.txt",
    col.names=F, quote=F, row.names=F)


deltaclus <- read.table("input/cluster_21A.Delta-2021-07-31.txt", header=F)
swissdelta <- grep("Switzerland", deltaclus$V1)
length(swissdelta)
# 142 - 1 Apr 22

delta21I <- read.table("input/cluster_21I.Delta-2021-07-31.txt")
swissDelta21I <- grep("Switzerland", delta21I$V1) 
length(swissDelta21I)
# 564 - 1 Apr 22

delta21J <- read.table("input/cluster_21J.Delta-2021-07-31.txt")
swissDelta21J <- grep("Switzerland", delta21J$V1) 
length(swissDelta21J)
# 4526 - 1 Apr 22

deltas <- c(deltaclus$V1[swissdelta], delta21I$V1[swissDelta21I], delta21J$V1[swissDelta21J]) 
#total is 5232 - 1 Apr 22

#of these sequences, not all will pass QC and make it into the final Nextstrain builds, we expect
# to loose a couple hundred or so.

write.table(deltas, "cluster_21A.Delta-swiss-2021-07-31.txt",
    col.names=F, quote=F, row.names=F)

