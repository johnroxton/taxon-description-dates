# TODO: Add comment
#
# Author: David Schellenberger Costa
###################################################################################################
# 1 Visualize data from LCVP and total numbers sent by Martin Freiberg
# 2 Try to reproduce the analyses from "Trends in botanical exploration in Nigeria forecast over
# 1000 yet undescribed vascular plant species", https://doi.org/10.1093/aob/mcad106
###################################################################################################

# 1 Visualize data from LCVP and total numbers sent by Martin Freiberg#############################
nextScript <- NULL

# load in libraries
library(data.table)
library(RColorBrewer)

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# read in data
lcvp1 <- fread(paste0(.brd, "taxonomy/LCVP/preprocessed/LCVP_2023-08-24_pp.gz")) # LCVP 2.0
lcvp2 <- fread(paste0(.brd, "taxon description dates/Tracheophyta beschriebene Arten pro Jahr.txt")) # data from Martin

table(lcvp1$Literature != "") # only very few entries (10%) have dates

# create dummy variable for year extraction
lcvp1[, temp := Literature]

# extract years from dummy variable
years <- strsplit(lcvp1$temp, ";")
years <- sapply(years, function(x) sapply(x, function(y) if (grepl("\\d{4}", y)) y else NA))
years <- sapply(years, function(x) sapply(x, function(y) regmatches(y, regexpr("\\d{4}", y))))
years <- sapply(years, function(x) if (length(unlist(x)) > 0) min(unlist(x), na.rm = TRUE) else NA)
lcvp1[, year := as.numeric(years)]

# check data

# years before publication of Linnés species descriptions
lcvp1[year < 1781] # problems with true year in parentheses
lcvp1[grepl("\\(\\d{4}\\)", Literature)] # show all entries with true year in parentheses
lcvp1[`global Id` == 1337197, year := 2019] # only one entry has an error

# years after data collection
lcvp1[year > 2022] # 0 has been read as 9
lcvp1[year > 2022, year := gsub("9", "0", year)] # repair all entries with errors
lcvp1[`global Id` == 1339189] # check first entry because of warning

# plot data
plot(table(lcvp1$year), type = "o", ylab = "revisions / descriptions",
		xlab = "year", xlim = c(1754, 2023), xaxt="n",
		col = brewer.pal(8, "Accent")[1])

# The data in LCVP are likely mostly no first descriptions but dates of taxonomic revisions.

colnames(lcvp2) <- c("year", "descriptions")
lcvp2 <- lcvp2[year > 1753]
lcvp2 <- lcvp2[year < 2023]

# add data from Martin
lines(lcvp2,type="o",col=brewer.pal(8,"Accent")[2],lwd=2)

# Martin's aggregated data does not allow for much analyses.

# 2 Try to reproduce the analyses from "Trends in botanical exploration in Nigeria forecast over
# 1000 yet undescribed vascular plant species", https://doi.org/10.1093/aob/mcad106
nextScript <- NULL

# check this URL for all authors from IPNI
# https://www.ipni.org/?perPage=500&q=author%20std%3A*


# load in libraries
library(data.table)
library(RColorBrewer)

# clear workspace
rm(list = ls())

# load help functions
source(paste0(.brd,"taxonomy/taxonomy author normalization.R"))

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# read in data
# use this version because primary author and year are given
wcvp <- fread(paste0(.brd, "taxonomy/WCVP/preprocessed/WCVP_2022_pp.gz"))
dist <- fread(paste0(.brd, "taxonomy/WCVP/raw/wcvp_220228/wcvp_distribution.txt"))
auth <- fread(paste0(.brd, "taxonomy/WCVP/raw/wcvp_220228/wcvp_reference.txt"))
load(paste0(paste0(.brd, "taxonomy/botanical author names IPNI_2023-02-13.RData")))

str(wcvp)
str(dist)
str(auth)

# get ids of plants from Nigeria
np <- dist[grepl("Nigeria",area)]$plant_name_id # Nigeria plants

# subset wcvp for those plants (including synonyms)
wcvpn <- wcvp[plant_name_id %in% np | accepted_plant_name_id %in% np]
str(wcvpn)

# test for problematic entries
wcvpn[taxon_status != "Accepted" & !(accepted_plant_name_id %in% np)]
# 10 unplaced names found

# reduce dataset to accepted names, as these should also be the oldest ones
wcvpn <- wcvpn[taxon_status == "Accepted"]

# convert first published column to year
wcvpn[, year := as.numeric(gsub("\\D","",first_published))]

# test for problematic entries
wcvpn[is.na(year) & first_published != ""]
# only found entries with unknown publication year

wcvpn <- wcvpn[!is.na(year)] # remove names with unknown publication year

# show primary authors
table(wcvpn$primary_author)

# remove "... ex" from primary author
wcvpn[, primary_author := sub(".*\\sex\\s","",primary_author)]

#repair primary authors
wcvpn[grepl("A J",primary_author),primary_author := sub("A J","A.J",primary_author)]
wcvpn[grepl("Ballard Hook.",primary_author),primary_author := sub("Ballard Hook.","Ballard",primary_author)]
wcvpn[grepl("Le Prieur",primary_author),primary_author := sub("Le Prieur","Leprieur",primary_author)]
wcvpn[ipni_id == "310558-1",primary_author := "Lorougnon & J.Raynal"]
wcvpn[,primary_author := sub("^Okafa$","Okafor",primary_author)]

wcvpn <- wcvpn[primary_author != ""] # remove names with unknown author

# check whether all primary authors can be found in IPNI
pa <- sort(unique(wcvpn$primary_author)) # primary authors
pa <- strsplit(pa,split=",|&")
pa <- unlist(pa)
pa <- gsub("\\s+"," ",pa)
pa <- sub("^\\s","",pa)
pa <- sub("\\s$","",pa)

aut <- pa[!(pa %in% aNI$authorAbbreviation)]
aut <- authorNormalization(aut)

