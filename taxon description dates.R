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

# years before publication of Linnaeus species descriptions
lcvp1[year < 1781] # problems with true year in parentheses
lcvp1[grepl("\\(\\d{4}\\)", Literature)] # show all entries with true year in parentheses
lcvp1[`global Id` == 1337197, year := 2019] # only one entry has an error

# years after data collection
lcvp1[year > 2022] # 0 has been read as 9
lcvp1[year > 2022, year := gsub("9", "0", year)] # repair all entries with errors
lcvp1[`global Id` == 1339189] # check first entry because of warning

# plot data
plot(table(lcvp1$year),
	type = "o", ylab = "revisions / descriptions",
	xlab = "year", xlim = c(1754, 2023), xaxt = "n",
	col = brewer.pal(8, "Accent")[1]
)

# The data in LCVP are likely mostly no first descriptions but dates of taxonomic revisions.

colnames(lcvp2) <- c("year", "descriptions")
lcvp2 <- lcvp2[year > 1753]
lcvp2 <- lcvp2[year < 2023]

# add data from Martin
lines(lcvp2, type = "o", col = brewer.pal(8, "Accent")[2], lwd = 2)

# Martin's aggregated data does not allow for much analyses.

# 2 Try to reproduce the analyses from "Trends in botanical exploration in Nigeria forecast over###
# 1000 yet undescribed vascular plant species", https://doi.org/10.1093/aob/mcad106################
nextScript <- NULL

# load in libraries
library(data.table)
library(RColorBrewer)
library(taxize)
library(RSelenium)

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxonomy"))

# load help functions
source("taxonomy help functions.R")
source("taxonomy author normalization.R")

# function to retrieve basionyms from POWO or year of first publication
POWOSearch <- function(ipniid, basionymSearch = TRUE, authCompThr = 0.4){
	resFull <- readLines(paste0("https://powo.science.kew.org/taxon/",ipniid))
	if (basionymSearch == TRUE){
		if (any(grepl("Homotypic Synonyms",resFull))){
			targetAuth <- sub("\\s*</small>.*","",sub(".*<small>\\s*","",resFull[grepl("c-summary__heading p-xl",resFull)]))
			targetAuth <- sub("\\(","",sub("\\).*","",targetAuth))
			resFull <- resFull[grep("Homotypic Synonyms",resFull):length(resFull)]
			resFull <- resFull[1:grep("</ul>$",resFull)[1]]
			resFull <- resFull[grepl("^\\s*<li",resFull)]
			synAuth <- sub("\\s*</a>.*","",sub(".*em>\\s*","",resFull))
			authComp <- colSums(sapply(synAuth,function(x)authorMatch(targetAuth,x)))
			if (any(authComp > authCompThr)){
				resFull <- resFull[order(authComp,decreasing=TRUE)[1]]
				oriIpniid <- sub("\".*","",sub(".*:names:","",resFull))
				oriName <- gsub("<[^<>]+>","",sub("</em>[^<>]+</a>.*","",sub(".*><em lang='la'>","",resFull)))
				oriAuthors <- synAuth[order(authComp,decreasing=TRUE)[1]]
				oriYear <- gsub("\\D","",regmatches(resFull,regexpr("\\(\\d{4}\\)",resFull)))
				return(c(oriIpniid,oriName,oriAuthors,oriYear))
			} else {
				return(rep(NA,4))
			}
		} else {
			return(rep(NA,4))
		}
	} else {
		res <- resFull[grepl("First published",resFull)]
		res <- as.numeric(gsub("\\D","",regmatches(res,regexpr("\\(\\d{4}\\)",res))))
		return(res)
	}
}

# function to retrieve basionyms from tropicos
tropicosBasionym <- function(name,authCompThr = 0.4){
	if (!("authorMatch" %in% ls(globalenv()))) {
		print("taxonomy help functions need to be sourced first.")
		return(c(NA,NA,NA))
	}
	# account for issues with infraspecies
	if (grepl("\\.",name)) tempName <- sub("\\s*\\S*\\..*","",name) else tempName <- name
	foundTro <- tp_search(tempName)
	if (tempName != name) foundTro <- foundTro[foundTro$scientificname == name,]
	if (nrow(foundTro) > 1){
		nw1 <- length(gregexpr("\\s",name))
		nw2 <- sapply(gregexpr("\\s",foundTro$scientificname),length)
		foundTro <- foundTro[nw2 == nw1,][1,]
	}
	if ("nameid" %in% colnames(foundTro)){
		res <- tp_synonyms(foundTro$nameid)
		if (res$accepted$nameid != "no syns found"){
			nchars1 <- nchar(res$accepted$scientificname)
			nchars2 <- nchar(res$accepted$scientificnamewithauthors)
			targetAuth <- substr(res$accepted$scientificnamewithauthors,nchars1 + 2,nchars2)
			targetAuth <- sub("\\(","",sub("\\).*","",targetAuth))
			nchars1 <- nchar(res$synonyms$scientificname)
			nchars2 <- nchar(res$synonyms$scientificnamewithauthors)
			synAuth <- substr(res$synonyms$scientificnamewithauthors,nchars1 + 2,nchars2)
			authComp <- colSums(sapply(synAuth,function(x)authorMatch(targetAuth,x)))
			if (any(authComp > authCompThr)){
				name <- res$synonyms$scientificname[order(authComp,decreasing=TRUE)[1]]
				# account for issues with infraspecies
				if (grepl("\\.",name)) tempName <- sub("\\s*\\S*\\..*","",name) else tempName <- name
				res <- tp_search(tempName)
				if (tempName != name) res <- res[res$scientificname == name,]
				if (nrow(res) > 1){
					nw1 <- length(gregexpr("\\s",name))
					nw2 <- sapply(gregexpr("\\s",res$scientificname),length)
					res <- res[nw2 == nw1,][1,]
				}
				return(c(foundTro$nameid,res$scientificname,res$author,res$displaydate))
			} else {
				return(c(foundTro$nameid,NA,NA,NA))
			}
		} else {
			return(c(foundTro$nameid,NA,NA,NA))
		}
	} else {
		return(rep(NA,4))
	}
}	

# function to extract the year from ipni table
ipniTableYear <- function(name, authors){
	if (!("remDr" %in% ls(globalenv()))) {
		print("Start an rsDriver instance first.")
		return(c(NA,NA))
	}
	remDr$navigate(paste0("https://www.ipni.org/?q=",gsub(" ","%20",name)))
	repeat {
		Sys.sleep(0.2)
		res <- remDr$findElements(using = "class", value = "result")
		if (length(res) > 0){
			break
		} else {
			res <- remDr$findElements(using = "class", value = "no-results")
			if (length(res) > 0){
				return(NA)
			}
		}
	}
	bas <- which(sapply(res,function(x) grepl(authors,x$getElementAttribute("outerHTML")[[1]])))[1]
	if (!is.na(bas)){
		year <- as.numeric(gsub("\\D","",regmatches(res[[bas]]$getElementAttribute("outerHTML")[[1]],
			regexpr("\\(\\d{4}\\)",res[[bas]]$getElementAttribute("outerHTML")[[1]]))))
		return(year)
	} else {
		return(NA)
	}
}

# function to extract the basionym from tropicos details
# very slow due to very slow tropicos website
tropicosDetailsBasionym <- function(nameid){
	if (!("remDr" %in% ls(globalenv()))) {
		print("Start an rsDriver instance first.")
		return(c(NA,NA))
	}
	remDr$navigate(paste0("https://tropicos.org/name/",nameid))
	repeat {
		Sys.sleep(2)
		res <- remDr$findElements(using = "class", value = "control-group")
		resYear <- remDr$findElements(using = "class", value = "section-group")
		if (length(res) > 0 && length(resYear) > 0) break
	}
	for (i in seq_along(res)){
		if (grepl("Basionym:",res[[i]]$getElementAttribute("outerHTML")[[1]])){
			res <- res[[i + 1]]
			break
		}
	}
	for (i in seq_along(resYear)){
		if (grepl("Published In:",resYear[[i]]$getElementAttribute("outerHTML")[[1]])){
			resYear <- resYear[[i]]
			break
		}
	}
	if (length(res) == 1 && length(resYear) == 1){
		author <- sub("</a>.*","",sub(".*</i>\\s*","",res$getElementAttribute("outerHTML")[[1]]))
		name <- gsub("</?i>","",sub("[^<>]+</a>.*","",sub(".*><i>\\s*","",res$getElementAttribute("outerHTML")[[1]])))
		year <- regmatches(resYear$getElementAttribute("outerHTML")[[1]],regexpr("\\d{4}",resYear$getElementAttribute("outerHTML")[[1]]))
		return(c(name,author,year))
	} else {
		return(rep(NA,3))
	}
}

# read in data
# use this version because primary author and year are given
wcvp <- fread(paste0(.brd, "taxonomy/WCVP/preprocessed/WCVP_2022_pp.gz"))
dist <- fread(paste0(.brd, "taxonomy/WCVP/raw/wcvp_220228/wcvp_distribution.txt"))
auth <- fread(paste0(.brd, "taxonomy/WCVP/raw/wcvp_220228/wcvp_reference.txt"))
versions <- list.files(pattern = "botanical author names IPNI")
aNI <- fread(versions[length(versions)])

# set working directory
setwd(paste0(.brd, "taxon description dates"))

str(wcvp)
str(dist)
str(auth)

# remove genera
table(wcvp$taxon_rank)
wcvp <- wcvp[!(taxon_rank %in% c("Genus"))]

# get ids of plants from Nigeria
np <- dist[grepl("Nigeria", area)]$plant_name_id # Nigeria plants

# subset wcvp for those plants (including synonyms)
wcvpn <- wcvp[plant_name_id %in% np | accepted_plant_name_id %in% np]
str(wcvpn)

# test for problematic entries
wcvpn[taxon_status != "Accepted" & !(accepted_plant_name_id %in% np)]
# 10 unplaced names found

# reduce dataset to accepted names, as these should also be the oldest ones
wcvpn <- wcvpn[taxon_status == "Accepted"]

# convert first published column to year
wcvpn[, year := as.numeric(gsub("\\D", "", sub(".*publ\\.\\s+","",first_published)))]

# test for problematic entries
wcvpn[is.na(year) & first_published != ""]
# only found entries with unknown publication year

wcvpn <- wcvpn[!is.na(year)] # remove names with unknown publication year

# show primary authors
# table(wcvpn$primary_author)

# remove "... ex" from primary author and taxon authors
wcvpn[, primary_author := sub(".*\\sex\\s+", "", primary_author)]
wcvpn[, taxon_authors := sub(".*\\sex\\s+", "", taxon_authors)]

# repair primary authors
wcvpn[grepl("A J", primary_author), primary_author := sub("A J", "A.J", primary_author)]
wcvpn[grepl("Ballard Hook.", primary_author), primary_author := sub("Ballard Hook.", "Ballard", primary_author)]
wcvpn[grepl("Le Prieur", primary_author), primary_author := sub("Le Prieur", "Leprieur", primary_author)]
wcvpn[ipni_id == "310558-1", primary_author := "Lorougnon & J.Raynal"]
wcvpn[, primary_author := sub("^Okafa$", "Okafor", primary_author)]

wcvpn <- wcvpn[primary_author != "" | taxon_authors != ""] # remove names with unknown author

# check whether all primary authors can be found in IPNI
pa <- sort(wcvpn$primary_author) # primary authors
pa <- strsplit(pa, split = ",|&")
pa <- unlist(pa)
pa <- gsub("\\s+", " ", pa)
pa <- sub("^\\s", "", pa)
pa <- unique(sub("\\s$", "", pa))

aut <- pa[!(pa %in% aNI$abbr)]
aut <- authorNormalization(aut) # normalize primary authors not found

# replace authors by normalized forms when necessary
for (i in seq_len(nrow(aut))){
	# print(wcvpn[grepl(aut$oriAuthor[i],primary_author)])
	# sometimes, replacements may be problematic, because correct and wrong versions of author
	# names may contain each other, therefore, safest is to just select entries that have wrong
	# but not correct versions of the author name in focus
	wcvpn[grepl(aut$oriAuthor[i],primary_author,fixed=TRUE) & !grepl(aut$newAuthor[i],primary_author,fixed=TRUE),
		primary_author := gsub(aut$oriAuthor[i],aut$newAuthor[i],primary_author,fixed=TRUE)]
}

# get flourishing estimate through checks of authors' appearances in taxon names
aNI[,pa := 0] # primary authoring
aNI[,aa := 0] # any authoring
aNI[,paStart := 0] # primary authoring start
aNI[,aaStart := 0] # ...
aNI[,paEnd := 0]
aNI[,aaEnd := 0]

for (i in seq_len(nrow(aNI))){
	temp <- wcvpn[grepl(paste0("(^|\\W)",aNI$abbr[i],"(\\W|$)"),taxon_authors)]
	aNI[i,aa := nrow(temp)]
	aNI[i,aaStart := min(temp$year,na.rm=TRUE)]
	aNI[i,aaEnd := max(temp$year,na.rm=TRUE)]
	temp <- temp[grepl(paste0("(^|\\W)",aNI$abbr[i],"(\\W|$)"),primary_author)]
	aNI[i,pa := nrow(temp)]
	aNI[i,paStart := min(temp$year,na.rm=TRUE)]
	aNI[i,paEnd := max(temp$year,na.rm=TRUE)]
}
checkCols <- c("pa","aa","paStart","aaStart","paEnd","aaEnd")
aNI[, (checkCols) := lapply(.SD, function(x) {
	x[!is.finite(x)] <- NA
	x[x == 0] <- NA
	x
}), .SDcols = (checkCols)]

# Think about what to extract. Remember that eventually, we need the year of the description.
# It may be enough to have the IPNI-ID of the record, but often, that will not come by itself.
# Probably, it is best to record the IPNI-ID if available and the name with authors otherwise.

# prepare results data frame

# idiosyncratic
wcvpn[plant_name_id == "1145127-az", basionym_plant_name_id := "1165611-az"]

nb <- wcvpn[grepl("\\(|\\)",taxon_authors) & basionym_plant_name_id == ""] # no basionym given
resAll <- data.table(
	accIpniid=character(),accName=character(),accAuthors=character(),
	oriIpniid=character(),oriName=character(),oriAuthors=character(),
	oriYear=numeric(),tropicosid=numeric(),authorErr=logical(),
	foundBas=logical(),foundRem=logical()
)
resAll <- rbind(resAll[0][seq_len(nrow(nb))])
resAll[,accIpniid := nb$ipni_id]
resAll[,accName := nb$taxon_name]
resAll[,accAuthors := nb$taxon_authors]
resAll[,oriAuthors := sub("\\(","",sub("\\).*","",accAuthors))]

resAll[plant_name_id == "1018043-az"]

# search for missing basionyms using POWO, IPNI and tropicos
for (i in seq_len(nrow(nb))){
	print(i)
	# try POWO
	resFull <- POWOSearch(resAll[i]$accIpniid)
	if (!all(is.na(resFull))){
		resAll[i,oriIpniid := resFull[1]]
		resAll[i,oriName := resFull[2]]
		if (resFull[3] == resAll[i]$oriAuthors) resAll[i,authorErr := FALSE] else resAll[i,authorErr := TRUE]
		resAll[i,oriAuthors := resFull[3]]
		resAll[i,oriYear := as.numeric(resFull[4])]
	} else {
		# try IPNI
		resFull <- readLines(paste0("https://www.ipni.org/n/",resAll[i]$accIpniid))
		# check whether correctly annotated basionym exists
		searchBas <- grepl("<dt>Basionym</dt>",resFull)
		resAll[i,foundBas := sum(searchBas) > 0]
		if (resAll[i]$foundBas){
			# enter data from IPNI
			res <- resFull[searchBas]
			resAll[i,oriIpniid := sub("/n/","",regmatches(res,regexpr("/n/\\d+-\\d",res)))]
			resAll[i,oriName := gsub("<[^<>]+>","",regmatches(res,regexpr("<i.*/i>",res)))]
			author <- sub(".*>\\s*","",sub(",?\\s*<.*","",regmatches(res,regexpr("/i>.*<",res))))
			if (author == resAll[i]$oriAuthors) resAll[i,authorErr := FALSE] else resAll[i,authorErr := TRUE]
			resAll[i,oriAuthors := author]
			resAll[i,oriYear := as.numeric(gsub("\\D","",regmatches(res,regexpr("\\(\\d{4}\\)",res))))]
		} else {
			# check whether badly annotated basionym exists
			searchRem <- grepl("<dt>Remarks</dt>",resFull)
			resAll[i,foundRem := sum(searchRem) > 0]
			if (resAll[i]$foundRem){
				res <- resFull[searchRem]
				res <- sub("\\(.*\\)","",res)
				res <- sub("^[[:upper:]]\\..*","",res)
				res <- sub(".*>","",sub("</dd>$","",res))
				epi1 <- sub("^\\s+","",regmatches(resAll[i]$accName,regexpr("\\s[[:lower:]]\\S+",resAll[i]$accName)))
				epi2 <- sub("^\\s+","",regmatches(res,regexpr("\\s[[:lower:]]\\S+",res))) # basionym is other genus
				epi3 <- sub("^\\.\\s+","",regmatches(res,regexpr("\\.\\s[[:lower:]]\\S+",res))) #basionym is infraspecies
				minChar <- min(nchar(epi1),nchar(epi2),nchar(epi3))
				epi1 <- substr(epi1,1,ceiling(minChar/2))
				epi2 <- substr(epi2,1,ceiling(minChar/2))
				epi3 <- substr(epi3,1,ceiling(minChar/2))
				if (length(epi2) < 1) epi2 <- ""
				if (length(epi3) < 1) epi3 <- ""
				if (epi1 == epi2 | epi1 == epi3){
					resAll[i,oriName := res] # year will be extracted later
				} else {
					# try tropicos
					res <- tropicosBasionym(resAll[i]$accName)
					if (!is.na(res[3])){
						if (res[3] == resAll[i]$oriAuthors) resAll[i,authorErr := FALSE] else resAll[i,authorErr := TRUE]
						resAll[i,oriAuthors := res[3]]
					}
					resAll[i,oriName := res[2]]
					resAll[i,tropicosid := res[1]]
					resAll[i,oriYear := as.numeric(res[4])]
					
				}
			} else {
				# try tropicos
				res <- tropicosBasionym(resAll[i]$accName)
				if (!is.na(res[3])){
					if (res[3] == resAll[i]$oriAuthors) resAll[i,authorErr := FALSE] else resAll[i,authorErr := TRUE]
					resAll[i,oriAuthors := res[3]]
				}
				resAll[i,oriName := res[2]]
				resAll[i,tropicosid := res[1]]
				resAll[i,oriYear := as.numeric(res[4])]
			}
		}
	}
}

# check whether there are problematic entries without data
resAll[is.na(oriName) & is.na(tropicosid)]

# perform details search for names with likely basionym in remarks
# POWO search
if (any(!is.na(resAll$oriName) & is.na(resAll$oriYear))){
	for (i in seq_len(nrow(resAll))){
		if (!is.na(resAll$oriName[i]) && is.na(resAll$oriYear[i])){
			print(i)
			res <- pow_search(resAll$oriName[i])$data
			res <- res[res$name == resAll$oriName[i],]
			authComp <- colSums(sapply(res$author,function(x)authorMatch(resAll[i]$oriAuthors,x)))
			if (any(authComp > 0.4)){
				res <- res[order(authComp,decreasing=TRUE)[1],]
				if (res$author == resAll[i]$oriAuthors) resAll[i,authorErr := FALSE] else resAll[i,authorErr := TRUE]
				resAll[i,oriAuthors := res$author]
				resAll[i,oriIpniid := sub(".*:","",res$url)]
				resAll[i,oriYear := POWOSearch(resAll$oriIpniid[i],basionymSearch=FALSE)]
			}
		}
	}
}

# IPNI search
if (any(!is.na(resAll$oriName) & is.na(resAll$oriYear))){
	# close all open ports
	try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)
	
	# necessary for tropicosDetailsBasionym
	rD <- rsDriver(browser = "firefox", verbose = FALSE, chromever = NULL)
	remDr <- rD[["client"]]
	remDr$navigate("https://www.ipni.org")
	
	for (i in seq_len(nrow(resAll))){
		if (!is.na(resAll$oriName[i]) && is.na(resAll$oriYear[i])){
			print(i)
			res <- ipniTableYear(resAll$oriName[i],resAll$oriAuthors[i])
			resAll[i,oriYear := res]
		}
	}
	# close all open ports
	try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)
}	

# perform details search for names with tropicosid but without synonyms
if (any(is.na(resAll$oriName) & !is.na(resAll[i]$tropicosid))){
	# close all open ports
	try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)
	
	# necessary for tropicosDetailsBasionym
	rD <- rsDriver(browser = "firefox", verbose = FALSE, chromever = NULL)
	remDr <- rD[["client"]]
	remDr$navigate("https://tropicos.org/home")
	
	for (i in seq_len(nrow(resAll))){
		if (is.na(resAll$oriName[i]) && !is.na(resAll$tropicosid[i])){
			print(i)
			res <- tropicosDetailsBasionym(resAll[i]$tropicosid)
			if (!all(is.na(res))){
				if (res[2] == resAll[i]$oriAuthors) resAll[i,authorErr := FALSE] else resAll[i,authorErr := TRUE]
			}
			resAll[i,oriName := res[1]]
			resAll[i,oriAuthors := res[2]]
			resAll[i,oriYear := res[3]]
		}
	}
	# close all open ports
	try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)
}

# check which entries remain
resAll[is.na(oriName) | is.na(oriYear)]

# HIER WEITER #
# save.image("temp.RData")
library(data.table)
library(RColorBrewer)
library(taxize)
library(RSelenium)
setwd(paste0(.brd, "taxon description dates"))
load("temp.RData")

# names not found in WCVP (irrelevant, because I just need the basionym publication years)
resAll[!(oriName %in% wcvp$nameIn)]
resAll[!(oriIpniid %in% wcvp$ipni_id)]

# add first publication year and authors based on basionyms for species who have an older basionym
wcvpn[,first_publication_date := as.numeric(gsub("\\D","",first_published))]
wcvp[,first_publication_date := as.numeric(gsub("\\D","",first_published))]
wcvpn[,first_authors := taxon_authors]
wcvp[,first_authors := taxon_authors]
# add entries with basionym ID where data can be retrieved from full dataset
bID <- wcvpn[basionym_plant_name_id != ""]$basionym_plant_name_id
setkey(wcvp,plant_name_id)
yearsTemp <- wcvp[bID]$first_publication_date
authorsTemp <- wcvp[bID]$first_authors
wcvpn[basionym_plant_name_id != "",first_publication_date := yearsTemp]
wcvpn[basionym_plant_name_id != "",first_authors := authorsTemp]
# add entries with missing basionym ID where data needed to be collected from POWO, IPNI, and tropicos
accIpniid <- wcvpn[grepl("\\(|\\)",taxon_authors) & basionym_plant_name_id == ""]$ipni_id
setkey(resAll,accIpniid)
wcvpn[grepl("\\(|\\)",taxon_authors) & basionym_plant_name_id == "", first_publication_date := resAll[accIpniid]$oriYear]
wcvpn[grepl("\\(|\\)",taxon_authors) & basionym_plant_name_id == "", first_authors := resAll[accIpniid]$oriAuthors]

# check taxon publication year for errors
# in general, the first date mentioned will be the correct, except for three dates mentioned, in which case its the last
wcvpn[(first_publication_date < 1753 | first_publication_date > 2023) & nchar(first_publication_date) == 12, first_publication_date :=
	as.numeric(sub("^\\d{8}","",first_publication_date))]
wcvpn[(first_publication_date < 1753 | first_publication_date > 2023) & nchar(first_publication_date) == 8, first_publication_date :=
	as.numeric(sub("\\d{4}$","",first_publication_date))]
wcvpn[(first_publication_date < 1753 | first_publication_date > 2023)]

# check authors for errors
wcvpn[grepl("\\(|\\)",first_authors)]
sort(unique(wcvpn$first_authors))

# remove author names before "ex"
wcvpn[, first_authors := sub (".*\\sex\\.?\\s+", "", first_authors)]

# plot descriptions over time
plot(table(wcvpn$first_publication_date),type="o",xlab="year",ylab="number of new descriptions / active authors")

# missing data in basionyms
wcvpn[is.na(first_publication_date)]
wcvpn[is.na(first_publication_date) & plant_name_id == "1055024-az",first_publication_date := 1843]
wcvpn[is.na(first_publication_date) & plant_name_id == "1159004-az",first_publication_date := 1836]

# check author contribution range (=first to last new description)
authors <- strsplit(wcvpn$first_authors,split=",|&|ex\\.?")
authors <- sort(unique(sub("^\\s+|\\s+$","",unlist(authors))))
pubRange <- data.table(author=authors,firstPub=numeric(),lastPub=numeric())
for (i in seq_along(authors)){
	pubDates <- wcvpn[grepl(authors[i],first_authors)]$first_publication_date
	pubRange[i,firstPub := min(pubDates)]
	pubRange[i,lastPub := max(pubDates)]
	pubDates <- c(wcvpn[grepl(authors[i],first_authors)]$first_publication_date,wcvpn[grepl(authors[i],taxon_authors)]$year)
	pubRange[i,firstRev := min(pubDates,na.rm=TRUE)]
	pubRange[i,lastRev := max(pubDates,na.rm=TRUE)]
}

# add number of authors active in a certain year 
# authors doing first descriptions (activity = after first and before last new description)
q1 <- as.numeric(names(table(wcvpn$first_publication_date)))
q2 <- sapply(q1,function(x)sum(pubRange$firstPub<=x & pubRange$lastPub>=x))
points(q1,q2,col="blue",type="o",lwd=2)
# authors doing first descriptions OR revisions
q1 <- as.numeric(names(table(wcvpn$first_publication_date)))
q3 <- sapply(q1,function(x)sum((pubRange$firstRev<=x|pubRange$firstPub<=x) & (pubRange$lastPub>=x|pubRange$lastRev>=x)))
points(q1,q3,col="red",type="o",lwd=2)

cor(table(wcvpn$first_publication_date),q2) # .37
cor(table(wcvpn$first_publication_date),q3) # no correlation

# the number of new descriptions is higher correlated with the number of authors
# working on first descriptions at the time than the number of authors
# working taxonomically (i.e. first descriptions and revisions)
# this is not suprising given that the high numbers of first descriptions require high numbers
# of first description authors in general

# the predictor variables in the Nigeria study are:
# 1) last year (autoregressive term)
# 2) year (linear function)
# 3) number of authors describing species in the actual year
# it is important to know that the number of species described per year is wrong in the Nigeria study,
# because they did not use the publication date of basionyms, but of the new names of certain species

# create variables
years <- seq(min(wcvpn$first_publication_date),max(wcvpn$first_publication_date))
firstAuthorsPerYear <- sapply(years,function(x) length(unique(unlist(strsplit(wcvpn[first_publication_date == x]$first_authors,split=",|&|ex\\.?")))))

points(years,firstAuthorsPerYear,col="orange")

# while I could try to redo the analysis from Bello et al., there is not much to learn here, as it suffers
# from the error that the long-term trend in species descriptions is likely driven by the descriptions from
# Linne in 1753, and this data point needs to be excluded.


