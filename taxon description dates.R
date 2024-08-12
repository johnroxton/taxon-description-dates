# These scripts use the taxon description dates that are part of the LifeGate data assembled by
# Martin Freiberg to visualize cumulative description curves and link them to potential predictors.
# Before this, data is compared between LCVP and LifeGate, and the estimation/extrapolation
# process of a paper about species numbers in Nigeria is examinated.
#
# Author: David Schellenberger Costa
###################################################################################################
# 1 Visualize data from LCVP and total numbers from LifeGate
# 2 Try to reproduce the analyses from "Trends in botanical exploration in Nigeria forecast over
#   1000 yet undescribed vascular plant species", https://doi.org/10.1093/aob/mcad106
# 3 Visualize LifeGate description dates
# 4 Approximate distribution of description dates using functions
# 5 Source distribution data from GBIF
# 6 Source body size data from Ulrich Brose's lab
# 7 Source data from BiL explorer
# 8 Source measure of public interest based on Wikipedia article lengths and Google hits
# 9 Analyse taxon description data
###################################################################################################

# 1 Visualize data from LCVP and total numbers sent by LifeGate####################################
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

# LifeGate's aggregated data does not allow for much analyses.

# 2 Try to reproduce the analyses from "Trends in botanical exploration in Nigeria forecast over###
#   1000 yet undescribed vascular plant species", https://doi.org/10.1093/aob/mcad106##############
nextScript <- NULL

# load in libraries
library(data.table)
library(RColorBrewer)
library(taxize)
library(RSelenium)
library(brms)

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxonomy"))

# load help functions
source("taxonomy help functions.R")
source("taxonomy author normalization.R")

# function to retrieve basionyms from POWO or year of first publication
POWOSearch <- function(ipniid, basionymSearch = TRUE, authCompThr = 0.4) {
	resFull <- readLines(paste0("https://powo.science.kew.org/taxon/", ipniid))
	if (basionymSearch == TRUE) {
		if (any(grepl("Homotypic Synonyms", resFull))) {
			targetAuth <- sub("\\s*</small>.*", "", sub(".*<small>\\s*", "", resFull[grepl("c-summary__heading p-xl", resFull)]))
			targetAuth <- sub("\\(", "", sub("\\).*", "", targetAuth))
			resFull <- resFull[grep("Homotypic Synonyms", resFull):length(resFull)]
			resFull <- resFull[1:grep("</ul>$", resFull)[1]]
			resFull <- resFull[grepl("^\\s*<li", resFull)]
			synAuth <- sub("\\s*</a>.*", "", sub(".*em>\\s*", "", resFull))
			authComp <- colSums(sapply(synAuth, function(x) authorMatch(targetAuth, x)))
			if (any(authComp > authCompThr)) {
				resFull <- resFull[order(authComp, decreasing = TRUE)[1]]
				oriIpniid <- sub("\".*", "", sub(".*:names:", "", resFull))
				oriName <- gsub("<[^<>]+>", "", sub("</em>[^<>]+</a>.*", "", sub(".*><em lang='la'>", "", resFull)))
				oriAuthors <- synAuth[order(authComp, decreasing = TRUE)[1]]
				oriYear <- gsub("\\D", "", regmatches(resFull, regexpr("\\(\\d{4}\\)", resFull)))
				return(c(oriIpniid, oriName, oriAuthors, oriYear))
			} else {
				return(rep(NA, 4))
			}
		} else {
			return(rep(NA, 4))
		}
	} else {
		res <- resFull[grepl("First published", resFull)]
		res <- as.numeric(gsub("\\D", "", regmatches(res, regexpr("\\(\\d{4}\\)", res))))
		return(res)
	}
}

# function to retrieve basionyms from tropicos
tropicosBasionym <- function(name, authCompThr = 0.4) {
	if (!("authorMatch" %in% ls(globalenv()))) {
		print("taxonomy help functions need to be sourced first.")
		return(c(NA, NA, NA))
	}
	# account for issues with infraspecies
	if (grepl("\\.", name)) tempName <- sub("\\s*\\S*\\..*", "", name) else tempName <- name
	foundTro <- tp_search(tempName)
	if (tempName != name) foundTro <- foundTro[foundTro$scientificname == name, ]
	if (nrow(foundTro) > 1) {
		nw1 <- length(gregexpr("\\s", name))
		nw2 <- sapply(gregexpr("\\s", foundTro$scientificname), length)
		foundTro <- foundTro[nw2 == nw1, ][1, ]
	}
	if ("nameid" %in% colnames(foundTro)) {
		res <- tp_synonyms(foundTro$nameid)
		if (res$accepted$nameid != "no syns found") {
			nchars1 <- nchar(res$accepted$scientificname)
			nchars2 <- nchar(res$accepted$scientificnamewithauthors)
			targetAuth <- substr(res$accepted$scientificnamewithauthors, nchars1 + 2, nchars2)
			targetAuth <- sub("\\(", "", sub("\\).*", "", targetAuth))
			nchars1 <- nchar(res$synonyms$scientificname)
			nchars2 <- nchar(res$synonyms$scientificnamewithauthors)
			synAuth <- substr(res$synonyms$scientificnamewithauthors, nchars1 + 2, nchars2)
			authComp <- colSums(sapply(synAuth, function(x) authorMatch(targetAuth, x)))
			if (any(authComp > authCompThr)) {
				name <- res$synonyms$scientificname[order(authComp, decreasing = TRUE)[1]]
				# account for issues with infraspecies
				if (grepl("\\.", name)) tempName <- sub("\\s*\\S*\\..*", "", name) else tempName <- name
				res <- tp_search(tempName)
				if (tempName != name) res <- res[res$scientificname == name, ]
				if (nrow(res) > 1) {
					nw1 <- length(gregexpr("\\s", name))
					nw2 <- sapply(gregexpr("\\s", res$scientificname), length)
					res <- res[nw2 == nw1, ][1, ]
				}
				return(c(foundTro$nameid, res$scientificname, res$author, res$displaydate))
			} else {
				return(c(foundTro$nameid, NA, NA, NA))
			}
		} else {
			return(c(foundTro$nameid, NA, NA, NA))
		}
	} else {
		return(rep(NA, 4))
	}
}

# function to extract the year from ipni table
ipniTableYear <- function(name, authors) {
	if (!("remDr" %in% ls(globalenv()))) {
		print("Start an rsDriver instance first.")
		return(c(NA, NA))
	}
	remDr$navigate(paste0("https://www.ipni.org/?q=", gsub(" ", "%20", name)))
	repeat {
		Sys.sleep(0.2)
		res <- remDr$findElements(using = "class", value = "result")
		if (length(res) > 0) {
			break
		} else {
			res <- remDr$findElements(using = "class", value = "no-results")
			if (length(res) > 0) {
				return(NA)
			}
		}
	}
	bas <- which(sapply(res, function(x) grepl(authors, x$getElementAttribute("outerHTML")[[1]])))[1]
	if (!is.na(bas)) {
		year <- as.numeric(gsub("\\D", "", regmatches(
			res[[bas]]$getElementAttribute("outerHTML")[[1]],
			regexpr("\\(\\d{4}\\)", res[[bas]]$getElementAttribute("outerHTML")[[1]])
		)))
		return(year)
	} else {
		return(NA)
	}
}

# function to extract the basionym from tropicos details
# very slow due to very slow tropicos website
tropicosDetailsBasionym <- function(nameid) {
	if (!("remDr" %in% ls(globalenv()))) {
		print("Start an rsDriver instance first.")
		return(c(NA, NA))
	}
	remDr$navigate(paste0("https://tropicos.org/name/", nameid))
	repeat {
		Sys.sleep(2)
		res <- remDr$findElements(using = "class", value = "control-group")
		resYear <- remDr$findElements(using = "class", value = "section-group")
		if (length(res) > 0 && length(resYear) > 0) break
	}
	for (i in seq_along(res)) {
		if (grepl("Basionym:", res[[i]]$getElementAttribute("outerHTML")[[1]])) {
			res <- res[[i + 1]]
			break
		}
	}
	for (i in seq_along(resYear)) {
		if (grepl("Published In:", resYear[[i]]$getElementAttribute("outerHTML")[[1]])) {
			resYear <- resYear[[i]]
			break
		}
	}
	if (length(res) == 1 && length(resYear) == 1) {
		author <- sub("</a>.*", "", sub(".*</i>\\s*", "", res$getElementAttribute("outerHTML")[[1]]))
		name <- gsub("</?i>", "", sub("[^<>]+</a>.*", "", sub(".*><i>\\s*", "", res$getElementAttribute("outerHTML")[[1]])))
		year <- regmatches(
			resYear$getElementAttribute("outerHTML")[[1]],
			regexpr("\\d{4}", resYear$getElementAttribute("outerHTML")[[1]])
		)
		return(c(name, author, year))
	} else {
		return(rep(NA, 3))
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
wcvpn[, year := as.numeric(gsub("\\D", "", sub(".*publ\\.\\s+", "", first_published)))]

# test for problematic entries
wcvpn[is.na(year) & first_published != ""]
# only found entries with unknown publication year

wcvpn <- wcvpn[!is.na(year)] # remove names with unknown publication year

# show primary authors
# table(wcvpn$primary_author)

# remove "... ex" from primary author and taxon authors
wcvpn[, primary_author := sub(".*\\se\\xs+", "", primary_author)]
wcvpn[, taxon_authors := sub(".*\\se\\xs+", "", taxon_authors)]

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
for (i in seq_len(nrow(aut))) {
	# print(wcvpn[grepl(aut$oriAuthor[i],primary_author)])
	# sometimes, replacements may be problematic, because correct and wrong versions of author
	# names may contain each other, therefore, safest is to just select entries that have wrong
	# but not correct versions of the author name in focus
	wcvpn[
		grepl(aut$oriAuthor[i], primary_author, fixed = TRUE) & !grepl(aut$newAuthor[i], primary_author, fixed = TRUE),
		primary_author := gsub(aut$oriAuthor[i], aut$newAuthor[i], primary_author, fixed = TRUE)
	]
}

# get flourishing estimate through checks of authors' appearances in taxon names
aNI[, pa := 0] # primary authoring
aNI[, aa := 0] # any authoring
aNI[, paStart := 0] # primary authoring start
aNI[, aaStart := 0] # ...
aNI[, paEnd := 0]
aNI[, aaEnd := 0]

for (i in seq_len(nrow(aNI))) {
	temp <- wcvpn[grepl(paste0("(^|\\W)", aNI$abbr[i], "(\\W|$)"), taxon_authors)]
	aNI[i, aa := nrow(temp)]
	aNI[i, aaStart := min(temp$year, na.rm = TRUE)]
	aNI[i, aaEnd := max(temp$year, na.rm = TRUE)]
	temp <- temp[grepl(paste0("(^|\\W)", aNI$abbr[i], "(\\W|$)"), primary_author)]
	aNI[i, pa := nrow(temp)]
	aNI[i, paStart := min(temp$year, na.rm = TRUE)]
	aNI[i, paEnd := max(temp$year, na.rm = TRUE)]
}
checkCols <- c("pa", "aa", "paStart", "aaStart", "paEnd", "aaEnd")
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

nb <- wcvpn[grepl("\\(|\\)", taxon_authors) & basionym_plant_name_id == ""] # no basionym given
resAll <- data.table(
	accIpniid = character(), accName = character(), accAuthors = character(),
	oriIpniid = character(), oriName = character(), oriAuthors = character(),
	oriYear = numeric(), tropicosid = numeric(), authorErr = logical(),
	foundBas = logical(), foundRem = logical()
)
resAll <- rbind(resAll[0][seq_len(nrow(nb))])
resAll[, accIpniid := nb$ipni_id]
resAll[, accName := nb$taxon_name]
resAll[, accAuthors := nb$taxon_authors]
resAll[, oriAuthors := sub("\\(", "", sub("\\).*", "", accAuthors))]

# search for missing basionyms using POWO, IPNI and tropicos
for (i in seq_len(nrow(nb))) {
	print(i)
	# try POWO
	resFull <- POWOSearch(resAll[i]$accIpniid)
	if (!all(is.na(resFull))) {
		resAll[i, oriIpniid := resFull[1]]
		resAll[i, oriName := resFull[2]]
		if (resFull[3] == resAll[i]$oriAuthors) resAll[i, authorErr := FALSE] else resAll[i, authorErr := TRUE]
		resAll[i, oriAuthors := resFull[3]]
		resAll[i, oriYear := as.numeric(resFull[4])]
	} else {
		# try IPNI
		resFull <- readLines(paste0("https://www.ipni.org/n/", resAll[i]$accIpniid))
		# check whether correctly annotated basionym exists
		searchBas <- grepl("<dt>Basionym</dt>", resFull)
		resAll[i, foundBas := sum(searchBas) > 0]
		if (resAll[i]$foundBas) {
			# enter data from IPNI
			res <- resFull[searchBas]
			resAll[i, oriIpniid := sub("/n/", "", regmatches(res, regexpr("/n/\\d+-\\d", res)))]
			resAll[i, oriName := gsub("<[^<>]+>", "", regmatches(res, regexpr("<i.*/i>", res)))]
			author <- sub(".*>\\s*", "", sub(",?\\s*<.*", "", regmatches(res, regexpr("/i>.*<", res))))
			if (author == resAll[i]$oriAuthors) resAll[i, authorErr := FALSE] else resAll[i, authorErr := TRUE]
			resAll[i, oriAuthors := author]
			resAll[i, oriYear := as.numeric(gsub("\\D", "", regmatches(res, regexpr("\\(\\d{4}\\)", res))))]
		} else {
			# check whether badly annotated basionym exists
			searchRem <- grepl("<dt>Remarks</dt>", resFull)
			resAll[i, foundRem := sum(searchRem) > 0]
			if (resAll[i]$foundRem) {
				res <- resFull[searchRem]
				res <- sub("\\(.*\\)", "", res)
				res <- sub("^[[:upper:]]\\..*", "", res)
				res <- sub(".*>", "", sub("</dd>$", "", res))
				epi1 <- sub("^\\s+", "", regmatches(resAll[i]$accName, regexpr("\\s[[:lower:]]\\S+", resAll[i]$accName)))
				epi2 <- sub("^\\s+", "", regmatches(res, regexpr("\\s[[:lower:]]\\S+", res))) # basionym is other genus
				epi3 <- sub("^\\.\\s+", "", regmatches(res, regexpr("\\.\\s[[:lower:]]\\S+", res))) # basionym is infraspecies
				minChar <- min(nchar(epi1), nchar(epi2), nchar(epi3))
				epi1 <- substr(epi1, 1, ceiling(minChar / 2))
				epi2 <- substr(epi2, 1, ceiling(minChar / 2))
				epi3 <- substr(epi3, 1, ceiling(minChar / 2))
				if (length(epi2) < 1) epi2 <- ""
				if (length(epi3) < 1) epi3 <- ""
				if (epi1 == epi2 || epi1 == epi3) {
					resAll[i, oriName := res] # year will be extracted later
				} else {
					# try tropicos
					res <- tropicosBasionym(resAll[i]$accName)
					if (!is.na(res[3])) {
						if (res[3] == resAll[i]$oriAuthors) resAll[i, authorErr := FALSE] else resAll[i, authorErr := TRUE]
						resAll[i, oriAuthors := res[3]]
					}
					resAll[i, oriName := res[2]]
					resAll[i, tropicosid := res[1]]
					resAll[i, oriYear := as.numeric(res[4])]
				}
			} else {
				# try tropicos
				res <- tropicosBasionym(resAll[i]$accName)
				if (!is.na(res[3])) {
					if (res[3] == resAll[i]$oriAuthors) resAll[i, authorErr := FALSE] else resAll[i, authorErr := TRUE]
					resAll[i, oriAuthors := res[3]]
				}
				resAll[i, oriName := res[2]]
				resAll[i, tropicosid := res[1]]
				resAll[i, oriYear := as.numeric(res[4])]
			}
		}
	}
}

# check whether there are problematic entries without data
resAll[is.na(oriName) & is.na(tropicosid)]

# perform details search for names with likely basionym in remarks
# POWO search
if (any(!is.na(resAll$oriName) & is.na(resAll$oriYear))) {
	for (i in seq_len(nrow(resAll))) {
		if (!is.na(resAll$oriName[i]) && is.na(resAll$oriYear[i])) {
			print(i)
			res <- pow_search(resAll$oriName[i])$data
			res <- res[res$name == resAll$oriName[i], ]
			authComp <- colSums(sapply(res$author, function(x) authorMatch(resAll[i]$oriAuthors, x)))
			if (any(authComp > 0.4)) {
				res <- res[order(authComp, decreasing = TRUE)[1], ]
				if (res$author == resAll[i]$oriAuthors) resAll[i, authorErr := FALSE] else resAll[i, authorErr := TRUE]
				resAll[i, oriAuthors := res$author]
				resAll[i, oriIpniid := sub(".*:", "", res$url)]
				resAll[i, oriYear := POWOSearch(resAll$oriIpniid[i], basionymSearch = FALSE)]
			}
		}
	}
}

# IPNI search
if (any(!is.na(resAll$oriName) & is.na(resAll$oriYear))) {
	# close all open ports
	try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)

	# necessary for tropicosDetailsBasionym
	rD <- rsDriver(browser = "firefox", verbose = FALSE, chromever = NULL)
	remDr <- rD[["client"]]
	remDr$navigate("https://www.ipni.org")

	for (i in seq_len(nrow(resAll))) {
		if (!is.na(resAll$oriName[i]) && is.na(resAll$oriYear[i])) {
			print(i)
			res <- ipniTableYear(resAll$oriName[i], resAll$oriAuthors[i])
			resAll[i, oriYear := res]
		}
	}
	# close all open ports
	try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)
}

# perform details search for names with tropicosid but without synonyms
if (any(is.na(resAll$oriName) & !is.na(resAll[i]$tropicosid))) {
	# close all open ports
	try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)

	# necessary for tropicosDetailsBasionym
	rD <- rsDriver(browser = "firefox", verbose = FALSE, chromever = NULL)
	remDr <- rD[["client"]]
	remDr$navigate("https://tropicos.org/home")

	for (i in seq_len(nrow(resAll))) {
		if (is.na(resAll$oriName[i]) && !is.na(resAll$tropicosid[i])) {
			print(i)
			res <- tropicosDetailsBasionym(resAll[i]$tropicosid)
			if (!all(is.na(res))) {
				if (res[2] == resAll[i]$oriAuthors) resAll[i, authorErr := FALSE] else resAll[i, authorErr := TRUE]
			}
			resAll[i, oriName := res[1]]
			resAll[i, oriAuthors := res[2]]
			resAll[i, oriYear := res[3]]
		}
	}
	# close all open ports
	try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)
}

# check which entries remain
resAll[is.na(oriName) | is.na(oriYear)]

# names not found in WCVP (irrelevant, because I just need the basionym publication years)
resAll[!(oriName %in% wcvp$nameIn)]
resAll[!(oriIpniid %in% wcvp$ipni_id)]

# add first publication year and authors based on basionyms for species who have an older basionym
wcvpn[, first_publication_date := as.numeric(gsub("\\D", "", first_published))]
wcvp[, first_publication_date := as.numeric(gsub("\\D", "", first_published))]
wcvpn[, first_authors := taxon_authors]
wcvp[, first_authors := taxon_authors]
# add entries with basionym ID where data can be retrieved from full dataset
bID <- wcvpn[basionym_plant_name_id != ""]$basionym_plant_name_id
setkey(wcvp, plant_name_id)
yearsTemp <- wcvp[bID]$first_publication_date
authorsTemp <- wcvp[bID]$first_authors
wcvpn[basionym_plant_name_id != "", first_publication_date := yearsTemp]
wcvpn[basionym_plant_name_id != "", first_authors := authorsTemp]
# add entries with missing basionym ID where data needed to be collected from POWO, IPNI, and tropicos
accIpniid <- wcvpn[grepl("\\(|\\)", taxon_authors) & basionym_plant_name_id == ""]$ipni_id
setkey(resAll, accIpniid)
wcvpn[
	grepl("\\(|\\)", taxon_authors) & basionym_plant_name_id == "",
	first_publication_date := resAll[accIpniid]$oriYear
]
wcvpn[grepl("\\(|\\)", taxon_authors) & basionym_plant_name_id == "", first_authors := resAll[accIpniid]$oriAuthors]

# check taxon publication year for errors
# in general, the first date mentioned will be the correct, except for three dates mentioned, in which case its the last
wcvpn[
	(first_publication_date < 1753 | first_publication_date > 2023) & nchar(first_publication_date) == 12,
	first_publication_date :=
		as.numeric(sub("^\\d{8}", "", first_publication_date))
]
wcvpn[
	(first_publication_date < 1753 | first_publication_date > 2023) & nchar(first_publication_date) == 8,
	first_publication_date := as.numeric(sub("\\d{4}$", "", first_publication_date))
]
wcvpn[(first_publication_date < 1753 | first_publication_date > 2023)]

# check authors for errors
wcvpn[grepl("\\(|\\)", first_authors)]
sort(unique(wcvpn$first_authors))

# remove author names before "ex"
wcvpn[, first_authors := sub(".*\\se\\x.?\\s+", "", first_authors)]

# plot descriptions over time
years <- range(wcvpn$first_publication_date, na.rm = TRUE)
dat <- data.table(year = seq(years[1], years[2]))
dat[, desc := sapply(year, function(x) sum(wcvpn$first_publication_date == x, na.rm = TRUE))]
plot(dat$year, dat$desc, type = "l", lwd = 3, xlab = "year", ylab = "number of new descriptions / active authors")

# missing data in basionyms
wcvpn[is.na(first_publication_date)]
wcvpn[is.na(first_publication_date) & plant_name_id == "1055024-az", first_publication_date := 1843]
wcvpn[is.na(first_publication_date) & plant_name_id == "1159004-az", first_publication_date := 1836]

# check author contribution range (=first to last new description)
authors <- strsplit(wcvpn$first_authors, split = ",|&|e\\x.?")
authors <- sort(unique(sub("^\\s+|\\s+$", "", unlist(authors))))
pubRange <- data.table(author = authors, firstPub = numeric(), lastPub = numeric())
for (i in seq_along(authors)) {
	pubDates <- wcvpn[grepl(authors[i], first_authors)]$first_publication_date
	pubRange[i, firstPub := min(pubDates)]
	pubRange[i, lastPub := max(pubDates)]
	pubDates <- c(
		wcvpn[grepl(authors[i], first_authors)]$first_publication_date,
		wcvpn[grepl(authors[i], taxon_authors)]$year
	)
	pubRange[i, firstRev := min(pubDates, na.rm = TRUE)]
	pubRange[i, lastRev := max(pubDates, na.rm = TRUE)]
}

# add number of active authors of descriptions (activity = after first and before last description)
dat[, descAuth := sapply(year, function(x) sum(pubRange$firstPub <= x & pubRange$lastPub >= x))]
lines(dat$year, dat$descAuth, col = "red", lwd = 2)

# add number of active authors of descriptions OR revisions
dat[, revAuth := sapply(year, function(x) {
	sum((pubRange$firstPub <= x | pubRange$firstRev <= x) & (pubRange$lastPub >= x | pubRange$lastRev >= x))
})]
lines(dat$year, dat$revAuth, col = "green", lwd = 2)

# the number of descriptions is higher correlated with the number of authors
# working on descriptions at the time than the number of authors
# working on descriptions and revisions
# this is not suprising given that the high numbers of first descriptions require high numbers
# of first description authors in general

# the predictor variables in the Nigeria study are:
# 1) last year (autoregressive term)
# 2) year (linear function)
# 3) number of authors describing species in the actual year
# it is important to know that the number of species described per year is wrong in the Nigeria study,
# because they did not use the publication date of basionyms, but of the new names of certain species

dat[, pubAuth := sapply(year, function(x) {
	length(unique(unlist(strsplit(wcvpn[first_publication_date == x]$first_authors, split = ",|&|e\\x.?"))))
})]
lines(dat$year, dat$pubAuth, col = "blue", lwd = 2)
dat[, descLastYear := c(0, desc[-length(desc)])]

round(cor(dat), 2)

# relationship between author total and mean number of authors with descriptions per year
length(unique(wcvpn$first_authors)) # 793
mean(dat$descAuth) # 34.7
mean(dat$pubAuth) # 7.5
793 / 34.7 # 22.85
793 / 7.5 # 105
# realistically, there are not publications every year, therefore, descAuth and pubAuth are quite different

# while I could try to redo the analysis from Bello et al., there is not much to learn here, as it suffers
# from the error that the long-term trend in species descriptions is likely driven by the descriptions from
# Linne in 1753, and this data point needs to be excluded.

# create models using the brms package
str(dat)
dat[, pubAuth := as.numeric(pubAuth)]
dat[pubAuth <= 0, pubAuth := 0.1]

dat[, desclastYear := c(0, dat$desc[-length(dat$desc)])]
# m1 <- brm(desc ~ 1, data = dat, family = negbinomial)
# m2 <- brm(desc ~ year, data = dat, family = negbinomial)
# m3 <- brm(desc ~ year + ar(p=1), data = dat, family = negbinomial, iter = 20000, control = list(adapt_delta = 0.95))
# m3.2 <- brm(desc ~ year + descLastYear, data = dat, family = negbinomial)
# m4 <- brm(desc ~ year + offset(log(pubAuth)), data = dat, family = negbinomial)
# m4.1 <- brm(desc ~ year + ar(p=1) + offset(log(pubAuth)), data = dat, family = negbinomial, iter = 20000, control = list(adapt_delta = 0.95)) # nolint: line_length_linter.
# m4.2 <- brm(desc ~ year + descLastYear + offset(log(pubAuth)), data = dat, family = negbinomial)
load("Nigeria descriptions models.RData")

# plot results
posteriorDraws <- 1000
models <- c("m1", "m2", "m3", "m3.2", "m4", "m4.1", "m4.2")
# save(list=models,file="models.RData")

pdf("Nigeria descriptions with 1753.pdf", height = 5, width = 15)
par(mar = c(0.1, 0.1, 0.1, 0.1))
par(oma = c(4, 6, 3, 1))
# plot descriptions per year
par(mfrow = c(2, 7))
responseFormula <- c(
	"rep(exp(pred[j, 1]), nrow(dat))",
	"exp(pred[j, 1] + pred[j, 2] * dat$year)",
	"exp(pred[j, 1] + pred[j, 2] * dat$year)",
	"exp(pred[j, 1] + pred[j, 2] * dat$year + pred[j, 3] * mean(dat$descLastYear))",
	"exp(pred[j, 1] + pred[j, 2] * dat$year) * mean(dat$pubAuth)",
	"exp(pred[j, 1] + pred[j, 2] * dat$year) * mean(dat$pubAuth)",
	"exp(pred[j, 1] + pred[j, 2] * dat$year + pred[j, 3] * mean(dat$descLastYear)) * mean(dat$pubAuth)"
)
for (i in seq_along(models)) {
	plot(NULL, xlim = range(dat$year), ylim = c(0, max(dat$desc)), xlab = "", xaxt = "n", ylab = "", yaxt = "n")
	mtext(models[i], 3, 1, font = 2)
	if (i < 2) {
		axis(2)
		mtext("descriptions per year", 2, 2.5)
	}
	abline(h = seq(0, max(dat$desc), 50), lty = 3)
	abline(v = seq(1750, 2000, 50), lty = 3)
	segments(dat$year, 0, dat$year, dat$desc, lwd = 2, col = "grey")
	# predictions
	pred <- fixef(eval(as.symbol(models[i])), summary = FALSE)
	for (j in seq_len(posteriorDraws)) {
		lines(dat$year, eval(parse(text = responseFormula[i])), col = "#add8e6")
	}
	lines(dat$year, eval(parse(text = gsub("pred\\[j,\\s*", paste0("fixef(", models[i], ")["), responseFormula[i]))),
		lwd = 2
	)
}
# plot total descriptions
responseFormula <- c(
	"cumsum(rep(exp(pred[j, 1]), nrow(dat)))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * dat$year))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * dat$year))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * dat$year + pred[j, 3] * mean(dat$descLastYear)))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * dat$year) * mean(dat$pubAuth))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * dat$year) * mean(dat$pubAuth))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * dat$year + pred[j, 3] * mean(dat$descLastYear)) * mean(dat$pubAuth))"
)
for (i in seq_along(models)) {
	plot(NULL, xlim = range(dat$year), ylim = c(0, 7000), xlab = "", ylab = "", yaxt = "n")
	if (i < 2) {
		axis(2)
		mtext("total descriptions", 2, 2.5)
	}
	mtext("year", 1, 2.5)
	abline(h = c(1:7) * 1e3, lty = 3)
	abline(v = seq(1750, 2000, 50), lty = 3)
	# predictions
	pred <- fixef(eval(as.symbol(models[i])), summary = FALSE)
	for (j in seq_len(posteriorDraws)) {
		lines(dat$year, eval(parse(text = responseFormula[i])), col = "#add8e6")
	}
	lines(dat$year, eval(parse(text = gsub("pred\\[j,\\s*", paste0("fixef(", models[i], ")["), responseFormula[i]))),
		lwd = 2
	)
	lines(dat$year, cumsum(dat$desc), lwd = 2)
}
dev.off()

# comparison of prediction methods and models
q1 <- predict(m3.2)
q2 <- posterior_predict(m3.2)
q2 <- colMeans(q2) # estimates
pred <- fixef(m3.2, summary = TRUE)
q3 <- exp(pred[1, 1] + pred[2, 1] * dat$year + pred[3, 1] * mean(dat$descLastYear))

q4 <- predict(m3)
q5 <- posterior_predict(m3)
q5 <- colMeans(q5) # estimates
pred <- fixef(m3, summary = TRUE)
q6 <- exp(pred[1, 1] + pred[2, 1] * dat$year)

par(mfrow = c(1, 1))
plot(NULL, xlim = range(dat$year), ylim = c(0, 100), xlab = "year", ylab = "desc")
abline(h = seq(0, max(dat$desc), 50), lty = 3)
abline(v = seq(1753, 2019, 50), lty = 3)
segments(dat$year, 0, dat$year, dat$desc, lwd = 2, col = "grey")
lines(dat$year, q1[, 1], lwd = 2)
lines(dat$year, q2, col = 2, lwd = 2)
lines(dat$year, q3, col = 3, lwd = 2)
lines(dat$year, q4[, 1], col = 4, lwd = 2)
lines(dat$year, q5, col = 5, lwd = 2)
lines(dat$year, q6, col = 6, lwd = 2)

# recalculate leaving out the first two years (1753, 1754) to avoid bias through
# Linnes work

datShort <- dat[year > 1754]
ms1 <- brm(desc ~ 1, data = datShort, family = negbinomial)
ms2 <- brm(desc ~ year, data = datShort, family = negbinomial)
ms3 <- brm(desc ~ year + ar(p = 1),
	data = datShort, family = negbinomial,
	iter = 20000, control = list(adapt_delta = 0.95)
)
ms3.2 <- brm(desc ~ year + descLastYear, data = datShort, family = negbinomial)
ms4 <- brm(desc ~ year + offset(log(pubAuth)), data = datShort, family = negbinomial)
ms4.1 <- brm(desc ~ year + ar(p = 1) + offset(log(pubAuth)),
	data = datShort,
	family = negbinomial, iter = 20000, control = list(adapt_delta = 0.95)
)
ms4.2 <- brm(desc ~ year + descLastYear + offset(log(pubAuth)),
	data = datShort,
	family = negbinomial
)

models <- c("ms1", "ms2", "ms3", "ms3.2", "ms4", "ms4.1", "ms4.2")
# save(list = models, file = "models.RData")

# plot results
posteriorDraws <- 1000
models <- c("m1", "m2", "m3", "m3.2", "m4", "m4.1", "m4.2")

pdf("Nigeria descriptions without 1753.pdf", height = 5, width = 15)
par(mar = c(0.1, 0.1, 0.1, 0.1))
par(oma = c(4, 6, 3, 1))
# plot descriptions per year
par(mfrow = c(2, 7))
responseFormula <- c(
	"rep(exp(pred[j, 1]), nrow(datShort))",
	"exp(pred[j, 1] + pred[j, 2] * datShort$year)",
	"exp(pred[j, 1] + pred[j, 2] * datShort$year)",
	"exp(pred[j, 1] + pred[j, 2] * datShort$year + pred[j, 3] * mean(datShort$descLastYear))",
	"exp(pred[j, 1] + pred[j, 2] * datShort$year) * mean(datShort$pubAuth)",
	"exp(pred[j, 1] + pred[j, 2] * datShort$year) * mean(datShort$pubAuth)",
	"exp(pred[j, 1] + pred[j, 2] * datShort$year + pred[j, 3] * mean(datShort$descLastYear)) * mean(datShort$pubAuth)"
)
for (i in seq_along(models)) {
	plot(NULL, xlim = range(datShort$year), ylim = c(0, max(datShort$desc)), xlab = "", xaxt = "n", ylab = "", yaxt = "n")
	mtext(models[i], 3, 1, font = 2)
	if (i < 2) {
		axis(2)
		mtext("descriptions per year", 2, 2.5)
	}
	abline(h = seq(0, max(datShort$desc), 50), lty = 3)
	abline(v = seq(1750, 2000, 50), lty = 3)
	segments(datShort$year, 0, datShort$year, datShort$desc, lwd = 2, col = "grey")
	# predictions
	pred <- fixef(eval(as.symbol(models[i])), summary = FALSE)
	for (j in seq_len(posteriorDraws)) {
		lines(datShort$year, eval(parse(text = responseFormula[i])), col = "#add8e6")
	}
	lines(datShort$year, eval(parse(text = gsub("pred\\[j,\\s*", paste0("fixef(", models[i], ")["), responseFormula[i]))),
		lwd = 2
	)
}
# plot total descriptions
responseFormula <- c(
	"cumsum(rep(exp(pred[j, 1]), nrow(datShort)))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * datShort$year))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * datShort$year))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * datShort$year + pred[j, 3] * mean(datShort$descLastYear)))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * datShort$year) * mean(datShort$pubAuth))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * datShort$year) * mean(datShort$pubAuth))",
	"cumsum(exp(pred[j, 1] + pred[j, 2] * datShort$year + pred[j, 3] * mean(datShort$descLastYear)) * mean(datShort$pubAuth))" # nolint: line_length_linter.
)
for (i in seq_along(models)) {
	plot(NULL, xlim = range(datShort$year), ylim = c(0, 7000), xlab = "", ylab = "", yaxt = "n")
	if (i < 2) {
		axis(2)
		mtext("total descriptions", 2, 2.5)
	}
	mtext("year", 1, 2.5)
	abline(h = c(1:7) * 1e3, lty = 3)
	abline(v = seq(1750, 2000, 50), lty = 3)
	# predictions
	pred <- fixef(eval(as.symbol(models[i])), summary = FALSE)
	for (j in seq_len(posteriorDraws)) {
		lines(datShort$year, eval(parse(text = responseFormula[i])), col = "#add8e6")
	}
	lines(datShort$year, eval(parse(text = gsub("pred\\[j,\\s*", paste0("fixef(", models[i], ")["), responseFormula[i]))),
		lwd = 2
	)
	lines(datShort$year, cumsum(datShort$desc), lwd = 2)
}
dev.off()

# 3 Visualize LifeGate description dates###########################################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets
library(rotl) # get open source phylogenies (should eventually be data from Martin)
library(ape) # plot phylogenies nicely
library(rphylopic) # get icons of taxonomic groups
library(png) # plot icons of taxonomic groups

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# read in metadata
groupsMeta <- data.table(file = list.files(path = "Arten pro Jahr beschrieben"))
groupsMeta[, name := sub(".* ", "", sub("\\s+beschriebene.*", "", file))]
# create regular expression to define groups receiving same colors (monophyletic groups of groups)
# Chordata is the only group that receives a distinguished color
groupRegs <- c(
	"[^O]omycota|sporidia", # fungi
	"Rhodophyta|Chlorophyta|Bryophyta|Marchantiophyta|Tracheophyta", # plants
	"ptera|Odo|Pso", # insects
	"Platy|Gastro|Ecto|Mollus|Nemer|Anne", # Spiralia
	# "Hapto|Myzo|Cilio|Foram|Oo|Ochro", # mostly unicellular, mostly uniflagellate
	"Chor" # chordates
)
groupCols <- c("brown", "darkgreen", "red", "lightblue", "blue") # colors of groups of groups
groupsMeta[, color := "black"] # standard color
for (i in seq_along(groupRegs)) {
	groupsMeta[grepl(groupRegs[i], name), color := groupCols[i]]
}
groupsVec <- seq_len(nrow(groupsMeta)) # create a convenience vector to process all groups
groupsMeta[, finalAuthorData := 1]
groupsMeta[grepl("Basidio|Asco|Micro|Chloro|Rhodo|Chytridio|Bryo|Cilio|Eugleno|Fora|Hapto|March|Myzo|Ochro|Oomyc", name), finalAuthorData := 0]

# read in data
groupsData <- fread(paste0("Arten pro Jahr beschrieben/", groupsMeta$file[1]))
for (i in groupsVec) {
	if (i > 1) {
		groupsData <- cbind(groupsData, fread(paste0(
			"Arten pro Jahr beschrieben/",
			groupsMeta$file[i]
		))$V2)
	}
}
colnames(groupsData) <- c("year", groupsMeta$name)
authorData <- fread(paste0(
	"AutorStat/",
	list.files(path = "AutorStat")[grepl(groupsMeta$name[1], list.files(path = "AutorStat"))]
))
for (i in groupsVec) {
	if (i > 1) {
		authorData <- cbind(authorData, fread(paste0("AutorStat/", list.files(path = "AutorStat")
		[grepl(groupsMeta$name[i], list.files(path = "AutorStat"))]))$V2)
	}
}
colnames(authorData) <- c("year", groupsMeta$name)

# reduce data to before 2018 (as it is not complete afterwards)
groupsData <- groupsData[year < 2018]
authorData <- authorData[year < 2018]

# retrieve phylogenetic information for the taxonomic groups

# retrieve OTT (Open Tree Taxonomy) IDs
# taxa <- data.table(tnrs_match_names(sub(".*\\s", "", groupsMeta$name))) # get OTT IDs from TNRS
# fwrite(taxa,file=paste0("taxa_",Sys.Date(),".txt")) # be independent of server availability
versions <- list.files(pattern = "taxa_.*\\.txt")
taxa <- fread(versions[length(versions)])
# repair some taxa
taxa[search_string == "trichoptera", ott_id := 457402] # Trichoptera, ambiguous
taxa[search_string == "myzozoa", ott_id := 435818] # Diplopsalopsis bomba
taxa[search_string == "ochrophyta", ott_id := 48614] # Phaeophyceae
taxa[search_string == "crustacea", ott_id := 212701] # Malacostraca
taxa[search_string == "chytridiomycota", ott_id := 535284] # Synchytrium endobioticum
taxa[search_string == "arachnida", ott_id := 614523] # Araneae
# make taxa search_string uppercase
taxa[, search_string := paste0(toupper(substr(search_string, 1, 1)), sub("^.", "", search_string))]
# check whether the IDs can be found in the Open Tree of Life
# ott_in_tree <- is_in_tree(taxa$ott_id)
# taxa[!ott_in_tree] # cannot be found
# the reason for taxa not to be found is that they are paraphyletic (= include other
# taxa), not monophyletic
# a workaround is to accept this and just use the node in the tree that contains all
# taxa from the group, and some others. however, this may lead to some taxa not being nodes
# anymore. it is therefore better to replace taxa with species belonging to them (or other lower, more
# specific taxa ranks)

# create phylogenetic tree
# pt <- tol_induced_subtree(ott_id = taxa$ott_id, label_format = "name")
# save(pt,file=paste0("pt_",Sys.Date(),".RData"))
versions <- list.files(pattern = "pt_.*\\.RData")
load(versions[length(versions)])
# dropping singleton nodes means that higher taxa that do not branch along the chosen
# phylogeny are not labelled
# check whether tree has correct number of nodes
# if not, one group is no tip, but just a node
if (length(pt$tip.label) < nrow(taxa)) {
	for (i in seq_along(taxa$ott_id[opt_in_tree])) {
		if (i > 1) {
			pt <- tol_induced_subtree(ott_id = taxa$ott_id[opt_in_tree][1:i])
			if (length(pt$tip.label) != i) break
		}
	}
}
# repair phylogenetic tree for correct fan/radial display
pt$edge <- pt$edge[order(pt$edge[, 2]), ]
# correct tip labels
# repair non-matching tip labels and unique_name values
pt$tip.label[!(pt$tip.label %in% taxa$unique_name)]
taxa[!(unique_name %in% pt$tip.label)]
pt$tip.label[!(pt$tip.label %in% taxa$unique_name)] <- c(
	"Ochrophyta (merged with Stramenopiles)",
	"Ciliophora (phylum in subkingdom SAR)",
	"Myzozoa",
	"Trichoptera (genus in order Diptera)",
	"Plecoptera (order in cohort Polyneoptera)",
	"Dictyoptera (genus in cohort Polyneoptera)",
	"Crustacea",
	"Arachnida",
	"Chytridiomycota"
)
# use search_string instead of unique_name as tip labels
pt$tip.label <- taxa$search_string[sapply(pt$tip.label, function(x) which(x == taxa$unique_name))]
# re-sort to show groups in order
setcolorder(groupsData, c(1, sapply(pt$tip.label, function(x) which(groupsMeta$name == x)) + 1))
setcolorder(authorData, c(1, sapply(pt$tip.label, function(x) which(groupsMeta$name == x)) + 1))
groupsMeta <- groupsMeta[sapply(pt$tip.label, function(x) which(sub(".* ", "", groupsMeta$name) == x))]

# write data for later use
fwrite(groupsMeta, file = "groupsMeta.txt")
fwrite(groupsData, file = "groupsData.txt")
fwrite(authorData, file = "authorData.txt")

# get icons of groups
# done manually from phylopic.org
# groupsMeta$name
# groupsMeta[1, icon :=
# 	"1b329337-8f8b-4380-8aae-23d50e0db14f",
# 	"4de6d5a6-bae9-432a-bb39-546791224857",
# 	"cc68fe27-0e9a-442d-a28d-0de60e30869e",
# 	"b31d2f98-7031-4bf3-8229-e84403a656df",
# 	"6a3bd364-ea78-4f1a-be98-0fe18e022b4e",
# 	"5a316192-4559-40b2-844d-47f622ae201e",
# 	"53b543aa-2a5e-4407-88cd-6acbb0b5c3f9",
# 	"0e450314-234c-40ae-bca2-c426ff655224",
# 	"456d3821-2308-4842-9bb3-d82620cdf101",
# 	"b21680aa-d2e0-4615-989a-76544f43c2e2",
# 	"07a536aa-4c8b-4089-b89b-3d69948a136f",
# 	"a862dac7-1e49-4f51-b67b-4c974ec29ab2",
# 	"a039e19e-a533-4ed9-92ad-ddba6da7801f",
# 	"dd7c7bc6-2c6d-48e4-860b-02bf37886b5b",
# 	"22f8d558-6763-426d-a651-f4087090fc41",
# 	"a3fb3d98-60a5-49bd-bf09-77cfdda1bfab",
# 	"191ad6ce-78f2-444d-b4ad-9c4162b1803f",
# 	"6c2e67f0-14e7-4ba0-ba73-2420cacfa9a3",
# 	"429429ed-d07d-463f-9316-62c7adde7e71",
# 	"5b80cf50-11f0-4e9a-b23b-452fac47615b",
# 	"660457de-6a6f-4f23-9611-20e2c8e950f8",
# 	"e1307c88-3e8f-4ba8-9f93-751df3deb739",
# 	"58eda733-f5dd-4a4d-8008-c85a333570b0",
# 	"8537a53a-541a-47af-a1db-b97e8a4bc81d",
# 	"7832cac2-f122-4112-961e-2d506789e1c1",
# 	"5053ea6b-faf3-438a-82b1-b9b742eed5a0",
# 	"b199a5f5-20c4-4cc9-9c54-1b51578c2487",
# 	"36a06b66-6eda-44a9-ae5f-cf28fcba7023",
# 	"37a0220e-de80-4ce2-aff2-ff3bf218ce52",
# 	"10742785-cc4f-4a08-bd7d-b62d1959f554",
# 	"a98546d9-e030-4501-85e6-e9af67f63b02",
# 	"3fd92960-bac6-4cd4-bcf6-584c27ee2758",
# 	"05cd564e-b3f4-45e1-95e6-5f4c14bad001",
# 	"e7ba011c-6172-47ca-84bf-0dfc87b47a32",
# 	"e642c1c1-76c1-4806-a024-aa9737d5bc41",
# 	"36543d34-f5d6-48a1-818f-eca72bfe1a19",
# 	"d96e18f2-c9d6-4b2e-914c-678ed3c72a28",
# 	"3ca7d0e6-4c15-48bb-b770-5da52d548197",
# 	"acc64610-abae-4bea-a9a5-46e2636145b9",
# 	"0957b8c0-a052-4f7d-a1ac-95f2179c6582",
# 	"b1146efb-589e-4046-8370-c712f7494377",
# 	"839b9df7-c97f-444a-a611-95609e950960",
# 	"9c529339-9a35-4c3b-9855-793fcd151ca1",
# 	"f722ba34-03d4-4b90-8409-e02829a3c5d4",
# 	"f96b5420-8b52-4b2e-aa79-9bb16d25adca",
# 	"8d07d74b-eafc-42bd-bd0a-120139d3cf36",
# 	"05def04a-f550-4b77-bdbe-d7c25a22f229"
# )]
# i <- 1
# get and save icons
# for (i in groupsVec){
# 	save_phylopic(get_phylopic(uuid = groupsMeta$icon[i]),path=paste0("icons/",groupsMeta$name[i],".svg"))
# }

parBackup <- par() # save graphics parameters for multiple trials

# plot overview phylogeny

# choose colors for monophyletic groups
# for coloring purposes change node labels into numbers (can be omitted after first run)
# pt$node.label <- length(pt$tip.label) + 1:length(pt$node.label)
xlim <- 1
ylim <- 1.6
cols <- rep("black", nrow(pt$edge))
groupCols <- c("brown", "darkgreen", "red", "lightblue", "green", "blue")
groupNode <- c(90, 52, 76, 66, 56, grep("Chor", pt$tip.label))
for (i in seq_along(groupCols)) {
	edgesOld <- c()
	edgesNew <- groupNode[i]
	while (length(edgesNew) != length(edgesOld) || any(edgesNew != edgesOld)) {
		edgesOld <- edgesNew
		edgesNew <- union(edgesOld, pt$edge[pt$edge[, 1] %in% edgesOld, 2])
	}
	cols[pt$edge[, 2] %in% edgesNew] <- groupCols[i]
}
cols[cols == "green"] <- "black" # remove unlabeled supergroup

# set show.node.label=FALSE to TRUE for coloring purposes
# pdf("groups phylogeny.pdf")
par(parBackup)
plot.phylo(pt,
	show.node.label = FALSE, cex = 0.7, font = 1, type = "fan", x.lim = c(-xlim, xlim),
	y.lim = c(-ylim, ylim), label.offset = 0.2, edge.color = cols, edge.width = 2
)
i <- 1
tipIncrease <- 1.1 # tip length increase for figure plotting
imageSize <- 0.07 # figure size factor
for (i in groupsVec) {
	segments(cos(2 * pi * (i - 1) / nrow(groupsMeta)), sin(2 * pi * (i - 1) / nrow(groupsMeta)),
		tipIncrease * cos(2 * pi * (i - 1) / nrow(groupsMeta)), tipIncrease * sin(2 * pi * (i - 1) / nrow(groupsMeta)),
		lwd = 2, col = groupsMeta$color[i]
	)
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	colRGB <- col2rgb(groupsMeta$color[i]) / 255
	for (j in 1:3) img[, , j] <- colRGB[[j]]
	rasterImage(img, tipIncrease * cos(2 * pi * (i - 1) / nrow(groupsMeta)) - imageSize,
		tipIncrease * sin(2 * pi * (i - 1) / nrow(groupsMeta)) - imageSize,
		tipIncrease * cos(2 * pi * (i - 1) / nrow(groupsMeta)) + imageSize,
		tipIncrease * sin(2 * pi * (i - 1) / nrow(groupsMeta)) + imageSize,
		col = groupsMeta$color[i]
	)
}
# dev.off()

# plot description history
polyYear <- c(groupsData$year, rev(groupsData$year))
zeros <- rep(0, nrow(groupsData))

# pdf("description history.pdf",width=11.7,height=8.3)
# plot all groups in one sheet
par(oma = c(5, 5, 1, 1))
par(mar = c(0, 0, 0, 0))
par(mfrow = c(5, 10))
for (i in groupsVec) {
	plot(NULL, xlim = range(groupsData$year), ylim = c(0, max(log(groupsData + 1))), xaxt = "n", yaxt = "n")
	abline(h = log(c(0, 10, 100, 1000) + 1), lty = 3)
	abline(v = c(1800, 1900, 2000), lty = 3)
	if ((i - 1) %% 10 < 1) axis(2, at = log(c(0, 10, 100, 1000) + 1), labels = c(0, 10, 100, 1000))
	if (i > 37) axis(1, at = c(1800, 1900, 2000))
	polygon(polyYear, c(log(groupsData[[i + 1]] + 1), zeros), border = NA, col = groupsMeta$col[i])
	text(1800, log(5001), labels = gsub("\\s", "\n", groupsMeta$name[i]), adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, max(log(groupsData + 1)))))
	rasterImage(img, 1760, log(2001), 1790, log(2001) + scale)
}
mtext("year", 1, 3, at = 1450)
mtext("descriptions", 2, 52.5, at = 23, xpd = TRUE)
# plot all groups in one sheet with cumulative relative values
par(oma = c(5, 5, 1, 1))
par(mar = c(0, 0, 0, 0))
par(mfrow = c(5, 10))
for (i in groupsVec) {
	plot(NULL, xlim = range(groupsData$year), ylim = c(0, 1), xaxt = "n", yaxt = "n")
	abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 3)
	abline(v = c(1800, 1900, 2000), lty = 3)
	if ((i - 1) %% 10 < 1) axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
	if (i > 37) axis(1, at = c(1800, 1900, 2000))
	polygon(polyYear, c(cumsum((groupsData[[i + 1]])) / sum(groupsData[[i + 1]]), zeros),
		border = NA, col = groupsMeta$col[i]
	)
	text(1800, log(5001) / max(log(groupsData + 1)), labels = groupsMeta$name[i], adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, 1)))
	rasterImage(img, 1760, log(2001) / max(log(groupsData + 1)), 1790, log(2001) / max(log(groupsData + 1)) + scale)
}
mtext("year", 1, 3, at = 1450)
mtext("cumulative descriptions", 2, 52.5, at = 2.7, xpd = TRUE)
# plot all groups separately
par(parBackup)
for (i in groupsVec) {
	plot(NULL,
		xlim = range(groupsData$year), ylim = c(0, max(log(groupsData + 1))), xaxt = "n", yaxt = "n",
		xlab = "year", ylab = "descriptions"
	)
	abline(h = log(c(0, 10, 100, 1000) + 1), lty = 3)
	abline(v = c(1800, 1900, 2000), lty = 3)
	axis(2, at = log(c(0, 10, 100, 1000) + 1), labels = c(0, 10, 100, 1000))
	axis(1, at = c(1800, 1900, 2000))
	polygon(polyYear, c(log(groupsData[[i + 1]] + 1), zeros), border = NA, col = groupsMeta$col[i])
	text(1800, log(5001), labels = groupsMeta$name[i], adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	scale <- diff(range(groupsData$year)) / diff(c(0, max(log(groupsData + 1)))) / 30
	rasterImage(img, 1760, log(2001), 1790, log(2001) + scale)
}
# plot all groups separately with cumulative relative values
par(parBackup)
i <- 1
for (i in groupsVec) {
	plot(NULL,
		xlim = range(groupsData$year), ylim = c(0, 1), xaxt = "n", yaxt = "n",
		xlab = "year", ylab = "cumulative descriptions"
	)
	abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 3)
	abline(v = c(1800, 1900, 2000), lty = 3)
	axis(2)
	axis(1, at = c(1800, 1900, 2000))
	polygon(polyYear, c(cumsum((groupsData[[i + 1]] + 1)) / sum(groupsData[[i + 1]] + 1), zeros),
		border = NA, col = groupsMeta$col[i]
	)
	text(1800, log(5001) / max(log(groupsData + 1)), labels = groupsMeta$name[i], adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, 1)))
	rasterImage(img, 1760, log(2001) / max(log(groupsData + 1)), 1790, log(2001) / max(log(groupsData + 1)) + scale)
}
# dev.off()

# 4 Approximate distribution of description dates using functions##################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets
library(rethinking) # approximate description times with functions

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# get data
groupsMeta <- fread("groupsMeta.txt")
groupsData <- fread("groupsData.txt")
groupsVec <- seq_len(nrow(groupsMeta)) # create a convenience vector to process all groups

# The data should be approximated by a function of two/three parameters, which could be:
# (1) technical ability to describe (i.e. some very small species), cutpoint with y-axis
# (2) estimated number of species in group, upper limit
# (3) estimated effort/(1/easiness) of discovery invested in group
# Effort or easiness of discovery cannot be distinguished from curves, they can only be
# inferred from external data, like number of researchers known to work on a specific group,
# or, as a proxy of interest, occurrence in popular literature.

parBackup <- par() # save graphics parameters for multiple trials

# function to draw cumulative descriptions
drawit <- function() {
	par(parBackup)
	plot(
		NULL,
		xlim = range(groupsData$year), ylim = c(0, 1), xaxt = "n", yaxt = "n",
		xlab = "year", ylab = "relative cumulative descriptions"
	)
	abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 3)
	abline(v = c(1800, 1900, 2000), lty = 3)
	axis(2)
	axis(1, at = c(1800, 1900, 2000))
	for (i in groupsVec) {
		lines(groupsData$year, cumsum((groupsData[[i + 1]])) / sum(groupsData[[i + 1]]), col = groupsMeta$col[i])
	}
}

# candidates for modelling

# normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
# vary standard deviation
drawit()
for (i in seq(50, 90, 10)) {
	q1 <- dnorm(seq(1753, 2022), 2022, i)
	q1 <- q1 / max(q1)
	lines(seq(1753, 2022), q1, lwd = 3, lty = 2)
}
# vary mean
for (i in seq(2022, 2072, 10)) {
	q1 <- dnorm(seq(1753, 2022), i, 80)
	q1 <- q1 / max(q1)
	lines(seq(1753, 2022), q1, lwd = 3)
}

# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
drawit()
for (i in seq(2.4, 3.8, l = 10)) {
	q1 <- seq(1753, 2022)
	q1 <- 3 * (1 - exp(-0.004 * (q1 - 1753)))^i
	lines(seq(1753, 2022), q1, lwd = 3)
}

# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
drawit()
for (i in seq(3, 16)) {
	q1 <- seq(1753, 2022)
	q1 <- exp(-i * exp(-.015 * (q1 - 1753)))
	lines(seq(1753, 2022), q1, lwd = 3, lty = 2)
}
for (i in seq(0.014, 0.02, l = 10)) {
	q1 <- seq(1753, 2022)
	q1 <- exp(-16 * exp(-i * (q1 - 1753)))
	lines(seq(1753, 2022), q1, lwd = 5)
}

# modelling

# test Bayesian models with dinosaurs data using the functions used with the dinosaur example

# prepare data
# dinosaurs
data(Dinosaurs)
d <- Dinosaurs
dd <- d[d$sp_id == 1, ]
dat4 <- list(
	A = dd$age,
	M = dd$mass / max(dd$mass) # normalize to [0,1]
)
A4seq <- seq(from = 0, to = 16, len = 50)
# descriptions
datn <- list(
	A = groupsData$year - min(groupsData$year),
	M = cumsum(groupsData$Tracheophyta) / max(cumsum(groupsData$Tracheophyta)) # normalize to [0,1]
)
A1seq <- seq(from = 0, to = max(groupsData$year - min(groupsData$year)), len = 50)

## linear function
#
## dinosaurs
# m4a <- ulam(
# 	alist(
# 		M ~ normal(mu, sigma),
# 		mu <- a + b * A,
# 		a ~ normal(0, 1),
# 		b ~ normal(0, 1),
# 		sigma ~ exponential(1)
# 	),
# 	data = dat4, chains = 4, log_lik = TRUE
# )
## descriptions
# m1a <- ulam(
# 	alist(
# 		M ~ normal(mu, sigma),
# 		mu <- a + b * A,
# 		a ~ normal(0, 1),
# 		b ~ normal(0, 1),
# 		sigma ~ exponential(1)
# 	),
# 	data = datn, chains = 4, log_lik = TRUE
# )
#
## negative exponential function
#
## dinosaurs
# m4b <- ulam(
# 	alist(
# 		M ~ normal(mu, sigma),
# 		mu <- k * (1 - exp(-b * A)),
# 		b ~ exponential(1),
# 		k ~ normal(1, 0.5),
# 		sigma ~ exponential(1)
# 	),
# 	data = dat4, chains = 4, log_lik = TRUE
# )
## descriptions
# m1b <- ulam(
# 	alist(
# 		M ~ normal(mu, sigma),
# 		mu <- k * (1 - exp(-b * A)),
# 		b ~ exponential(1),
# 		k ~ normal(1, 0.5),
# 		sigma ~ exponential(1)
# 	),
# 	data = datn, chains = 4, log_lik = TRUE
# )
#
## custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
#
## dinosaurs
# m4c <- ulam(
# 	alist(
# 		M ~ normal(mu, sigma),
# 		mu <- k * (1 - exp(-b * A))^a,
# 		a ~ exponential(0.1),
# 		b ~ exponential(1),
# 		k ~ normal(1, 0.5),
# 		sigma ~ exponential(1)
# 	),
# 	data = dat4, chains = 4, log_lik = TRUE
# )
## descriptions
# m1c <- ulam(
# 	alist(
# 		M ~ normal(mu, sigma),
# 		mu <- k * (1 - exp(-b * A))^a,
# 		a ~ exponential(0.1),
# 		b ~ exponential(2),
# 		k ~ normal(2, 0.5),
# 		sigma ~ exponential(1)
# 	),
# 	data = datn, chains = 4, log_lik = TRUE
# )
load("models1.RData")

# compare the three approximations for dinosaurs and description data
# pdf("Bayesian modelling test.pdf",width=11.7,height=8.3)
par(mfrow = c(3, 2))
for (i in 1:3) {
	plot(
		dat4$A, dat4$M,
		xlab = "age", ylab = "mass (normalized)",
		lwd = 3, xlim = c(0, 16), col = "green"
	)
	mu <- link(eval(as.symbol(c("m4a", "m4b", "m4c")[i])), data = list(A = A4seq))
	lines(A4seq, apply(mu, 2, mean), lwd = 3)
	shade(apply(mu, 2, PI), A4seq)
	plot(
		datn$A, datn$M,
		xlab = "year", ylab = "cumulative descriptions", type = "l",
		lwd = 3, xaxt = "n", xlim = c(0, max(groupsData$year - min(groupsData$year))), col = "green"
	)
	axis(1, at = 47 + c(0, 100, 200), labels = c(1800, 1900, 2000))
	mu <- link(eval(as.symbol(c("m1a", "m1b", "m1c")[i])), data = list(A = A1seq))
	lines(A1seq, apply(mu, 2, mean), lwd = 3)
	shade(apply(mu, 2, PI), A1seq)
}
# dev.off()
# save(list=c("m1a","m1b","m1c","m4a","m4b","m4c"),file="models1.RData")

# The algorithm works. There may be some issues wit details that have to be addressed, but overall,
# this method can be implemented. The third model matches the data very well. It only predicts
# a large number of species to be described over the next hundred years, but this may not be a
# big issue at the moment. It may well be impossible to model this without extra assumptions.

# try different functions for prediction

# function draw cumulative descriptions of Tracheophyta
drawit <- function() {
	par(parBackup)
	plot(
		datn$A, datn$M,
		xlab = "year", ylab = "relative cumulative descriptions", type = "l",
		lwd = 3, xlim = c(0, 500), ylim = c(0, 2), col = "green", xaxt = "n"
	)
	abline(h = (1:5) * 2 / 5, lty = 3)
	abline(v = (0:5) * 100 + 47, lty = 3)
	axis(1, at = c(47, 147, 247), labels = c(1800, 1900, 2000))
}

# normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)

# predictive prior simulation
drawit()
a <- 75
b <- 280
k <- 200
mu <- k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((A1seq - b) / a)^2)
lines(A1seq, mu, col = "red", lwd = 3)

# run model
# m2c <- ulam(
# 	alist(
# 		M ~ normal(mu, sigma),
# 		mu <- k * 1 / ((2 * 3.141593)^0.5 * a) * exp(-0.5 * ((A - b) / a)^2),
# 		a ~ exponential(0.1),
# 		b ~ exponential(0.1),
# 		k ~ exponential(0.1),
# 		sigma ~ exponential(1)
# 	),
# 	data = datn, chains = 4, log_lik = TRUE
# )

# Gompertz function: f(x) = k * exp(-a * exp(-b * x))

# predictive prior simulation
drawit()
k <- 1
a <- 3
b <- 1 / 100
mu <- k * exp(-a * exp(-b * A1seq))
lines(A1seq, mu, col = "red")

# run model
# m3c <- ulam(
# 	alist(
# 		M ~ normal(mu, sigma),
# 		mu <- k * exp(-a * exp(-b * A)),
# 		k ~ normal(1.5, 2),
# 		a ~ exponential(0.1),
# 		b ~ exponential(1),
# 		sigma ~ exponential(1)
# 	),
# 	data = datn, chains = 4, log_lik = TRUE
# )
load("models2.RData")

# compare the three approximation functions
# pdf("approximation functions test.pdf",width=11.7,height=8.3)
plot(
	datn$A, datn$M,
	xlab = "year", ylab = "cumulative descriptions", type = "l",
	lwd = 3, xaxt = "n", xlim = c(0, 500), ylim = c(0, 2)
)
axis(1, at = 47 + c(0, 100, 200, 300, 400), labels = c(1800, 1900, 2000, 2100, 2200))
A1seq <- seq(from = 0, to = 500, len = 50)
for (i in 1:3) {
	mu <- link(eval(as.symbol(c("m1c", "m2c", "m3c")[i])), data = list(A = A1seq))
	lines(A1seq, apply(mu, 2, mean), lwd = 3, col = i + 1)
	shade(apply(mu, 2, PI), A1seq)
}
# dev.off()
# save(list=c("m1c","m2c","m3c"),file="models2.RData")

# All three approximations work nicely with this first trial. However, there are issues
# with the the values of the priors. Do not forget that exponential(1) specifies the
# exponential distribution, smaller values flatten it, larger ones make it tighter.

# create Bayesian models for all groups

# plot prior distributions
plot(NULL, xlim = c(0, 10), ylim = c(0, 1))
xseq <- seq(0, 10, l = 100)
counter <- 0
for (i in seq_len(5)) {
	counter <- counter + 1
	yseq <- dnorm(xseq, 0, i)
	lines(xseq, yseq / max(yseq), col = counter + 1)
}
counter <- 0
for (i in seq(0.1, 0.5, l = 5)) {
	counter <- counter + 1
	yseq <- dexp(xseq, i)
	lines(xseq, yseq / max(yseq), col = counter + 1, lty = 2)
}

# m <- list()
# coefs <- list()
#
# i <- 5
## custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
# m$m1 <- list()
# for (i in groupsVec) {
# for (i in which(coefs[[1]][1,]>6)) {
# 	print(paste0("Bertalanffy model ", i, "/", nrow(groupsMeta)))
# 	datn <- list(
# 		Y = groupsData$year - min(groupsData$year),
# 		D = cumsum(groupsData[[i + 1]]) / max(cumsum(groupsData[[i + 1]]))
# 	)
# 	m$m1[[i]] <- ulam(
# 		alist(
# 			D ~ normal(mu, sigma),
# 			mu <- k * (1 - exp(-b * Y))^a,
# 			k ~ exponential(5),
# 			a ~ exponential(1),
# 			b ~ exponential(1),
# 			sigma ~ exponential(1)
# 		),
# 		data = datn, chains = 4, log_lik = TRUE
# 	)
# 	print(precis(m$m1[[i]]))
# 	# traceplot(q1,pars=c("k","a","b","sigma"),n_cols=2)
# }
# coefs[[1]] <- sapply(m$m1, coef)
## normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
# m$m2 <- list()
# for (i in groupsVec) {
# 	print(paste0("Normal model ", i, "/", nrow(groupsMeta)))
# 	datn <- list(
# 		Y = groupsData$year - min(groupsData$year),
# 		D = cumsum(groupsData[[i + 1]]) / max(cumsum(groupsData[[i + 1]]))
# 	)
# 	datn$Y <- standardize(datn$Y) # get more stable results
# 	m$m2[[i]] <- ulam(
# 		alist(
# 			D ~ normal(mu, sigma),
# 			mu <- k * 1 / ((2 * 3.141593)^0.5 * a) * exp(-0.5 * ((Y - b) / a)^2),
# 			k ~ exponential(1),
# 			a ~ exponential(1),
# 			b ~ exponential(1),
# 			sigma ~ exponential(1)
# 		),
# 		data = datn, chains = 4, log_lik = TRUE
# 	)
# 	#precis(m$m2[[i]])
# 	#precis(test)
# 	#
# 	#traceplot(m$m2[[i]],pars=c("k","a","b","sigma"),n_cols=2)
# 	#traceplot(test,pars=c("k","a","b","sigma"),n_cols=2)
# }
# coefs[[2]] <- sapply(m$m2, coef)
#
## Gompertz function: f(x) = k * exp(-a * exp(-b * x))
# m$m3 <- list()
# for (i in groupsVec) {
# 	print(paste0("Gompertz model ", i, "/", nrow(groupsMeta)))
# 	datn <- list(
# 		Y = groupsData$year - min(groupsData$year),
# 		D = cumsum(groupsData[[i + 1]]) / max(cumsum(groupsData[[i + 1]]))
# 	)
# 	m$m3[[i]] <- ulam(
# 		alist(
# 			D ~ normal(mu, sigma),
# 			mu <- k * exp(-a * exp(-b * Y)),
# 			k ~ exponential(1),
# 			a ~ exponential(0.1),
# 			b ~ exponential(1),
# 			sigma ~ exponential(1)
# 		),
# 		data = datn, chains = 4, log_lik = TRUE
# 	)
# 	# precis(m$m3[[i]])
# 	# traceplot(m$m3[[i]],pars=c("k","a","b","sigma"),n_cols=2)
# }
# coefs[[3]] <- sapply(m$m3, coef)
load("models3.RData")

# plot all groups separately with cumulative relative values and model approximations
# added information:
# - predicted asymptote/maximum
# - year of predicted asymptote/maximum
# - approximation quality

# plot description history
polyYear <- c(groupsData$year, rev(groupsData$year))
zeros <- rep(0, nrow(groupsData))

i <- 7
par(parBackup)
pdf("description history fits.pdf", width = 11.7, height = 8.3)
# plot all groups separately with all approximation functions
xseq <- seq(from = 0, to = max(groupsData$year) - min(groupsData$year), len = 50)
for (i in groupsVec) {
	plot(NULL,
		xlim = range(groupsData$year), ylim = c(0, 1), xaxt = "n", yaxt = "n",
		xlab = "year", ylab = "cumulative descriptions"
	)
	abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 3)
	abline(v = c(1800, 1900, 2000), lty = 3)
	axis(2)
	axis(1, at = c(1800, 1900, 2000))
	polygon(polyYear, c(cumsum((groupsData[[i + 1]])) / sum(groupsData[[i + 1]]), zeros),
		border = NA, col = groupsMeta$col[i]
	)
	text(1800, log(5001) / max(log(groupsData + 1)), labels = groupsMeta$name[i], adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, 1)))
	rasterImage(img, 1760, log(2001) / max(log(groupsData + 1)), 1790, log(2001) / max(log(groupsData + 1)) + scale)
	for (j in 1:3) {
		if (j != 2) {
			mu <- link(m[[j]][[i]], data = list(Y = xseq))
		} else {
			# for normal distributions, it was necessary to scale the year variable to get good results
			mu <- link(m[[j]][[i]], data = list(Y = standardize(xseq)))
		}
		lines(xseq + min(groupsData$year), apply(mu, 2, mean), lwd = 3, col = j + 1)
		shade(apply(mu, 2, PI), xseq, col = j + 1)
	}
	legend("bottomright",
		legend = c("custom Bertalanffy", "normal distribution", "Gompertz"),
		lwd = 3, col = 2:4, bg = "white"
	)
}
# plot all on one page only showing the normal distribution approximation
par(oma = c(5, 5, 1, 1))
par(mar = c(0, 0, 0, 0))
par(mfrow = c(5, 10))
for (i in groupsVec) {
	plot(NULL, xlim = range(groupsData$year), ylim = c(0, 1), xaxt = "n", yaxt = "n")
	abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 3)
	abline(v = c(1800, 1900, 2000), lty = 3)
	if ((i - 1) %% 10 < 1) axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
	if (i > 37) axis(1, at = c(1800, 1900, 2000))
	polygon(polyYear, c(cumsum((groupsData[[i + 1]])) / sum(groupsData[[i + 1]]), zeros),
		border = NA, col = groupsMeta$col[i]
	)
	text(1800, log(5001) / max(log(groupsData + 1)), labels = groupsMeta$name[i], adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, 1)))
	rasterImage(img, 1760, log(2001) / max(log(groupsData + 1)), 1790, log(2001) / max(log(groupsData + 1)) + scale)
	# for normal distributions, it was necessary to scale the year variable to get good results
	mu <- link(m[[2]][[i]], data = list(Y = standardize(xseq)))
	if (groupsMeta$col[i] %in% c("brown", "red")) {
		lines(xseq + min(groupsData$year), apply(mu, 2, mean), lwd = 3, col = "black")
	} else {
		lines(xseq + min(groupsData$year), apply(mu, 2, mean), lwd = 3, col = "red")
	}
	shade(apply(mu, 2, PI), xseq, col = j + 1)
}
mtext("year", 1, 3, at = 1450)
mtext("cumulative descriptions", 2, 52.5, at = 2.7, xpd = TRUE)
dev.off()

# interpret coefficients

# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
# a = difference in start of descriptions
# ~ (proportional)
# end only relevant if maximum close to be reached
# b = description pace
# ~
# par(mfrow = c(2, 2))
# range(coefs[[1]][2, ])
# xseq <- seq(0, 10, l = 100)
# a <- seq(0, 150, l = 20)
# b <- 1
# plot(NULL, xlim = range(xseq), ylim = c(0, 1), xlab = "", ylab = "")
# for (i in a) {
# 	yseq <- (1 - exp(-b * xseq))^i
# 	lines(xseq, yseq, lwd = 3)
# 	abline(v = log(1 / i) / -b)
# 	Sys.sleep(0.25)
# }
# range(coefs[[1]][3, ])
# xseq <- seq(0, 1000, l = 100)
# a <- 1
# b <- seq(0, 0.03, l = 20)
# plot(NULL, xlim = range(xseq), ylim = c(0, 1), xlab = "", ylab = "")
# for (i in b) {
# 	yseq <- (1 - exp(-i * xseq))^a
# 	lines(xseq, yseq)
# 	abline(v = log(1 / a) / -i)
# 	Sys.sleep(0.25)
# }

# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
# a = difference in start of descriptions
# ~ (proportional)
# end only relevant if maximum close to be reached
# b = description pace
# ~
# range(coefs[[3]][2, ])
# xseq <- seq(0, 10, l = 100)
# a <- seq(1, 200, l = 20)
# b <- 1
# plot(NULL, xlim = range(xseq), ylim = c(0, 1), xlab = "", ylab = "")
# for (i in a) {
# 	yseq <- exp(-i * exp(-b * xseq))
# 	lines(xseq, yseq)
# 	abline(v = log(1 / i) / -b)
# 	Sys.sleep(0.25)
# }
# range(coefs[[3]][3, ])
# xseq <- seq(0, 1000, l = 100)
# a <- 1
# b <- seq(0.003, 0.05, l = 20)
# plot(NULL, xlim = range(xseq), ylim = c(0, 1), xlab = "", ylab = "")
# for (i in b) {
# 	yseq <- exp(-a * exp(-i * xseq))
# 	lines(xseq, yseq)
# 	abline(v = log(1 / a) / -i)
# 	Sys.sleep(0.25)
# }

# functions to calculate values for the respective model functions
calcB <- function(k, a, b, x) {
	k * (1 - exp(-b * x))^a
}
calcN <- function(k, a, b, x) {
	k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
}
calcG <- function(k, a, b, x) {
	k * exp(-a * exp(-b * x))
}

# extract response variables
funNames <- c("Betalanffy", "normal", "Gompertz")
funCoefs <- data.table(group = groupsMeta$name)
vars <- c("estFutureDesc", "tenPercDesc", "maxDescPace", "descVar")
# estimated future descriptions relative to current descriptions
for (i in seq_along(funNames)) {
	funCoefs[, new := numeric()]
	if (i < 2) {
		# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
		funCoefs[, new := coefs[[i]][1, ]]
	} else if (i < 3) {
		# normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
		funCoefs[, new := coefs[[i]][1, ] * 1 / ((2 * pi)^0.5 * coefs[[i]][2, ])]
	} else {
		# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
		funCoefs[, new := coefs[[i]][1, ]]
	}
	colnames(funCoefs)[ncol(funCoefs)] <- paste0(vars[1], "_", funNames[i])
	funCoefs[, new := numeric()]
	if (i < 2) {
		# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
		funCoefs[, new := log(1 - (0.1 / coefs[[i]][1, ])^(1 / coefs[[i]][2, ])) / -coefs[[i]][3, ]]
	} else if (i < 3) {
		# normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
		yVals <- 2 * coefs[[i]][3, ] - (coefs[[i]][2, ] * (-2 * log(((2 * pi)^0.5 * coefs[[i]][2, ]) /
			(10 * coefs[[i]][1, ])))^0.5 + coefs[[i]][3, ])
		st <- standardize(groupsData$year - min(groupsData$year))
		funCoefs[, new := yVals * attr(st, "scaled:scale") + attr(st, "scaled:center")]
	} else {
		# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
		funCoefs[, new := log((log(0.1 / coefs[[i]][1, ]) / -coefs[[i]][2, ])) / -coefs[[i]][3, ]]
	}
	colnames(funCoefs)[ncol(funCoefs)] <- paste0(vars[2], "_", funNames[i])
	funCoefs[, new := numeric()]
	if (i < 2) {
		# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
		yVals <- log(1 / coefs[[i]][2, ]) / -coefs[[i]][3, ]
		yVals[yVals > 2017 - 1753] <- 2017 - 1753
		yVals <- coefs[[i]][1, ] * (1 - exp(-coefs[[i]][3, ] * yVals))^coefs[[i]][2, ]
		funCoefs[, new := yVals * colSums(groupsData)[-1]]
	} else if (i < 3) {
		# normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
		yVals <- coefs[[i]][3, ] - coefs[[i]][2, ]
		st <- standardize(groupsData$year - min(groupsData$year))
		yVals[yVals > max(st)] <- max(st)
		yVals <- coefs[[i]][1, ] * 1 / ((2 * pi)^0.5 * coefs[[i]][2, ]) *
			exp(-0.5 * ((yVals - coefs[[i]][3, ]) / coefs[[i]][2, ])^2)
		funCoefs[, new := yVals * colSums(groupsData)[-1]]
	} else {
		# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
		yVals <- log(1 / coefs[[i]][2, ]) / -coefs[[i]][3, ]
		yVals[yVals > 2017 - 1753] <- 2017 - 1753
		yVals <- coefs[[i]][1, ] * exp(-coefs[[i]][2, ] * exp(-coefs[[i]][3, ] * yVals))
		funCoefs[, new := yVals * colSums(groupsData)[-1]]
	}
	colnames(funCoefs)[ncol(funCoefs)] <- paste0(vars[3], "_", funNames[i])
	funCoefs[, new := numeric()]
	if (i < 2) {
		# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
		# estimates
		ests <- sapply(seq_len(ncol(coefs[[i]])), function(x) {
			calcB(coefs[[i]][1, x], coefs[[i]][2, x], coefs[[i]][3, x], groupsData$year - min(groupsData$year))
		})
		# measurements
		meas <- groupsData[, -"year"]
		meas <- apply(meas, 2, function(x) cumsum(x))
		meas <- apply(meas, 2, function(x) x / max(x))
		funCoefs[, new := colSums((ests - meas)^2) / nrow(meas)]
	} else if (i < 3) {
		# estimates
		ests <- sapply(seq_len(ncol(coefs[[i]])), function(x) {
			calcN(coefs[[i]][1, x], coefs[[i]][2, x], coefs[[i]][3, x], standardize(groupsData$year - min(groupsData$year)))
		})
		# measurements
		meas <- groupsData[, -"year"]
		meas <- apply(meas, 2, function(x) cumsum(x))
		meas <- apply(meas, 2, function(x) x / max(x))
		funCoefs[, new := colSums((ests - meas)^2) / nrow(meas)]
	} else {
		# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
		# estimates
		ests <- sapply(seq_len(ncol(coefs[[i]])), function(x) {
			calcG(coefs[[i]][1, x], coefs[[i]][2, x], coefs[[i]][3, x], groupsData$year - min(groupsData$year))
		})
		# measurements
		meas <- groupsData[, -"year"]
		meas <- apply(meas, 2, function(x) cumsum(x))
		meas <- apply(meas, 2, function(x) x / max(x))
		funCoefs[, new := colSums((ests - meas)^2) / nrow(meas)]
	}
	colnames(funCoefs)[ncol(funCoefs)] <- paste0(vars[4], "_", funNames[i])
}
selCols <- colnames(funCoefs)[grepl("normal", colnames(funCoefs))]
groupsResponses <- cbind(groupsMeta, funCoefs[, ..selCols])
colnames(groupsResponses) <- sub("_normal", "", colnames(groupsResponses))
rmCols <- c("file", "color", "finalAuthorData")
groupsResponses[, (rmCols) := NULL]
fwrite(groupsResponses, file = "groupsResponses.txt")

# show estimates of response variables extracted from curves

# pdf("response variables.pdf", width = 8.3, height = 11.7)
# estimated future descriptions relative to current descriptions
par(mfrow = c(4, 3))
par(oma = c(0, 1, 4, 1))
par(mar = c(5.1, 2.1, 2.1, 1.1))
i <- 1
for (i in seq_along(funNames)) {
	plot(colSums(groupsData)[-1], funCoefs[[which(colnames(funCoefs) == paste0(vars[1], "_", funNames[i]))]],
		col = "white", xlab = "current descriptions",
		ylab = "", xaxt = "n", main = ""
	)
	axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
	abline(v = seq(0, 500000, 50000), lty = 3)
	if (max(yVals) < 10) abline(h = seq(0, max(yVals + 1)), lty = 3) else abline(h = seq(0, max(yVals + 1), 5), lty = 3)
	text(colSums(groupsData)[-1], funCoefs[[which(colnames(funCoefs) == paste0(vars[1], "_", funNames[i]))]],
		label = 1:47, col = groupsMeta$color, cex = 0.9
	)
	mtext(funNames[i], 3, 3, font = 2)
	if (i == 2) mtext("estimated future descriptions / current descriptions", 3, 1, font = 2, cex = 0.9)
}
# differences in start of descriptions
# identify moment in time when 10% of descriptions where made
for (i in seq_along(funNames)) {
	plot(colSums(groupsData)[-1], funCoefs[[which(colnames(funCoefs) == paste0(vars[2], "_", funNames[i]))]],
		col = "white", xlab = "current descriptions",
		ylab = "", xaxt = "n", main = ""
	)
	axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
	abline(v = seq(0, 500000, 50000), lty = 3)
	abline(h = seq(0, 300, 20), lty = 3)
	text(colSums(groupsData)[-1], funCoefs[[which(colnames(funCoefs) == paste0(vars[2], "_", funNames[i]))]],
		label = 1:47, col = groupsMeta$color, cex = 0.9
	)
	if (i == 2) mtext("years until 10% of current descriptions made", 3, 1, font = 2, cex = 0.9)
}
# description pace
# identify maximum description velocity and multiply by actual description number to get
# comparable values
for (i in seq_along(funNames)) {
	plot(colSums(groupsData)[-1], funCoefs[[which(colnames(funCoefs) == paste0(vars[3], "_", funNames[i]))]],
		col = "white", xlab = "current descriptions",
		ylab = "", xaxt = "n", main = ""
	)
	axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
	abline(v = seq(0, 500000, 50000), lty = 3)
	if (max(yVals) > 150000) abline(h = seq(0, 500000, 50000), lty = 3) else abline(h = seq(0, 500000, 20000), lty = 3)
	text(colSums(groupsData)[-1], funCoefs[[which(colnames(funCoefs) == paste0(vars[3], "_", funNames[i]))]],
		label = 1:47, col = groupsMeta$color, cex = 0.9
	)
	if (i == 2) mtext("maximum descriptions pace", 3, 1, font = 2, cex = 0.9)
}
# variability
for (i in seq_along(funNames)) {
	plot(colSums(groupsData)[-1], funCoefs[[which(colnames(funCoefs) == paste0(vars[4], "_", funNames[i]))]],
		col = "white", xlab = "current descriptions",
		ylab = "", xaxt = "n", main = ""
	)
	axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
	abline(v = seq(0, 500000, 50000), lty = 3)
	abline(h = seq(0, 2, 0.1), lty = 3)
	text(colSums(groupsData)[-1], funCoefs[[which(colnames(funCoefs) == paste0(vars[4], "_", funNames[i]))]],
		label = 1:47, col = groupsMeta$color, cex = 0.9
	)
	if (i == 2) mtext("sum of squares (estimates - measurements)", 3, 1, font = 2, cex = 0.9)
}
# dev.off()

# nicer display of the above
# pdf("testDisplay.pdf")
# b1 <- barplot(coefs[[1]][1, ], col = groupsMeta$color, border = NA, horiz = TRUE, space = 0.5)
# i <- 1
# for (i in seq_len(ncol(coefs[[1]]))) {
# 	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
# 	scale <- max(coefs[[1]][1, ]) / max(b1)
# 	colRGB <- col2rgb(groupsMeta$color[i]) / 255
# 	for (j in 1:3) img[, , j] <- colRGB[[j]]
# 	rasterImage(img, coefs[[1]][1, i] + 0.5, b1[i] - 1, scale * 2 + coefs[[1]][1, i] + 0.5, b1[i] + 1, xpd = TRUE)
# }
# dev.off()

# plot numbers and groups
par(parBackup)
# pdf("numbers and groups.pdf")
plot(NULL, xlim = c(0, 5), ylim = c(0, 47), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
for (i in groupsVec) {
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	colRGB <- col2rgb(groupsMeta$color[i]) / 255
	for (j in 1:3) img[, , j] <- colRGB[[j]]
	rasterImage(img, 1.452, i - 1, 1.55, i)
	text(2, i - 0.5, groupsMeta$name[i], adj = 0, col = groupsMeta$color[i], cex = 0.7)
	text(1, i - 0.5, i, col = groupsMeta$color[i], cex = 0.7)
}
# dev.off()

# 5 Source distribution data from GBIF#############################################################
nextScript <- NULL

# load in libraries
library(data.table)
library(RJSONIO)
library(rgbif)

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# get data
groupsData <- fread("groupsData.txt")

# groupsGBIFTaxonKeys <- rep(NA,ncol(groupsData)-1)
# for (i in seq_along(groupsGBIFTaxonKeys)){
# 	res <- fromJSON(paste0("https://api.gbif.org/v1/species/match?name=",
# colnames(groupsData)[i + 1],"&strict=TRUE&verbose=TRUE"))
# 	if ("usageKey" %in% names(res) && res$matchType=="EXACT"){
# 		groupsGBIFTaxonKeys[i] <- res$usageKey
# 	} else {
# 		if ("alternatives" %in% names(res)){
# 			res <- res$alternatives
# 			for (j in seq_along(res)) print(unlist(res[[j]]))
# 			myNum <- as.numeric(readline(prompt = "Which one?"))
# 			groupsGBIFTaxonKeys[i] <- res[[myNum]]$usageKey
# 		}
# 	}
# }
## four problematic keys: Plecoptera gets Insecta key, Ectoprocta not found, Myriapoda not found, Crustacea is insect taxon
# groupsGBIFTaxonKeys[which(colnames(groupsData) == "Crustacea") - 1] <- NA
# groupsGBIFTaxonKeys[which(colnames(groupsData) == "Ectoprocta") - 1] <- 53 # Bryozoa == Ectoprocta
# groupsGBIFTaxonKeys[which(colnames(groupsData) == "Plecoptera") - 1] <- 787 # alternative
## Myriapoda is not acknowledged as a group
myriapods <- c("Chilopoda", "Symphyla", "Pauropoda", "Diplopoda")
## Crustacea not acknowledged as a group
crustaceans <- c("Malacostraca", "Copepoda", "Ostracoda", "Branchiopoda", "Maxillopoda", "Cephalocarida", "Remipedia")
# newNames <- c(myriapods,crustaceans)
# newKeys <- rep(NA,length(newNames))
# for (i in seq_along(newNames)){
# 	newKeys[i] <- fromJSON(paste0("https://api.gbif.org/v1/species/match?name=",newNames[i],"&strict=TRUE&verbose=TRUE"))$usageKey
# }
# groupsGBIFTaxonKeys <- c(groupsGBIFTaxonKeys,newKeys)
# names(groupsGBIFTaxonKeys) <- c(colnames(groupsData)[-1],myriapods,crustaceans)
# save(groupsGBIFTaxonKeys, file ="groupsGBIFTaxonKeys.RData")
load("groupsGBIFTaxonKeys.RData")

continents <- c("africa", "antarctica", "asia", "europe", "north_america", "oceania", "south_america")

## GBIF occurrences - NOT SPECIES numbers
# resOccs <- matrix(NA,length(GBIFTaxonKeys) + length(MyriapodTaxonKeys),length(continents))
#
# for (i in seq_along(GBIFTaxonKeys)){
# 	if (!is.na(GBIFTaxonKeys)[i]){
# 		for (j in seq_along(continents)){
# 			resOccs[i,j] <- occ_count(taxonKey=GBIFTaxonKeys[i],continent=continents[j])
# 		}
# 	}
# }
# for (i in seq_along(MyriapodTaxonKeys)){
# 	if (!is.na(MyriapodTaxonKeys)[i]){
# 		for (j in seq_along(continents)){
# 			resOccs[i+length(GBIFTaxonKeys),j] <- occ_count(taxonKey=MyriapodTaxonKeys[i],continent=continents[j])
# 		}
# 	}
# }
# resOccs <- data.table(resOccs)
# resOccs[,group := c(colnames(groupsData)[-1],c("Chilopoda","Symphyla", "Pauropoda","Diplopoda"))]
# colnames(resOccs) <- c(continents,"group")
# resOccs[, oldWorld := north_america + europe]
# resOccs[, newWorld := africa + antarctica + asia + oceania + south_america]
# resOccs[, ratio := oldWorld/newWorld]

## THERE IS A DIFFERENCE IN THE NUMBER OF SPECIES RETURNED THROUGH OCCURRENCES AND THROUGH THE SPECIES API
#
## get usage key for genus
# eKey <- fromJSON(paste0("https://api.gbif.org/v1/species/match?name=Entandrophragma&strict=TRUE&verbose=TRUE"))$usageKey
## get children
# eChilds <- fromJSON(paste0("https://api.gbif.org/v1/species/",eKey,"/childrenAll"))
# sapply(eChilds,function(x) x$key)
# sapply(eChilds,function(x) x$name)
# sapply(eChilds,function(x) x$rank)
## get occurrences
# eOccs <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=",eKey,"&limit=0&facet=speciesKey&facetLimit=10000"))
# sapply(eOccs$facets[[1]]$counts,function(x) x$name)
# sapply(eOccs$facets[[1]]$counts,function(x) x$count)
# sapply(eChilds,function(x) x$name)[!(sapply(eChilds,function(x) x$key)) %in% sapply(eOccs$facets[[1]]$counts,function(x) x$name)]
#
## The occurrence API seems to correctly merge the two occurrences of Entandrophragma delevoyi to one record.
#
## get usage key for genus
# eKey <- fromJSON(paste0("https://api.gbif.org/v1/species/match?name=Bellis&strict=TRUE&verbose=TRUE"))$usageKey
## get children
# eChilds <- fromJSON(paste0("https://api.gbif.org/v1/species/",eKey,"/childrenAll"))
# sapply(eChilds,function(x) x$key)
# sapply(eChilds,function(x) x$name)
# sapply(eChilds,function(x) x$rank)
## get occurrences
# eOccs <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=",eKey,"&limit=0&facet=speciesKey&facetLimit=10000"))
# sapply(eOccs$facets[[1]]$counts,function(x) x$name)
# sapply(eOccs$facets[[1]]$counts,function(x) x$count)
# sapply(eChilds,function(x) x$name)[!(sapply(eChilds,function(x) x$key)) %in% sapply(eOccs$facets[[1]]$counts,function(x) x$name)]
#
## The occurrence API does not have records for every species. This is likely the main reason for different species numbers.

# Retrieve species numbers of all taxa per continent

# GBIF occurrences with facetted species numbers - yes, species keys correspond to accepted species (in theory)
resSpecs <- list()
for (i in seq_along(GBIFTaxonKeys)) {
	if (i > 0) {
		print(names(GBIFTaxonKeys)[i])
		resSpecs[[i]] <- list()
		if (!is.na(GBIFTaxonKeys)[i]) {
			for (j in seq_along(continents)) {
				resSpecs[[i]][[j]] <- NA
				offset <- 0
				repeat {
					res <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=", GBIFTaxonKeys[i], "&continent=", continents[j], "&limit=0&facet=speciesKey&facetLimit=500000&facetOffset=", offset))
					resSpecs[[i]][[j]] <- c(resSpecs[[i]][[j]], as.numeric(sapply(res$facets[[1]]$counts, function(x) x$name)))
					if (length(res$facets[[1]]$counts) < 500000) break else offset <- offset + 500000
				}
			}
		}
	}
}
# calculate numbers for Crustaceae and Myriapoda (warning: has several NA values in it)
for (i in seq_along(continents)) {
	resSpecs[[which(names(GBIFTaxonKeys) == "Crustacea")]][[i]] <- unlist(sapply(which(names(GBIFTaxonKeys) %in% crustaceans), function(x) {
		resSpecs[[x]][[i]]
	}))
	resSpecs[[which(names(GBIFTaxonKeys) == "Myriapoda")]][[i]] <- unlist(sapply(which(names(GBIFTaxonKeys) %in% myriapods), function(x) {
		resSpecs[[x]][[i]]
	}))
}
resSpecsTable <- data.table(GBIFTaxonKey = GBIFTaxonKeys, name = names(GBIFTaxonKeys))
# sum up europe and north america
temp <- sapply(seq_along(resSpecs), function(x) length(unique(unlist(resSpecs[[x]][4:5]))) - 1)
resSpecsTable[, europe_north_america := temp]
temp <- sapply(seq_along(resSpecs), function(x) length(unique(unlist(resSpecs[[x]][c(1:3, 6:7)]))) - 1)
resSpecsTable[, africa_antarctica_asia_oceania_south_america := temp]
temp <- sapply(seq_along(resSpecs), function(x) length(unique(unlist(resSpecs[[x]][c(1:7)]))) - 1)
resSpecsTable[, world := temp]
# calculate ratio
resSpecsTable <- resSpecsTable[!name %in% c(myriapods, crustaceans)]
resSpecsTable[, ratio := europe_north_america / world]

# check ratios across groups
setorder(resSpecsTable, "ratio")
resSpecsTable[, c("name", "ratio")]
hist(resSpecsTable$ratio)

# reorder table
resSpecsTable[, oriOrder := sapply(resSpecsTable$name, function(x) which(colnames(groupsData) == x)) - 1]
setorder(resSpecsTable, "oriOrder")

# check whether total numbers match LifeGate
# could be higher, if taxa occur in more than one continent, and lower,
# if taxa have no occurrence data
plot(resSpecsTable$world ~ colSums(groupsData[, -"year"]), xlab = "LifeGate total species numbers", ylab = "GBIF total species numbers", col = "white")
names(GBIFTaxonKeys)
abline(0, 1)
text(colSums(groupsData[, -"year"]), resSpecsTable$world, seq_len(nrow(resSpecsTable)))

# write data
fwrite(resSpecsTable, "groupsOccurrences.txt")

# 6 Source body size data from Ulrich Brose's lab##################################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large data sets
library(RJSONIO) # read json

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxonomy"))

# load help functions
source("taxonomy help functions.R")

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# read in data
groupsMeta <- fread("groupsMeta.txt")

b1 <- fread("bodysizes.txt", sep = "\t", fill = TRUE)
b2 <- fread("bodysizes_2008.txt")
# repair CP1252 encoding
# make sure later regex use works as expected
selCols <- colnames(b1)[sapply(b1, is.character)]
b1[, (selCols) := lapply(.SD, function(x) {
	x[!validUTF8(x)] <- iconv(x[!validUTF8(x)], from = "CP1252", to = "UTF-8")
	Encoding(x) <- "unknown"
	x
}), .SDcols = selCols]
selCols <- colnames(b2)[sapply(b2, is.character)]
b2[, (selCols) := lapply(.SD, function(x) {
	x[!validUTF8(x)] <- iconv(x[!validUTF8(x)], from = "CP1252", to = "UTF-8")
	Encoding(x) <- "unknown"
	x
}), .SDcols = selCols]

rmCols <- colnames(b1)[sapply(b1, function(x) all(is.na(x)))]
b1[, (rmCols) := NULL]
b1[b1 == -999] <- NA
b2[b2 == -999] <- NA
rmRows <- apply(b1, 1, function(x) all(is.na(x) | x == ""))
b1 <- b1[!rmRows]
str(b1)
str(b2)
which(b1 != b2[seq_len(nrow(b1))], arr.ind = TRUE) # many
which(b1 != b2[seq_len(nrow(b1))], arr.ind = TRUE)[1]
rbind(b1[348], b2[348])
rbind(b1[348], b2[grepl("Procladius sagittalis", `Taxonomy consumer`)][3])
# use b2 and hope it's better

apply(b2, 2, function(x) sum(!is.na(x)))
table(b2$`Common name(s) consumer`)
b2[is.na(`Mean length (m) consumer`) & (!is.na(`Minimum length (m) consumer`) | !is.na(`Maximum length (m) consumer`))]
b2[is.na(`Mean length (m) resource`) & (!is.na(`Minimum length (m) resource`) | !is.na(`Maximum length (m) resource`))]
str(b2)

# create datatable with relevant information
datn <- data.table(
	name1 = b2$`Taxonomy consumer`, name2 = b2$`Common name(s) consumer`,
	lengthMin = b2$`Minimum length (m) consumer`, lengthMean = b2$`Mean length (m) consumer`,
	lengthMax = b2$`Maximum length (m) consumer`, massMin = b2$`Minimum mass (g) consumer`,
	massMean = b2$`Mean mass (g) consumer`, massMax = b2$`Maximum mass (g) consumer`
)
dat2 <- data.table(
	name1 = b2$`Taxonomy resource`, name2 = b2$`Common name(s) resource`,
	lengthMin = b2$`Minimum length (m) resource`, lengthMean = b2$`Mean length (m) resource`,
	lengthMax = b2$`Maximum length (m) resource`, massMin = b2$`Minimum mass (g) resource`,
	massMean = b2$`Mean mass (g) resource`, massMax = b2$`Maximum mass (g) resource`
)
dat <- rbind(datn, dat2)
# create mean from min and max
dat[is.na(lengthMean) & !is.na(lengthMin), lengthMean := (lengthMin + lengthMax) / 2]
dat[is.na(massMean) & !is.na(massMin), massMean := (massMin + massMax) / 2]
# create mean from either min or max
mins <- dat[is.na(lengthMean)]$lengthMin
maxs <- dat[is.na(lengthMean)]$lengthMax
mins[is.na(mins)] <- 0
maxs[is.na(maxs)] <- 0
means <- mins + maxs
means[means == 0] <- NA
dat[is.na(lengthMean), lengthMean := means]
mins <- dat[is.na(massMean)]$massMin
maxs <- dat[is.na(massMean)]$massMax
mins[is.na(mins)] <- 0
maxs[is.na(maxs)] <- 0
means <- mins + maxs
means[means == 0] <- NA
dat[is.na(massMean), massMean := means]
# remove not needed cols
rmCols <- colnames(dat)[grepl("Min|Max", colnames(dat))]
dat[, (rmCols) := NULL]

# taxonomic classification of names using GBIF
namesDict <- unique(dat[, c("name1", "name2")])
namesDict[, order := character()]
namesDict[, class := character()]
namesDict[, phylum := character()]
# create searchName that avoids problematic characters
namesDict[, searchName := sub("\\s+sp(ec)?p?\\.?(ies)?( 1)?(\\s.*|$)", "", name1)]
namesDict[, searchName := gsub("\"|\\\\", "", searchName)]
namesDict[, searchName := sub("^Other\\s+", "", searchName)]
namesDict[, searchName := sub("\\s([[:upper:]]|\\().*", "", searchName)]
namesDict[!validUTF8(searchName), searchName := iconv(searchName[!validUTF8(searchName)], from = "CP1252", to = "UTF8")]
Encoding(namesDict$searchName) <- "unknown"
for (i in seq_len(nrow(substAcc))) {
	namesDict[, searchName := gsub(substAcc[i, 1], substAcc[i, 2], searchName)]
}
for (i in seq_len(nrow(namesDict))) {
	if (i %% 100 == 0) print(i)
	if (i > 0) {
		if (is.na(namesDict$order[i])) {
			for (j in seq_len(trials)) {
				res <- tryCatch(fromJSON(paste0(
					"https://api.gbif.org/v1/species/match?name=", gsub(" ", "%20", namesDict[i]$searchName),
					"&verbose=true"
				)), error = function(e) NA)
				if (!is.atomic(res) || j > trials - 1) break else Sys.sleep(1)
			}
			if ((res$matchType == "NONE" || res$matchType == "HIGHERRANK" && res$rank != "SPECIES") &&
				"alternatives" %in% names(res)) {
				names <- sapply(res$alternatives, function(x) if ("canonicalName" %in% names(x)) x$canonicalName else "")
				conf <- sapply(res$alternatives, function(x) if ("confidence" %in% names(x)) x$confidence else "")
				order <- sapply(res$alternatives, function(x) if ("order" %in% names(x)) x$order else "")
				class <- sapply(res$alternatives, function(x) if ("class" %in% names(x)) x$class else "")
				phylum <- sapply(res$alternatives, function(x) if ("phylum" %in% names(x)) x$phylum else "")
				resTable <- data.table(num = seq_along(names), names, conf, order, class, phylum)
				resTable <- resTable[conf == max(conf)]
				if (length(unique(resTable$order[resTable$order != ""])) == 1) {
					namesDict[i, order := resTable$order[resTable$order != ""][1]]
					namesDict[i, class := resTable$class[resTable$class != ""][1]]
					namesDict[i, phylum := resTable$phylum[resTable$phylum != ""][1]]
				} else if (length(unique(resTable$class[resTable$class != ""])) == 1) {
					namesDict[i, class := resTable$class[resTable$class != ""][1]]
					namesDict[i, phylum := resTable$phylum[resTable$phylum != ""][1]]
				} else if (length(unique(resTable$phylum[resTable$phylum != ""])) == 1) {
					namesDict[i, phylum := resTable$phylum[resTable$phylum != ""][1]]
				}
			} else {
				if ("order" %in% names(res)) namesDict[i, order := res$order]
				if ("class" %in% names(res)) namesDict[i, class := res$class]
				if ("phylum" %in% names(res)) namesDict[i, phylum := res$phylum]
			}
		}
	}
}
namesDict[is.na(order) & is.na(class) & is.na(phylum)]$searchName
# create searchName that avoids problematic characters
namesDict[, searchName := sub("(/|:|,).*", "", name2)]
namesDict[, searchName := gsub("\"|\\\\", "", searchName)]
namesDict[, searchName := sub(".*(Benthic|[[:lower:]]+plankton|Sediment|Terrestrial|Other|Large|Small)(\\s.*|$)", "", name2)]
namesDict[, searchName := gsub("-sp$", "", searchName)]
namesDict[, searchName := gsub("\\s+and\\s.*$", "", searchName)]
namesDict[!validUTF8(searchName), searchName := iconv(searchName[!validUTF8(searchName)], from = "CP1252", to = "UTF8")]
Encoding(namesDict$searchName) <- "unknown"
for (i in seq_len(nrow(substAcc))) {
	namesDict[, searchName := gsub(substAcc[i, 1], substAcc[i, 2], searchName)]
}
for (i in seq_len(nrow(namesDict))) {
	if (i %% 100 == 0) print(i)
	if (is.na(namesDict$order[i]) && is.na(namesDict$class[i]) && is.na(namesDict$phylum[i])) {
		if (!is.na(namesDict[i]$searchName) && namesDict[i]$searchName != "") {
			for (j in seq_len(trials)) {
				res <- tryCatch(fromJSON(paste0(
					"https://api.gbif.org/v1/species/search?q=", gsub(" ", "%20", namesDict[i]$searchName),
					"&qField=VERNACULAR"
				)), error = function(e) NA)
				if (!is.atomic(res) || j > trials - 1) break else Sys.sleep(1)
			}
			if (length(res$results) > 0) {
				groupNames <- list()
				groupNames[[1]] <- sapply(res$results, function(x) if ("order" %in% names(x)) tolower(x$order) else "")
				groupNames[[2]] <- sapply(res$results, function(x) if ("class" %in% names(x)) tolower(x$class) else "")
				groupNames[[3]] <- sapply(res$results, function(x) if ("phylum" %in% names(x)) tolower(x$phylum) else "")
				for (j in seq_along(groupNames)) {
					groupNames[[j]] <- groupNames[[j]][groupNames[[j]] != ""]
					groupNames[[j]] <- table(groupNames[[j]])
					if (length(groupNames[[j]]) > 0) groupNames[[j]] <- groupNames[[j]][groupNames[[j]] == max(groupNames[[j]])]
					if (length(groupNames[[j]]) != 1) {
						groupNames[[j]] <- NA
						names(groupNames[[j]]) <- ""
					}
				}
				namesDict[i, order := names(groupNames[[1]])]
				namesDict[i, class := names(groupNames[[2]])]
				namesDict[i, phylum := names(groupNames[[3]])]
			} else {
				for (j in seq_len(trials)) {
					res <- tryCatch(fromJSON(paste0(
						"https://api.gbif.org/v1/species/match?name=", gsub(" ", "%20", namesDict[i]$searchName),
						"&verbose=true"
					)), error = function(e) NA)
					if (!is.atomic(res) || j > trials - 1) break else Sys.sleep(1)
				}
				if ((res$matchType == "NONE" || res$matchType == "HIGHERRANK" && res$rank != "SPECIES") &&
					"alternatives" %in% names(res)) {
					names <- sapply(res$alternatives, function(x) if ("canonicalName" %in% names(x)) x$canonicalName else "")
					conf <- sapply(res$alternatives, function(x) if ("confidence" %in% names(x)) x$confidence else "")
					order <- sapply(res$alternatives, function(x) if ("order" %in% names(x)) x$order else "")
					class <- sapply(res$alternatives, function(x) if ("class" %in% names(x)) x$class else "")
					phylum <- sapply(res$alternatives, function(x) if ("phylum" %in% names(x)) x$phylum else "")
					resTable <- data.table(num = seq_along(names), names, conf, order, class, phylum)
					resTable <- resTable[conf == max(conf)]
					if (length(unique(resTable$order[resTable$order != ""])) == 1) {
						namesDict[i, order := resTable$order[resTable$order != ""][1]]
						namesDict[i, class := resTable$class[resTable$class != ""][1]]
						namesDict[i, phylum := resTable$phylum[resTable$phylum != ""][1]]
					} else if (length(unique(resTable$class[resTable$class != ""])) == 1) {
						namesDict[i, class := resTable$class[resTable$class != ""][1]]
						namesDict[i, phylum := resTable$phylum[resTable$phylum != ""][1]]
					} else if (length(unique(resTable$phylum[resTable$phylum != ""])) == 1) {
						namesDict[i, phylum := resTable$phylum[resTable$phylum != ""][1]]
					}
				} else {
					if ("order" %in% names(res)) namesDict[i, order := res$order]
					if ("class" %in% names(res)) namesDict[i, class := res$class]
					if ("phylum" %in% names(res)) namesDict[i, phylum := res$phylum]
				}
			}
		}
	}
}
namesDict[namesDict == ""] <- NA
namesDict[is.na(order) & is.na(class) & is.na(phylum)]$searchName
namesDict[is.na(phylum) & (!is.na(order) | !is.na(class))]
namesDict[is.na(phylum) & is.na(class) & !is.na(order)]
namesDict[is.na(phylum), order := NA]
namesDict[is.na(phylum), class := NA]
selCols <- c("order", "class", "phylum")
namesDict[, (selCols) := lapply(.SD, function(x) {
	x[!is.na(x)] <- paste0(toupper(substr(x[!is.na(x)], 1, 1)), sub("^.", "", x[!is.na(x)]))
}), .SDcols = (selCols)]
shortNames <- groupsMeta$name
myriapods <- c("Chilopoda", "Symphyla", "Pauropoda", "Diplopoda")
crustaceans <- c("Malacostraca", "Copepoda", "Ostracoda", "Branchiopoda", "Maxillopoda", "Cephalocarida", "Remipedia")
namesDict[, nameShort := character()]
namesDict[order %in% shortNames, nameShort := order]
namesDict[class %in% myriapods, nameShort := "Myriapoda"]
namesDict[class %in% crustaceans, nameShort := "Crustacea"]
namesDict[class %in% shortNames, nameShort := class]
namesDict[phylum == "Bryozoa", nameShort := "Ectoprocta"]
namesDict[phylum %in% shortNames, nameShort := phylum]
namesDict[name1 == "Anisoptera", nameShort := "Odonata"]
namesDict[grepl("Oligochaet", name1), nameShort := "Annelida"]
namesDict[grepl("Distephanus", name1), nameShort := "Ochrophyta"]
namesDict[grepl("Aran", name1), nameShort := "Arachnida"]
namesDict[grepl("Homoptera", name1) | grepl("Homoptera", name2), nameShort := "Hemiptera"]
namesDict[grepl("Chaetognatha", name1), nameShort := NA]
table(namesDict[(!is.na(order) | !is.na(class) | !is.na(phylum)) & is.na(nameShort)]$order)
table(namesDict[(!is.na(order) | !is.na(class) | !is.na(phylum)) & is.na(nameShort)]$class)
table(namesDict[(!is.na(order) | !is.na(class) | !is.na(phylum)) & is.na(nameShort)]$phylum)
table(namesDict$nameShort)
setdiff(shortNames, names(table(namesDict$nameShort)))

namesDict[, key := paste(name1, name2)]
setkey(namesDict, key)
dat[, name := namesDict[paste(dat$name1, dat$name2)]$nameShort]
dat <- dat[!is.na(nameShort)]

datSummary <- data.table(
	nameShort = shortNames, minLength = numeric(), meanLength = numeric(), maxLength = numeric(), nLength = numeric(),
	minMass = numeric(), meanMass = numeric(), maxMass = numeric(), nMass = numeric()
)
i <- 4
for (i in seq_len(nrow(datSummary))) {
	print(datSummary$nameShort[i])
	temp <- dat[nameShort == datSummary$nameShort[i]]
	print(unique(temp, by = c("name1", "name2"))[, c("name1", "name2")])
	print("--------------------------------")
	if (nrow(temp) > 0) {
		datSummary[i, nLength := sum(!is.na(temp$lengthMean))]
		if (datSummary[i]$nLength > 0) {
			datSummary[i, minLength := min(temp$lengthMean, na.rm = TRUE)]
			datSummary[i, meanLength := mean(temp$lengthMean, na.rm = TRUE)]
			datSummary[i, maxLength := max(temp$lengthMean, na.rm = TRUE)]
		}
		datSummary[i, nMass := sum(!is.na(temp$massMean))]
		if (datSummary[i]$nMass > 0) {
			datSummary[i, minMass := min(temp$massMean, na.rm = TRUE)]
			datSummary[i, meanMass := mean(temp$massMean, na.rm = TRUE)]
			datSummary[i, maxMass := max(temp$massMean, na.rm = TRUE)]
		}
	}
}
fwrite(datSummary, file = "groupsLengthSize.txt")

# 7 Source data from BiL explorer##################################################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets
library(RSelenium) # for web scraping
library(doSNOW)

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxonomy"))

# load help functions
source("taxonomy help functions.R")

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# get groups
groupsMeta <- fread("groupsMeta.txt")

# get target names
names <- fread("unique-taxons.csv")

# function to retrieve data
webScraper <- function(i) {
	tryCatch(
		{
			# open the website
			remDr$navigate("http://ch01.informatik.uni-leipzig.de:5100/bil-explorer/")
			Sys.sleep(5)
			repeat {
				inp <- tryCatch(remDr$findElements(using = "id", value = "react-select-2-input"), error = function(e) NA)
				if (is.na(inp)) Sys.sleep(1) else break
			}
			Sys.sleep(1)
			inp[[1]]$sendKeysToElement(list(names$occId[i], key = "enter")) # press button
			Sys.sleep(1)
			repeat {
				waitSymbol <- remDr$findElements(using = "class", value = "animate-spin")[[1]]
				if (!grepl("aria-hidden=.true.", waitSymbol$getElementAttribute("outerHTML")[[1]])) {
					Sys.sleep(1)
				} else {
					break
				}
			}
			Sys.sleep(1)
			occLine <- remDr$findElements(using = "id", value = "viz-container")[[1]]
			winSize <- occLine$getElementSize()
			html_content <- occLine$getElementAttribute("outerHTML")[[1]]
			values <- regmatches(html_content, regexpr("<path.*></path>", html_content))
			values <- sub("\".*", "", sub(".*d=\"M", "", values))
			values <- strsplit(values, "L")[[1]]
			# catch special cases with 0 or 1 data point
			if (length(values) == 1) {
				if (grepl("Z", values)) values <- sub("Z", "", values) else values <- NA
			}
			if (!all(is.na(values))) {
				values <- data.table(ori = values, x = as.numeric(sub(",.*", "", values)), y = as.numeric(sub(".*,", "", values)))
				# plot(-values$y~values$x,type="l", ylab = "relative importance", xlab = "relative time", main = names$occId[i])

				values[, xRel := x - winSize$width / 2]
				values[, yRel := y - winSize$height / 2]
				values[, year := numeric()]
				values[, count := numeric()]
				for (j in seq_len(nrow(values))) {
					test <- tryCatch(remDr$mouseMoveToLocation(x = values$xRel[j], y = values$yRel[j], webElement = occLine),
						error = function(e) {
							print(e)
							NA
						}
					)
					if (is.null(test)) {
						remDr$click()
						repeat {
							waitSymbol <- remDr$findElements(using = "class", value = "animate-spin")
							if (length(waitSymbol) < 1) break else Sys.sleep(0.5)
						}
						curYear <- remDr$findElements(using = "css", value = "h5")[[1]]
						curYear <- as.numeric(sub(":.*", "", sub(".*year ", "", curYear$getElementAttribute("outerHTML")[[1]])))
						curCount <- remDr$findElements(using = "css", value = "th")[[1]]
						values[j, year := curYear]
						values[j, count := as.numeric(gsub("\\D", "", curCount$getElementAttribute("outerHTML")[[1]]))]
						print(values[j])
						backButton <- remDr$findElements(using = "css", value = "svg")
						backButton[[6]]$clickElement()
					}
				}
			} else {
				# empty data.table due to missing data
				values <- data.table(ori = character(), x = numeric(), y = numeric(), xRel = numeric(), yRel = numeric(), year = numeric(), count = numeric())
			}
			save(values, file = paste0("temp/", i, ".RData"))
			return(values)
		}, # error somewhere within the function
		error = function(e) NA
	)
}

# create temporary directory to store results in
if (!("temp" %in% list.files())) dir.create("temp")

# close all open ports
try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)

# rD <- rsDriver(browser = "firefox", verbose = FALSE, chromever = NULL)
# remDr <- rD[["client"]]

resAll <- list()
for (i in seq_len(nrow(names))) resAll[[i]] <- NA

# load gathered data
if ("resAllTemp.RData" %in% list.files()) load("resAllTemp.RData")
table(!sapply(resAll, is.atomic)) # number of done searches

# gather information from last runs to store in one variable
resFiles <- sort(as.numeric(sub(".RData", "", list.files("temp"))))
if (length(resFiles) > 0) {
	for (i in seq_along(resFiles)) {
		load(paste0("temp/", resFiles[i], ".RData"))
		if (!is.atomic(resAll[[resFiles[i]]])) {
			if (!all(resAll[[resFiles[i]]] == values)) {
				print(paste("Discrepancy detected between stored data and new data for", resFiles[i]))
				break
			}
		} else {
			resAll[[resFiles[i]]] <- values
		}
	}
	save(resAll, file = "resAllTemp.RData")
	invisible(file.remove(paste0("temp/", list.files("temp"))))
}

# cluster preparation
clustNumAll <- 8 # number of cores to use
cl <- parallel::makeCluster(clustNumAll)
registerDoSNOW(cl)

# clearly assign jobs to workers so that they make best use of the firefox instances opened
nJobs <- rep(sum(sapply(resAll, is.atomic)) %/% clustNumAll, clustNumAll)
nJobs[1] <- nJobs[1] + (sum(!sapply(resAll, is.atomic))) %% clustNumAll

resAll <- foreach(ii = seq_len(clustNumAll), .combine = c, .packages = c("data.table", "RSelenium")) %dopar% {
	rD <- rsDriver(browser = "firefox", verbose = FALSE, port = as.integer(4567L + ii), chromever = NULL)
	remDr <- rD[["client"]]
	for (i in which(sapply(resAll, is.atomic))[1:nJobs[ii] + sum(nJobs[0:(ii - 1)])]) {
		temp <- webScraper(i)
	}
	ii
}
stopCluster(cl)

# close all open ports
try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)

# create results table
resTable <- data.table(year = 1700:2000)
for (i in seq_len(nrow(names))) {
	resTable[, new := numeric()]
	colnames(resTable)[ncol(resTable)] <- names$occId[i]
}
# fill with data
for (i in seq_along(resAll)) {
	nameCol <- names$occId[i]
	if (nrow(resAll[[i]]) > 0) {
		resTable[sapply(resAll[[i]]$year, function(x) which(resTable$year == x)), (nameCol) := resAll[[i]]$count]
	}
}
# remove empty years
resRange <- range(which(rowSums(resTable[, -"year"], na.rm = TRUE) > 0))
resTable <- resTable[resRange[1]:resRange[2]]
range(resTable$year)

fwrite(resTable, file = "groupsLitOccurrences.txt")

trials <- 3
# do taxonomic classification using term variable
# manual selection is necessary in some cases due to ambiguity of scientific names
# provided
resVars <- c("canonicalName", "authorship", "family", "order", "class", "kingdom", "phylum", "rank")
resTax <- data.table(term = names$term)
resTax[, (resVars) := character()]

i <- 15
for (i in which(is.na(resTax$canonicalName))) {
	if (i %% 50 == 0) print(i)
	for (j in seq_len(trials)) {
		res <- tryCatch(fromJSON(paste0(
			"https://api.gbif.org/v1/species/match?name=", gsub(" ", "%20", names[i]$term),
			"&strict=false&verbose=true"
		)), error = function(e) NA)
		if (!is.atomic(res) || j > trials - 1) break else Sys.sleep(1)
	}
	if (!all(c("scientificName", "canonicalName") %in% names(res)) && "alternatives" %in% names(res)) {
		res <- res$alternatives
		resTable <- data.table(confidence = sapply(res, function(x) x$confidence))
		for (j in seq_along(resVars)) {
			resTable[, new := sapply(res, function(x) if (resVars[j] %in% names(x)) x[[which(names(x) == resVars[j])]] else NA)]
			colnames(resTable)[ncol(resTable)] <- resVars[j]
		}
		print(names[i])
		print(resTable)
		myNum <- as.numeric(readline(prompt = "Enter number: "))
		res <- resTable[myNum]
	}
	if (!is.null(dim(res)) && nrow(res) > 0) {
		if (all(c("scientificName", "canonicalName") %in% names(res))) {
			res$authorship <- trimws(sub(res$canonicalName, "", res$scientificName))
		}
		for (j in seq_along(resVars)) {
			if (resVars[j] %in% names(res)) resTax[i, (resVars[j]) := res[[which(names(res) == resVars[j])]]]
		}
	}
}

# do taxonomic classification using vernacular names
resVars <- c("canonicalName", "authorship", "family", "order", "class", "kingdom", "phylum", "rank", "vernacularNames")
resVern <- data.table(occId = names$occId)
resVern[, (resVars) := character()]

# check kingdoms and ranks
kingdoms <- c() # should be "Animalia","Plantae","Fungi"
ranks <- c() # should be "GENUS","SPECIES","SUBSPECIES","VARIETY","FAMILY","SUBCLASS","CLASS","PHYLUM", "ORDER"
for (i in which(is.na(resVern$canonicalName))) {
	if (i %% 50 == 0) print(i)
	for (k in seq_len(trials)) {
		res <- tryCatch(fromJSON(gsub(" ", "%20", paste0(
			"https://api.gbif.org/v1/species/search?q=",
			names[i]$occId, "&qField=VERNACULAR&limit=20"
		))), error = function(e) NULL)
		if (!is.null(res) || k > 2) break else Sys.sleep(1)
	}
	if (length(res$results) > 0) {
		# reformat results into table
		res <- res$results
		rows <- sapply(res, function(x) length(x$vernacularNames))
		for (k in seq_along(resVars)) {
			if (resVars[k] == "canonicalName") {
				resTable <- data.table(unlist(rep(sapply(res, function(x) {
					if (resVars[k] %in% names(x)) x[names(x) == resVars[k]] else NA
				}), rows)))
			} else if (resVars[k] != "vernacularNames") {
				resTable <- cbind(resTable, unlist(rep(sapply(res, function(x) {
					if (resVars[k] %in% names(x)) x[names(x) == resVars[k]] else NA
				}), rows)))
			} else if (resVars[k] == "vernacularNames") {
				resTable <- cbind(resTable, unlist(as.vector(sapply(res, function(x) {
					temp <- unlist(x[names(x) == resVars[k]])
					temp <- temp[!grepl("language", names(temp))]
				}))))
			}
		}
		colnames(resTable) <- resVars
		# substitute non-ASCII characters
		for (k in seq_len(nrow(substAcc))) {
			resTable[, vernacularNames :=
				gsub(substAcc[k, 1], substAcc[k, 2], vernacularNames)]
		}
		# remove hyphens
		resTable[, vernacularNames := gsub("-", " ", tolower(vernacularNames))]
		resTable[, nameFound := FALSE]
		resTable[tolower(names[i]$occId) == vernacularNames, nameFound := TRUE]
		if (all(resTable$nameFound == FALSE)) {
			resTable[tolower(names[i]$occId) == gsub(" ", "", vernacularNames), nameFound := TRUE]
		}
		if (all(resTable$nameFound == FALSE)) {
			resTable[tolower(names[i]$occId) == gsub(".* ", "", vernacularNames), nameFound := TRUE]
		}
		# gather data for check whether some data values are not covered
		kingdoms <- union(kingdoms, resTable$kingdom)
		ranks <- union(ranks, resTable$rank)
		# pre-select data based on kingdom, rank, etc.
		resTable <- resTable[nameFound == TRUE]
		resTable[, species := ""]
		resTable[!is.na(canonicalName) & grepl("species|variety|infraspecific|form", tolower(rank)) & grepl("^\\S+\\s+\\S+", canonicalName), species :=
			regmatches(canonicalName, regexpr("^\\S+\\s+\\S+", canonicalName))]
		resTable[, genus := ""]
		# species includes subspecies
		resTable[grepl("species|variety|infraspecific|form", tolower(rank)), genus := sub("\\s.*", "", species)]
		resTable[tolower(rank) == "genus", genus := canonicalName]
		resTable[, c("family", "order", "class", "kingdom") :=
			list(tolower(family), tolower(order), tolower(class), tolower(kingdom))]
		# select data
		if (nrow(resTable) > 0) {
			taxLevels <- c("species", "genus", "family", "order", "class", "kingdom")
			for (k in seq_along(taxLevels)) {
				if (taxLevels[k] %in% colnames(resTable)) {
					datTable <- table(resTable[[which(colnames(resTable) == taxLevels[k])]])
					datTable <- datTable / sum(datTable)
					datTable <- datTable[names(datTable) != ""]
					if (length(datTable) > 0 && max(datTable) >= 0.75) {
						resTable <- resTable[resTable[[which(colnames(resTable) == taxLevels[k])]] == names(datTable[datTable == max(datTable)])][1]
						break
					} else {
						resTable[, (taxLevels[k]) := NA]
					}
				}
			}
			if (nrow(resTable) < 2) {
				for (j in seq_along(resVars)) {
					if (resVars[j] %in% colnames(resTable)) {
						resVern[i, (resVars[j]) := resTable[[which(colnames(resTable) == resVars[j])]]]
					}
				}
			}
		}
	}
}
kingdoms
ranks

# combine data and search results
colnames(resTax) <- paste0("resTax_", colnames(resTax))
colnames(resVern) <- paste0("resVern_", colnames(resVern))
names <- cbind(names, resTax)
names <- cbind(names, resVern)
colnames(resTax <- sub("resTax_", "", colnames(resTax)))
colnames(resVern <- sub("resVern_", "", colnames(resVern)))

# check which group names can be found in the rows of the data
res <- matrix(NA, nrow = nrow(names), ncol = nrow(groupsMeta))
for (i in seq_len(nrow(groupsMeta))) {
	for (j in seq_len(nrow(names))) {
		res[j, i] <- tolower(groupsMeta$name[i]) %in% tolower(names[j])
	}
}
res <- data.table(res)
colnames(res) <- groupsMeta$name
res[, occId := names$occId]
checkSums <- rowSums(res[, -"occId"])
names[
	checkSums > 1,
	c(
		"occId", "term", "resTax_order", "resTax_class", "resTax_phylum",
		"resTax_kingdom", "resVern_order", "resVern_class", "resVern_phylum", "resVern_kingdom"
	)
]
names[checkSums < 1]

# create special column for ambiguous cases that contains the group the name
# is classified into
names[, group := character()]
probs <- which(checkSums != 1)
print(paste0(length(probs), " cases to check."))
for (i in seq_along(probs)) {
	if (is.na(names[probs[i]]$group)) {
		dCols <- colnames(names)[!is.na(names[probs[i]])]
		print(names[probs[i], ..dCols])
		if (checkSums[probs[i]] < 1) {
			print(groupsMeta$name)
			browseURL(paste0("https://en.wikipedia.org/wiki/", names$occId[probs[i]]), browser = NULL)
			myNum <- as.numeric(readline(prompt = "Select number: "))
			if (myNum %in% seq_along(groups$nameShort)) {
				names[probs[i], group := groupsMeta$name[myNum]]
			} else if (myNum == 0) {
				names[probs[i], group := ""]
			}
		} else {
			dCols <- colnames(res[, -"occId"])[res[probs[i], -"occId"] != FALSE]
			print(names(res[probs[i], ..dCols]))
			browseURL(paste0("https://en.wikipedia.org/wiki/", names$occId[probs[i]]), browser = NULL)
			myNum <- as.numeric(readline(prompt = "Select number: "))
			if (myNum %in% seq_along(names(res[probs[i], ..dCols]))) {
				names[probs[i], group := names(res[probs[i], ..dCols])[myNum]]
			} else if (myNum == 0) {
				names[probs[i], group := ""]
			}
		}
	}
}
names[is.na(group) & checkSums != 1]

# add non-ambiguous data
resGroup <- sapply(seq_len(nrow(res)), function(x) {
	if (rowSums(res[x, -"occId"]) == 1) colnames(res[, -"occId"])[res[x, -"occId"] == 1] else NA
})
names[!is.na(resGroup), group := resGroup[!is.na(resGroup)]]

# check unclassified records
names[group == ""]

fwrite(names, file = "groupsLitTaxonomy.txt")

# 8 Source measure of public interest based on Wikipedia article lengths and Google hits###########
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets
library(RSelenium) # web scraping

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# read in data

# groups metadata and response variables
groupsMeta <- fread("groupsMeta.txt")

# close all open ports
try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)

rD <- rsDriver(browser = "firefox", verbose = FALSE, chromever = NULL)
remDr <- rD[["client"]]

res <- data.table(name = groupsMeta$name, gLines = numeric(), wLines = numeric())
i <- 1
for (i in seq_len(nrow(groupsMeta))) {
	remDr$navigate(paste0("https://www.google.com/search?client=firefox-b-d&q=", groupsMeta$name[i]))
	gHit <- remDr$findElement(using = "id", value = "result-stats")
	gHit <- as.numeric(gsub("\\D", "", sub(" Ergebnisse.*", "", gHit$getElementAttribute("innerHTML")[[1]])))
	res[i, gHits := gHit]
	remDr$navigate(paste0("https://en.wikipedia.org/wiki/", groupsMeta$name[i]))
	wLine <- remDr$findElement(using = "css", value = "html")
	wLine <- wLine$getElementAttribute("innerHTML")[[1]]
	wLine <- length(strsplit(wLine, split = "\n")[[1]])
	res[i, wLines := wLine]
	print(res[i])
}

# close all open ports
try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)

fwrite(res, "groupsWebOccurrences.txt")

# 9 Analyse taxon description data#################################################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets
library(openxlsx) # handle excel files
library(rethinking) # Bayesian inference
library(blavaan)
library(lavaan)

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxonomy"))

# authors for mosses and vascular plants from IPNI
ipniAuthors <- fread("botanical author names IPNI_2024-01-25.gz")

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# read in data

# groups metadata and response variables
groupsMeta <- fread("groupsMeta.txt")
groupsResponses <- fread("groupsResponses.txt")
groupsData <- fread("groupsData.txt")
groupsTotal <- colSums(groupsData[, -"year"])

# groups size classification, soil/endoparasitic vs free/non-parasitic, aquatic vs terrestrial
size <- data.table(read.xlsx("groupsPredictors.xlsx", sep.names = "_"))

# alternative group sizes from Brose et al.
sizeB <- fread("groupsLengthSize.txt")

# groups occurrence records
regions <- fread("groupsOccurrences.txt")

# public interest data from BiL
interestOcc <- fread("groupsLitOccurrences.txt")
interestTax <- fread("groupsLitTaxonomy.txt")

# public interest alternative sources
interest <- fread("groupsWebOccurrences.txt") # will include other interest data later on

# author numbers (although over time, no totals)
authors <- fread("authorData.txt")

# compare size data
size$name == sizeB$name
par(mfrow = c(1, 2))
plot(sizeB$`minLength` ~ I(10^size$`min_body_length_(mm)` / 1000))
abline(0, 1)
plot(sizeB$`maxLength` ~ I(10^size$`max_body_length_(mm)` / 1000))
abline(0, 1)
rm(sizeB)
colnames(size) <- gsub("\\([^\\(\\)]+\\)", "", colnames(size))
colnames(size) <- gsub("^_|_$", "", gsub("_{2,}", "_", colnames(size)))

# use data from Wikipedia because it is more complete and not (as) prone to sampling bias
# as Brose's data

# compare public interest data
str(interestOcc)
str(interestTax)
litOccs <- colSums(interestOcc[, -"year"], na.rm = TRUE)
litOccs <- data.table(sum = litOccs, name = interestTax$group)
litOccs <- tapply(litOccs$sum, litOccs$name, sum)
litOccs <- litOccs[names(litOccs) != ""]
interest[, litOcc := 0]
interest[sapply(names(litOccs), function(x) which(interest$name == x)), litOcc := litOccs]
plot(interest[, -"name"])
plot(interest$gHits, interest$litOcc)
plot(interest$gHits, log(interest$litOcc))

# Google hits and log(literature occurrences) show similar patterns,
# while wikipedia lines do not
interest[, wLines := NULL]

# check author numbers from LifeGate
colSums(authors[, -"year"])
par(mfrow = c(1, 1))
plot(NULL, xlim = range(authors$year), ylim = range(authors[, -"year"]), xlab = "year", ylab = "authors")
for (i in seq_len(ncol(authors))) {
	if (i > 1) {
		lines(authors$year, authors[[i]])
	}
}
# add means
abline(h = colMeans(authors[, -"year"]), lty = 2)
authorMeans <- colMeans(authors[, -"year"])
authorVar <- sapply(authors[, -"year"], var)
authorCv <- authorVar^.5 / authorMeans
colSums(authors[, -"year"])

# calculate relationship between hypothetical total author number and mean per year for
# hypothetical author numbers
# reps <- 100
# auths <- c(1,100000)
# years <- c(5,30)
# starts <- c(1753,2017)
# res <- data.table(auths=rep(0,reps),authorMeans=rep(0,reps))
# i <- 100
# for (i in seq_len(reps)){
# 	auth <- round(runif(1,auths[1],auths[2]))
# 	year <- round(runif(auth,years[1],years[2]))
# 	#start <- runif(auth,0,1)
# 	start <- rexp(auth,rate=10)
# 	start[start > 1] <- 1
# 	start <- (1 - start) * c(2017 - 1753) + 1753
# 	#hist(start)
# 	res[i, auths := auth]
# 	yearSum <- rep(0,length(seq(starts[1],starts[2])))
# 	for (j in seq_len(starts[2]-starts[1])){
# 		yearSum[j] <- sum(start <= starts[1] + j & start + year >= starts[1] + j)
# 	}
# 	res[i, authorMeans := mean(yearSum)]
# 	print(res[i])
# }
# cor(res)
## (reg <- lm(res))
# mean(res$auths/res$authorMeans)
# plot(res$auths~res$authorMeans)
# abline(reg)
# total author numbers and mean authors per year are always highly correlated,
# independent of the actual distribution of authors across the time period
# the actual publication number per author is 1/5 per year (see script 2),
# therefore, res$authorMeans should be divided by 5
# however, as the relationship between total species numbers and mean authors is very different
# across taxa, it is likely that this number (1/5) differs across groups
# mean(res$auths / (res$authorMeans / 5))
# nrow(ipniAuthors) # total number of authors for terrestrial plants (Tracheophyta, Marchantiophyta, Bryophyta)
# nrow(ipniAuthors) / sum(authorMeans[grepl("Trache|March|Bryo", names(authorMeans))])
# the relationship does not hold perfectly here, but it is comparable
# IPNI also has some bad data included and may therefore be too high, and author numbers for mosses are too low
# anyway, as this would only result in a linear transformation of the data, it is not needed

# create a common data.table
dat <- data.table(
	name = groupsResponses$name,
	currDesc = groupsTotal,
	meanAuthors = authorMeans,
	varAuthors = authorVar,
	cvAuthors = authorCv,
	authorFinal = groupsMeta$finalAuthorData,
	litOcc = interest$litOcc,
	gHits = interest$gHits,
	lengthMin = size$min_body_length,
	lengthMax = size$max_body_length,
	lengthMedian = size$median_body_length,
	lengthRange = size$max_body_length - size$min_body_length,
	soilEndo = size$`not_soil-dwelling_or_endoparasitic_-_soil-dwelling_or_endoparasitic`,
	aqua = size$`terrestrial_-_aquatic`,
	occIn = regions$ratio,
	estFutureDesc = groupsResponses$estFutureDesc,
	tenPercDesc = groupsResponses$tenPercDesc,
	maxDescPace = groupsResponses$maxDescPace,
	descVar = groupsResponses$descVar
)
dat[, c("name", "currDesc")]

fwrite(dat, file = "groupsVariables.txt")

# compare mean author numbers and total number of species
# pdf("mean authors per year vs total description number.pdf", height=8.3,width=11.7)
par(mfrow = c(1, 2))
plot(NULL,
	xlim = range(dat$currDesc), ylim = range(dat$meanAuthors),
	xlab = "current descriptions", ylab = "mean authors w descriptions per year"
)
text(dat$currDesc, dat$meanAuthors, dat$name,
	cex = .75,
	col = c("red", "blue")[dat$authorFinal + 1]
)
abline(v = 20000, h = 25, col = "grey")
abline(lm(dat$meanAuthors ~ dat$currDesc + 0), col = "grey")
plot(NULL,
	xlim = c(0, 20000), ylim = c(0, 25), xlab = "current descriptions",
	ylab = "mean authors w descriptions per year"
)
text(dat$currDesc, dat$meanAuthors, dat$name,
	cex = .75,
	col = c("red", "blue")[dat$authorFinal + 1]
)
abline(lm(dat$meanAuthors ~ dat$currDesc + 0), col = "grey")
# dev.off()
# Plants have much more authors per species than other groups,
# i.e. less descriptions per author.

# Data analysis

# response variables
# - estimated future descriptions
# - ten percent of current descriptions
# - maximum description pace
# - description variance

# predictor variables
# - total author number/mean authors per year
# - google hits / literature occurrences
# - body length
# - endoparasitic/soil-dwelling or not
# - aquatic or not
# - occurrences inside Europe + North America / occurrences worldwide

# Hypotheses
#
# estimated future descriptions ~ ??
# reason -> ten percent of what we have now and maximum pace are no deterministic predictors,
# only complete species numbers
#
# ten percent of current descriptions ~ author number + body length + occurrences +
# endoparasitic/soil-dwelling or not + aquatic or not
# reason -> depends on technical limitations, manpower, hardness to find species
#
# maximum description pace ~ author number + endoparasitic/soil-dwelling or not + aquatic or not
# reason -> body length + occurrences are no longer a limitation these days
#
# description variance ~ author number +ok (bad)+
# reason -> with few authors, individual authors have more impact on the description pace
#
# author number ~ public interest +ok (good)+
# reason -> more public interest should lead to more funding and scientists working on the topic

cors <- round(cor(dat[, -c("name", "authorFinal")]), 2)
cors[abs(cors) < .3] <- 0
cors

# pdf("variable pairs.pdf",width=20,height=9)
par(mfrow = c(3, 6))
par(mar = c(4, 4, 2, 3))
for (i in seq_len(ncol(dat))) {
	if (is.numeric(dat[[i]])) {
		for (j in seq_len(ncol(dat))) {
			if (is.numeric(dat[[j]]) && i != j) {
				plot(NULL, xlim = range(dat[[j]]), ylim = range(dat[[i]]), xlab = colnames(dat)[j], ylab = colnames(dat)[i])
				text(dat[[j]], dat[[i]], labels = dat$name, cex = .75, col = groupsMeta$color)
				points(dat[[j]], dat[[i]], cex = .75)
			}
		}
	}
}
# dev.off()

# Playground - test assumed relationships with OLS models first

# author number vs public interest
par(mfrow = c(2, 2))
reg <- lm((meanAuthors) ~ (gHits) + 0, data = dat)
plot((meanAuthors) ~ (gHits), data = dat, main = paste0("Google hits - R\xc2\xb2 = ", round(summary(reg)$adj.r.squared, 2)))
abline(reg, col = "red", lwd = 2)
reg <- lm((meanAuthors) ~ (litOcc) + 0, data = dat)
plot((meanAuthors) ~ (litOcc), data = dat, main = paste0("LitOcc - R\xc2\xb2 = ", round(summary(reg)$adj.r.squared, 2)))
abline(reg, col = "red", lwd = 2)
reg <- lm(log(meanAuthors) ~ log(gHits) + 0, data = dat)
plot(log(meanAuthors) ~ log(gHits), data = dat, main = paste0("Log Google hits - R\xc2\xb2 = ", round(summary(reg)$adj.r.squared, 2)))
abline(reg, col = "red", lwd = 2)
reg <- lm(log(meanAuthors) ~ log(litOcc + 1) + 0, data = dat)
plot(log(meanAuthors) ~ log(litOcc + 1),
	data = dat,
	main = paste0("Log LitOcc - R\xc2\xb2 = ", round(summary(reg)$adj.r.squared, 2))
)
abline(reg, col = "red", lwd = 2)
# highly correlated, best seen in log-log, literature occurrences

# description variability vs meanAuthors
par(mfrow = c(1, 2))
reg <- lm(descVar ~ meanAuthors, data = dat)
summary(reg) # no relationship
plot(descVar ~ meanAuthors, data = dat, main = paste0("R\xc2\xb2 = ", round(summary(reg)$adj.r.squared, 2)))
abline(reg, col = "red", lwd = 2)
abline(v = seq(0, 300, 10), col = "grey", lty = 2)
# no correlation, but a clear trend of high variability in description variance
# at small mean author numbers and lower variability at large mean author numbers
plot(abs(resid(reg)) ~ fitted(reg))
# try weighted regression to account for heteroskedasticity
weights <- 1 / resid(reg)^2
regW <- lm(descVar ~ meanAuthors, weight = weights, data = dat)
summary(regW) # looks good
par(mfrow = c(2, 2))
plot(descVar ~ meanAuthors, data = dat, main = paste0("Unweighted R\xc2\xb2 = ", round(summary(reg)$adj.r.squared, 2)))
abline(reg, col = "red", lwd = 2)
plot(descVar ~ meanAuthors, data = dat, main = paste0("Weighted R\xc2\xb2 = ", round(summary(regW)$adj.r.squared, 2)))
abline(regW, col = "red", lwd = 2)
plot(abs(resid(reg)) ~ fitted(reg))
plot(abs(I(resid(regW) * weights)) ~ fitted(regW))

# maximum description pace ~ author number + endoparasitic/soil-dwelling or not + aquatic or not
# reason -> body length + occurrences are no longer a limitation these days
par(mfrow = c(1, 1))
dat[, maxDescPace := maxDescPace / colSums(groupsData)[-1]]
plot(maxDescPace ~ I(maxDescPace * colSums(groupsData)[-1]), data = dat)
text(I(dat$maxDescPace * colSums(groupsData)[-1]), dat$maxDescPace, labels = dat$name, cex = .75, col = groups$color)
# it might be more interesting to investigate the relative description pace
# instead of the absolute description pace. the absolute one is just related
# to the number of authors
plot(colSums(groupsData)[-1] ~ dat$meanAuthors)
text(dat$meanAuthors, colSums(groupsData)[-1], labels = dat$name, cex = .75, col = groups$color)
reg <- lm(maxDescPace ~ tenPercDesc, data = dat)
summary(reg)
reg <- lm(maxDescPace ~ tenPercDesc + soilEndo + aqua, data = dat)
summary(reg)
reg <- lm(maxDescPace ~ tenPercDesc + aqua, data = dat)
summary(reg)
# it looks as if nothing can be learned from soildwelling/endoparasitic
par(mfrow = c(2, 3))
plot(tenPercDesc ~ meanAuthors, data = dat)
plot(maxDescPace ~ aqua, data = dat)
plot(maxDescPace ~ lengthMedian, data = dat)
plot(maxDescPace ~ lengthMin, data = dat)
plot(maxDescPace ~ lengthMax, data = dat)
summary(lm(maxDescPace ~ meanAuthors, data = dat))
summary(lm(tenPercDesc ~ meanAuthors + lengthMedian, data = dat))
# there seems to hardly be an answer on the maxDescPace,
# if it is not taken to be the absolute value

# ten percent of current descriptions ~ author number + body length + occurrences +
# endoparasitic/soil-dwelling or not + aquatic or not
# reason -> depends on technical limitations, manpower, hardness to find species
summary(lm(tenPercDesc ~ meanAuthors + lengthMedian + occIn + soilEndo + aqua, data = dat))
summary(lm(tenPercDesc ~ meanAuthors + lengthMedian + occIn, data = dat))
# this seems to be ok
summary(lm(tenPercDesc ~ meanAuthors + lengthMedian + occIn + I(log(aqua)), data = dat))
# study soil endo and aqua for relationships with something (also check outliers!)

summary(lm(estFutureDesc ~ maxDescPace, data = dat))

# Concept figure (true diversity is unknown, could be modelled as latent variable)
#
# See manuscripts/taxon description dates/path diagram.pptx

# Analyse in a Bayesian framework
str(dat)

# create scaled dataset
dats <- copy(dat)
numCols <- colnames(dats)[sapply(dats, is.numeric)]
dats[, (numCols) := lapply(.SD, scale), .SDcols = (numCols)]
str(dats)
round(cor(dats[, -c("name")]), 2)

# create short name dataset
datn <- data.table(
	C = dats$currDesc,
	A = dats$meanAuthors,
	DV = dats$descVar,
	L = dats$litOcc,
	B = dats$lengthMedian, # body size
	S = dats$soilEndo,
	AQ = dats$aqua,
	O = dats$occIn,
	TT = dats$tenPercDesc, # time until ten percent
	MP = dats$maxDescPace, # maximum description pace
	EF = dats$estFutureDesc
)

# rethinking implementation (no SEM)

# reg1 <- ulam(
# 	alist(	# relationships between variables
# 		L ~ normal(alpha, rho),
# 		alpha <- a * C,
# 		A ~ normal(beta, sigma),
# 		beta <-  b * L + c * C,
# 		DV ~ normal(gamma, tau),
# 		tau <-  d * exp(-A) + e * exp(-O),
# 		TT ~ normal(delta, ypsilon),
# 		delta <- f * B + g * A + h * O + ii * S + j * AQ + k * DV,
# 		MP ~ normal(epsilon, phi),
# 		epsilon <- l * B + m * A + n * O + o * S + p * AQ + q * TT,
# 		EF ~ normal(zeta, chi),
# 		zeta <- r * TT + s * MP,
# 		# regression priors
# 		c(a, b, c, f, g, h, ii, j, k, l, m, n, o, p, q, r, s) ~ dnorm(0,1),
# 		c(d, e) ~ dexp(0.1),
# 		# variance priors
# 		c(rho, sigma, gamma, ypsilon, phi, chi) ~ dexp(1)
# 	),
# 	data = datn, chains = 4, log_lik = TRUE
# )

# (res1 <- precis(reg1,depth=2))
# mean   sd  5.5% 94.5% rhat ess_bulk
# s        0.70 0.12  0.51  0.89    1  2534.55
# r        0.04 0.12 -0.13  0.23    1  2267.37
# q        0.84 0.18  0.53  1.14    1  2081.66
# p        0.14 0.13 -0.09  0.35    1  2611.31
# o       -0.12 0.14 -0.35  0.11    1  2825.66
# n       -0.18 0.16 -0.43  0.09    1  2266.67
# m        0.14 0.14 -0.09  0.37    1  2941.90
# l        0.42 0.17  0.14  0.68    1  1924.44
# k       -0.38 0.09 -0.52 -0.24    1  3574.56
# j       -0.04 0.09 -0.19  0.10    1  2791.82
# ii       0.09 0.10 -0.07  0.26    1  2466.16
# h        0.33 0.10  0.17  0.49    1  2292.19
# g       -0.37 0.10 -0.52 -0.21    1  2951.02
# f       -0.42 0.11 -0.59 -0.25    1  2287.30
# c        0.74 0.06  0.65  0.83    1  3287.11
# b        0.36 0.06  0.27  0.45    1  3105.01
# a        0.35 0.14  0.13  0.57    1  3439.43
# e        0.24 0.09  0.12  0.38    1  2825.21
# d        0.75 0.12  0.57  0.96    1  2785.12
# chi      0.71 0.07  0.60  0.83    1  3649.51
# phi      0.85 0.09  0.72  1.01    1  3120.19
# ypsilon  0.59 0.07  0.49  0.70    1  2377.10
# gamma    0.04 0.04  0.00  0.11    1  2072.67
# sigma    0.36 0.04  0.30  0.43    1  3460.39
# rho      0.95 0.11  0.80  1.13    1  4173.36

# len1 <- length(res1@.Data[[1]])
# par(mfrow=c(1,1))
# plot(NULL,xlim=c(-1.5,1.5),ylim=c(0,7),xlab="",ylab="")
# xseq <- seq(-1.5,1.5,l=500)
# for (i in seq_len(len1)){
# 	if (nchar(res1@row.names[i]) < 3){
# 		yseq <- dnorm(xseq,res1@.Data[[1]][i],res1@.Data[[2]][i])
# 		if (res1@.Data[[3]][i] * res1@.Data[[4]][i] < 0){
# 			polygon(c(xseq,rev(xseq)),c(rep(0,length(xseq)),rev(yseq)),col=paste0(rainbow(len1)[i],"33"), border = "red")
# 		} else {
# 			polygon(c(xseq,rev(xseq)),c(rep(0,length(xseq)),rev(yseq)),col=paste0(rainbow(len1)[i],"33"), border = "dark grey")
# 		}
# 		text(res1@.Data[[1]][i], 1 / (res1@.Data[[2]][i] * sqrt(2 * pi)),res1@row.names[i])
# 	}
# }

# reg2 <- ulam(
# 	alist(	# relationships between variables
# 		L ~ normal(alpha, rho),
# 		alpha <- a * C,
# 		A ~ normal(beta, sigma),
# 		beta <-  b * L + c * C,
# 		DV ~ normal(gamma, tau),
# 		tau <-  d * exp(-A) + e * exp(-O),
# 		TT ~ normal(delta, ypsilon),
# 		delta <- f * B + g * A + h * O + k * DV,
# 		MP ~ normal(epsilon, phi),
# 		epsilon <- l * B + q * TT,
# 		EF ~ normal(zeta, chi),
# 		zeta <- s * MP,
# 		# regression priors
# 		c(a, b, c, f, g, h, k, l, q, s) ~ dnorm(0,1),
# 		c(d, e) ~ dexp(0.1),
# 		# variance priors
# 		c(rho, sigma, gamma, ypsilon, phi, chi) ~ dexp(1)
# 	),
# 	data = datn, chains = 4, log_lik = TRUE
# )

# (res2 <- precis(reg2,depth=2))
# mean   sd  5.5% 94.5% rhat ess_bulk
# s        0.72 0.11  0.56  0.89 1.01  2657.38
# q        0.63 0.15  0.38  0.87 1.00  2242.07
# l        0.32 0.15  0.09  0.56 1.00  1969.44
# k       -0.38 0.09 -0.52 -0.24 1.00  1897.78
# h        0.37 0.09  0.22  0.52 1.00  2208.08
# g       -0.34 0.10 -0.49 -0.20 1.00  1920.52
# f       -0.45 0.10 -0.61 -0.30 1.01  2267.86
# c        0.74 0.06  0.65  0.83 1.00  2443.34
# b        0.36 0.06  0.27  0.45 1.00  2162.38
# a        0.35 0.14  0.14  0.57 1.00  2631.49
# e        0.24 0.09  0.11  0.40 1.00  2169.46
# d        0.74 0.13  0.56  0.96 1.00  2257.95
# chi      0.70 0.08  0.59  0.83 1.00  2600.21
# phi      0.86 0.09  0.72  1.03 1.00  2310.79
# ypsilon  0.58 0.06  0.49  0.69 1.00  2160.27
# gamma    0.04 0.04  0.00  0.12 1.00  1351.30
# sigma    0.36 0.04  0.30  0.43 1.00  2304.37
# rho      0.95 0.10  0.81  1.12 1.00  2476.98

# len2 <- length(res2@.Data[[1]])
# par(mfrow=c(1,1))
# plot(NULL,xlim=c(-1.5,1.5),ylim=c(0,7),xlab="",ylab="")
# xseq <- seq(-1.5,1.5,l=500)
# for (i in seq_len(len2)){
# 	if (nchar(res2@row.names[i]) < 3){
# 		yseq <- dnorm(xseq,res2@.Data[[1]][i],res2@.Data[[2]][i])
# 		if (res2@.Data[[3]][i] * res2@.Data[[4]][i] < 0){
# 			polygon(c(xseq,rev(xseq)),c(rep(0,length(xseq)),rev(yseq)),col=paste0(rainbow(len2)[i],"33"), border = "red")
# 		} else {
# 			polygon(c(xseq,rev(xseq)),c(rep(0,length(xseq)),rev(yseq)),col=paste0(rainbow(len2)[i],"33"), border = "dark grey")
# 		}
# 		text(res2@.Data[[1]][i], 1 / (res2@.Data[[2]][i] * sqrt(2 * pi)),res2@row.names[i])
# 	}
# }

# way to present data
# precis(reg1,depth=2)
# post <- extract.samples(reg1) #get samples from the posterior distribution
# par(mfrow=c(1,2))
# plot(datn$P,datn$A,col=2,lwd=3,xlab="log(literature occurences)",ylab="log(mean author number per year)")
# for (i in 1:100) abline(a=0, b=post$a[i],lwd=1)
# dens(post$a,lwd=2,col=2)

# Blavaan implementation

# Blavaan cannot model the effect of a variable on the variance of another.
# Therefore, a workaround would be to estimate the variance and regress on
# the variance instead of the variable itself

## test estimation of the variance
# v1 <- seq(1,100,l=1000)
# v2 <- rnorm(1000,0, sd = 1.5 * v1)
# plot(v2 ~ v1)
# ints <- 20
# sdInts <- length(v2)/ints
# sdVals <- rep(NA,ints)
# for (i in seq_along(sdVals)) sdVals[i] <- sd(v2[seq_len(sdInts) + sdInts * (i - 1)])
# sdPrds <- rep(NA,ints)
# for (i in seq_along(sdPrds)) sdPrds[i] <- mean(v1[seq_len(sdInts) + sdInts * (i - 1)])
# dat <- data.table(sdPrds,sdVals)
# plot(sdVals~sdPrds,data=dat)
# model <- '
# 		# relationships between measured variables
# 		sdVals ~ sdPrds
# 		'
## fit the model to the data
# fit <- bsem(model, data = dat, sample=2000)
# summary(fit)
# This approach works, but there is of course some additional uncertainty associated,
# especially considering that the real data is very sparse

# do the estimation for the real data

# I want to predict description variability with two predictors, which are A (the number
# of author) and O (the occurrence of species in Europe)
# I need to calculate mean variability for each combination of A and O, or for certain intervals
# then, I can use this variable to regress it on A and O

# get means within intervals
datn[, DVVar := numeric()]
range(datn$A)
range(datn$O)
for (i in seq(-1, 6, l = 6)) {
	for (j in seq(-2, 2.5, l = 6)) {
		datn[A > i & A < i + 1 & O > j & O < j + 1, DVVar := mean(DV, na.rm = TRUE)]
	}
}

plot(datn$A, datn$O, cex = 10 * (datn$DVVar - min(datn$DVVar, na.rm = TRUE)) / (max(datn$DVVar, na.rm = TRUE) - datn$DVVar - min(datn$DVVar, na.rm = TRUE)))
points(datn$A, datn$O, pch = 16)
abline(v = seq(-1, 6, l = 6), lty = 2)
abline(h = seq(-2, 2.5, l = 6), lty = 2)
par(mfrow = c(1, 2))
plot(datn$DVVar ~ datn$A, pch = 16)
plot(datn$DVVar ~ datn$O, pch = 16)

model1 <- "
	# latent variable definitions
	# DIV =~ C + TT + MP
	# relationships between measured variables
	L ~ a * C
	A ~  b * L + c * C
	DVVar ~  d * A + e * O
	TT ~ f * B + g * A + h * O + ii * S + j * AQ + k * DVVar
	MP ~ l * B + m * A + n * O + o * S + p * AQ + q * TT
	EF ~ r * TT + s * MP
"

fit1 <- bsem(model1, data = datn)
summary(fit1)
# fitMeasures(fit)
# plot(fit1,plot.type = "trace")
# good overall fit, but several problematic coefficients

model2 <- "
	# latent variable definitions
	# DIV =~ C + TT + MP
	# relationships between measured variables
	L ~ a * C
	A ~  b * L + c * C
	DVVar ~  d * A + e * O
	TT ~ f * B + g * A + h * O + k * DVVar
	MP ~ l * B + q * TT
	EF ~ s * MP
"

fit2 <- bsem(model2, data = datn)
summary(fit2)
fitMeasures(fit2)
# plot(fit2,plot.type = "trace")
# good overall fit, could try to include latent variable "true diversity"

# model3 <- '
# 	# latent variable definitions
# 	DIV =~ C + TT + MP
# 	# relationships between measured variables
# 	L ~ a * C
# 	A ~  b * L + c * C
# 	DVVar ~  d * A + e * O
# 	TT ~ f * B + g * A + h * O + k * DVVar
# 	MP ~ l * B + q * TT
# 	EF ~ s * MP
#'
#
# fit3 <- bsem(model3, data = datn)
# summary(fit3)
#fitMeasures(fit3)
#plot(fit3,plot.type = "trace")
# problem with effective sample size adding latent variable, and not
# needed for this analysis

# Entry point ###################################

# load in libraries
library(data.table) # handle large datasets
library(openxlsx) # handle excel files
library(rethinking) # Bayesian inference
library(blavaan)
library(lavaan)

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# save.image("temp.RData")
load("temp.RData")

#################################################
