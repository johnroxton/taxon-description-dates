# These scripts use the taxon description dates that are part of the LifeGate data assembled by
# Martin Freiberg to visualize cumulative description curves and link them to potential predictors.
# Before this, data is compared between LCVP and LifeGate, and the estimation/extrapolation
# process of a paper about species numbers in Nigeria is examined.
#
# Author: David Schellenberger Costa
###################################################################################################
# 1 Visualize data from LCVP and total numbers from LifeGate
# 2 Reproduce the analyses from "Trends in botanical exploration in Nigeria forecast over
#   1000 yet undescribed vascular plant species", https://doi.org/10.1093/aob/mcad106
# 3 Visualize LifeGate description dates
# 4 Approximate distribution of description dates using functions
# 5 Source distribution data from GBIF
# 6 Source body size data from Ulrich Brose's lab
# 7 Source data from BiL explorer
# 8 Source measure of public interest based on Wikipedia article lengths and Google hits
# 9 Analyze taxon description data
# 10 Compare species numbers between LifeGate, GBIF, and CoL
# 11 Compare checked species numbers for Lepidoptera and Bryophyta
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

# 2 Reproduce the analyses from "Trends in botanical exploration in Nigeria forecast over##########
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
pubRange <- data.table(author = authors, firstPub = NA_real_, lastPub = NA_real_)
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

# While I could try to redo the analysis from Bello et al., there is not much to learn here, as it suffers
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
library(RColorBrewer) # color palettes
library(treemap)

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# read in metadata
groupsMeta <- data.table(file = list.files(path = "Arten pro Jahr beschrieben"))
groupsMeta <- groupsMeta[!grepl("_CoL_GBIF", groupsMeta$file)]
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
groupCols <- brewer.pal(12, "Paired")[c(12, 4, 6, 1, 2)] # colors of groups of groups
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
		if (grepl("Lepidoptera|Bryophyta", groupsMeta$file[i])) {
			groupsData <- cbind(groupsData, fread(paste0(
				"Arten pro Jahr beschrieben/",
				sub("Jahr", "Jahr_CoL_GBIF", groupsMeta$file[i])
			))$V2)
		} else {
			groupsData <- cbind(groupsData, fread(paste0(
				"Arten pro Jahr beschrieben/",
				groupsMeta$file[i]
			))$V2)
		}
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

# get icons of groups from phylopic.org
# IDs were retrieved manually
# groupsMeta[, icon := c(
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
# # get and save icons
# for (i in groupsVec){
# 	save_phylopic(get_phylopic(uuid = groupsMeta$icon[i]),path=paste0("icons2/",groupsMeta$name[i],".svg"))
# }

parBackup <- par() # save graphics parameters for multiple trials

# plot overview phylogeny

# choose colors for monophyletic groups
# for coloring purposes change node labels into numbers (can be omitted after first run)
# pt$node.label <- length(pt$tip.label) + 1:length(pt$node.label)
xlim <- 1
ylim <- 1.6
cols <- rep("black", nrow(pt$edge))
groupCols <- brewer.pal(12, "Paired")[c(12, 4, 6, 1, 2)] # colors of groups of groups
groupNode <- c(90, 52, 76, 66, grep("Chor", pt$tip.label))
for (i in seq_along(groupCols)) {
	edgesOld <- c()
	edgesNew <- groupNode[i]
	while (length(edgesNew) != length(edgesOld) || any(edgesNew != edgesOld)) {
		edgesOld <- edgesNew
		edgesNew <- union(edgesOld, pt$edge[pt$edge[, 1] %in% edgesOld, 2])
	}
	cols[pt$edge[, 2] %in% edgesNew] <- groupCols[i]
}

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
	abline(h = log(c(0, 10, 100, 1000) + 1), lty = 1, col = "grey")
	abline(v = c(1800, 1900, 2000), lty = 1, col = "grey")
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
	abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 1, col = "grey")
	abline(v = c(1800, 1900, 2000), lty = 1, col = "grey")
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
# dev.off()

# plot groups together based on phylogeny
plotGroups <- list()
plotGroups[[1]] <- c(1:5, 44:47) # plants and fungi
plotGroups[[2]] <- c(6:13) # unicellular species and "lower" plants
plotGroups[[3]] <- c(37:43) # "lower" animals, crustaceans, myriapods, arachnids, tardigrads
plotGroups[[4]] <- c(14:20) # Chordata, Spiralia
plotGroups[[5]] <- c(21:28) # insects 1
plotGroups[[6]] <- c(29:36) # insects 2

# pdf("description history comparative.pdf",width=10,height=5)
cols <- c(brewer.pal(8, "Dark2"), "black")
par(mfrow = c(2, 3))
par(oma = c(5, 5, 0.1, 0.1))
par(mar = c(0, 0, 0, .45))
i <- 1
groupsMeta[color == "black", color := "#000000"]
for (i in seq_along(plotGroups)) {
	plot(NULL, xlim = range(groupsData$year), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
	abline(h = c(0.2, 0.4, 0.6, 0.8, 1), lty = 3)
	abline(v = c(1800, 1850, 1900, 1950, 2000), lty = 3)
	if ((i - 1) %% 3 == 0) axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
	if (i > 3) axis(1, at = c(1800, 1900, 2000))
	for (j in groupsVec) {
		lines(groupsData$year, cumsum((groupsData[[j + 1]])) / sum(groupsData[[j + 1]]), col = "lightgrey")
	}
	k <- 1
	for (j in plotGroups[[i]]) {
		lines(groupsData$year, cumsum((groupsData[[j + 1]])) / sum(groupsData[[j + 1]]), col = cols[k], lwd = 2)
		text(1790, 1 - k * 0.07, groupsMeta$name[j], adj = 0, col = cols[k], font = 2, cex = 1.1)
		img <- readPNG(paste0("group icons/", groupsMeta$name[j], ".png"))
		scale <- 30 / diff(range(groupsData$year)) * 0.6
		colRGB <- col2rgb(cols[k]) / 255
		for (l in 1:3) img[, , l] <- colRGB[[l]]
		rasterImage(img, 1760, 1 - k * 0.07 - 0.03, 1775, 1 - k * 0.07 - 0.03 + scale)
		k <- k + 1
	}
	col1 <- groupsMeta$color[plotGroups[[i]][1]]
	col2 <- groupsMeta$color[plotGroups[[i]][length(plotGroups[[i]])]]
	colParts1 <- c(substr(col1, 2, 3), substr(col1, 4, 5), substr(col1, 6, 7))
	colParts2 <- c(substr(col2, 2, 3), substr(col2, 4, 5), substr(col2, 6, 7))
	colRange <- list()
	for (k in 1:3) colRange[[k]] <- sprintf("%02x", as.hexmode(round(seq(paste0("0x", colParts2)[k], paste0("0x", colParts1)[k], l = 100))))
	for (k in seq_len(100)) {
		currCol <- paste0("#", colRange[[1]][k], colRange[[2]][k], colRange[[3]][k])
		rect(2028, -.1 + 1.2 / 100 * (k - 1), 2033, -.1 + 1.2 / 100 * (k) + 0.01, col = currCol, xpd = TRUE, border = NA)
	}
}
# dev.off()

# plot groups based on curve shape
groupsDataRel <- groupsData[, -"year"]
groupsDataRel <- data.table(sapply(groupsDataRel, function(x) x / sum(x)))
groupsDataRel[, year := groupsData$year]
groupsDataRelSpeed <- as.matrix(groupsDataRel[, -"year"])
for (j in seq_len(nrow(groupsDataRel) - 4)) {
	groupsDataRelSpeed[j, ] <- colSums(groupsDataRel[j:(j + 0), -"year", drop = FALSE])
}
plotGroups <- list()
# groups that were over-sampled before 1900
i <- 1
res <- colSums(groupsDataRel[year <= 1900, -"year"])
plotGroups[[i]] <- c(order(res, decreasing = TRUE)[1:6])
# groups that were over-sampled after 1975
i <- i + 1
res <- colSums(groupsDataRel[year <= 1975, -"year"])
plotGroups[[i]] <- setdiff(order(res, decreasing = TRUE), plotGroups[[i - 1]])[1:9]
# groups that were under-sampled before 1975
i <- i + 1
res <- colSums(groupsDataRel[year <= 1975, -"year"])
plotGroups[[i]] <- c(order(res)[1:8])
# groups that strongly increased in specific time intervals
i <- i + 1
maxRate <- apply(groupsDataRelSpeed[1:which(groupsData$year == 2010), ], 2, max)
plotGroups[[i]] <- setdiff(order(maxRate, decreasing = TRUE), unlist(plotGroups[1:(i - 1)]))[1:8]
maxYears <- apply(groupsDataRelSpeed[1:which(groupsData$year == 2010), plotGroups[[i]]], 2, function(x) which(x == max(x))) - 0.5
# plot(cumsum(groupsDataRel[[7]]),type="l")
# abline(v=258-0.5)
# groups that slowly increased after 1900 in specific time intervals
i <- i + 1
minRate <- apply(groupsDataRelSpeed[which(groupsData$year == 1900):which(groupsData$year == 2010), ], 2, min)
plotGroups[[i]] <- setdiff(order(minRate, decreasing = FALSE), unlist(plotGroups[1:(i - 1)]))[1:8]
minYears <- apply(groupsDataRelSpeed[which(groupsData$year == 1900):which(groupsData$year == 2010), plotGroups[[i]]], 2, function(x) which(x == min(x))[1]) - 0.5
minYears <- minYears - 1753 + 1900
# rest
i <- i + 1
rest <- setdiff(1:47, unlist(plotGroups[1:(i - 1)]))
plotGroups[[i]] <- rest

pdf("description history comparative_functional.pdf", width = 12.3, height = 11.7)
cols <- c(brewer.pal(8, "Dark2"), rep("black", 4))
par(mfrow = c(3, 2))
par(oma = c(5, 5, 0.1, 0.1))
par(mar = c(0, 0, 0, 0))
i <- 5
groupsMeta[color == "black", color := "#000000"]
for (i in seq_along(plotGroups)) {
	plot(NULL, xlim = range(groupsData$year), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
	abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 3)
	abline(v = c(1800, 1850, 1900, 1950, 2000), lty = 3)
	if ((i - 1) %% 2 == 0) axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = 2)
	if (i > 4) axis(1, at = c(1800, 1900, 2000), cex.axis = 2, mgp = c(3, 1.75, 0))
	for (j in groupsVec) {
		lines(groupsData$year, cumsum((groupsData[[j + 1]])) / sum(groupsData[[j + 1]]), col = "lightgrey")
	}
	text(1760, 0.97, paste0("(", letters[i], ")"), adj = 0, font = 2, cex = 2)
	k <- 1
	for (j in plotGroups[[i]]) {
		lines(groupsData$year, cumsum((groupsData[[j + 1]])) / sum(groupsData[[j + 1]]), col = cols[k], lwd = 2)
		text(1790, 1 - (k + 1) * 0.07, groupsMeta$name[j], adj = 0, col = cols[k], font = 2, cex = 1.1)
		img <- readPNG(paste0("group icons/", groupsMeta$name[j], ".png"))
		scale <- 30 / diff(range(groupsData$year)) * 0.6
		colRGB <- col2rgb(cols[k]) / 255
		for (l in 1:3) img[, , l] <- colRGB[[l]]
		rasterImage(img, 1760, 1 - (k + 1) * 0.07 - 0.03, 1775, 1 - (k + 1) * 0.07 - 0.03 + scale)
		k <- k + 1
	}
	if (i == 1) abline(v = 1900, lty = 2, col = "black")
	if (i %in% c(2, 3)) abline(v = 1975, lty = 2, col = "black")
	if (i %in% c(4, 5)) {
		if (i == 4) selYears <- maxYears else selYears <- minYears
		k <- 1
		for (j in plotGroups[[i]]) {
			points(
				unlist(selYears)[k] + 1752,
				(cumsum((groupsData[[j + 1]])) / sum(groupsData[[j + 1]]))[unlist(selYears)[k] - 0.5], ,
				col = "black", cex = 1.6
			)
			segments(
				unlist(selYears)[k] + 1752,
				0,
				unlist(selYears)[k] + 1752,
				0.05, ,
				col = cols[k], lwd = 2
			)
			k <- k + 1
		}
	}
}
dev.off()

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

# plot number of descriptions per group
dtf <- data.table(index = colnames(groupsData)[-1], vSize = colSums(groupsData[, -"year"]), color = groupsMeta$color, fontsize = 10)
setorder(dtf, vSize)
dtf[vSize > 75000, fontsize := 10000]
dtf[vSize > 75000, fontsize := 1000]

# plot all groups, remove labels of smaller groups manually
# create a larger plot to be able to copy paste a fraction of it, whose labels are not displayed otherwise
pdf("number of descriptions per group.pdf", width = 117, height = 83)
treemap(dtf, index = "index", vSize = "vSize", title = "", vColor = "color", type = "color", border.col = "grey", border.lwds = 25, fontsize.labels = 200, lowerbound.cex.labels = 0)
dev.off()

# 4 Approximate distribution of description dates using functions##################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets
library(rethinking) # approximate description times with functions
library(png) # plot icons of taxonomic groups

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
# 	data = datq, chains = 4, log_lik = TRUE
# )
#
# precis(m2c, depth=2)

# modification of normal distribution with constant tail after maximum
# mu <- k * 1 / ((2 * 3.141593)^0.5 * a) * (exp(-0.5 * ((A - b) / a)^2) * (A <= b) + (A > b))

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
#  	print(paste0("Bertalanffy model ", i, "/", nrow(groupsMeta)))
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
# # normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
# m$m2 <- list()
# for (i in groupsVec) {
# 	print(paste0("Normal model ", i, "/", nrow(groupsMeta)))
# 	datn <- list(
# 		Y = scale(groupsData$year),
# 		D = cumsum(groupsData[[i + 1]]) / max(cumsum(groupsData[[i + 1]])),
# 		# minY = min(scale(groupsData$year)) # make sure max is not reached before max(datn$Y) - DEPRECATED
# 	)
# 	#datn$minY <- 0 # try to reproduce former results
# 	m$m2[[i]] <- ulam(
# 		alist(
# 			D ~ normal(mu, sigma),
# 			# note that + minY ensures the maximum is not reached before present, as in this data,
#           # minY = - max(Y), meaning max(Y) is added to b to ensure it has at least this size
# 			#mu <- (k / ((2 * 3.141593)^0.5 * a) + 1) * exp(-0.5 * ((Y - b + minY) / a)^2),
# 			# model specification to account for possible constant slope (no minY needed)
# 			mu <- (k / ((2 * 3.141593)^0.5 * a) + 1) * (exp(-0.5 * ((Y - b) / a)^2) * (Y <= b) + (Y > b)),
# 			k ~ exponential(1),
# 			a ~ exponential(1),
# 			b ~ exponential(1),
# 			sigma ~ exponential(1)
# 		),
# 		data = datn, chains = 4, log_lik = TRUE
# 	)
# 	#	print(coefs[[2]][,i])
# 	#	print(coef(m$m2[[i]]))
# 	#	print(precis(m$m2[[i]])
# 	# 	traceplot(m$m2[[i]],pars=c("k","a","b","sigma"),n_cols=2)
# }
# coefs[[2]] <- sapply(m$m2, coef)
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

# account for coefficient constraints in normal distribution fit - DEPRECATED due to constant function after maximum
# maximum not before 2017 - add minY to b of coefs[[2]]
# coefs[[2]][rownames(coefs[[2]]) == "b",] <- coefs[[2]][rownames(coefs[[2]]) == "b",] - min(scale(groupsData$year))
# maximum not smaller than data - add sqrt(2*pi)*a to k of coefs[[2]]
# coefs[[2]][rownames(coefs[[2]]) == "k",] <- coefs[[2]][rownames(coefs[[2]]) == "k",] + sqrt(2*pi)*coefs[[2]][rownames(coefs[[2]]) == "a",]
load("models3.RData")

# plot all groups separately with cumulative relative values and model approximations
# added information:
# - predicted asymptote/maximum
# - year of predicted asymptote/maximum
# - approximation quality

# plot description history
polyYear <- c(groupsData$year, rev(groupsData$year))
zeros <- rep(0, nrow(groupsData))

i <- 24
par(parBackup)
# pdf("description history fits_new.pdf", width = 11.7, height = 8.3)
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
	rasterImage(img, 1765, log(2001) / max(log(groupsData + 1)), 1785, log(2001) / max(log(groupsData + 1)) + scale)
	for (j in 1:3) {
		# shade(apply(mu, 2, PI), xseq, col = j + 1)
		if (j != 2) {
			mu <- link(m[[j]][[i]], data = list(Y = xseq))
		} else {
			# for normal distributions, it was necessary to scale the year variable to get good results
			mu <- link(m[[j]][[i]], data = list(Y = standardize(xseq), minY = min(scale(xseq))))
		}
		lines(xseq + min(groupsData$year), apply(mu, 2, mean), lwd = 3, col = j + 1)
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
	abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 1, col = "grey")
	abline(v = c(1800, 1900, 2000), lty = 1, col = "grey")
	if ((i - 1) %% 10 < 1) axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
	if (i > 37) axis(1, at = c(1800, 1900, 2000))
	polygon(polyYear, c(cumsum((groupsData[[i + 1]])) / sum(groupsData[[i + 1]]), zeros),
		border = NA, col = groupsMeta$col[i]
	)
	# for normal distributions, it was necessary to scale the year variable to get good results
	mu <- link(m[[2]][[i]], data = list(Y = standardize(xseq), minY = min(scale(xseq))))
	lines(xseq + min(groupsData$year), apply(mu, 2, mean), lwd = 3, col = "#999999", lty = 1)
	text(1800, log(5001) / max(log(groupsData + 1)), labels = groupsMeta$name[i], adj = c(0, 1))
	text(1800, log(1501) / max(log(groupsData + 1)), labels = colSums(groupsData)[i + 1], adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, 1)))
	rasterImage(img, 1755, log(2001) / max(log(groupsData + 1)), 1795, log(2001) / max(log(groupsData + 1)) + scale)
}
mtext("year", 1, 3, at = 1450)
mtext("cumulative descriptions", 2, 52.5, at = 2.7, xpd = TRUE)
# dev.off()
# save(list=c("coefs","m"),file="models3.RData")

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


for (i in seq_along(funNames)) {
	# estimated future descriptions relative to current descriptions
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
	# time until 10% found
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
	# maximum description pace
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
	# description variability
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

# plot an example on one page to explain the fitting process
# pdf("description fit example.pdf", width = 11.7, height = 8.3)
i <- 9
par(parBackup)
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(NULL,
	xlim = range(groupsData$year), ylim = c(0, 1), xaxt = "n", yaxt = "n",
	xlab = "year", ylab = "cumulative descriptions", cex.lab = 1.5
)
abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 3)
abline(v = c(1800, 1850, 1900, 1950, 2000), lty = 3)
axis(2, cex.axis = 1.5)
axis(1, at = c(1800, 1900, 2000), cex.axis = 1.5)
polygon(polyYear, c(cumsum((groupsData[[i + 1]])) / sum(groupsData[[i + 1]]), zeros),
	border = NA, col = "grey"
)
img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, 1))) * 11.7 / 8.3
j <- 2
# shade(apply(mu, 2, PI), xseq, col = j + 1)
# for normal distributions, it was necessary to scale the year variable to get good results
mu <- link(m[[j]][[i]], data = list(Y = standardize(xseq), minY = min(scale(xseq))))
lines(xseq + min(groupsData$year), apply(mu, 2, mean), lwd = 3, col = "black")
# add time until 10% found
abline(h = 0.1, col = "red", lwd = 2)
# add inflection point
xMax <- coefs[[2]][3, i] - coefs[[2]][2, i]
st <- standardize(groupsData$year - min(groupsData$year))
xMax[xMax > max(st)] <- max(st)
xMax <- xMax * attr(st, "scaled:scale") + attr(st, "scaled:center") + 1753
abline(v = xMax, col = "red", lwd = 2)
# add maximum
abline(h = funCoefs$estFutureDesc_normal[i], col = "red", lwd = 2)
rasterImage(img, 1755, .85, 1785, .85 + scale)
text(1790, .86 + scale / 2, labels = groupsMeta$name[i], adj = c(0, 0.5), cex = 1.5)
# dev.off()

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
par(parBackup)
# pdf("testDisplay.pdf")
b1 <- barplot(coefs[[1]][1, ],
	col = groupsMeta$color, border = NA, horiz = TRUE, xlim = c(0, 11), space = 0.6,
	main = "Bertalanffy function k"
)
textCols <- rep("black", nrow(groupsMeta))
textCols[groupsMeta$color == "black"] <- "white"
for (i in seq_len(ncol(coefs[[1]]))) {
	text(0.025, b1[i] - 1 + 1, groupsMeta$name[i], cex = 0.3, adj = 0, col = textCols[i])
	img <- readPNG(paste0("group icons/", groupsMeta$name[i], ".png"))
	scale <- max(coefs[[1]][1, ]) / max(b1)
	colRGB <- col2rgb(groupsMeta$color[i]) / 255
	for (j in 1:3) img[, , j] <- colRGB[[j]]
	rasterImage(img, coefs[[1]][1, i] + 0.075, b1[i] - 1, scale * 2 + coefs[[1]][1, i] + 0.075, b1[i] + 1, xpd = TRUE)
}
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

# get a measure of decrease in velocity due to idiosyncratic events (wars, other crisis)
# pdf("description pace anomalies.pdf",width=11.7,height=8.3)
plot(NULL, ylim = c(0, 1.9), xlim = range(groupsData$year), type = "l", xlab = "year", ylab = "", yaxt = "n")
numMax <- max(rowSums(speeds == -1), rowSums(speeds == 1))
axis(2, c(0, 0.5, 1, 1.5 + c(0, 1 / 6, -1 / 6, 1 / 3, -1 / 3)), c(0, 0.5, 1, 0, numMax / 2, numMax / 2, numMax, numMax))
text(1724, 0.5, "descriptions", srt = 90, cex = 1.2, xpd = TRUE)
text(1724, 1.5, "description pace anomalies", srt = 90, cex = 1.2, xpd = TRUE)
abline(h = 1.1)
for (i in seq_len(ncol(groupsData) - .8)) {
	lines(groupsData$year, cumsum(groupsData[[i + 1]]) / sum(groupsData[[i + 1]]), col = "grey")
}
lines(groupsData$year, cumsum(rowSums(groupsData[, -1])) / sum(groupsData[, -1]), lwd = 2)
# Franco-German war
polygon(c(c(1870:1871), rev(c(1870:1871))), c(rep(0, 2), rep(1.9, 2)), col = "#99555555", border = NA)
# first world war
polygon(c(c(1914:1918), rev(c(1914:1918))), c(rep(0, 5), rep(1.9, 5)), col = "#99555555", border = NA)
# second world war
polygon(c(c(1939:1945), rev(c(1939:1945))), c(rep(0, 7), rep(1.9, 7)), col = "#99555555", border = NA)
# first high-quality microscope from Carl Zeiss # no visible effect
# polygon(c(c(1872:1873),rev(c(1872:1873))),c(rep(0,2),rep(1.9,2)),col="#55559955", border=NA)
# first electrone microscope from Ernst Ruska and Max Knoll # no visible effect
# polygon(c(c(1931:1932),rev(c(1931:1932))),c(rep(0,2),rep(1.9,2)),col="#55559955", border=NA)
# publication of Grundzüge einer Theorie der pyhlogenetischen Systematik
polygon(c(c(1950:1951), rev(c(1950:1951))), c(rep(0, 2), rep(1.9, 2)), col = "#55559955", border = NA)
# publication of Illustrations of the Zoology of South Africa
polygon(c(c(1838:1839), rev(c(1838:1839))), c(rep(0, 2), rep(1.9, 2)), col = "#55559955", border = NA)
speeds <- matrix(0, nrow(groupsData), ncol(groupsData) - 1)
for (i in seq_len(ncol(groupsData) - 1)) {
	for (j in seq_along(groupsData[[i + 1]])) {
		if (j > 5 & j < nrow(groupsData)) {
			if (mean(groupsData[[i + 1]][j + c(0:4)], na.rm = TRUE) < 1 / 2 * mean(groupsData[[i + 1]][j - (1:5)], na.rm = TRUE)
			# & any(groupsData[[i + 1]][j+c(-5:4)]/sum(groupsData[[i + 1]]) > .001, na.rm=TRUE)
			# & groupsData[[i + 1]][j] < mean(groupsData[[i + 1]][groupsData[[i + 1]] > 0])
			& (mean(groupsData[[i + 1]][j + c(0:4)], na.rm = TRUE) - mean(groupsData[[i + 1]][j - (1:5)], na.rm = TRUE)) / sum(groupsData[[i + 1]]) < -0.001
			) {
				speeds[j, i] <- -1
			}
			if (mean(groupsData[[i + 1]][j + c(0:4)], na.rm = TRUE) > 2 * mean(groupsData[[i + 1]][j - (1:5)], na.rm = TRUE)
			# & any(groupsData[[i + 1]][j+c(-5:4)]/sum(groupsData[[i + 1]]) > .001, na.rm=TRUE)
			# & groupsData[[i + 1]][j] < mean(groupsData[[i + 1]][groupsData[[i + 1]] > 0])
			& (mean(groupsData[[i + 1]][j + c(0:4)], na.rm = TRUE) - mean(groupsData[[i + 1]][j - (1:5)], na.rm = TRUE)) / sum(groupsData[[i + 1]]) > 0.001
			) {
				speeds[j, i] <- 1
			}
		}
	}
}
abline(h = 1.5, lty = 2)
for (i in seq_len(nrow(speeds) - 10)) {
	if (i > 10) {
		segments(groupsData$year[i], 1.5, groupsData$year[i], 1.5 - sum(speeds[i, ] == -1) / max(rowSums(speeds == -1), rowSums(speeds == 1)) / 3, lwd = 1.1)
		segments(groupsData$year[i], 1.5, groupsData$year[i], 1.5 + sum(speeds[i, ] == 1) / max(rowSums(speeds == -1), rowSums(speeds == 1)) / 3, lwd = 1.1)
		if (sum(speeds[i, ] == -1) > 0.4 * max(rowSums(speeds == -1), rowSums(speeds == 1)) &
			sum(speeds[i - 1, ] == -1) <= 0.4 * max(rowSums(speeds == -1), rowSums(speeds == 1))) {
			text(groupsData$year[i] - 1, 1.5 - sum(speeds[i, ] == -1) / max(rowSums(speeds == -1), rowSums(speeds == 1)) / 3 - 0.14, groupsData$year[i], adj = 0, srt = 90, cex = 0.9)
		}
		if (sum(speeds[i, ] == 1) > 0.65 * max(rowSums(speeds == -1), rowSums(speeds == 1)) &
			sum(speeds[i - 1, ] == 1) <= 0.65 * max(rowSums(speeds == -1), rowSums(speeds == 1))) {
			text(groupsData$year[i] - 1, 1.5 + sum(speeds[i, ] == 1) / max(rowSums(speeds == -1), rowSums(speeds == 1)) / 3 + 0.02, groupsData$year[i], adj = 0, srt = 90, cex = 0.9)
		}
	}
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
# resOccs <- matrix(NA,length(groupsGBIFTaxonKeys) + length(MyriapodTaxonKeys),length(continents))
#
# for (i in seq_along(groupsGBIFTaxonKeys)){
# 	if (!is.na(groupsGBIFTaxonKeys)[i]){
# 		for (j in seq_along(continents)){
# 			resOccs[i,j] <- occ_count(taxonKey=groupsGBIFTaxonKeys[i],continent=continents[j])
# 		}
# 	}
# }
# for (i in seq_along(MyriapodTaxonKeys)){
# 	if (!is.na(MyriapodTaxonKeys)[i]){
# 		for (j in seq_along(continents)){
# 			resOccs[i+length(groupsGBIFTaxonKeys),j] <- occ_count(taxonKey=MyriapodTaxonKeys[i],continent=continents[j])
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

# There is a problem with oceans: Species found there are not assigned to any continent. Therefore, our data
# is biased, as terrestrial groups are treated differently from such the also live in the seas.
# I can calculate overall species numbers, but this has names with no occurrences and such with occurrences
# outside continents. A possibility is to us entries that also have coordinates, as this excludes the no occurrences.
# However, this also excludes some occurrences where locations were given like "Brazil", or "Patagonia".
# To make things comparable and calculate the proportion of (ignored) marine species, I therefore calculate everyting
# with coordinates. This eliminates about 8% of the data apparently (read in blog).

# GBIF occurrences with facetted species numbers - yes, species keys correspond to accepted species (in theory)
resSpecs <- list() # one entry per taxonomic group
for (i in seq_along(groupsGBIFTaxonKeys)) {
	if (i > 0) {
		print(names(groupsGBIFTaxonKeys)[i])
		resSpecs[[i]] <- list() # seven entries, one per continent
		if (!is.na(groupsGBIFTaxonKeys)[i]) {
			for (j in seq_along(continents)) {
				resSpecs[[i]][[j]] <- list() # all and with coordinates only
				resSpecs[[i]][[j]][[1]] <- NA # species keys for each species
				resSpecs[[i]][[j]][[2]] <- NA
				# continent with or without coordinates
				offset <- 0
				repeat {
					res <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=", groupsGBIFTaxonKeys[i], "&continent=", continents[j], "&limit=0&facet=speciesKey&facetLimit=500000&facetOffset=", offset))
					resSpecs[[i]][[j]][[1]] <- c(resSpecs[[i]][[j]][[1]], as.numeric(sapply(res$facets[[1]]$counts, function(x) x$name)))
					if (length(res$facets[[1]]$counts) < 500000) break else offset <- offset + 500000
				}
				# continent with coordinates
				offset <- 0
				repeat {
					res <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=", groupsGBIFTaxonKeys[i], "&hasCoordinate=True&continent=", continents[j], "&limit=0&facet=speciesKey&facetLimit=500000&facetOffset=", offset))
					resSpecs[[i]][[j]][[2]] <- c(resSpecs[[i]][[j]][[2]], as.numeric(sapply(res$facets[[1]]$counts, function(x) x$name)))
					if (length(res$facets[[1]]$counts) < 500000) break else offset <- offset + 500000
				}
			}
			# whole world with coordinates
			resSpecs[[i]][[j + 1]] <- NA
			offset <- 0
			repeat {
				res <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=", groupsGBIFTaxonKeys[i], "&hasCoordinate=TRUE&limit=0&facet=speciesKey&facetLimit=500000&facetOffset=", offset))
				resSpecs[[i]][[j + 1]] <- c(resSpecs[[i]][[j + 1]], as.numeric(sapply(res$facets[[1]]$counts, function(x) x$name)))
				if (length(res$facets[[1]]$counts) < 500000) break else offset <- offset + 500000
			}
		}
	}
}

# check whether more information could be gained by using country instead of continent
# countries<-c("AF","AX","AL","DZ","AS","AD","AO","AI","AQ","AG","AR","AM","AW","AU","AT","AZ","BS",
# 	"BH","BD","BB","BY","BE","BZ","BJ","BM","BT","BO","BQ","BA","BW","BV","BR","IO","BN","BG","BF",
# 	"BI","KH","CM","CA","CV","KY","CF","TD","CL","CN","CX","CC","CO","KM","CD","CG","CK","CR","CI",
# 	"HR","CU","CW","CY","CZ","DK","DJ","DM","DO","EC","EG","SV","GQ","ER","EE","ET","FK","FO","FJ",
# 	"FI","FR","GF","PF","TF","GA","GM","GE","DE","GH","GI","GR","GL","GD","GP","GU","GT","GG","GN",
# 	"GW","GY","HT","HM","VA","HN","HK","HU","IS","IN","ID","IR","IQ","IE","IM","IL","IT","JM","JP",
# 	"JE","JO","KZ","KE","KI","KP","KR","KW","KG","LA","LV","LB","LS","LR","LY","LI","LT","LU","MO",
# 	"MK","MG","MW","MY","MV","ML","MT","MH","MQ","MR","MU","YT","MX","FM","MD","MC","MN","ME","MS",
# 	"MA","MZ","MM","NA","NR","NP","NL","NC","NZ","NI","NE","NG","NU","NF","MP","NO","OM","PK","PW",
# 	"PS","PA","PG","PY","PE","PH","PN","PL","PT","PR","QA","RE","RO","RU","RW","BL","SH","KN","LC",
# 	"MF","PM","VC","WS","SM","ST","SA","SN","RS","SC","SL","SG","SX","SK","SI","SB","SO","ZA","GS",
# 	"SS","ES","LK","SD","SR","SJ","SZ","SE","CH","SY","TW","TJ","TZ","TH","TL","TG","TK","TO","TT",
# 	"TN","TR","TM","TC","TV","UG","UA","AE","GB","US","UM","UY","UZ","VU","VE","VN","VG","VI","WF",
# 	"EH","YE","ZM","ZW","AA","XK","XZ","ZZ")
# resTestCountry <- list()
# for (i in seq_along(countries)){
# 	resTestCountry[[i]] <- NA
# 	offset <- 0
# 	repeat {
# 		res <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=1470&country=", countries[i], "&limit=0&facet=speciesKey&facetLimit=500000&facetOffset=", offset))
# 		resTestCountry[[i]] <- c(resTestCountry[[i]], as.numeric(sapply(res$facets[[1]]$counts, function(x) x$name)))
# 		if (length(res$facets[[1]]$counts) < 500000) break else offset <- offset + 500000
# 	}
# }
# resTestContinent <- list()
# for (i in seq_along(continents)){
# 	resTestContinent[[i]] <- NA
# 	offset <- 0
# 	repeat {
# 		res <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=1470&continent=", continents[i], "&limit=0&facet=speciesKey&facetLimit=500000&facetOffset=", offset))
# 		resTestContinent[[i]] <- c(resTestContinent[[i]], as.numeric(sapply(res$facets[[1]]$counts, function(x) x$name)))
# 		if (length(res$facets[[1]]$counts) < 500000) break else offset <- offset + 500000
# 	}
# }
# sort(setdiff(unique(unlist(resTestCountry)),unique(unlist(resTestContinent))))
# sort(setdiff(unique(unlist(resTestContinent)),unique(unlist(resTestCountry))))
# more can be gained by using country data. apparently, there is more information on countries than on
# continents. however, i will not run this again for now because we also lose occurrences. additionally, we
# would have to assign continents to countries and calculate everything for 253 instead of 7 areas.
# for vascular plants, we would win 7314 species and lose 133. for haptophyta, we would win 133 and lose 2.
# for beetles, we would gain 8569 and lose 315.

## check effect of being marine with the blue whale as an example
# res <- c()
# for (i in seq_along(continents)){
# 	res <- c(res,fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=2440735&continent=", continents[i], "&limit=0"))$count)
# }
## found per continent
# sum(res)
## found in the whole world, including museums with unknown provenance
# resWorld <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=2440735&limit=0"))$count
# resWorld
## found in the whole world without coordinates
# resNoCoords <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=2440735&limit=0&hasCoordinate=FALSE"))$count
## subtract no coords from world
# resWorld - resNoCoords
# resWaterbody <- fromJSON(paste0("https://api.gbif.org/v1/occurrence/search?taxon_key=2440735&limit=0&facet=waterBody"))
# resWaterbody
# it is possible to retrieve waterbody information, however, for our approach, this is not helpful, as not all species occur in water
# and we cannot use nested facets.

# calculate numbers for Crustaceae and Myriapoda (warning: has several NA values in it)
for (i in seq_along(continents)) {
	resSpecs[[which(names(groupsGBIFTaxonKeys) == "Crustacea")]][[i]] <- list()
	resSpecs[[which(names(groupsGBIFTaxonKeys) == "Myriapoda")]][[i]] <- list()
	for (j in 1:2) {
		resSpecs[[which(names(groupsGBIFTaxonKeys) == "Crustacea")]][[i]][[j]] <-
			unique(unlist(sapply(which(names(groupsGBIFTaxonKeys) %in% crustaceans), function(x) {
				resSpecs[[x]][[i]][[j]]
			})))
		resSpecs[[which(names(groupsGBIFTaxonKeys) == "Myriapoda")]][[i]][[j]] <-
			unique(unlist(sapply(which(names(groupsGBIFTaxonKeys) %in% myriapods), function(x) {
				resSpecs[[x]][[i]][[j]]
			})))
	}
}
resSpecs[[which(names(groupsGBIFTaxonKeys) == "Crustacea")]][[8]] <-
	unique(unlist(sapply(which(names(groupsGBIFTaxonKeys) %in% crustaceans), function(x) {
		resSpecs[[x]][[8]]
	})))
resSpecs[[which(names(groupsGBIFTaxonKeys) == "Myriapoda")]][[8]] <-
	unique(unlist(sapply(which(names(groupsGBIFTaxonKeys) %in% myriapods), function(x) {
		resSpecs[[x]][[8]]
	})))

# create summary table
resSpecsTable <- data.table(GBIFTaxonKey = groupsGBIFTaxonKeys, name = names(groupsGBIFTaxonKeys))
resSpecsTable[, c(
	"europe_north_america", "africa_antarctica_asia_oceania_south_america", "world",
	"europe_north_america_coords", "africa_antarctica_asia_oceania_south_america_coords", "world_coords"
) := NA_real_]
for (i in seq_len(nrow(resSpecsTable))) {
	en <- enc <- aaaos <- aaaosc <- w <- wc <- NA
	for (j in seq_along(continents)) {
		if (continents[j] %in% c("europe", "north_america")) {
			en <- union(en, unlist(resSpecs[[i]][[j]][[1]]))
			enc <- union(enc, unlist(resSpecs[[i]][[j]][[2]]))
		} else {
			aaaos <- union(aaaos, unlist(resSpecs[[i]][[j]][[1]]))
			aaaosc <- union(aaaosc, unlist(resSpecs[[i]][[j]][[2]]))
		}
		w <- union(w, unlist(resSpecs[[i]][[j]][[1]]))
		wc <- union(wc, unlist(resSpecs[[i]][[j]][[2]]))
	}
	resSpecsTable[i, europe_north_america := length(en) - 1]
	resSpecsTable[i, europe_north_america_coords := length(enc) - 1]
	resSpecsTable[i, africa_antarctica_asia_oceania_south_america := length(aaaos) - 1]
	resSpecsTable[i, africa_antarctica_asia_oceania_south_america_coords := length(aaaosc) - 1]
	resSpecsTable[i, world := length(w) - 1]
	resSpecsTable[i, world_coords := length(wc) - 1]
	resSpecsTable[i, world_coords_land_and_sea := length(resSpecs[[i]][[8]]) - 1]
	resSpecsTable[i, seaborne := round(length(setdiff(resSpecs[[i]][[8]], wc)) / world_coords, 1)]
}
# remove myriapod and crustaceen sub-groups
resSpecsTable <- resSpecsTable[!name %in% c(myriapods, crustaceans)]

# calculate ratio
resSpecsTable[, ratio := round(europe_north_america / world, 2)]
resSpecsTable[, ratio_coords := round(europe_north_america_coords / world_coords, 2)]

# check taxa with large values for seaborne
resSpecsTable[seaborne > 0.1, c("name", "seaborne", "ratio", "ratio_coords")]

# reorder table
resSpecsTable[, oriOrder := sapply(resSpecsTable$name, function(x) which(colnames(groupsData) == x)) - 1]
setorder(resSpecsTable, "oriOrder")

# check whether total numbers match LifeGate
# could be higher, if taxa occur in more than one continent, and lower,
# if taxa have no occurrence data
plot(resSpecsTable$world ~ colSums(groupsData[, -"year"]), xlab = "LifeGate total species numbers", ylab = "GBIF total species numbers", col = "white")
names(groupsGBIFTaxonKeys)
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
	namesDict[, searchName := gsub(substAcc$UTF8[i], substAcc$ASCII[i], searchName, fixed = TRUE)]
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
	namesDict[, searchName := gsub(substAcc$UTF8[i], substAcc$ASCII[i], searchName, fixed = TRUE)]
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
	nameShort = shortNames, minLength = NA_real_, meanLength = NA_real_, maxLength = NA_real_, nLength = NA_real_,
	minMass = NA_real_, meanMass = NA_real_, maxMass = NA_real_, nMass = NA_real_
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
if ("resBiL.RData" %in% list.files()) load("resBiL.RData")
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
	save(resAll, file = "resBiL.RData")
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
# do taxonomic classification using term variable, i.e. scientific names
# manual selection is necessary in some cases due to ambiguity of scientific names
# provided
resVars <- c("canonicalName", "authorship", "family", "order", "class", "kingdom", "phylum", "rank")
resTax <- data.table(term = names$term)
resTax[, (resVars) := character()]

i <- 2
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
	if (is.list(res) || (!is.null(dim(res)) && nrow(res) > 0)) {
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
				gsub(substAcc$UTF8[k], substAcc$ASCII[k], vernacularNames, fixed = TRUE)]
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

# create special column for ambiguous cases that contain the group the name
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

rD <- rsDriver(browser = "firefox", verbose = FALSE, chromever = NULL, phantomver = NULL)
remDr <- rD[["client"]]

res <- data.table(name = groupsMeta$name, gHits = NA_real_, wLines = NA_real_, wViews = NA_real_)
i <- 1
for (i in seq_len(nrow(groupsMeta))) {
	remDr$navigate(paste0("https://www.google.com/search?client=firefox-b-d&q=", groupsMeta$name[i]))
	gHit <- NULL
	while (length(gHit) < 1) {
		try(gHit <- remDr$findElement(using = "id", value = "result-stats"), silent = TRUE)
		Sys.sleep(0.5)
	}
	gHit <- as.numeric(gsub("\\D", "", sub(" Ergebnisse.*", "", gHit$getElementAttribute("innerHTML")[[1]])))
	res[i, gHits := gHit]
	remDr$navigate(paste0("https://en.wikipedia.org/wiki/", groupsMeta$name[i]))
	wLine <- remDr$findElement(using = "css", value = "html")
	wLine <- wLine$getElementAttribute("innerHTML")[[1]]
	wLine <- length(strsplit(wLine, split = "\n")[[1]])
	res[i, wLines := wLine]
	remDr$navigate(paste0(
		"https://pageviews.wmcloud.org/?project=en.wikipedia.org&platform=all-access&agent=user&redirects=0&range=all-time&pages=",
		groupsMeta$name[i]
	))
	wView <- NULL
	while (length(wView) < 1) {
		wView <- remDr$findElements(using = "class", value = "linear-legend--counts")
		Sys.sleep(0.5)
	}
	wView <- wView[[1]]$getElementAttribute("innerHTML")[[1]]
	wView <- strsplit(wView, split = "\n")[[1]]
	wView <- as.numeric(gsub("\\.", "", wView[4]))
	res[i, wViews := wView]
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
library(lavaan) # help functions for Bayesian SEM
library(blavaan) # Bayesian SEM
library(RColorBrewer) # color palettes

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
# all(size$name == sizeB$name)
par(mfrow = c(1, 2))
plot(sizeB$`minLength` ~ I(10^size$`min_body_length` / 1000), ylab = "size Brose", xlab = "size Wikipedia")
abline(0, 1)
plot(sizeB$`maxLength` ~ I(10^size$`max_body_length` / 1000), ylab = "size Brose", xlab = "size Wikipedia")
abline(0, 1)
rm(sizeB)
colnames(size) <- gsub("\\([^\\(\\)]+\\)", "", colnames(size))
colnames(size) <- gsub("^_|_$", "", gsub("_{2,}", "_", colnames(size)))
# use data from Wikipedia because it is more complete and not (as) prone to sampling bias
# as Brose's data

# compare public interest data
# str(interestOcc)
# str(interestTax)
litOccs <- colSums(interestOcc[, -"year"], na.rm = TRUE)
litOccs <- data.table(sum = litOccs, name = interestTax$group)
litOccs <- tapply(litOccs$sum, litOccs$name, sum)
litOccs <- litOccs[names(litOccs) != ""]
interest[, litOcc := 0]
interest[sapply(names(litOccs), function(x) which(interest$name == x)), litOcc := litOccs]
plot(interest[, -"name"])
plot(interest$gHits, interest$litOcc, ylab = "Google Hits", xlab = "Public interest")
plot(interest$gHits, log(interest$litOcc), ylab = "log(Google Hits)", xlab = "Public interest")
interest[, wLines := NULL]
# Google hits and log(literature occurrences) show similar patterns,
# while wikipedia lines do not

# check author numbers from LifeGate
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
#
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

# create a joint data.table
dat <- data.table(
	name = groupsResponses$name,
	C = groupsTotal,
	A = authorMeans,
	AVar = authorVar,
	ACv = authorCv,
	L = interest$litOcc,
	LG = interest$gHits,
	B = size$median_body_length,
	S = size$`not_soil-dwelling_or_endoparasitic_-_soil-dwelling_or_endoparasitic`,
	AQ = size$`terrestrial_-_aquatic`,
	O = regions$ratio,
	EF = groupsResponses$estFutureDesc,
	TT = groupsResponses$tenPercDesc,
	MP = groupsResponses$maxDescPace,
	DV = groupsResponses$descVar
)

# use relative description pace (as maxDescPace was multiplied by groupsTotal on creation)
dat[, MP := MP / groupsTotal]

# It might be more interesting to investigate the relative description pace
# instead of the absolute description pace. The absolute one is just related
# to the number of authors, therefore the change.

# fwrite(dat, file = "groupsVariables.txt")

# compare mean author numbers and total number of species
# pdf("mean authors per year vs total description number.pdf", height=8.3,width=11.7)
# par(mfrow = c(1, 2))
# plot(NULL,
# 	xlim = range(dat$C), ylim = range(dat$A),
# 	xlab = "current descriptions", ylab = "mean authors w descriptions per year"
# )
# text(dat$C, dat$A, dat$name,
# 	cex = .75,
# 	col = "red"
# )
# abline(v = 20000, h = 25, col = "grey")
# abline(lm(dat$A ~ dat$C + 0), col = "grey")
## zoom into small groups
# plot(NULL,
# 	xlim = c(0, 20000), ylim = c(0, 25), xlab = "current descriptions",
# 	ylab = "mean authors w descriptions per year"
# )
# text(dat$C, dat$A, dat$name,
# 	cex = .75,
# 	col = "red"
# )
# abline(lm(dat$A ~ dat$C + 0), col = "grey")
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
# - mean authors per year
# - google hits / literature occurrences
# - body length
# - endoparasitic/soil-dwelling or not
# - aquatic or not
# - occurrences inside Europe + North America / occurrences worldwide

# Hypotheses
#
# estimated future descriptions ~ ??
# reason -> ten percent of what we have now and maximum pace are no deterministic predictors,
# however, calculation-wise, maximum pace will be strongly linked
#
summary(lm(EF ~ MP + TT, data = dat)) # mostly MP

# ten percent of current descriptions ~ author number + body length + occurrences +
# endoparasitic/soil-dwelling or not + aquatic or not
# reason -> depends on technical limitations, manpower, hardness to find species
#
summary(lm(TT ~ A + B + O + S + AQ, data = dat)) # A and B

# maximum description pace ~ author number + endoparasitic/soil-dwelling or not + aquatic or not
# reason -> body length + occurrences are no longer a limitation these days
#
summary(lm(MP ~ A + B + O + S + AQ + TT, data = dat)) # B and TT

# description variance ~ author number + public interest + occurrence
# reason -> with few authors, individual authors have more impact on the description pace, public interest governs money flow,
# occurrence in Europe reduces variability
#
summary(lm(DV ~ A + L + O, data = dat)) # nothing

# author number ~ public interest
# reason -> more public interest should lead to more funding and scientists working on the topic
#
summary(lm(A ~ L, data = dat)) # L

# pdf("variable pairs.pdf",width=20,height=9)
# par(mfrow = c(3, 6))
# par(mar = c(4, 4, 2, 3))
# for (i in seq_len(ncol(dat))) {
# 	if (is.numeric(dat[[i]])) {
# 		for (j in seq_len(ncol(dat))) {
# 			if (is.numeric(dat[[j]]) && i != j) {
# 				plot(NULL, xlim = range(dat[[j]]), ylim = range(dat[[i]]), xlab = colnames(dat)[j], ylab = colnames(dat)[i])
# 				text(dat[[j]], dat[[i]], labels = dat$name, cex = .75, col = groupsMeta$color)
# 				points(dat[[j]], dat[[i]], cex = .75)
# 			}
# 		}
# 	}
# }
# dev.off()

# Concept figure
#
# See manuscripts/taxon description dates/path diagram.pptx

# create scaled dataset
dats <- copy(dat)
rmCols <- c("name", "AVar", "ACv", "LG")
dats[, (rmCols) := NULL]
dats[, (colnames(dats)) := lapply(.SD, scale)]
round(cor(dats), 2)

# rethinking implementation (no SEM)

rmCols <- colnames(dats)[grepl("DVVar", colnames(dats))]
if (length(rmCols) > 0) dats[, (rmCols) := NULL]

# reg1 <- ulam(
# 	alist(	# relationships between variables
# 		L ~ normal(alpha, rho),
# 		alpha <- a * C,
# 		A ~ normal(beta, sigma),
# 		beta <-  b * L + c * C,
# 		DV ~ normal(gamma, tau),
# 		tau <-  d * exp(-A) + e * exp(-L),
# 		TT ~ normal(delta, ypsilon),
# 		delta <- f * B + g * A + h * O + ii * S + j * AQ,
# 		MP ~ normal(epsilon, phi),
# 		epsilon <- k * B + l * A + m * O + n * S + o * AQ + p * TT,
# 		EF ~ normal(zeta, chi),
# 		zeta <- q * TT + r * MP,
# 		# regression priors
# 		c(a, b, c, f, g, h, ii, j, k, l, m, n, o, p, q, r) ~ dnorm(0,1),
# 		c(d, e) ~ dexp(0.1),
# 		# variance priors
# 		c(rho, sigma, gamma, ypsilon, phi, chi) ~ dexp(1)
# 	),
# 	data = dats, chains = 4, log_lik = TRUE
# )
#
# (res1 <- precis(reg1,depth=2))
# mean   sd  5.5% 94.5% rhat ess_bulk
# r        0.70 0.11  0.51  0.88    1  2247.52
# q        0.04 0.12 -0.15  0.22    1  2331.34
# p        0.79 0.16  0.52  1.05    1  2047.67
# o        0.17 0.16 -0.09  0.41    1  2042.28
# n       -0.20 0.15 -0.45  0.04    1  2365.81
# m        0.08 0.21 -0.26  0.42    1  1424.54
# l        0.10 0.15 -0.14  0.32    1  2298.79
# k        0.50 0.18  0.21  0.78    1  2121.51
# j       -0.13 0.15 -0.36  0.11    1  1657.26
# ii       0.22 0.14  0.00  0.44    1  1878.45
# h        0.03 0.20 -0.29  0.35    1  1450.92
# g       -0.35 0.13 -0.55 -0.15    1  2737.67
# f       -0.36 0.16 -0.62 -0.10    1  2094.20
# c        0.76 0.05  0.68  0.85    1  2723.81
# b        0.34 0.05  0.26  0.42    1  2394.90
# a        0.37 0.14  0.15  0.60    1  2922.39
# e        4.84 0.58  4.00  5.81    1  3500.12
# d        0.80 0.26  0.48  1.27    1  3908.61
# chi      0.71 0.08  0.60  0.84    1  3427.81
# phi      0.82 0.09  0.68  0.97    1  2514.57
# ypsilon  0.79 0.09  0.66  0.95    1  2660.67
# gamma    0.02 0.02  0.00  0.05    1  1899.20
# sigma    0.33 0.04  0.28  0.40    1  3092.41
# rho      0.94 0.10  0.80  1.11    1  3360.44

# len1 <- length(res1@.Data[[1]])
# par(mfrow=c(1,1))
# plot(NULL,xlim=c(-1.5,1.5),ylim=c(0,8),xlab="",ylab="", main = "Model 1 coefficients")
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

# Blavaan implementation

# Blavaan cannot model the effect of a variable on the variance of another.
# Therefore, a workaround would be to estimate the variance and regress on
# the variance instead of the variable itself

# test estimation of the variance

# # first variable, just a linear sequence
# v1 <- seq(1,100,l=1000)
# # second variable, random number around zero with sd defined by first variable
# v2 <- rnorm(1000,0, sd = 1.5 * v1)
# plot(v2 ~ v1)
# # create intervals
# ints <- 20
# sdInts <- length(v2)/ints # number of samples per interval
# sdVals <- rep(NA,ints) # empty sd results
# # calculate variance per interval
# for (i in seq_along(sdVals)) sdVals[i] <- sd(v2[seq_len(sdInts) + sdInts * (i - 1)])
# sdPrds <- rep(NA,ints)
# # calculate mean x value per interval
# for (i in seq_along(sdPrds)) sdPrds[i] <- mean(v1[seq_len(sdInts) + sdInts * (i - 1)])
# dat <- data.table(sdPrds,sdVals)
# plot(sdVals~sdPrds,data=dat)
# # try to retrieve relationship
# model <- '
# 		# relationships between measured variables
# 		sdVals ~ sdPrds
# 		'
# # fit the model to the data
# fit <- bsem(model, data = dat, sample=2000)
# summary(fit)
# # This approach works, but there is of course some additional uncertainty associated,
# # especially considering that the real data is very sparse

# do the estimation for the real data

# I want to predict description variability with two or three predictors,
# which are
# - A (the number of authors)
# - O (the occurrence of species in Europe and North America)
# - L (public interest)
# I need to calculate mean variability of description variability for each combination
# of variables or for certain intervals then, I can use this variable to regress it
# A, O, or L, respectivly

## get means within intervals
# rmCols <- colnames(dats)[grepl("DVVar",colnames(dats))]
# if (length(rmCols) > 0) dats[, (rmCols) := NULL]
# newCols <- c("DVVarA","DVVarO","DVVarL","DVVarAO","DVVarAL","DVVarOL","DVVarAOL")
# dats[, (newCols) := numeric()]
# range(dats$A)
# range(dats$O)
# range(dats$L)
# ints <- 10
# seqi <- seq(-1, 5.5, l = ints)
# seqj <- seq(-2, 2.5, l = ints)
# seqk <- seq(-0.5,5.5, l = ints)
# k <- 1
# for (i in seq_len(ints-1)){
# 	dats[A > seqi[i] & A < seqi[i+1], DVVarA := sd(DV, na.rm = TRUE)]
# 	for (j in seq_len(ints-1)){
# 		dats[O > seqj[j] & O < seqj[j + 1], DVVarO := sd(DV, na.rm = TRUE)]
# 		dats[A > seqi[i] & A < seqi[i + 1] & O > seqj[j] & O < seqj[j + 1], DVVarAO := sd(DV, na.rm = TRUE)]
# 		for (k in seq_len(ints-1)){
# 			dats[L > seqk[k] & L < seqk[k + 1], DVVarL := sd(DV, na.rm = TRUE)]
# 			dats[A > seqi[i] & A < seqi[i + 1] & O > seqj[j] & O < seqj[j + 1] & L > seqk[k] & L < seqk[k + 1], DVVarAOL := sd(DV, na.rm = TRUE)]
# 			dats[A > seqi[i] & A < seqi[i + 1] & L > seqk[k] & L < seqk[k + 1], DVVarAL := sd(DV, na.rm = TRUE)]
# 			dats[L > seqk[k] & L < seqk[k + 1] & O > seqj[j] & O < seqj[j + 1], DVVarOL := sd(DV, na.rm = TRUE)]
# 		}
# 	}
# }
# dats[, (colnames(dats)) := lapply(.SD, scale)]

## data imputation for missing values
# for (i in which(is.na(dats$DVVar))) {
# 	# calculate distance between points
# 	dists <- sqrt((dats$A[i] - dats$A[-i])^2 + (dats$O[i] - dats$O[-i])^2)
# 	# calculate meann weighted by distance
# 	dats$DVVar[i] <- sum(dats$DVVar[-i] * dists,na.rm=TRUE) / sum(dists)
# }

model1 <- "
		# relationships between measured variables
		L ~ a * C
		A ~ b * L + c * C
		DV ~ d * B + e * A + f * O + g * S + h * AQ
		TT ~ ii * B + j * A + k * O + l * S + m * AQ + n * DV
		# MP ~ o * B + p * A + q * O + r * S + s * AQ + t * DV
		EF ~ u * B + v * A + w * O + x * S + y * AQ + z * DV
		"
fit1 <- bsem(model1, data = dats, sample = 2000)
summary(fit1)

# Statistic                                 MargLogLik         PPP
# Value                                       -341.229       0.500
#
# Parameter Estimates:
#
#
# 		Regressions:
# 		Estimate  Post.SD pi.lower pi.upper     Rhat    Prior
# L ~
# 		C          (a)    0.378    0.139    0.100    0.656    1.000    normal(0,10)
# A ~
# 		L          (b)    0.340    0.054    0.235    0.447    1.000    normal(0,10)
# C          (c)    0.765    0.054    0.658    0.871    1.000    normal(0,10)
# DV ~
# 		B          (d)   -0.059    0.213   -0.483    0.357    1.000    normal(0,10)
# A          (e)   -0.122    0.168   -0.461    0.206    1.000    normal(0,10)
# O          (f)   -0.496    0.269   -1.030    0.032    1.000    normal(0,10)
# S          (g)    0.115    0.188   -0.260    0.490    1.000    normal(0,10)
# AQ         (h)    0.285    0.199   -0.111    0.685    1.000    normal(0,10)
# TT ~
# 		B         (ii)   -0.400    0.146   -0.684   -0.121    1.000    normal(0,10)
# A          (j)   -0.404    0.115   -0.628   -0.172    1.000    normal(0,10)
# O          (k)   -0.197    0.194   -0.586    0.176    1.000    normal(0,10)
# S          (l)    0.274    0.125    0.027    0.520    1.000    normal(0,10)
# AQ         (m)   -0.002    0.141   -0.277    0.275    1.001    normal(0,10)
# DV         (n)   -0.442    0.109   -0.657   -0.227    1.001    normal(0,10)
# EF ~
# 		B          (o)    0.123    0.221   -0.316    0.555    1.000    normal(0,10)
# A          (p)   -0.246    0.168   -0.567    0.090    1.000    normal(0,10)
# O          (q)   -0.012    0.294   -0.578    0.576    1.000    normal(0,10)
# S          (r)   -0.128    0.192   -0.504    0.247    1.001    normal(0,10)
# AQ         (s)   -0.060    0.210   -0.464    0.352    1.001    normal(0,10)
# DV         (t)   -0.339    0.163   -0.669   -0.015    1.000    normal(0,10)
#
# Covariances:
# 		Estimate  Post.SD pi.lower pi.upper     Rhat    Prior
# .TT ~~
# 		.EF                0.251    0.120    0.041    0.515    1.001       beta(1,1)

# fitMeasures(fit1)
# plot(fit1,plot.type = "trace")
# good sampling and fit

# get samples from the posterior
postSamples1 <- standardizedPosterior(fit1, type = "std.lv")
len1 <- sum(fit1@ParTable$op == "~")

# pdf("Model 1 coefficients.pdf", width=11.7,height=8.3)
par(mfrow = c(1, 1))
plot(NULL, xlim = c(-1.5, 1.5), ylim = c(0, 7.75), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5)
xseq <- seq(-1.5, 1.5, l = 500)
for (i in seq_len(len1)) {
	# normal distribution values
	yseq <- dnorm(xseq, fit1@ParTable$est[i], fit1@ParTable$se[i])
	# density of sampled values
	dens <- density(postSamples1[, i])
	# plot sampled values
	polygon(c(dens$x, rev(dens$x)), c(rep(0, length(dens$x)), rev(dens$y)), col = paste0(rep(brewer.pal(12, "Paired"), 2)[i], "55"), border = "dark grey")
	# plot normal distribution
	lines(xseq, yseq)
	# show variable name
	text(fit1@ParTable$est[i], 1 / (fit1@ParTable$se[i] * sqrt(2 * pi)) * 1.05, fit1@ParTable$label[i], font = 2)
}
# dev.off()

# remove those coefficients that have zero within 90% PI
testPI <- sapply(seq_len(sum(fit1@ParTable$label != "")), function(x) qnorm(c(.1, .9), fit1@ParTable$est[x], fit1@ParTable$se[x]))
colnames(testPI) <- fit1@ParTable$label[fit1@ParTable$label != ""]
testPI # those that include zero in this interval need to be removed
# d, e, g, k, m, o, q, r, s, u, w, x, y

model2 <- "
		# relationships between measured variables
		L ~ a * C
		A ~  b * L + c * C
		DV ~  e * A + f * O + h * AQ
		TT ~ ii * B + j * A + l * S + n * DV
		MP ~  p * A + t * DV
		EF ~  v * A + z * DV
		"
fit2 <- bsem(model2, data = dats)
summary(fit2)

# Statistic                                 MargLogLik         PPP
# Value                                       -370.003       0.635
#
# Parameter Estimates:
#
#
# 		Regressions:
# 		Estimate  Post.SD pi.lower pi.upper     Rhat    Prior
# L ~
# 		C          (a)    0.382    0.138    0.116    0.656    1.000    normal(0,10)
# A ~
# 		L          (b)    0.339    0.054    0.232    0.443    0.999    normal(0,10)
# C          (c)    0.763    0.054    0.658    0.872    0.999    normal(0,10)
# DV ~
# 		A          (e)   -0.135    0.155   -0.455    0.158    1.000    normal(0,10)
# O          (f)   -0.375    0.178   -0.725   -0.019    0.999    normal(0,10)
# AQ         (h)    0.228    0.178   -0.118    0.567    0.999    normal(0,10)
# TT ~
# 		B         (ii)   -0.364    0.106   -0.571   -0.160    1.000    normal(0,10)
# A          (j)   -0.363    0.108   -0.573   -0.145    0.999    normal(0,10)
# S          (l)    0.229    0.094    0.043    0.409    1.000    normal(0,10)
# DV         (n)   -0.407    0.108   -0.615   -0.196    0.999    normal(0,10)
# MP ~
# 		A          (p)   -0.200    0.138   -0.474    0.072    1.000    normal(0,10)
# DV         (t)   -0.449    0.139   -0.724   -0.173    1.000    normal(0,10)
# EF ~
# 		A          (v)   -0.167    0.148   -0.456    0.129    1.000    normal(0,10)
# DV         (z)   -0.306    0.146   -0.598   -0.025    1.000    normal(0,10)
#
# Covariances:
# 		Estimate  Post.SD pi.lower pi.upper     Rhat    Prior
# .TT ~~
# 		.MP                0.267    0.109    0.084    0.519    1.001     lkj_corr(1)
# .EF                0.233    0.111    0.044    0.471    1.000     lkj_corr(1)
# .MP ~~
# 		.EF                0.580    0.165    0.323    0.964    0.999     lkj_corr(1)

# fitMeasures(fit2)
# plot(fit2,plot.type = "trace")
# good overall fit

# extract results
res1 <- as.data.table(summary(fit1))
fwrite(res1, file = "SEM results_fit1.txt")
res2 <- as.data.table(summary(fit2))
fwrite(res2, file = "SEM results_fit2.txt")

# Entry point ###################################

# load in libraries
library(data.table) # handle large datasets
library(openxlsx) # handle excel files
library(rethinking) # Bayesian inference
library(lavaan) # help functions for Bayesian SEM
library(blavaan) # Bayesian SEM
library(RColorBrewer) # color palettes

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# save.image("temp.RData")
# load("temp.RData")

#################################################

# 10 Compare species numbers between LifeGate, GBIF, and CoL#######################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets
library(RSelenium) # for web scraping
library(RJSONIO)

# clean workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# read in data
groupsData <- fread("groupsData.txt")
groupsTotal <- colSums(groupsData[, -"year"])

## close all open ports
# try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)
#
# rD <- rsDriver(browser = "firefox", verbose = FALSE, chromever = NULL, phantomver = NULL)
# remDr <- rD[["client"]]
#
# get CoL species numbers
# colNums <- rep(NA,length(groupsTotal))
# for (i in seq_along(groupsTotal)){
# 	if (is.na(colNums[i])){
# 		remDr$navigate("https://www.catalogueoflife.org/")
# 		searchField <- remDr$findElement(using = "id", value = "rc_select_0")
# 		searchField$sendKeysToElement(list(names(groupsTotal)[i]))
# 		Sys.sleep(1)
# 		searchField <- remDr$findElements(using = "class",  value ="ant-select-item")
# 		if (length(searchField) == 0){
# 			print(paste0(names(groupsTotal)[i]," not found."))
# 		} else {
# 			for (j in seq_along(searchField)){
# 				dat <- searchField[[j]]$getElementAttribute("outerHTML")[[1]]
# 				if (grepl(paste0(">",names(groupsTotal)[i]), dat)){
# 					searchField <- searchField[[j]]
# 					break
# 				}
# 			}
# 			if (length(searchField) > 1){
# 				print(paste0("Problem encountered with ",names(groupsTotal)[i],"."))
# 			} else {
# 				searchField$clickElement()
# 				taxonKey <- sub(".*=","",remDr$getCurrentUrl()[[1]])
# 				remDr$navigate(paste0("https://www.catalogueoflife.org/data/taxon/",taxonKey))
# 				sleeper <- 0
# 				while (length(remDr$findElements(using = "class",  value ="highcharts-series-group")) == 0){
# 					Sys.sleep(1)
# 					sleeper <- sleeper + 1
# 					if (sleeper > 10) break
# 				}
# 				numField <- remDr$findElements(using = "xpath", value = "//a[@href]")
# 				if (length(numField) == 0){
# 					print(paste0("Problem encountered with ",names(groupsTotal)[i],"."))
# 				} else {
# 					for (j in seq_along(numField)) {
# 						dat <- numField[[j]]$getElementAttribute("outerHTML")[[1]]
# 						if (grepl("rank=species&", dat)){
# 							dat <- strsplit(dat, split = "rank")[[1]]
## 							dat <- dat[grepl("=species&",dat)]
# 							colNums[i] <- as.numeric(sub(".*>","",sub("<.*","",dat)))
# 							break
# 						}
# 					}
# 					if (is.na(colNums[i])){
# 						print(paste0("No species number found for ",names(groupsTotal)[i],"."))
# 					} else {
# 						print(paste0("Found ",colNums[i]," species for ",names(groupsTotal)[i],"."))
# 					}
# 				}
# 			}
# 		}
# 	}
# }
# fwrite(as.list(colNums),file="CoL species numbers.txt")
colNums <- fread("CoL species numbers.txt")
colNums <- unlist(colNums)
names(colNums) <- names(groupsTotal)

# get GBIF species numbers
load("groupsGBIFTaxonKeys.RData")
groupsGBIFTaxonKeys
gbifNums <- rep(NA, length(groupsGBIFTaxonKeys))
for (i in seq_along(groupsGBIFTaxonKeys)) {
	if (!is.na(groupsGBIFTaxonKeys)[i]) {
		res <- fromJSON(paste0("https://api.gbif.org/v1/species/", groupsGBIFTaxonKeys[i], "/metrics"))
		if ("numSpecies" %in% names(res)) {
			gbifNums[i] <- res[names(res) == "numSpecies"]
		}
	}
}
gbifNums[37] <- sum(gbifNums[52:58])
gbifNums[38] <- sum(gbifNums[48:51])
gbifNums <- gbifNums[1:47]

dat <- data.table(group = names(groupsTotal), LifeGate = groupsTotal, GBIF = gbifNums, CoL = colNums)
par()$mar
pdf(paste0("Comparison of LifeGate, GBIF, and CoL species numbers_", Sys.Date(), ".pdf"), width = 8.3, height = 11.7)
par(mar = c(5.1, 8.1, 4.1, 2.1))
barplot(t(dat[, -"group"]), names.arg = dat$group, las = 2, beside = TRUE, legend.text = TRUE, horiz = TRUE, las = 1)
dev.off()

# compare insect species numbers (number taken from CoL directly)
names(groupsTotal)
sum(groupsTotal[21:36]) # 1184055 in LifeGate, 1114071 in GBIF, 994767 in CoL

# ratio
dat[, ratioGBIF := round(dat$GBIF / dat$LifeGate, 2)]
dat[, ratioCoL := round(dat$CoL / dat$LifeGate, 2)]
dat[(!is.na(ratioGBIF) & (ratioGBIF > 1.25 | ratioGBIF < 0.75)) & (!is.na(ratioCoL) & (ratioCoL > 1.25 | ratioCoL < 0.75)), ]
# group LifeGate   GBIF    CoL ratioGBIF ratioCoL
# 1:    Bryophyta    25570  13285  12233      0.52     0.48
# 2: Foraminifera    30026  51490  49654      1.71     1.65
# 3:     Mollusca    98947 182230 137768      1.84     1.39
# 4:  Lepidoptera   310708 184453 166393      0.59     0.54
# 5:     Cnidaria    12785  24359  19494      1.91     1.52
# Several taxonomic groups are much less in LifeGate (Foraminifera, Mollusca, Cnidaria).
# Several taxonomic groups are much more in LifeGate (Bryophyta, Lepidoptera).

# 10 Get Lepidoptera and Bryophyta names from GBIF and CoL#########################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets
library(RSelenium) # for web scraping

# clean workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# get author name comparison function
source(paste0(.brd, "taxonomy/taxonomy help functions.R"))

# close all open ports
try(system("taskkill /im java.exe /f", intern = FALSE, ignore.stdout = FALSE), silent = TRUE)

rD <- rsDriver(browser = "firefox", verbose = FALSE, chromever = NULL)
remDr <- rD[["client"]]

# function to extract names data from CoL
CoLGetLowerRankNames <- function(taxon, taxonID, targetRank) {
	print(taxon)
	# open search page
	remDr$navigate(paste0(
		"https://www.catalogueoflife.org/data/search?TAXON_ID=",
		taxonID, "&extinct=false&extinct=&facet=rank&facet=issue&facet=status&facet=nomStatus&facet=nameType&facet=field&facet=authorship&facet=extinct&facet=environment&limit=10&offset=0&rank=",
		targetRank, "&status=accepted"
	))
	Sys.sleep(4) # let the page build
	# extract number of names
	resNum <- remDr$findElements(using = "class", value = "ant-col-12")
	for (j in seq_along(resNum)) {
		namesNum <- resNum[[j]]$getElementAttribute("innerHTML")
		if (grepl("results: ", namesNum)) {
			namesNum <- as.numeric(gsub("\\D", "", namesNum))
		}
	}
	# calculate number of iterations for extraction
	resPages <- ceiling(namesNum / 100)
	# create data.table for data
	res <- data.table(id = seq_len(namesNum), name = NA_character_, authors = NA_character_, year = NA_real_, authorsInBrackets = NA_character_, yearInBrackets = NA_real_, taxonID = NA_character_)
	# collect species names
	options(scipen = 999) # make sure large numbers are shown as numbers
	for (j in seq_len(resPages)) {
		remDr$navigate(paste0(
			"https://www.catalogueoflife.org/data/search?TAXON_ID=",
			taxonID, "&extinct=false&extinct=&facet=rank&facet=issue&facet=status&facet=nomStatus&facet=nameType&facet=field&facet=authorship&facet=extinct&facet=environment&limit=100&offset=",
			(j - 1) * 100, "&rank=",
			targetRank, "&status=accepted"
		))
		print(paste0(j, "/", resPages, " - ", Sys.time()))
		repeat{
			Sys.sleep(1) # let the page build
			nameLinks <- remDr$findElements(using = "tag", value = "a")
			if (j < resPages) {
				if (length(nameLinks) > 300) break
			} else {
				if (length(nameLinks) > 3 * namesNum %% 100) break
			}
		}
		l <- 0
		for (k in seq_along(nameLinks)) {
			# get names only
			potName <- nameLinks[[k]]$getElementAttribute("outerHTML")[[1]]
			if (grepl("<i>", potName)) {
				l <- l + 1
				# get taxon ID
				res[(j - 1) * 100 + l, taxonID := sub("\".*", "", sub(".*taxon/", "", potName))]
				potName <- sub("</a>", "", sub(".*(?=<i>)", "", potName, perl = TRUE))
				res[(j - 1) * 100 + l, name := sub("<i>", "", sub("</i>.*", "", potName))]
				potName <- gsub("\\]|\\[", "", potName)
				if (grepl("\\(", potName)) {
					temp <- regmatches(potName, regexpr("\\(.*\\)", potName))
					temp <- gsub("\\(|\\)", "", temp)
					res[(j - 1) * 100 + l, authorsInBrackets := sub(",\\s\\d+\\s?$", "", temp)]
					res[(j - 1) * 100 + l, yearInBrackets := as.numeric(gsub("\\D", "", temp))]
					potName <- gsub("\\(.*\\)", "", potName)
				}
				res[(j - 1) * 100 + l, authors := sub(",\\s\\d+\\s?$", "", sub(".*</i>", "", potName))]
				res[(j - 1) * 100 + l, year := as.numeric(gsub("\\D", "", potName))]
			}
		}
	}
	return(res)
}
# Problem: CoL has a maximum offset of 100000, which means that not all Lepidoptera species can be retrieved
# using the species parameter. Therefore, genera have to be retrieved first.

# resBryo <- CoLGetLowerRankNames("Bryophyta","BJ5TM","species")
# fwrite(resBryo,file="Bryophytes.gz")

# resLepiGen <- CoLGetLowerRankNames("Lepidoptera genera","CB2MR","genus")
# fwrite(resLepiGen,file="Lepidoptera genera.gz")

# get Lepidoptera species from CoL
# genera <- fread("Lepidoptera genera.gz") # results from old trial (i.e. "Lepidoptera","CB2MR","species")
# missingGenera <- genera[name > "Monoctenia"]
#
# for (i in seq_len(nrow(missingGenera))){
# 	resTemp <- CoLGetLowerRankNames(missingGenera[i]$name,missingGenera[i]$taxonID,"species")
# 	fwrite(resTemp,file=paste0(missingGenera[i]$name,".gz"))
# }

# get everything from GBIF
# While the below code should work, it does not, because GBIF has an offset limit of 10000.
# However, I was able to download some data here https://www.checklistbank.org/.
# Also the CoL data can conveniently be downloaded there :P.

# usageKey <- readLines("https://api.gbif.org/v1/species/match?name=Lepidoptera")
# usageKey <- sub(".*:","",sub(",.*","",usageKey))
# namesNum <- readLines(paste0("https://api.gbif.org/v1/species/search?datasetKey=d7dddbf4-2cf0-4f39-9b2a-bb099caae36c&rank=SPECIES&higherTaxonKey=",
# 	usageKey,"&status=ACCEPTED&limit=1&offset=0"))
# namesNum <- as.numeric(sub(",.*","",sub(".*count\":","",namesNum)))
## calculate number of iterations for extraction
# resPages <- ceiling(namesNum/1000)
## create data.table for data
# resTable <- data.table(id=seq_len(namesNum),name=NA_character_,author=NA_character_,year=NA_real_,authorsInBrackets=NA_character_,yearInBrackets=NA_real_,nubKey=NA_real_)
## extract data
# for (i in seq_len(resPages)){
# 	print(i)
# 	res <- readLines(paste0("https://api.gbif.org/v1/species/search?datasetKey=d7dddbf4-2cf0-4f39-9b2a-bb099caae36c&rank=SPECIES&higherTaxonKey=",
# 		usageKey,"&status=ACCEPTED&limit=1000&offset=",
# 		(i-1)*1000))
# 	res <- strsplit(res,split="\"key\":\\d+")[[1]][-1]
# 	nubKeys <- as.numeric(gsub("\\D","",regmatches(res,regexpr("\"nubKey\":\\d+",res))))
# 	resTable[(i-1)*1000+1:1000, nubKey := nubKeys]
# 	species <- sub("\"$","",sub("\"species\":\"","",regmatches(res,regexpr("\"species\":\"[^\"]+\"",res))))
# 	resTable[(i-1)*1000+1:1000, name := species]
# 	authors <- sub("\"$","",sub("\"authorship\":\"","",regmatches(res,regexpr("\"authorship\":\"[^\"]*\"",res))))
# 	resTable[(i-1)*1000+1:1000, author := authors]
# 	years <- as.numeric(gsub("\\D","",sub("\"$","",sub("\"authorship\":\"","",regmatches(res,regexpr("\"authorship\":\"[^\"]*\"",res))))))
# 	resTable[(i-1)*1000+1:1000, year := years]
# }

# Check downloads

# Lepidoptera
LepGBIF <- fread("Lepidoptera GBIF/nameUsage.tsv")
table(LepGBIF$`col:rank`)
table(LepGBIF$`col:status`)
LepGBIF <- LepGBIF[`col:rank` == "species"]
LepGBIF[, name := `col:scientificName`]
LepGBIF[, authors := `col:authorship`]
LepGBIF[, fullName := trimws(paste(name, authors))]
LepGBIF[, match := 0]
LepGBIF[, yearMatch := numeric()]
table(LepGBIF$match)

LepCoL <- fread("Lepidoptera CoL/nameUsage.tsv")
table(LepCoL$`col:rank`)
table(LepCoL$`col:status`)
LepCoL <- LepCoL[`col:rank` == "species"]
LepCoL[, name := `col:scientificName`]
LepCoL[, authors := `col:authorship`]
LepCoL[, fullName := trimws(paste(name, authors))]
LepCoL[, match := 0]
LepCoL[, yearMatch := numeric()]
table(LepCoL$match)

LepCoL[fullName %in% LepGBIF$fullName, match := 4]
LepGBIF[fullName %in% LepCoL$fullName, match := 4]

# check differences
for (i in seq_len(nrow(LepGBIF))) {
	if (LepGBIF$match[i] < 1) {
		if (LepGBIF$name[i] %in% LepCoL$name) {
			LepGBIF[i, match := 1]
			temp <- LepCoL[name == LepGBIF$name[i]]
			tempYear1 <- min(as.numeric(regmatches(LepGBIF$authors[i], gregexpr("\\d{4}", LepGBIF$authors[i]))[[1]]))
			tempYear2 <- regmatches(temp$authors, gregexpr("\\d+", temp$authors))
			tempYear2 <- sapply(tempYear2, function(x) min(as.numeric(x)))
			author1 <- gsub("\\s*,?\\s*\\(?\\[?\\d{4}\\]?\\)?\\s*,?\\s*", "", LepGBIF$authors[i])
			author2 <- gsub("\\s*,?\\s*\\(?\\[?\\d{4}\\]?\\)?\\s*,?\\s*", "", temp$authors)
			# print(paste(author1,author2))
			resAuthors <- sapply(author2, function(x) sum(authorMatch(author1, x)))
			# print(resAuthors)
			if (author1 == "" || all(author2 == "")) {
				LepGBIF[i, match := 4]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || any(tempYear2 == tempYear1)) {
					LepGBIF[i, yearMatch := 1]
				} else {
					LepGBIF[i, yearMatch := 0]
				}
			} else if (max(resAuthors) > 4 / 3) {
				LepGBIF[i, match := 3]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || tempYear2[order(resAuthors, decreasing = TRUE)[1]] == tempYear1) {
					LepGBIF[i, yearMatch := 1]
				} else {
					LepGBIF[i, yearMatch := 0]
				}
			} else if (max(resAuthors) > 2 / 3) {
				LepGBIF[i, match := 2]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || tempYear2[order(resAuthors, decreasing = TRUE)[1]] == tempYear1) {
					LepGBIF[i, yearMatch := 1]
				} else {
					LepGBIF[i, yearMatch := 0]
				}
			} else {
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || any(tempYear2 == tempYear1)) {
					LepGBIF[i, yearMatch := 1]
				} else {
					LepGBIF[i, yearMatch := 0]
				}
			}
		}
	}
}
for (i in seq_len(nrow(LepCoL))) {
	if (LepCoL$match[i] < 1) {
		if (LepCoL$name[i] %in% LepGBIF$name) {
			LepCoL[i, match := 1]
			temp <- LepGBIF[name == LepCoL$name[i]]
			tempYear1 <- min(as.numeric(regmatches(LepCoL$authors[i], gregexpr("\\d{4}", LepCoL$authors[i]))[[1]]))
			tempYear2 <- regmatches(temp$authors, gregexpr("\\d+", temp$authors))
			tempYear2 <- sapply(tempYear2, function(x) min(as.numeric(x)))
			author1 <- gsub("\\s*,?\\s*\\(?\\[?\\d{4}\\]?\\)?\\s*,?\\s*", "", LepCoL$authors[i])
			author2 <- gsub("\\s*,?\\s*\\(?\\[?\\d{4}\\]?\\)?\\s*,?\\s*", "", temp$authors)
			# print(paste(author1,author2))
			resAuthors <- sapply(author2, function(x) sum(authorMatch(author1, x)))
			# print(resAuthors)
			if (author1 == "" || all(author2 == "")) {
				LepCoL[i, match := 4]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || any(tempYear2 == tempYear1)) {
					LepCoL[i, yearMatch := 1]
				} else {
					LepCoL[i, yearMatch := 0]
				}
			} else if (max(resAuthors, na.rm = TRUE) > 4 / 3) {
				LepCoL[i, match := 3]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || tempYear2[order(resAuthors, decreasing = TRUE)[1]] == tempYear1) {
					LepCoL[i, yearMatch := 1]
				} else {
					LepCoL[i, yearMatch := 0]
				}
			} else if (max(resAuthors, na.rm = TRUE) > 2 / 3) {
				LepCoL[i, match := 2]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || tempYear2[order(resAuthors, decreasing = TRUE)[1]] == tempYear1) {
					LepCoL[i, yearMatch := 1]
				} else {
					LepCoL[i, yearMatch := 0]
				}
			} else {
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || any(tempYear2 == tempYear1)) {
					LepCoL[i, yearMatch := 1]
				} else {
					LepCoL[i, yearMatch := 0]
				}
			}
		}
	}
}
table(LepGBIF$match)
table(LepCoL$match)
table(LepGBIF$match, LepGBIF$yearMatch)
table(LepCoL$match, LepCoL$yearMatch)

# create combined table
str(LepGBIF)
table(LepGBIF[LepGBIF$authors == ""]$match)
table(LepCoL[LepCoL$authors == ""]$match)

LepCoL[, database := "CoL"]
LepGBIF[, database := "GBIF"]
Lep <- rbind(LepGBIF, LepCoL, fill = TRUE)

Lep[, authorsNoYear := gsub("\\d{4}", "", authors)]
Lep[, authorsNoYear := gsub(", (\\[-\\])?\\)", ")", authorsNoYear)]
Lep[, authorsNoYear := sub(", \\[\\]", "", authorsNoYear)]
Lep[, authorsNoYear := gsub("\\s\\)", ")", authorsNoYear)]
Lep[, authorsNoYear := gsub("\\(\\s", "(", authorsNoYear)]
Lep[, authorsNoYear := gsub("\\(\\[", "(", authorsNoYear)]
Lep[, authorsNoYear := gsub("\\]\\)", ")", authorsNoYear)]
Lep[, authorsNoYear := gsub(", -\\)", ")", authorsNoYear)]
Lep[, authorsNoYear := gsub("\\s&\\)", ")", authorsNoYear)]
Lep[, authorsNoYear := gsub("(\\s*,\\s*)*$", "", authorsNoYear)]
Lep[, authorsNoYear := gsub("\\s{2,}", " ", authorsNoYear)]
Lep[, year := as.numeric(gsub("\\D", "", authors))]
Lep[year > 2024, year := as.numeric(substr(year, 1, 4))]
Lep <- Lep[is.na(year) | year <= 2024]

table(Lep[grepl("\\([^A-Za-z]", authorsNoYear, perl = TRUE)]$authorsNoYear)
table(Lep[grepl("\\[[^A-Za-z]", authorsNoYear, perl = TRUE)]$authorsNoYear)
table(Lep[grepl("[^A-Za-z\\.]\\)", authorsNoYear, perl = TRUE)]$authorsNoYear)
table(Lep[grepl("[^A-Za-z\\.]\\]", authorsNoYear, perl = TRUE)]$authorsNoYear)

setorder(Lep, name, authors)
colnames(Lep)
setcolorder(Lep, c(colnames(Lep)[1:8], c("authorsNoYear", "year", "fullName", "database", "match", "yearMatch")))
Lep[1:3]
fwrite(Lep, file = "Lepidoptera GBIF and CoL.csv")

# Bryophyta
BryoGBIF <- fread("Bryophyta GBIF/nameUsage.tsv")
table(BryoGBIF$`col:rank`)
table(BryoGBIF$`col:status`)
BryoGBIF <- BryoGBIF[`col:rank` == "species"]
BryoGBIF[, name := `col:scientificName`]
BryoGBIF[, authors := `col:authorship`]
BryoGBIF[, fullName := trimws(paste(name, authors))]
BryoGBIF[, match := 0]
BryoGBIF[, yearMatch := numeric()]

BryoCoL <- fread("Bryophyta CoL/nameUsage.tsv")
table(BryoCoL$`col:rank`)
table(BryoCoL$`col:status`)
BryoCoL <- BryoCoL[`col:rank` == "species"]
BryoCoL[, name := `col:scientificName`]
BryoCoL[, authors := `col:authorship`]
BryoCoL[, fullName := trimws(paste(name, authors))]
BryoCoL[, match := 0]
BryoCoL[, yearMatch := numeric()]

BryoCoL[fullName %in% BryoGBIF$fullName, match := 4]
BryoGBIF[fullName %in% BryoCoL$fullName, match := 4]
table(BryoGBIF$match)
table(BryoCoL$match)

# check differences
for (i in seq_len(nrow(BryoGBIF))) {
	if (BryoGBIF$match[i] < 1) {
		if (BryoGBIF$name[i] %in% BryoCoL$name) {
			BryoGBIF[i, match := 1]
			temp <- BryoCoL[name == BryoGBIF$name[i]]
			tempYear1 <- min(as.numeric(regmatches(BryoGBIF$authors[i], gregexpr("\\d{4}", BryoGBIF$authors[i]))[[1]]))
			tempYear2 <- regmatches(temp$authors, gregexpr("\\d+", temp$authors))
			tempYear2 <- sapply(tempYear2, function(x) min(as.numeric(x)))
			author1 <- gsub("\\s*,?\\s*\\(?\\[?\\d{4}\\]?\\)?\\s*,?\\s*", "", BryoGBIF$authors[i])
			author2 <- gsub("\\s*,?\\s*\\(?\\[?\\d{4}\\]?\\)?\\s*,?\\s*", "", temp$authors)
			# print(paste(author1,author2))
			resAuthors <- sapply(author2, function(x) sum(authorMatch(author1, x)))
			# print(resAuthors)
			if (author1 == "" || all(author2 == "")) {
				BryoGBIF[i, match := 4]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || any(tempYear2 == tempYear1)) {
					BryoGBIF[i, yearMatch := 1]
				} else {
					BryoGBIF[i, yearMatch := 0]
				}
			} else if (max(resAuthors) > 4 / 3) {
				BryoGBIF[i, match := 3]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || tempYear2[order(resAuthors, decreasing = TRUE)[1]] == tempYear1) {
					BryoGBIF[i, yearMatch := 1]
				} else {
					BryoGBIF[i, yearMatch := 0]
				}
			} else if (max(resAuthors) > 2 / 3) {
				BryoGBIF[i, match := 2]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || tempYear2[order(resAuthors, decreasing = TRUE)[1]] == tempYear1) {
					BryoGBIF[i, yearMatch := 1]
				} else {
					BryoGBIF[i, yearMatch := 0]
				}
			} else {
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || any(tempYear2 == tempYear1)) {
					BryoGBIF[i, yearMatch := 1]
				} else {
					BryoGBIF[i, yearMatch := 0]
				}
			}
		}
	}
}
for (i in seq_len(nrow(BryoCoL))) {
	if (BryoCoL$match[i] < 1) {
		if (BryoCoL$name[i] %in% BryoGBIF$name) {
			BryoCoL[i, match := 1]
			temp <- BryoGBIF[name == BryoCoL$name[i]]
			tempYear1 <- min(as.numeric(regmatches(BryoCoL$authors[i], gregexpr("\\d{4}", BryoCoL$authors[i]))[[1]]))
			tempYear2 <- regmatches(temp$authors, gregexpr("\\d+", temp$authors))
			tempYear2 <- sapply(tempYear2, function(x) min(as.numeric(x)))
			author1 <- gsub("\\s*,?\\s*\\(?\\[?\\d{4}\\]?\\)?\\s*,?\\s*", "", BryoCoL$authors[i])
			author2 <- gsub("\\s*,?\\s*\\(?\\[?\\d{4}\\]?\\)?\\s*,?\\s*", "", temp$authors)
			# print(paste(author1,author2))
			resAuthors <- sapply(author2, function(x) sum(authorMatch(author1, x)))
			# print(resAuthors)
			if (author1 == "" || all(author2 == "")) {
				BryoCoL[i, match := 4]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || any(tempYear2 == tempYear1)) {
					BryoCoL[i, yearMatch := 1]
				} else {
					BryoCoL[i, yearMatch := 0]
				}
			} else if (max(resAuthors, na.rm = TRUE) > 4 / 3) {
				BryoCoL[i, match := 3]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || tempYear2[order(resAuthors, decreasing = TRUE)[1]] == tempYear1) {
					BryoCoL[i, yearMatch := 1]
				} else {
					BryoCoL[i, yearMatch := 0]
				}
			} else if (max(resAuthors, na.rm = TRUE) > 2 / 3) {
				BryoCoL[i, match := 2]
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || tempYear2[order(resAuthors, decreasing = TRUE)[1]] == tempYear1) {
					BryoCoL[i, yearMatch := 1]
				} else {
					BryoCoL[i, yearMatch := 0]
				}
			} else {
				if (length(tempYear1) < 1 || all(sapply(tempYear2, length) < 1) || any(tempYear2 == tempYear1)) {
					BryoCoL[i, yearMatch := 1]
				} else {
					BryoCoL[i, yearMatch := 0]
				}
			}
		}
	}
}
table(BryoGBIF$match)
table(BryoCoL$match)
table(BryoGBIF$match, BryoGBIF$yearMatch)
table(BryoCoL$match, BryoCoL$yearMatch)

# create combined table
str(BryoGBIF)
table(BryoGBIF[BryoGBIF$authors == ""]$match)
table(BryoCoL[BryoCoL$authors == ""]$match)

BryoCoL[, database := "CoL"]
BryoGBIF[, database := "GBIF"]
Bryo <- rbind(BryoGBIF, BryoCoL, fill = TRUE)

Bryo[, authorsNoYear := gsub("\\d{4}", "", authors)]
Bryo[, authorsNoYear := gsub(", (\\[-\\])?\\)", ")", authorsNoYear)]
Bryo[, authorsNoYear := sub(", \\[\\]", "", authorsNoYear)]
Bryo[, authorsNoYear := gsub("\\s\\)", ")", authorsNoYear)]
Bryo[, authorsNoYear := gsub("\\(\\s", "(", authorsNoYear)]
Bryo[, authorsNoYear := gsub("\\(\\[", "(", authorsNoYear)]
Bryo[, authorsNoYear := gsub("\\]\\)", ")", authorsNoYear)]
Bryo[, authorsNoYear := gsub(", -\\)", ")", authorsNoYear)]
Bryo[, authorsNoYear := gsub("\\s&\\)", ")", authorsNoYear)]
Bryo[, authorsNoYear := gsub("(\\s*,\\s*)*$", "", authorsNoYear)]
Bryo[, authorsNoYear := gsub("\\s{2,}", " ", authorsNoYear)]
Bryo[, year := as.numeric(gsub("\\D", "", authors))]
Bryo[year > 2024, year := as.numeric(substr(year, 1, 4))]
Bryo <- Bryo[is.na(year) | year <= 2024]

table(Bryo[grepl("\\([^A-Za-z]", authorsNoYear, perl = TRUE)]$authorsNoYear)
table(Bryo[grepl("\\[[^A-Za-z]", authorsNoYear, perl = TRUE)]$authorsNoYear)
table(Bryo[grepl("[^A-Za-z\\.]\\)", authorsNoYear, perl = TRUE)]$authorsNoYear)
table(Bryo[grepl("[^A-Za-z\\.]\\]", authorsNoYear, perl = TRUE)]$authorsNoYear)

setorder(Bryo, name, authors)
colnames(Bryo)
setcolorder(Bryo, c(colnames(Bryo)[1:8], c("authorsNoYear", "year", "fullName", "database", "match", "yearMatch")))
Bryo[1:10]
fwrite(Bryo, file = "Bryophyta GBIF and CoL.csv")

# 11 Compare checked species numbers for Lepidoptera and Bryophyta#################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates/Arten pro Jahr beschrieben"))

# get data
# Lepidoptera
lep1 <- fread("Arthropoda Insecta Lepidoptera beschriebene Arten pro Jahr.txt")
lep2 <- fread("Arthropoda Insecta Lepidoptera beschriebene Arten pro Jahr_CoL_GBIF.txt")
# Bryophyta
bryo1 <- fread("Bryophyta beschriebene Arten pro Jahr.txt")
bryo2 <- fread("Bryophyta beschriebene Arten pro Jahr_CoL_GBIF.txt")

# have a look
plot(NULL, xlim = range(lep1$V1), ylim = c(0, max(lep1$V2, lep2$V2)), xlab = "year", ylab = "descriptions")
lines(lep1$V1, lep1$V2, lwd = 2)
lines(lep2$V1, lep2$V2, lwd = 2, col = "red")
plot(NULL, xlim = range(bryo1$V1), ylim = c(0, max(bryo1$V2, bryo2$V2)), xlab = "year", ylab = "descriptions")
lines(bryo1$V1, bryo1$V2, lwd = 2)
lines(bryo2$V1, bryo2$V2, lwd = 2, col = "red")

# looks good
