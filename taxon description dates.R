# TODO: Add comment
#
# Author: David Schellenberger Costa
###################################################################################################
# 1 Visualize data from LCVP and total numbers sent by Martin Freiberg
# 2 Try to reproduce the analyses from "Trends in botanical exploration in Nigeria forecast over
# 1000 yet undescribed vascular plant species", https://doi.org/10.1093/aob/mcad106
# 3 Work on Martin's data of all groups
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
wcvpn[, first_authors := sub(".*\\sex\\.?\\s+", "", first_authors)]

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
authors <- strsplit(wcvpn$first_authors, split = ",|&|ex\\.?")
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
	length(unique(unlist(strsplit(wcvpn[first_publication_date == x]$first_authors, split = ",|&|ex\\.?"))))
})]
lines(dat$year, dat$pubAuth, col = "blue", lwd = 2)
dat[, descLastYear := c(0, desc[-length(desc)])]

round(cor(dat), 2)

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
save(list = models, file = "models.RData")

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

# 3 Work on Martin's data of all groups############################################################
nextScript <- NULL

# load in libraries
library(data.table) # handle large datasets
library(rotl) # get open source phylogenies (should eventually be data from Martin)
library(ape) # plot phylogenies nicely
library(rphylopic) # get icons of taxonomic groups
library(png) # plot icons of taxonomic groups
library(rethinking) # approximate description times with functions

# clear workspace
rm(list = ls())

# set working directory
setwd(paste0(.brd, "taxon description dates"))

# read in metadata
groups <- data.table(file = list.files(path = "Arten pro Jahr beschrieben"))
groups[, name := sub("\\s+beschriebene.*", "", file)]
# create regular expression to define groups receiving same colors (monophyletic groups of groups)
# Chordata is the only group that receives a distinguished color
groupRegs <- c(
	"[^O]omycota|sporidia", # fungi
	"Rhodophyta|Chlorophyta|Bryophyta|Marchantiophyta|Tracheophyta", # plants
	"ptera|Odo|Pso", # insects
	"Platy|Gastro|Ecto|Mollus|Nemer|Anne", # Spiralia
	"Hapto|Myzo|Cilio|Foram|Oo|Ochro", # mostly unicellular, mostly uniflagellate
	"Chor" # chordates
)
groupCols <- c("brown", "darkgreen", "red", "lightblue", "green", "blue") # colors of groups of groups
groups[, color := "black"] # standard color
for (i in seq_along(groupRegs)) {
	groups[grepl(groupRegs[i], name), color := groupCols[i]]
}
groupsVec <- seq_len(nrow(groups)) # create a convenience vector to process all groups

# read in data
groupsData <- fread(paste0("Arten pro Jahr beschrieben/", groups$file[1]))
for (i in groupsVec) {
	if (i > 1) {
		groupsData <- cbind(groupsData, fread(paste0(
			"Arten pro Jahr beschrieben/",
			groups$file[i]
		))$V2)
	}
}
colnames(groupsData) <- c("year", paste0("g", seq_len(ncol(groupsData) - 1)))
# reduce data to before 2019 (as it is not complete afterwards)
groupsData <- groupsData[year < 2018]

# retrieve phylogenetic information for the taxonomic groups

# retrieve OTT (Open Tree Taxonomy) IDs
# taxa <- data.table(tnrs_match_names(sub(".*\\s", "", groups$name))) # get OTT IDs from TNRS
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
# anymore. it is therefore better to replace taxa with species belonging to them (or other lower
# taxa ranks)

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
setcolorder(groupsData, c(1, sapply(pt$tip.label, function(x) which(sub(".* ", "", groups$name) == x)) + 1))
groups <- groups[sapply(pt$tip.label, function(x) which(sub(".* ", "", groups$name) == x))]

# get icons of groups
# done manually from phylopic.org
# groups$name
# groups[1, icon :=
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
# 	save_phylopic(get_phylopic(uuid = groups$icon[i]),path=paste0("icons/",groups$name[i],".svg"))
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
# set show.node.label=FALSE to TRUE for coloring purposes
# pdf("group phylogeny.pdf")
par(parBackup)
plot.phylo(pt,
	show.node.label = FALSE, cex = 0.7, font = 1, type = "fan", x.lim = c(-xlim, xlim),
	y.lim = c(-ylim, ylim), label.offset = 0.2, edge.color = cols, edge.width = 2
)
i <- 1
tipIncrease <- 1.1 # tip length increase for figure plotting
imageSize <- 0.07 # figure size factor
for (i in seq_len(nrow(groups))) {
	segments(cos(2 * pi * (i - 1) / nrow(groups)), sin(2 * pi * (i - 1) / nrow(groups)),
		tipIncrease * cos(2 * pi * (i - 1) / nrow(groups)), tipIncrease * sin(2 * pi * (i - 1) / nrow(groups)),
		lwd = 2, col = groups$color[i]
	)
	img <- readPNG(paste0("group icons/", groups$name[i], ".png"))
	colRGB <- col2rgb(groups$color[i]) / 255
	for (j in 1:3) img[, , j] <- colRGB[[j]]
	rasterImage(img, tipIncrease * cos(2 * pi * (i - 1) / nrow(groups)) - imageSize,
		tipIncrease * sin(2 * pi * (i - 1) / nrow(groups)) - imageSize,
		tipIncrease * cos(2 * pi * (i - 1) / nrow(groups)) + imageSize,
		tipIncrease * sin(2 * pi * (i - 1) / nrow(groups)) + imageSize,
		col = groups$color[i]
	)
}
# dev.off()

# plot description history
polyYear <- c(groupsData$year, rev(groupsData$year))
zeros <- rep(0, nrow(groupsData))

# pdf("description history2.pdf",width=11.7,height=8.3)
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
	polygon(polyYear, c(log(groupsData[[i + 1]] + 1), zeros), border = NA, col = groups$col[i])
	text(1800, log(5001), labels = gsub("\\s", "\n", groups$name[i]), adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groups$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, max(log(groupsData + 1)))))
	rasterImage(img, 1760, log(2001), 1790, log(2001) + scale)
}
mtext("year", 1, 3, at = 1450)
mtext("descriptions", 2, 41.5, at = 23, xpd = TRUE)
# plot all groups in one sheet with cumulative relative values
par(oma = c(5, 5, 1, 1))
par(mar = c(0, 0, 0, 0))
par(mfrow = c(5, 10))
for (i in groupsVec) {
	plot(NULL, xlim = range(groupsData$year), ylim = c(0, 1), xaxt = "n", yaxt = "n")
	abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 3)
	abline(v = c(1800, 1900, 2000), lty = 3)
	if ((i - 1) %% 10 < 1) axis(2)
	if (i > 37) axis(1, at = c(1800, 1900, 2000))
	polygon(polyYear, c(cumsum((groupsData[[i + 1]])) / sum(groupsData[[i + 1]]), zeros),
		border = NA, col = groups$col[i]
	)
	text(1800, log(5001) / max(log(groupsData + 1)), labels = gsub("\\s", "\n", groups$name[i]), adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groups$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, 1)))
	rasterImage(img, 1760, log(2001) / max(log(groupsData + 1)), 1790, log(2001) / max(log(groupsData + 1)) + scale)
}
mtext("year", 1, 3, at = 1450)
mtext("cumulative descriptions", 2, 38.5, at = 2.7, xpd = TRUE)
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
	polygon(polyYear, c(log(groupsData[[i + 1]] + 1), zeros), border = NA, col = groups$col[i])
	text(1800, log(5001), labels = gsub("\\s", "\n", groups$name[i]), adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groups$name[i], ".png"))
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
		border = NA, col = groups$col[i]
	)
	text(1800, log(5001) / max(log(groupsData + 1)), labels = gsub("\\s", "\n", groups$name[i]), adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groups$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, 1)))
	rasterImage(img, 1760, log(2001) / max(log(groupsData + 1)), 1790, log(2001) / max(log(groupsData + 1)) + scale)
}
# dev.off()

# approximate data by function of two/three parameters
# a) technical ability to describe, cutpoint with y-axis
# b) estimated number of species in group, upper limit
# c) estimated effort/1/easiness of discovery invested in group ->
# effort or easiness of discovery cannot be distinguished from curves,
# they can only be inferred from external data, like number of researchers known to work on
# a specific group, or, as a proxy of interest, occurrence in popular literature

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
		lines(groupsData$year, cumsum((groupsData[[i + 1]])) / sum(groupsData[[i + 1]]), col = groups$col[i])
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

# implement modelling

# test Bayesian models with dinosaurs data

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
dat1 <- list(
	A = groupsData$year - min(groupsData$year),
	M = cumsum(groupsData$g47) / max(cumsum(groupsData$g47)) # normalize to [0,1]
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
# 	data = dat1, chains = 4, log_lik = TRUE
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
# 	data = dat1, chains = 4, log_lik = TRUE
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
# 	data = dat1, chains = 4, log_lik = TRUE
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
		dat1$A, dat1$M,
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

# the algorithm works. there may be some issues wit details that have to be addressed, but overall,
# this method can be implemented.
# the third model matches the data very well. it only predicts a large number of species
# to be described over the next hundred years, but this may not be a big issue at the moment
# it may well be impossible to model this without extra assumptions.

# try the different functions for prediction

drawit <- function() {
	par(parBackup)
	plot(
		dat1$A, dat1$M,
		xlab = "year", ylab = "cumulative descriptions", type = "l",
		lwd = 3, xlim = c(0, 500), ylim = c(0, 2), col = "green"
	)
	abline(h = (1:5) * 2 / 5, lty = 3)
	abline(v = (1:5) * 100, lty = 3)
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
# 	data = dat1, chains = 4, log_lik = TRUE
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
# 	data = dat1, chains = 4, log_lik = TRUE
# )
load("models2.RData")

# compare the three approximation functions
# pdf("approximation functions test.pdf",width=11.7,height=8.3)
plot(
	dat1$A, dat1$M,
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

# all three approximations work nicely with this first trial
# however, there are issues with the the values of the priors
# do not forget that exponential(1) specifies the exponential distribution,
# smaller values flatten it, larger ones make it tighter

# try the whole dataset

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
# for (i in seq_len(nrow(groups))) {
# 	print(paste0("Bertalanffy model ", i, "/", nrow(groups)))
# 	dat1 <- list(
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
# 		data = dat1, chains = 4, log_lik = TRUE
# 	)
# 	# precis(q1)
# 	# traceplot(q1,pars=c("k","a","b","sigma"),n_cols=2)
# }
# coefs[[1]] <- sapply(m$m1, coef)
## normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
# m$m2 <- list()
# for (i in seq_len(nrow(groups))) {
# 	print(paste0("Normal model ", i, "/", nrow(groups)))
# 	dat1 <- list(
# 		Y = groupsData$year - min(groupsData$year),
# 		D = cumsum(groupsData[[i + 1]]) / max(cumsum(groupsData[[i + 1]]))
# 	)
# 	dat1$Y <- standardize(dat1$Y) # get more stable results
# 	m$m2[[i]] <- ulam(
# 		alist(
# 			D ~ normal(mu, sigma),
# 			mu <- k * 1 / ((2 * 3.141593)^0.5 * a) * exp(-0.5 * ((Y - b) / a)^2),
# 			k ~ exponential(0.1),
# 			a ~ exponential(0.1),
# 			b ~ exponential(1),
# 			sigma ~ exponential(1)
# 		),
# 		data = dat1, chains = 4, log_lik = TRUE
# 	)
# 	# precis(m$m2[[i]])
# 	# traceplot(m$m2[[i]],pars=c("k","a","b","sigma"),n_cols=2)
# }
# coefs[[2]] <- sapply(m$m2, coef)
#
## Gompertz function: f(x) = k * exp(-a * exp(-b * x))
# m$m3 <- list()
# for (i in seq_len(nrow(groups))) {
# 	print(paste0("Gompertz model ", i, "/", nrow(groups)))
# 	dat1 <- list(
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
# 		data = dat1, chains = 4, log_lik = TRUE
# 	)
# 	# precis(m$m3[[i]])
# 	# traceplot(m$m3[[i]],pars=c("k","a","b","sigma"),n_cols=2)
# }
# coefs[[3]] <- sapply(m$m3, coef)
load("models3.RData")

# compare model runs
mo <- m
load("models3New.RData")

# extract coefficients
coefsm <- list()
coefsmo <- list()
for (i in seq_along(m)) {
	coefsm[[i]] <- sapply(m[[i]], coef)
	coefsmo[[i]] <- sapply(mo[[i]], coef)
}

save(list = c("coefsm", "coefsmo"), file = "compareCoef.RData")

# coefsOld <- coefs
# load("coefs models3.RData")
#
# str(coefs)
# str(coefsOld)
#
# plot(coefs[[3]][2,]~coefsOld[[3]][2,])
# abline(0,1)



# plot all groups separately with cumulative relative values and model approximations
# added information:
# - predicted asymptote/maximum
# - year of predicted asymptote/maximum
# - approximation quality

i <- 7
par(parBackup)
# pdf("description history fits2.pdf", width = 11.7, height = 8.3)
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
		border = NA, col = groups$col[i]
	)
	text(1800, log(5001) / max(log(groupsData + 1)), labels = gsub("\\s", "\n", groups$name[i]), adj = c(0, 1))
	img <- readPNG(paste0("group icons/", groups$name[i], ".png"))
	scale <- 30 / (diff(range(groupsData$year)) / diff(c(0, 1)))
	rasterImage(img, 1760, log(2001) / max(log(groupsData + 1)), 1790, log(2001) / max(log(groupsData + 1)) + scale)
	xseq <- seq(from = 0, to = max(groupsData$year) - min(groupsData$year), len = 50)
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
# dev.off()
# save(list = c("coefs", "m"), file = "models3New.RData")
# save("coefs", file = "coefs models3.RData")

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

# show parameter estimates

pdf("parameters.pdf", width = 8.3, height = 11.7)
# estimated future descriptions relative to current descriptions
par(mfrow = c(4, 3))
par(oma = c(0, 1, 4, 1))
par(mar = c(5.1, 2.1, 2.1, 1.1))
# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
plot(colSums(groupsData)[-1], coefs[[1]][1, ],
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(1, max(coefs[[1]][1, ]) + 1), lty = 3)
text(colSums(groupsData)[-1], coefs[[1]][1, ], label = 1:47, col = groups$color, cex = 0.9)
mtext("Bertalanffy", 3, 3, font = 2)
# normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
plot(colSums(groupsData)[-1], coefs[[2]][1, ] * 1 / ((2 * pi)^0.5 * coefs[[2]][2, ]),
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(1, max(coefs[[2]][1, ] * 1 / ((2 * pi)^0.5 * coefs[[2]][2, ])) + 1), lty = 3)
text(colSums(groupsData)[-1], coefs[[2]][1, ] * 1 / ((2 * pi)^0.5 * coefs[[2]][2, ]),
	label = 1:47,
	col = groups$color, cex = 0.9
)
mtext("estimated future descriptions / current descriptions", 3, 1, font = 2, cex = 0.9)
mtext("normal", 3, 3, font = 2)
# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
plot(colSums(groupsData)[-1], coefs[[3]][1, ],
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, max(coefs[[3]][1, ]) + 10, 10), lty = 3)
text(colSums(groupsData)[-1], coefs[[3]][1, ], label = 1:47, col = groups$color, cex = 0.9)
mtext("Gompertz", 3, 3, font = 2)
# differences in start of descriptions
# identify moment in time when 10% of descriptions where made
# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
vals <- log(1 - (0.1 / coefs[[1]][1, ])^(1 / coefs[[1]][2, ])) / -coefs[[1]][3, ]
plot(colSums(groupsData)[-1], vals,
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, 300, 20), lty = 3)
text(colSums(groupsData)[-1], vals, label = 1:47, col = groups$color, cex = 0.9)
# normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
vals <- 2 * coefs[[2]][3, ] - (coefs[[2]][2, ] * (-2 * log(((2 * pi)^0.5 * coefs[[2]][2, ]) /
	(10 * coefs[[2]][1, ])))^0.5 + coefs[[2]][3, ])
st <- standardize(groupsData$year - min(groupsData$year))
vals <- vals * attr(st, "scaled:scale") + attr(st, "scaled:center")
plot(colSums(groupsData)[-1], vals,
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, 300, 20), lty = 3)
text(colSums(groupsData)[-1], vals, label = 1:47, col = groups$color, cex = 0.9)
mtext("years until 10% of current descriptions made", 3, 1, font = 2, cex = 0.9)
# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
vals <- log((log(0.1 / coefs[[3]][1, ]) / -coefs[[3]][2, ])) / -coefs[[3]][3, ]
plot(colSums(groupsData)[-1], vals,
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, 300, 20), lty = 3)
text(colSums(groupsData)[-1], vals, label = 1:47, col = groups$color, cex = 0.9)
# description pace
# identify maximum description velocity and multiply by actual description number to get
# comparable values
# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
vals <- log(1 / coefs[[1]][2, ]) / -coefs[[1]][3, ]
vals[vals > 2017 - 1753] <- 2017 - 1753
vals <- coefs[[1]][1, ] * (1 - exp(-coefs[[1]][3, ] * vals))^coefs[[1]][2, ]
vals <- vals * colSums(groupsData)[-1]
plot(colSums(groupsData)[-1], vals,
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, 2, 0.2), lty = 3)
text(colSums(groupsData)[-1], vals, label = 1:47, col = groups$color, cex = 0.9)
# normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
vals <- coefs[[2]][3, ] - coefs[[2]][2, ]
st <- standardize(groupsData$year - min(groupsData$year))
vals[vals > max(st)] <- max(st)
vals <- coefs[[2]][1, ] * 1 / ((2 * pi)^0.5 * coefs[[2]][2, ]) *
	exp(-0.5 * ((vals - coefs[[2]][3, ]) / coefs[[2]][2, ])^2)
vals <- vals * colSums(groupsData)[-1]
plot(colSums(groupsData)[-1], vals,
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, 5, 0.1), lty = 3)
text(colSums(groupsData)[-1], vals, label = 1:47, col = groups$color, cex = 0.9)
mtext("maximum descriptions pace", 3, 1, font = 2, cex = 0.9)
# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
vals <- log(1 / coefs[[3]][2, ]) / -coefs[[3]][3, ]
vals[vals > 2017 - 1753] <- 2017 - 1753
vals <- coefs[[3]][1, ] * exp(-coefs[[3]][2, ] * exp(-coefs[[1]][3, ] * vals))
vals <- vals * colSums(groupsData)[-1]
plot(colSums(groupsData)[-1], vals,
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, 30, 2), lty = 3)
text(colSums(groupsData)[-1], vals, label = 1:47, col = groups$color, cex = 0.9)
# variability
# custom Bertalanffy function: f(x) = k * (1 - exp(-b * x))^a
# estimates
ests <- sapply(seq_len(ncol(coefs[[1]])), function(x) {
	calcB(coefs[[1]][1, x], coefs[[1]][2, x], coefs[[1]][3, x], groupsData$year - min(groupsData$year))
})
# measurements
meas <- groupsData[, -"year"]
meas <- apply(meas, 2, function(x) cumsum(x))
meas <- apply(meas, 2, function(x) x / max(x))
vals <- colSums((ests - meas)^2)
plot(colSums(groupsData)[-1], vals,
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, 2, 0.2), lty = 3)
text(colSums(groupsData)[-1], vals, label = 1:47, col = groups$color, cex = 0.9)
# normal distribution function: f(x) = k * 1 / ((2 * pi)^0.5 * a) * exp(-0.5 * ((x - b) / a)^2)
# estimates
ests <- sapply(seq_len(ncol(coefs[[2]])), function(x) {
	calcN(coefs[[2]][1, x], coefs[[2]][2, x], coefs[[2]][3, x], standardize(groupsData$year - min(groupsData$year)))
})
# measurements
meas <- groupsData[, -"year"]
meas <- apply(meas, 2, function(x) cumsum(x))
meas <- apply(meas, 2, function(x) x / max(x))
vals <- colSums((ests - meas)^2)
plot(colSums(groupsData)[-1], vals,
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, 2, 0.2), lty = 3)
text(colSums(groupsData)[-1], vals, label = 1:47, col = groups$color, cex = 0.9)
mtext("sum of squares (estimates - measurements)", 3, 1, font = 2, cex = 0.9)
# Gompertz function: f(x) = k * exp(-a * exp(-b * x))
# estimates
ests <- sapply(seq_len(ncol(coefs[[3]])), function(x) {
	calcG(coefs[[3]][1, x], coefs[[3]][2, x], coefs[[3]][3, x], groupsData$year - min(groupsData$year))
})
# measurements
meas <- groupsData[, -"year"]
meas <- apply(meas, 2, function(x) cumsum(x))
meas <- apply(meas, 2, function(x) x / max(x))
vals <- colSums((ests - meas)^2)
plot(colSums(groupsData)[-1], vals,
	col = "white", xlab = "current descriptions",
	ylab = "", xaxt = "n", main = ""
)
axis(1, at = seq(0, 500000, 50000), labels = sub("^\\s*", "", format(seq(0, 500000, 50000), scientific = FALSE)))
abline(v = seq(0, 500000, 50000), lty = 3)
abline(h = seq(0, 2, 0.2), lty = 3)
text(colSums(groupsData)[-1], vals, label = 1:47, col = groups$color, cex = 0.9)
dev.off()

# nicer display of the above
pdf("test2.pdf")
b1 <- barplot(coefs[[1]][1, ], col = groups$color, border = NA, horiz = TRUE, space = 0.5)
i <- 1
for (i in seq_len(ncol(coefs[[1]]))) {
	img <- readPNG(paste0("group icons/", groups$name[i], ".png"))
	scale <- max(coefs[[1]][1, ]) / max(b1)
	colRGB <- col2rgb(groups$color[i]) / 255
	for (j in 1:3) img[, , j] <- colRGB[[j]]
	rasterImage(img, coefs[[1]][1, i] + 0.5, b1[i] - 1, scale * 2 + coefs[[1]][1, i] + 0.5, b1[i] + 1, xpd = TRUE)
}
dev.off()

# plot numbers and groups
pdf("numbers and groups.pdf", width = 8.3, height = 11.7)
plot(NULL, xlim = c(0, 5), ylim = c(0, 47), xlab = "", ylab = "")
for (i in seq_len(nrow(groups))) {
	img <- readPNG(paste0("group icons/", groups$name[i], ".png"))
	scale <- max(coefs[[1]][1, ]) / max(b1)
	colRGB <- col2rgb(groups$color[i]) / 255
	for (j in 1:3) img[, , j] <- colRGB[[j]]
	rasterImage(img, 1.33, i - 1, 1.66, i)
	text(2, i - 0.5, groups$name[i], adj = 0, col = groups$color[i])
	text(1, i - 0.5, i, col = groups$color[i])
}
dev.off()

# it depends a lot on the models used
# need to make them more stable by re-running and optimizing model parameters
