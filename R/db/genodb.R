## functions for interacting with array genotype databases

MMDB.PATH = "/db/arrays/megamuga/megamuga.db"
MDADB.PATH = "/Volumes/apm_passport/db/arrays/mda/MDA.db"

## internal utilities 
.types.quotes <- function(col, ...) {
	
	col.type <- "int"
	col.quotes <- ""
	
	if (col == "name") {
		col.type <- "varchar"
		col.quotes <- "'"
	}
	
	return( list(type = col.type, quotes = col.quotes) )
	
}

.insert.samples <- function(ids, db, temp.table = "_mysamples", by = c("name","id"), ...) {
	
	require(RSQLite)
	stopifnot(!is.null(ids))
	stopifnot(inherits(db, "SQLiteConnection"))
	
	.by <- match.arg(by)
	tnq <- .types.quotes(.by)
	
	sql <- paste("CREATE TEMP TABLE", temp.table, "(", .by, tnq$type,");")
	dbGetQuery(db, sql)
	
	for (s in ids) {
		sql <- paste0("INSERT INTO ", temp.table, " (", .by, ") VALUES (", tnq$quotes, s, tnq$quotes, ");")
		dbGetQuery(db, sql)
	}
	
}

.insert.markers <- function(ids, db, temp.table = "_mymarkers", by = c("name","id"), ...) {
	
	require(RSQLite)
	stopifnot(!is.null(ids))
	stopifnot(inherits(db, "SQLiteConnection"))
	
	.by <- match.arg(by)
	tnq <- .types.quotes(.by)
	
	sql <- paste("CREATE TEMP TABLE", temp.table, "(", .by, tnq$type,");")
	dbGetQuery(db, sql)
	
	for (s in ids) {
		sql <- paste0("INSERT INTO ", temp.table, " (", .by, ") VALUES (", tnq$quotes, s, tnq$quotes, ");")
		dbGetQuery(db, sql)
	}
	
}

.chunk.query <- function(db, sql, batch.size = -1, ...) {
	
	require(RSQLite)
	stopifnot(inherits(db, "SQLiteConnection"))
	
	rez <- dbSendQuery(db, sql)
	df <- data.frame()
	while (!dbHasCompleted(rez)) {
		df <- rbind(df, fetch(rez, n = batch.size))
	}
	dbClearResult(rez)
	dbDisconnect(db)
	
	return(df)
	
}

.dbnotfound <- function() {
	stop("Specify a database in which to search ('mda' = MDA, 'mm' = MegaMuga). GigaMUGA not yet supported.")
}

## entry point for external calls
fetch.samples <- function(..., db = c("mm","mda")) {
	
	if (db == "mm")
		.fetch.samples.mm(...)
	else if (db == "mda")
		.fetch.samples.mda(...)
	else
		.dbnotfound()
	
}

fetch.intensities <- function(..., db = c("mm","mda")) {
	
	if (db == "mm")
		.fetch.intensities.mm(...)
	else if (db == "mda")
		.fetch.intensities.mda(...)
	else
		.dbnotfound()
	
}

fetch.calls <- function(..., db = c("mm","mda")) {
	
	if (db == "mm")
		## NB: for MM, intensities and calls always returned together
		.fetch.intensities.mm(...)
	else if (db == "mda")
		.fetch.calls.mda(...)
	else
		.dbnotfound()
	
}

fetch.controls <- function(..., db = c("mm","mda")) {
	
	if (db == "mm")
		## NB: for MM, intensities and calls always returned together
		.fetch.controls.mm(...)
	else
		stop("Control samples only defined for MegaMUGA so far.")
	
}

## array-specific functions
.fetch.samples.mm <- function(ids = NULL, group = NULL, db = MMDB.PATH, by = c("name","id"), exact = TRUE, strict.case = TRUE, verbose = TRUE, ...) {
	
	require(RSQLite)
	stopifnot( !all(!is.null(ids), !is.null(group)) )
	
	db <- dbConnect(SQLite(), dbname = db)
	
	cols <- paste("s", c("id", "name", "well", "batch", "sex", "flags", "timestamp"), sep = ".", collapse = ", ")
	sql <- paste("SELECT", cols)
	if (!is.null(ids)) {
		if (exact) {
			.ids <- na.omit(ids)
			.insert.samples(.ids, db, by = match.arg(by))
			sql <- paste0(sql, " FROM samples s ",
										"INNER JOIN _mysamples as sg ON s.", by, " = sg.", by)
		}
		else {
			sql <- paste0(sql, " FROM SAMPLES s ",
										"WHERE s.name LIKE '", ids[1], "'")
		}
	}
	else if (!is.null(group)) {
		sql <- paste0(sql, ", g.name as gname FROM samples s ",
									"INNER JOIN samples_groups sg ON sg.sid = s.id ",
									"INNER JOIN groups g ON sg.gid = g.id ",
									"WHERE g.name LIKE '", group[1], "'")
	}
	if (!strict.case)
		sql <- paste0(sql, " COLLATE NOCASE")
	sql <- paste0(sql, ";")
	
	if (verbose)
		cat(sql, "\n")
	
	.chunk.query(db, sql, -1)
	
}

.fetch.samples.mda <- function(ids = NULL, db = MDADB.PATH, by = c("name","id"), exact = TRUE, ...) {
	
	require(RSQLite)
	
	db <- dbConnect(SQLite(), dbname = db)
	
	cols <- paste("s", c("id", "name", "cel_file"), sep = ".", collapse = ", ")
	sql <- paste("SELECT", cols)
	if (!is.null(ids)) {
		if (exact) {
			.ids <- na.omit(ids)
			.insert.samples(.ids, db, by = match.arg(by))
			sql <- paste0(sql, " FROM samples s ",
										"INNER JOIN _mysamples as sg ON s.", by, " = sg.", by)
		}
		else {
			sql <- paste0(sql, " FROM samples s ",
										"WHERE s.name LIKE '", ids[1], "'")
		}
	}
	sql <- paste0(sql, ";")
	cat(sql, "\n")
	
	.chunk.query(db, sql, -1)
	
}

.fetch.intensities.mm <- function(ids, markers = NULL, chr = NULL, start = NULL, end = NULL, by = c("name","id"), db = MMDB.PATH, batch.size = 1000, verbose = TRUE, ...) {
	
	require(RSQLite)
	stopifnot(!is.null(ids))
	
	db <- dbConnect(SQLite(), dbname = db)
	
	.insert.samples(ids, db, by = by)
	
	sql <- paste("SELECT m.name as marker, s.id as sid, s.name as id, m.chromosome as chr, m.position as pos, g.x, g.y, g.call",
							 "FROM samples as s",
							 paste0("INNER JOIN _mysamples as sg ON s.", by, " = sg.", by),
							 "INNER JOIN genotypes as g ON g.sampleID = s.id",
							 "INNER JOIN snps as m on g.snpID = m.id", sep = "\n")
	if (!is.null(markers)) {
		.insert.markers(markers, db, by ="name")
		sql <- paste0(sql, "\nINNER JOIN _mymarkers as mym ON m.name = mym.name")
	}
	if (!is.null(chr))
		sql <- paste0(sql, "\nWHERE chr = '", gsub("chr","",chr), "'")
	if (!is.null(start))
		sql <- paste(sql, "AND\npos >=", formatC(start, format = "d"))
	if (!is.null(end))
		sql <- paste(sql, "AND\npos <=", formatC(end, format = "d"))
	sql <- paste0(sql, ";")
	
	if (verbose)
		cat(sql, "\n")
	
	.chunk.query(db, sql, batch.size)
	
}

.fetch.intensities.mda <- function(ids, markers = NULL, chr = NULL, start = NULL, end = NULL, by = c("name","id"), db = MDADB.PATH, batch.size = 1000, verbose = TRUE, ...) {
	
	require(RSQLite)
	stopifnot(!is.null(ids))
	
	db <- dbConnect(SQLite(), dbname = db)
	
	.insert.samples(ids, db, by = by)
	
	sql <- paste("SELECT s.id as sid, s.name as name,",
							 "m.snpID as marker, m.chrID as chr, m.positionBp as pos, m.positioncM as cM, m.flagged as flag, m.mm10positionBp as mm10pos,",
							 "i.average as avg, i.contrast as contr",
							 "FROM samples as s",
							 paste0("INNER JOIN _mysamples as sg ON s.", by, " = sg.", by),
							 #"INNER JOIN genotypes as g ON g.sampleID = s.id",
							 "INNER JOIN intensity as i on i.sampleID = s.id",
							 "INNER JOIN snpInfo as m on i.snpID = m.snpID", sep = "\n")
	if (!is.null(markers)) {
		.insert.markers(markers, db, by ="name")
		sql <- paste0(sql, "\nINNER JOIN _mymarkers as mym ON m.snpID = mym.name")
	}
	if (!is.null(chr))
		sql <- paste0(sql, "\nWHERE m.chrID = '", gsub("chr","",chr), "'")
	if (!is.null(start))
		sql <- paste(sql, "AND\n m.positionBp >=", formatC(start, format = "d"))
	if (!is.null(end))
		sql <- paste(sql, "AND\n m.positionBp <=", formatC(end, format = "d"))
	sql <- paste0(sql, ";")
	
	if (verbose)
		cat(sql, "\n")
	
	.chunk.query(db, sql, batch.size)
	
}

.fetch.calls.mda <- function(ids, markers = NULL, chr = NULL, start = NULL, end = NULL, by = c("name","id"), db = MDADB.PATH, batch.size = 1000, verbose = TRUE, ...) {
	
	require(RSQLite)
	stopifnot(!is.null(ids))
	
	db <- dbConnect(SQLite(), dbname = db)
	
	.insert.samples(ids, db, by = by, verbose = TRUE)
	
	sql <- paste("SELECT s.id as sid, s.name as name,",
							 "m.snpID as marker, m.chrID as chr, m.positionBp as pos, m.positioncM as cM, m.flagged as flag, m.mm10positionBp as mm10pos,",
							 "m.alleleA, m.alleleB, g.call as call",
							 "FROM samples as s",
							 paste0("INNER JOIN _mysamples as sg ON s.", by, " = sg.", by),
							 "INNER JOIN genotypes as g ON g.sampleID = s.id",
							 #"INNER JOIN intensity as i on i.snpID = m.snpID",
							 "INNER JOIN snpInfo as m on g.snpID = m.snpID", sep = "\n")
	if (!is.null(markers)) {
		.insert.markers(markers, db, by ="name")
		sql <- paste0(sql, "\nINNER JOIN _mymarkers as mym ON m.snpID = mym.name")
	}
	if (!is.null(chr))
		sql <- paste0(sql, "\nWHERE m.chrID = '", gsub("chr","",chr), "'")
	if (!is.null(start))
		sql <- paste(sql, "AND\n m.positionBp >=", formatC(start, format = "d"))
	if (!is.null(end))
		sql <- paste(sql, "AND\n m.positionBp <=", formatC(end, format = "d"))
	sql <- paste0(sql, ";")
	
	if (verbose)
		cat(sql, "\n")
	
	rez <- .chunk.query(db, sql, batch.size)
	rez$allele <- with(rez, ifelse(call == "A", as.character(alleleA),
																 ifelse(call == "B", as.character(alleleB),
																 			 ifelse(call == "H", "H", "N"))))
	return(rez)
	
}

.fetch.controls.mm <- function(type = c("all","classical","wild"), ...) {
	
	require(plyr)
	
	tt <- match.arg(type)
	filters <- list(classical = c(A = "A/J%",B = "C57BL/6J%", C = "129S1%", D = "NOD/ShiLtJ%",E = "NZO/HILtJ%"),
									wild = c(F = "CAST/EiJ%", G = "PWK/PhJ%", H = "WSB/EiJ%"))
	ff <- character()
	if (tt == "classical") {
		ff <- c(ff, filters[[1]])
	}
	else if (tt == "wild") {
		ff <- c(ff, filters[[2]])
	}
	else {
		ff <- unlist(unname(filters))
	}
	
	rez <- ldply(ff, fetch.samples, by = "name", exact = FALSE, ...)
	colnames(rez)[1] <- "strain"
	return(rez)
	
}

## auxiliary functions to summarize results
summarize.intensity <- function(..., db = c("mm","mda")) {
	
	if (db == "mm")
		.summarize.intensity.mm(...)
	else if (db == "mda")
		.summarize.intensity.mda(...)
	else
		.dbnotfound()
	
}


.summarize.intensity.mm <- function(df, markers = NULL, by = .(sid), ...) {
	
	require(plyr)
	if (!is.null(markers))
		df <- subset(df, marker %in% markers)
	
	ddply(df, by, summarize,
				si = sum(x+y), n = length(x),
				theta = atan2(sum(y), sum(x)), hypo = sum(sqrt(x^2 + y^2)))	
	
}

.summarize.intensity.mda <- function(df, markers = NULL, by = .(sid), ...) {
	
	require(plyr)
	if (!is.null(markers))
		df <- subset(df, marker %in% markers)
	
	ddply(df, by, summarize,
				si = sum(average), n = length(x))

}