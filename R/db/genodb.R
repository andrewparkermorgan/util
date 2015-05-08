## functions for interacting with array genotype databases

MUGADB.PATH = "/db/arrays/muga/geneseek.db"
MMDB.PATH = "/db/arrays/megamuga/megamuga.db"
#MDADB.PATH = "/Volumes/apm_passport/db/arrays/mda/MDA.db"
MDADB.PATH = "/Volumes/apm_passport/db/arrays/mda/snp.db"
GIGADB.PATH = "/db/arrays/gigamuga/gigamuga.db"

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
		print(sql)
		dbGetQuery(db, sql)
	}
	print( dbGetQuery(db, "select * from _mymarkers;") )
}

.chunk.query <- function(db, sql, batch.size = -1, one.piece = FALSE, ...) {

	require(RSQLite)
	stopifnot(inherits(db, "SQLiteConnection"))

	i = 0
	if (!one.piece) {
		rez <- dbSendQuery(db, sql)
		df <- data.frame()
		while (!dbHasCompleted(rez)) {
			next.rows <- fetch(rez, n = batch.size)
			df <- rbind(df, next.rows)
			i <- i + nrow(next.rows)
			cat(i, " records ...\n")
		}
		dbClearResult(rez)
		dbDisconnect(db)
	}
	else {
		df <- dbGetQuery(db, sql)
	}

	return(df)

}

.dbnotfound <- function() {
	stop("Specify a database in which to search ('mda' = MDA, 'mm' = MegaMuga). GigaMUGA not yet supported.")
}

## entry point for external calls
fetch.samples <- function(..., db = c("mm","mda","giga")) {

	if (db == "mm")
		.fetch.samples.mm(..., bad.flags = "u")
	else if (db == "muga")
		.fetch.samples.mm(..., db = MUGADB.PATH)
	else if (db == "mda")
		.fetch.samples.mda(...)
	else if (db == "giga")
		.fetch.samples.mm(..., db = GIGADB.PATH, bad.flags = "B$")
	else
		.dbnotfound()

}

fetch.intensities <- function(..., db = c("mm","muga","mda","giga")) {

	if (db == "mm")
		#.fetch.intensities.mm(...)
		.fetch.intensities.python(..., pos.col = "position", db = MMDB.PATH)
	else if (db == "muga")
		#.fetch.intensities.mm(...)
		.fetch.intensities.python(..., pos.col = "position", db = MUGADB.PATH)
	else if (db == "mda")
		#.fetch.intensities.mda(...)
		.fetch.intensities.python(..., db = MDADB.PATH, mda = TRUE)
	else if (db == "giga")
		#.fetch.intensities.giga(..., pos.col = "position38", db = GIGADB.PATH)
		.fetch.intensities.python(..., pos.col = "position38", db = GIGADB.PATH)
	else
		.dbnotfound()

}

fetch.calls <- function(..., db = c("mm","muga","mda","giga")) {
	
	if (db == "mm")
		## NB: for MM, intensities and calls always returned together
		#.fetch.intensities.mm(...)
		.fetch.intensities.python(..., pos.col = "position", db = MMDB.PATH)
	else if (db == "muga")
		## NB: for MM, intensities and calls always returned together
		#.fetch.intensities.mm(...)
		.fetch.intensities.python(..., pos.col = "position", db = MUGADB.PATH)
	else if (db == "giga")
		## NB: for MM, intensities and calls always returned together
		#.fetch.intensities.mm(..., pos.col = "position38", db = GIGADB.PATH, one.piece = TRUE)
		.fetch.intensities.python(..., pos.col = "37", db = GIGADB.PATH)
	else if (db == "mda")
		.fetch.calls.mda(...)
	else
		.dbnotfound()

}

fetch.controls <- function(..., db = c("mm","mda","giga")) {

	if (db == "mm")
		.fetch.controls.mm(...)
	if (db == "giga")
		.fetch.controls.giga(...)
	else
		stop("Control samples only defined for MegaMUGA, GigaMUGA so far.")

}

## array-specific functions
.fetch.samples.mm <- function(ids = NULL, flags = NULL, bad.flags = "", group = NULL, db = MMDB.PATH, by = c("name","id"),
							  exact = TRUE, strict.case = TRUE, verbose = TRUE, ...) {

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
	else {
		sql <- paste0(sql, " FROM SAMPLES s")
	}
	if (!strict.case)
		sql <- paste0(sql, " COLLATE NOCASE")
	sql <- paste0(sql, ";")

	if (verbose)
		cat(sql, "\n")

	rez <- .chunk.query(db, sql, -1)
	rez$bad <- grepl(bad.flags[1], rez$flags)
	if (!is.null(flags)) {
		.flags <- flags
		rez <- subset(rez, grepl(.flags, flags))
		message(paste0(nrow(rez), " records with requested flags."))
	}
	return(rez)

}

.fetch.samples.mda <- function(ids = NULL, db = MDADB.PATH, by = c("name","id"), exact = TRUE, ...) {

	require(RSQLite)

	db <- dbConnect(SQLite(), dbname = db)

	#cols <- paste("s", c("id", "name", "cel_file"), sep = ".", collapse = ", ")
	#sql <- paste("SELECT", cols)
	#if (!is.null(ids)) {
	#	if (exact) {
	#		.ids <- na.omit(ids)
	#		.insert.samples(.ids, db, by = match.arg(by))
	#		sql <- paste0(sql, " FROM samples s ",
	#					  "INNER JOIN _mysamples as sg ON s.", by, " = sg.", by)
	#	}
	#	else {
	#		sql <- paste0(sql, " FROM samples s ",
	#					  "WHERE s.name LIKE '", ids[1], "'")
	#	}
	#}
	sql <- paste0("SELECT s.celid as id, s.name as name, s.species, s.subspecies, s.arrayrunat as runat ",
				  "FROM sample_info as s ")
	if (!is.null(ids)) {
		if (exact) {
			.ids <- na.omit(ids)
			.insert.samples(.ids, db, by = match.arg(by))
			sql <- paste0(sql, "INNER JOIN _mysamples as sg ON s.", by, " = sg.", by)
		}
		else {
			sql <- paste0(sql,"WHERE name LIKE '", ids[1], "'")
		}
	}
	
	sql <- paste0(sql, ";")
	cat(sql, "\n")

	.chunk.query(db, sql, -1)

}

# .fetch.intensities.mm <- function(ids, markers = NULL, chr = NULL, start = NULL, end = NULL, pos.col = "position",
# 								  by = c("name","id"), db = MMDB.PATH, batch.size = -1, verbose = TRUE,
# 								  one.piece = FALSE, bluff = FALSE, ...) {

# 	require(RSQLite)
# 	stopifnot(!is.null(ids))

# 	db <- dbConnect(SQLite(), dbname = db)

# 	.insert.samples(ids, db, by = by)

# 	sql <- paste("SELECT m.name as marker, s.id as sid, s.name as id, m.chromosome as chr, ",
# 				 paste0("m.", pos.col)," as pos, g.x, g.y, g.call",
# 				 "FROM samples as s",
# 				 paste0("INNER JOIN _mysamples as sg ON s.", by, " = sg.", by),
# 				 "INNER JOIN genotypes as g ON g.sampleID = s.id",
# 				 "INNER JOIN snps as m on g.snpID = m.id", sep = "\n")
# 	if (!is.null(markers)) {
# 		if (!is.numeric(markers)) {
# 			.insert.markers(markers, db, by = "name")
# 			sql <- paste0(sql, "\nINNER JOIN _mymarkers as mym ON m.name = mym.name")
# 		}
# 		else {
# 			.insert.markers(markers, db, by = "id")
# 			sql <- paste0(sql, "\nINNER JOIN _mymarkers as mym ON m.id = mym.id")
# 		}
# 	}
# 	if (!is.null(chr))
# 		sql <- paste0(sql, "\nWHERE chr = '", gsub("chr","",chr), "'")
# 	if (!is.null(start))
# 		sql <- paste(sql, "AND\npos >=", formatC(start, format = "d"))
# 	if (!is.null(end))
# 		sql <- paste(sql, "AND\npos <=", formatC(end, format = "d"))
# 	sql <- paste0(sql, ";")

# 	if (verbose)
# 		cat(sql, "\n")

# 	if (!bluff) {
# 		rez <- .chunk.query(db, sql, batch.size, one.piece = one.piece)
# 		rez$marker <- make.names(as.character(rez$marker))
# 		return(rez)
# 	}
# 	else {
# 		return(NULL)
# 	}

# }

.fetch.intensities.python <- function(ids, markers = NULL, chr = NULL, start = NULL, end = NULL, pos.col = "position38",
									by = c("id","name"), operation = "intensities", db = MMDB.PATH, mda = FALSE, script = "~/lib/util/db/genodb.py", ... ) {

	require(data.table)

	ff = tempfile()
	ff.id = tempfile()
	write.table(ids, ff.id, row.names = FALSE, col.names = FALSE, quote = FALSE)
	mm.id = tempfile()
	write.table(markers, mm.id, row.names = FALSE, col.names = FALSE, quote = FALSE)

	cmd <- paste0(script, " -o ", operation," --db ", db, " --samples ", ff.id, " --by ", by[1], " --out ", ff)
	if (mda) {
		cmd <- paste0(cmd, " --mda ")
	}
	if (!is.null(markers))
		cmd <- paste0(cmd, " --markers ", mm.id)
	if (!is.null(chr))
		cmd <- paste0(cmd, " --chr ", chr)
	if (!is.null(start))
		cmd <- paste0(cmd, " --from-bp ", formatC(start, format = "d"))
	if (!is.null(end))
		cmd <- paste0(cmd, " --to-bp ", formatC(end, format = "d"))
	cmd <- paste0(cmd, " --pos ", pos.col)

	cat(cmd)
	system(cmd)
	rez <- fread(ff)
	if (nrow(rez)) {
		message("Deleting temporary file...")
		system(paste("rm", ff), intern = FALSE)
	}
	return(rez)

}

.fetch.intensities.giga <- function(...) {
	.fetch.intensities.python(..., db = GIGADB.PATH)
}

.fetch.intensities.mm <- function(...) {
	.fetch.intensities.python(..., db = MMDB.PATH)
}

# .fetch.intensities.mda <- function(ids, markers = NULL, chr = NULL, start = NULL, end = NULL, by = c("name","id"), db = MDADB.PATH, batch.size = 1000, verbose = TRUE, ...) {
# 
# 	require(RSQLite)
# 	stopifnot(!is.null(ids))
# 
# 	db <- dbConnect(SQLite(), dbname = db)
# 
# 	.insert.samples(ids, db, by = by)
# 
# 	sql <- paste("SELECT s.id as sid, s.name as name,",
# 				 "m.snpID as marker, m.chrID as chr, m.positionBp as pos, m.positioncM as cM, m.flagged as flag, m.mm10positionBp as mm10pos,",
# 				 "i.average as avg, i.contrast as contr",
# 				 "FROM samples as s",
# 				 paste0("INNER JOIN _mysamples as sg ON s.", by, " = sg.", by),
# 				 #"INNER JOIN genotypes as g ON g.sampleID = s.id",
# 				 "INNER JOIN intensity as i on i.sampleID = s.id",
# 				 "INNER JOIN snpInfo as m on i.snpID = m.snpID", sep = "\n")
# 	if (!is.null(markers)) {
# 		.insert.markers(markers, db, by ="name")
# 		sql <- paste0(sql, "\nINNER JOIN _mymarkers as mym ON m.snpID = mym.name")
# 	}
# 	if (!is.null(chr))
# 		sql <- paste0(sql, "\nWHERE m.chrID = '", gsub("chr","",chr), "'")
# 	if (!is.null(start))
# 		sql <- paste(sql, "AND\n m.positionBp >=", formatC(start, format = "d"))
# 	if (!is.null(end))
# 		sql <- paste(sql, "AND\n m.positionBp <=", formatC(end, format = "d"))
# 	sql <- paste0(sql, ";")
# 
# 	if (verbose)
# 		cat(sql, "\n")
# 
# 	.chunk.query(db, sql, batch.size)
# 
# }

.fetch.calls.mda <- function(..., db = MDADB.PATH) {

	# require(RSQLite)
	# stopifnot(!is.null(ids))
	#
	# db <- dbConnect(SQLite(), dbname = db)
	#
	# .insert.samples(ids, db, by = by, verbose = TRUE)
	#
	# sql <- paste("SELECT s.id as sid, s.name as name,",
	# 			 "m.snpID as marker, m.chrID as chr, m.positionBp as pos, m.positioncM as cM, m.flagged as flag, m.mm10positionBp as mm10pos,",
	# 			 "m.alleleA, m.alleleB, g.call as call",
	# 			 "FROM samples as s",
	# 			 paste0("INNER JOIN _mysamples as sg ON s.", by, " = sg.", by),
	# 			 "INNER JOIN genotypes as g ON g.sampleID = s.id",
	# 			 #"INNER JOIN intensity as i on i.snpID = m.snpID",
	# 			 "INNER JOIN snpInfo as m on g.snpID = m.snpID", sep = "\n")
	# if (!is.null(markers)) {
	# 	.insert.markers(markers, db, by ="name")
	# 	sql <- paste0(sql, "\nINNER JOIN _mymarkers as mym ON m.snpID = mym.name")
	# }
	# if (!is.null(chr))
	# 	sql <- paste0(sql, "\nWHERE m.chrID = '", gsub("chr","",chr), "'")
	# if (!is.null(start))
	# 	sql <- paste(sql, "AND\n m.positionBp >=", formatC(start, format = "d"))
	# if (!is.null(end))
	# 	sql <- paste(sql, "AND\n m.positionBp <=", formatC(end, format = "d"))
	# sql <- paste0(sql, ";")
	#
	# if (verbose)
	# 	cat(sql, "\n")

	rez <- .fetch.intensities.python(..., operation = "genotypes", db = db, mda = TRUE)

	#rez$allele <- with(rez, ifelse(call == "A", as.character(alleleA),
	#							   ifelse(call == "B", as.character(alleleB),
	#							   	   ifelse(call == "H", "H", "N"))))
	return(rez)

}

.fetch.controls.mm <- function(type = c("all","classical","wild"), f1 = FALSE, ...) {

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

	if (f1) {
		pairs <- apply( combn(seq_along(ff), 2), 2, function(i) paste0("%", ff[ i[1] ], "x", ff[ i[2] ]) )
		names(pairs) <- apply( combn(seq_along(ff), 2), 2, function(i) paste0(names(ff)[ i[1] ], names(ff[ i[2] ])) )
		print(pairs)
		ff <- c(ff, pairs)
	}
	
	rez <- ldply(ff, fetch.samples, by = "name", exact = FALSE, db = "mm", ...)
	colnames(rez)[1] <- "strain"
	rez <- transform(rez, type = ifelse(nchar(strain) == 1, "inbred","F1"))
	return(rez)

}

.fetch.controls.giga <- function(type = c("all","classical","wild"), f1 = FALSE, ...) {
	
	require(plyr)
	
	tt <- match.arg(type)
	filters <- list(classical = c(A = "AA_%",B = "BB_%", C = "CC_%", D = "NOD/ShiLtJ%",E = "NZO/HILtJ%"),
					wild = c(F = "CAST/EiJ%", G = "PWK/PhJ%", H = "WSB/EiJ%"))
	ff <- character()
	if (tt == "classical") {
		ff <- sapply(toupper(letters[1:5]), function(x) paste0(x, x))
	}
	else if (tt == "wild") {
		ff <- sapply(toupper(letters[6:8]), function(x) paste0(x, x))
	}
	else {
		ff <- sapply(toupper(letters[1:8]), function(x) paste0(x, x))
	}
	
	if (f1) {
		pairs <- apply( combn(ff, 2), 2, function(i) paste0(substr(i[1], 1, 1), substr(i[2], 2, 2)) )
		names(pairs) <- pairs
		print(pairs)
		ff <- c(ff, pairs)
	}
	ff <- setNames( paste0(ff, "_%"), names(ff) )
	rez <- ldply(ff, fetch.samples, by = "name", exact = FALSE, db = "giga", ...)
	colnames(rez)[1] <- "strain"
	rez <- transform(rez, type = ifelse(nchar(strain) == 1, "inbred","F1"))
	
	## dump a lot of DO samples with names starting with 'AA'
	rez <- subset(rez, !grepl("AA\\-DO", name))
	
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
