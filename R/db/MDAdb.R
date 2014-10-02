## Most logic can be ported over from functions written for MM database.
## Functions to extract samples and intensity have been modified to reflect MDA database schema.

source("~/Dropbox/pmdvlab/util/MMdb.R")
MDADB.PATH <- "/db/arrays/mda/MDA.db"

fetch.samples <- function(ids = NULL, db = MDADB.PATH, by = c("name","id"), exact = TRUE, ...) {
	
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

fetch.intensities <- function(ids, markers = NULL, chr = NULL, start = NULL, end = NULL, by = c("name","id"), db = MDADB.PATH, batch.size = 1000, verbose = TRUE, ...) {
	
	require(RSQLite)
	stopifnot(!is.null(ids))
	
	db <- dbConnect(SQLite(), dbname = db)
	
	.insert.samples(ids, db, by = by)
	
	sql <- paste("SELECT s.id as sid, s.name as name,",
							 "m.snpID as marker, m.chrID as chr, m.positionBp as pos, m.positioncM as cM, m.flagged as flag, m.mm10positionBp as mm10pos,",
							 "i.average as avg, i.contrast as contr, g.call as call",
							 "FROM samples as s",
							 paste0("INNER JOIN _mysamples as sg ON s.", by, " = sg.", by),
							 "INNER JOIN genotypes as g ON g.sampleID = s.id",
							 "INNER JOIN intensity as i on i.snpID = m.snpID",
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
	
	.chunk.query(db, sql, batch.size)
	
}