## useful stuff about mouse genome for interval calculations

library(GenomicRanges)

## chromosome lengths as seqinfo object
mm9 <- Seqinfo(seqnames = paste0("chr", c(1:19, "X","Y","M")),
							 seqlengths = c(197195432, 181748087, 159599783, 155630120, 152537259,
							 							 149517037, 152524553, 131738871, 124076172, 129993255,
							 							 121843856, 121257530, 120284312, 125194864, 103494974,
							 							 98319150, 95272651, 90772031, 61342430,
							 							 166650296, 91744698, 16299),
							 genome = "mm9")