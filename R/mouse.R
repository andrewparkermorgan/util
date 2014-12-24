## useful stuff about mouse genome for interval calculations

library(GenomicRanges)

## chromosome lengths as seqinfo object
mm9 <- Seqinfo(seqnames = paste0("chr", c(1:19, "X","Y","M")),
							 seqlengths = c(197195432, 181748087, 159599783, 155630120, 152537259,
							 							 149517037, 152524553, 131738871, 124076172, 129993255,
							 							 121843856, 121257530, 120284312, 125194864, 103494974,
							 							 98319150, 95272651, 90772031, 61342430,
							 							 166650296, 91744698, 16299),
							 isCircular = c(rep(FALSE, 21),TRUE),
							 genome = "mm9")

mm10 <- Seqinfo(seqnames = paste0("chr", c(1:19, "X","Y","M")),
								seqlengths = c(195471971,182113224,160039680,156508116,151834684,
															 149736546,145441459,129401213,124595110,130694993,
															 122082543,120129022,120421639,124902244,104043685,
															 98207768,94987271,90702639,61431566,
															 171031299,91744698,16299),
								isCircular = c(rep(FALSE, 21),TRUE),
								genome = "mm10")

## standard colours for Collaborative Cross strains
CC.STRAINS <- toupper(letters[1:8])
cc.strains <- c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ","NZO/HILtJ",
								"CAST/EiJ","PWK/PhJ","WSB/EiJ")
CC.COLORS <- c("#DAA520","#404040","#F08080","#1010F0","#00A0F0","#00A000","#F00000","#9000E0")
mouse.colors <- setNames( c(CC.COLORS, "grey60","black"), c(cc.strains, "FVB/NJ","SPRET/EiJ") )

## standard colors for Mus taxa
MUS.TAXA <- c("mus","dom","cas","musculus","domesticus","castaneus",
			  "molossinus",
			  "spretus","spicilegus","cypriacus","macedonicus",
			  "famulus","caroli","pahari",
			  "cookii","fragilicauda")
MUS.COLORS <- c("#e41a1c", "#377eb8","#4daf4a","#e41a1c", "#377eb8","#4daf4a",
				"brown",
				"grey40","cadetblue4","darkgoldenrod3","darkolivegreen",
				"darkorange4","burlywood4","aquamarine4",
				"black","black")
names(MUS.COLORS) <- c(MUS.TAXA)