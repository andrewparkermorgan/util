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

## pseudoautosomal regions
# first coordinate is CAST/EiJ 'extended PAR', second is normal PAR
PAR.mm9 <- c(165980013,166.41e6)
PAR.mm10 <- c(169542082,169969759)

## standard colours for Collaborative Cross strains
CC.STRAINS <- toupper(letters[1:8])
cc.strains <- c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ","NZO/HILtJ",
								"CAST/EiJ","PWK/PhJ","WSB/EiJ")
CC.COLORS <- setNames( c("#DAA520","#404040","#F08080","#1010F0","#00A0F0","#00A000","#F00000","#9000E0"), CC.STRAINS )
CC.COLORS <- c( CC.COLORS, setNames( rep("grey60",28), apply(combn(CC.STRAINS, 2), 2, paste, collapse = "") ) )
mouse.colors <- c(CC.COLORS, setNames( c("grey40","grey20"), c("FVB/NJ","SPRET/EiJ") ))
cc.colors <- setNames( c(CC.COLORS[1:8], CC.COLORS[c(5,1)]),  c(cc.strains, "NZO/HlLtJ","129S1SvlmJ") )

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

## Sequence Ontology terms used by Sanger MGP
SO.TERMS <- c(
	"3_prime_UTR_variant",
	"5_prime_UTR_variant",
	"coding_sequence_variant",
	"downstream_gene_variant",
	"feature_elongation",
	"feature_truncation",
	"frameshift_variant",
	"incomplete_terminal_codon_variant",
	"inframe_deletion",
	"inframe_insertion",
	"initiator_codon_variant",
	"intergenic_variant",
	"intron_variant",
	"mature_miRNA_variant",
	"missense_variant",
	"NMD_transcript_variant",
	"nc_transcript_variant",
	"non_coding_exon_variant",
	"regulatory_region_ablation",
	"regulatory_region_amplification",
	"regulatory_region_variant",
	"splice_acceptor_variant",
	"splice_donor_variant",
	"splice_region_variant",
	"stop_gained",
	"stop_lost",
	"stop_retained_variant",
	"synonymous_variant",
	"TF_binding_site_variant",
	"TFBS_ablation",
	"TFBS_amplification",
	"transcript_ablation",
	"transcript_amplification",
	"upstream_gene_variant" )