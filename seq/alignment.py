from Bio import AlignIO
from collections import Counter

## get variable columns from a multiple sequence alignment
def get_variable_sites(aln):

    nsites = len(aln[0])
    variable = []
    for i in range(0, nsites):
        if len(set(aln[:,i])) > 1:
            variable.append(i)

    return variable

## get consensus character from each column of a multiple sequence alignment
def get_consensus(aln):
    nsites = len(aln[0])
    return [ Counter(aln[:,i]).most_common(1)[0][0] for i in range(0, nsites) ]

## given list of variable positions in a multiple sequence alignment, clump them by proximity
def cluster_sites(pos, radius = 40, maxwidth = 40):
    clusters = []
    this_cluster = []
    pos = sorted(pos)
    if len(pos):
        last_p = pos.pop(0)
        this_cluster.append(last_p)
        while len(pos):
            p = pos.pop(0)
            if (p - last_p) > radius or (p - min(this_cluster)) > maxwidth:
                clusters.append(this_cluster)
                this_cluster = [p]
            else:
                this_cluster.append(p)
            last_p = p
        clusters.append(this_cluster)
        return clusters
    else:
        return None

## extract unique haplotypes for given intervals in a multiple sequence alignment
def get_haplotypes(aln, intervals, flank = 5):
    haps = []
    nsites = len(aln[0])
    for i in intervals:
        s = max( min(i)-flank, 0 )
        e = min( max(i)+flank, nsites )
        this_hap = set( [ str(aln[ j,s:e ].seq) for j in range(0, len(aln)) ] )
        haps.append(this_hap)
    return haps

## for each cluster of variable sites in a multiple sequence alignment, return k-mers representing all unique haplotypes
def get_variable_kmers(aln, k = 40):
    pos = get_variable_sites(aln)
    intervals = cluster_sites(pos, radius = k, maxwidth = k)
    return get_haplotypes(aln, intervals, 0)
