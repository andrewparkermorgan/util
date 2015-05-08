#!/usr/bin/python

import AffyAPI


def getSnpIntensity(row, cel):
    '''returns a list of intensities in the form of [senseA, senseA, antisenseA, antisenseA, senseB, senseB, antisenseB, antisenseB]'''
    cols=cel.info["Cols"]
    rows=cel.info["Rows"]
    xList = [row[9], row[13], row[25], row [29], row[17], row[21], row[33], row[37]] #x coordinates for 2 sense A probes, 2 antisense A, 2 sense B, and 2 antisense B
    yList = [row[10], row[14], row[26], row [30], row[18], row[22], row[34], row[38]] #y coordinates for 2 sense A probes, 2 antisense A, 2 sense B, and 2 antisense B
    intensity=[cel.getIntensity(x, rows-y-1) for (x, y) in zip(xList, yList)] #row coordinates are flipped in Affy cel data so use (x, rows-y-1) for coordinates (x,y)
    return intensity

def getOtherIntensity(row, cel):
    '''For exon1probes, exon2probes, and otherprobes tables. Input: row in sqlite3 cursor, cel file object.
    returns a list of intensities in the form of [sense, antisense]'''
    cols=cel.info["Cols"]
    rows=cel.info["Rows"]
    xList = [row[11], row[14]] #x coordinates for sense, antisense
    yList = [row[12], row[15]] #y coordinates for sense, antisense
    intensity=[cel.getIntensity(x, rows-y-1) for (x, y) in zip(xList, yList)]
    return intensity
    