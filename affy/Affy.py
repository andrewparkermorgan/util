#!/usr/bin/python

import sys
import struct

def printCELHeaderList(header, depth=0):
    if (type(header) == type(tuple())):
        parens = "()"
    else:
        parens = "[]"
    print "%s%s" % (depth*'    ', parens[0])
    depth += 1
    for item in header:
        if (type(item) == type(dict())):
            if (len(item) == 0):
                print "%s{}," % (depth*'    ')
            else:
                print "%s{" % (depth*'    ')
                for key, value in item.items():
                    print "%s  %s = %s %s" % (depth*'    ', key, value, str(type(value)))
                print "%s}," % (depth*'    ')
        elif (type(item) == type(list())) or (type(item) == type(tuple())):
            printCELHeaderList(item, depth)
        elif (type(item) == type(u'')):
            print '%su"%s",' % (depth*'    ', item)
        elif (type(item) == type('')):
            print '%s"%s",' % (depth*'    ', item)
        else:
            print '%s%s,' % (depth*'    ', str(item))
    print "%s%s" % ((depth-1)*'    ', parens[-1])

def printQCProbe(probeData):
    print "%8d:" % probeData[0]
    for i in xrange(len(probeData[1])):
        print "  Cell %2d: %s" % (i, str(probeData[1][i]))

def printProbe(probeData):
    print "%8d, %1d, %1d:" % probeData[0]
    for i in xrange(len(probeData[1])):
        print "  Block %2d: %s" % (i, str(probeData[1][i]))
    for i in xrange(len(probeData[2])):
        print "  Cell %3d: %s" % (i, str(probeData[2][i]))
     
def printCELHeirarchy(cel, header=None, depth=0):
    print "%s%2d Parameters" % (depth*'  ', len(cel.parameters(header)))
    parList = cel.parents(header)
    if (len(parList) > 0):
        print "%s%2d Parents" % (depth*'  ', len(parList))
        for parent in parList:
            printCELHeirarchy(cel, parent, depth+1)
    
class CelData:
    def __init__(self, filename):
        self.fp = open(filename, 'rb')
        sbuffer = self.fp.read(1)
        magic = ord(sbuffer[0]) 
        if (magic == 59):
            # Command Console generic data format
            sbuffer += self.fp.read(9)
            self.FileHeader = struct.unpack(">BBiI", sbuffer)
            self.GenericDataHeader = self.readGenericDataHeader()
            self.DataGroups = self.readDataGroups()
        elif (magic == 64):
            sbuffer += self.fp.read(19)
            self.FileHeader = struct.unpack("<iiiii", sbuffer)
            self.GenericDataHeader = self.readVer4Header()
        else:
            print "Invalid Magic Number: %d" % magic
            return
    
    def typeId(self, header = None):
        if (header is None):
            header = self.GenericDataHeader
        try:
            return header[0]
        except IndexError:
            return ""

    def fileId(self, header = None):
        if (header is None):
            header = self.GenericDataHeader
        try:
            return header[1]
        except IndexError:
            return ""

    def dateTime(self, headerList = None):
        if (header is None):
            header = self.GenericDataHeader
        try:
            return header[2]
        except IndexError:
            return ""

    def locale(self, header = None):
        if (header is None):
            header = self.GenericDataHeader
        try:
            return header[3]
        except IndexError:
            return ""

    def parameters(self, header = None):
        if (header is None):
            header = self.GenericDataHeader
        try:
            return header[4]
        except IndexError:
            return dict()

    def parents(self, header = None):
        if (header is None):
            header = self.GenericDataHeader
        try:
            return header[5:]
        except IndexError:
            return list()    
    
    def getValue(self, groupSet, rowCol):
        validTypes = [">b", ">B", ">h", ">H", ">i", ">I", ">f"]
        dataSet = self.DataGroups[groupSet[0]][1][groupSet[1]]
        dataOffset = dataSet[-1] + dataSet[-2]*rowCol[0] + dataSet[-3][rowCol[1]]
        colInfo = dataSet[4][rowCol[1]]
        self.fp.seek(dataOffset)
        sbuffer = self.fp.read(colInfo[2])
        if (len(sbuffer) != colInfo[2]):
            print rowCol, dataOffset, colInfo
        return struct.unpack(validTypes[colInfo[1]], sbuffer)[0]    
    
    def readVer4Header(self, header = None):
        if (header is None):
            header = list()
            self.fp.seek(20)
        header.append("affymetrix-GCOS-version-4")
        header.append("")
        header.append(u"")
        header.append(u"en-US")
        # read header
        fieldSize = struct.unpack("<i", self.fp.read(4))[0]
        line = self.fp.read(fieldSize).split('\n')
        parameters = dict()
        for pString in line:
            separator = pString.find('=')
            if (separator < 0):
                continue
            pName = pString[:separator]
            pValue = pString[separator+1:]
            try:
                parameters[pName] = int(pValue)
            except ValueError:
                parameters[pName] = pValue
        fieldSize = struct.unpack("<i", self.fp.read(4))[0]
        line = self.fp.read(fieldSize)
        parameters["Algorithm"] = line
        fieldSize = struct.unpack("<i", self.fp.read(4))[0]
        line = self.fp.read(fieldSize)
        parameters["AlgorithmParameters"] = line
        header.append(parameters)
        margin, outliers, masked, subgrids = struct.unpack("<iIIi", self.fp.read(16))
        parameters["Margin"] = margin
        # print margin, outliers, masked, subgrids
        return header

    def readGenericDataHeader(self, header = None):
        if (header is None):
            header = list()
            self.fp.seek(10)
        header.append(self.readString())                # data type identifier
        header.append(self.readString())                # file identifier
        header.append(self.readWString())               # Date Time        
        locale = self.readWString()                     # Locale
        header.append(locale)
        parameterCount = struct.unpack(">i", self.fp.read(4))[0]
        parameters = dict()
        while (parameterCount > 0):
            pName = self.readWString()
            pValue = self.readString()
            pType = self.readWString()
            parameters[pName] = self.decodeMime(pValue, pType)
            parameterCount -= 1
        header.append(parameters)
        parentCount = struct.unpack(">i", self.fp.read(4))[0]
        while (parentCount > 0):
            header.append(self.readGenericDataHeader([]))
            parentCount -= 1
        return header

    def readDataGroups(self):
        dataGroup = list()
        nextGroup = self.FileHeader[3]
        while (nextGroup):
            self.fp.seek(nextGroup)
            groupInfo = struct.unpack(">IIi", self.fp.read(12))
            groupName = self.readWString()
            dataSet = list()
            dataOffset = groupInfo[1]
            for i in xrange(groupInfo[2]):
                self.fp.seek(dataOffset)
                dataInfo = struct.unpack(">II", self.fp.read(8))
                dataName = self.readWString()
                parameterCount = struct.unpack(">i", self.fp.read(4))[0]
                parameters = dict()
                while (parameterCount > 0):
                    pName = self.readWString()
                    pValue = self.readString()
                    pType = self.readWString()
                    parameters[pName] = self.decodeMime(pValue, pType)
                    parameterCount -= 1
                dataColumns = struct.unpack(">I", self.fp.read(4))[0]
                colLabel = list()
                colOffsets = list()
                rowSize = 0
                for j in xrange(dataColumns):
                    colName = self.readWString()
                    colInfo = struct.unpack(">bi", self.fp.read(5))
                    colLabel.append([colName, colInfo[0], colInfo[1]])
                    colOffsets.append(rowSize)
                    rowSize += colInfo[1]
                dataRows = struct.unpack(">I", self.fp.read(4))[0]
                dataSet.append((dataName, parameters, dataColumns, dataRows, colLabel, colOffsets, rowSize, dataInfo[0]))
                dataOffset = dataInfo[1]
            dataGroup.append((groupName, dataSet))
            nextGroup = groupInfo[0]
        return dataGroup
    
    def readString(self):
        strLen = struct.unpack(">i", self.fp.read(4))[0]
        return self.fp.read(strLen)

    def readWString(self):
        strLen = struct.unpack(">i", self.fp.read(4))[0]
        return self.fp.read(2*strLen).decode("utf-16-be")
    
    def decodeMime(self, value, type):
        if (type == "text/ascii"):
            return value
        elif (type == "text/plain"):
            return value.decode("utf-16-be")
        elif (type == "text/x-calvin-float"):
            return struct.unpack(">f", value[:4])[0]
        elif (type == "text/x-calvin-integer-32"):
            return struct.unpack(">i", value[:4])[0]
        elif (type == "text/x-calvin-integer-16"):
            return struct.unpack(">h", value[2:4])[0]
        elif (type == "text/x-calvin-integer-8"):
            return struct.unpack(">b", value[3:4])[0]
        elif (type == "text/x-calvin-unsigned-integer-32"):
            return struct.unpack(">I", value[:4])[0]
        elif (type == "text/x-calvin-unsigned-integer-16"):
            return struct.unpack(">H", value[2:4])[0]
        elif (type == "text/x-calvin-unsigned-integer-8"):
            return struct.unpack(">B", value[3:4])[0]
        else:
            raise TypeError, 'Unknown MIME type: "%s"' % type

class CdfData:
    def __init__(self, filename):
        # currently only handles the XDA binary format
        self.fp = open(filename, 'rb')
        self.FileHeader = struct.unpack("<iiHH", self.fp.read(12))
        self.columns = self.FileHeader[2]
        self.rows = self.FileHeader[3]
        units, qcUnits, seqLen = struct.unpack("<iii", self.fp.read(12))
        self.customSequence = self.fp.read(seqLen)
        probeList = [self.fp.read(64).strip('\0') for i in xrange(units)]
        self.qcOffset = list(struct.unpack("<%di" % qcUnits, self.fp.read(4*qcUnits)))
        probeOffset = list(struct.unpack("<%di" % units, self.fp.read(4*units)))
        self.probeTable = dict(zip(probeList, probeOffset))
    
    def getProbeData(self, probe):
        data = self.probeTable[probe]
        if (type(data) == type(int())):
            data = self.readProbe(data)
            self.probeTable[probe] = data
        return data
    
    def getQCData(self, index):
        data = self.qcOffset[index]
        if (type(data) == type(int())):
            data = self.readQC(data)
            self.qcOffset[index] = data
        return data        
            
    def readProbe(self, offset):
        data = list()
        self.fp.seek(offset)
        type, direction, numAtoms, numBlocks, numCells, probeID, cellsPerAtom = struct.unpack("<HBiiiiB", self.fp.read(20))
        data.append((probeID, type, direction))
        # read blocks and cells
        blockList = list()
        cellList = list()
        for blk in xrange(numBlocks):
            blkAtoms, blkCells, blkCellsPerAtom, blkDirection, blkCellPos, blkUnused = struct.unpack("<iiBBii", self.fp.read(18))
            blkName = self.fp.read(64).strip('\0')
            if (self.FileHeader[1] == 2):
                blkWobble, blkAllele = struct.unpack("<HH", self.fp.read(4))
                blockList.append((blkName, blkDirection, blkWobble, blkAllele))
            else:
                blockList.append((blkName, blkDirection))
            for cel in xrange(blkCells):
                atomNum, cellX, cellY, cellIndex, baseProbe, baseTarget = struct.unpack("<iHHicc", self.fp.read(14))
                if (self.FileHeader[1] == 2):
                    seqLength, probeGroup = struct.unpack("<HH", self.fp.read(4))
                    cellList.append((cellX, cellY, cellIndex, baseProbe, baseTarget, seqLength, probeGroup, blk))
                else:
                    cellList.append((cellX, cellY, cellIndex, baseProbe, baseTarget, blk))
            numAtoms -= blkAtoms
            numCells -= blkCells
        data.append(blockList)
        data.append(cellList)
        # Assert numCells == 0 and numAtoms == 0
        if (numCells != 0) or (numAtoms != 0):
            print "Inconsistent cell count", numCells, numAtoms, cellsPerAtom, blkAtoms, blkCellsPerAtom
            printProbe(data)
        return data

    def readQC(self, offset):
        data = list()
        self.fp.seek(offset)
        type, numProbes = struct.unpack("<Hi", self.fp.read(6))
        data.append(type)
        probeList = list()
        for probe in xrange(numProbes):
            probeX, probeY, probeLen, matchFlag, backFlag = struct.unpack("<HHBBB", self.fp.read(7))
            probeList.append((probeX, probeY, probeLen, matchFlag, backFlag))
        data.append(probeList)
        return data
   
if __name__ == "__main__":
    celTest = CelData(sys.argv[1])
    printCELHeaderList(celTest.GenericDataHeader)
    print celTest.DataGroups
    parameter = celTest.parameters()
    rows = parameter["affymetrix-cel-rows"]
    cols = parameter["affymetrix-cel-cols"]
    print rows, cols
        
