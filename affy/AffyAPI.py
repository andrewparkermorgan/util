#!/usr/bin/python

import sys
import struct
from pylab import *
from AffyParam import *
import numpy

maxVal=10000

class CEL:
    """
    CEL class
    Public methods are named without underscore (except for the constructor).
    Protected methods are named with one heading underscore.        
    """
    # Public Interface -------------------------------------------
    def __init__(self, celFileName):        
        """ 
        CEL class constructor
        """
        # Initialize flags.
        self._intensityInMem = self._stdDevInMem = self._pixCntInMem = 0 
        self._mkcellInMem = self._outlierInMem = self._subgridInMem = 0
        self.info = dict() # Initialize the dictionary       
        self.info_exception = dict() # Initialize the exception dictionary 
        self.size = (0, 0)
        self._fileHeader = list()        
        self._intensityOffset = self._stdDevOffset = self._pixCntOffset = -1
        self._mkcellOffset = self._outlierOffset = self._subgridOffset = -1                
        self._intensityCol = self._stdDevCol = self._pixCntCol = -1
        self._mkcellCol = self._outlierCol = self._subgridCol = -1
        self.debug = False
        # Read the file header.
        self._fp = open(celFileName, 'rb')
        sbuffer = self._fp.read(1)
        magic = ord(sbuffer[0]) 
        self.info["MagicNumber"] = magic
        if (magic == 59): # Command Console generic data format            
            sbuffer += self._fp.read(9)
            self._fileHeader = struct.unpack(">BBiI", sbuffer)
            #self._genericDataHeader = self._readGenericDataHeader()
            self._readGenericDataHeader()
            self.version = self._fileHeader[1]
            self.info["FileVersion"] = self.version      
            self.info["DataGroups"] = self._fileHeader[2]
            self.info["FilePositionForFirstDataGroup"] = self._fileHeader[3]            
            self._readDataGroups() # Parse the group data
##            if (self.debug):
##                print self._fileHeader           
            
        elif (magic == 64): # Version 4 format            
            sbuffer += self._fp.read(19)            
            self._fileHeader = struct.unpack("<iiiii", sbuffer)
            self.version = self._fileHeader[1]
            self.info["FileVersion"] = self.version
            self.info["Cols"] = self._fileHeader[3] # The Cols and Rows are flipped for the purpose of image display 
            self.info["Rows"] = self._fileHeader[2]     
            self.info["Cells"] = self._fileHeader[4]      
            self.size = (self._fileHeader[2], self._fileHeader[3])
            self.nCell = self._fileHeader[4]
            #self._genericDataHeader = self._readVer4Header()            
            self._readVer4Header()
##            if (self.debug):        
##                print self._fileHeader
                
        else:
            #print "Unknown file format!! Invalid Magic Number: %d" % magic
            return
        
        return
    
    def __del__(self):
        """ 
        Finalizer 
        """    
##        if (self.debug):
##            print "Delete CEL object."
##        close(self._fp)
    
    def _parseDatHeader(self, header):
        """ 
        Parse parameters in the DAT file into the dictionary. 
        """        
        delim = header.find(":")
        self.info["CLS"] = int(header[delim+5:delim+10])               
        self.info["RWS"] = int(header[delim+14:delim+19])         
        self.info["XIN"] = int(header[delim+23:delim+26])        
        self.info["YIN"] = int(header[delim+30:delim+33])         
        self.info["VE"] = float(header[delim+36:delim+38])         
        self.info["Temperature"] = 0
        if (not header[delim+38:delim+45]):
            self.info["Temperature"] = float(header[delim+38:delim+45])         
        self.info["LaserPower"] = float(header[delim+45:delim+49])         
        self.info["ScanDate"] = header[delim+50:delim+58]
        self.info["ScanTime"] = header[delim+59:delim+67]
        self.info["DAT_Undefined"] = header[delim+68:]        
        
    def _stringParser(self, input, delimiter, string_delim=None):
        """ 
        Parse the given input using the delimiter. Use "string_delim" if the second delimiter is required.
        """
        result = list()        
        
        if (string_delim != None):
            input = input.rsplit(string_delim)
        
        for pString in input:
            separator = pString.find(delimiter)            
            if (separator < 0):
                continue
            pName = pString[:separator]
            pValue = pString[separator+1:]
            result.append((pName, pValue))
        return result
    
    def _safeAssign(self, paramName, value, mode):
        """
        Prevent two same parameters with different values from being overwritten.             
        """
        if (mode == 43): # Version 4 overwrites version 3: Copy the version 3 value into exception.
            try:
                if (self.info[paramName] != value):                
                    self.info_exception[paramName] = self.info[paramName]
                    self.info[paramName] = value
                    #print "[Warning] The parameter %s is not consistent between version 3 and 4!" % paramName
            except KeyError:
                #print "[Warning] The parameter %s doesn't exist!" % paramName 
                pass
        elif (mode == 34): # Version 3 wants to overwrite version 4: Copy the version 3 value into exception.
            try:
                if (self.info[paramName] != value):                
                    self.info_exception[paramName] = value                    
                    #print "[Warning] The parameter %s is not consistent between version 3 and 4!" % paramName
            except KeyError:
##                print "[Warning] The parameter %s doesn't exist!" % paramName
                pass
        elif (mode == 1): # Parent overwrite the child in Command Console version 1
            try:
                if (self.info[paramName] != value):                               
                    self.info_exception[paramName] = value                    
                    #print "[Warning] The parameter %s is not consistent between parent and child (Command Console version 1)!" % paramName
            except KeyError:
                self.info[paramName] = value                    
    
    def _readVer4Header(self):
        """
        Read version 4 file header, parameters, and embedded version 3 headers (including its file header, parameters, and the DAT header).
        """
        self._fp.seek(20)
        # Parse version 4 file header ---------------------------------------------------------------------------------
        fieldSize = struct.unpack("<i", self._fp.read(4))[0] # header length       
        line = self._fp.read(fieldSize).split('\n') # Version 3 file header & its parameters                
        # Parse version 3 file header, DAT header, and algorithm parameters ----------------------------------            
        self.info_exception = dict() # Initialize the exception dictionary 
        pString = self._stringParser(line, "=");        
        for i in pString:
            pName, pValue = i            
            try:                
                self._safeAssign(v3_to_global(pName), int(pValue), 34)                                 
            except ValueError:
                if (pName.find("GridCorner") != -1):
                    sep = pValue.find(' ')                    
                    self.info[v3_to_global(pName)] = (float(pValue[:sep]), float(pValue[sep+1:]))
                elif (pName == "DatHeader"):
                    self._parseDatHeader(pValue)
                elif (pName == "AlgorithmParameters"): 
                    algString = self._stringParser(pValue, ":", ";") # v3 algorithm parameters                    
                    for j in algString:
##                        if (self.debug):              
##                            print "v3 algo param: ", j[0], j[1]          
                        try:
                            self.info[v3_to_global(j[0])] = int(j[1])
                        except ValueError: # Catch a non-integer value as a string.                           
                            self.info[v3_to_global(j[0])] = j[1]                            
                else:    
                    self.info[v3_to_global(pName)] = pValue                               
        # Parse version 4 file header ---------------------------------------------------------------------------------
        fieldSize = struct.unpack("<i", self._fp.read(4))[0]
        line = self._fp.read(fieldSize)        
        self._safeAssign("AlgorithmName", line, 43)        
        # Parse version 4 algorithm parameters ---------------------------------------------------------------------
        fieldSize = struct.unpack("<i", self._fp.read(4))[0]
        line = self._fp.read(fieldSize)        
        pString = self._stringParser(line, ":", ";");  # v4 algorithm parameters      
        for i in pString:            
##            if (self.debug):    
##                print "v4 algo param: ", i[0], i[1]            
            try:
                self._safeAssign(i[0], int(i[1]), 43)
            except ValueError:
                self._safeAssign(i[0], i[1], 43)                               
        # Parse version 4 file header -----------------------------------------------------------------------------------
        margin, self.nOutlierCell, self.nMarkedCell, self.nSubgridCell = struct.unpack("<iIIi", self._fp.read(16))         
        self.info["CellMargin"] = margin        
        self.info["Outliers"] = self.nOutlierCell
        self.info["Masked"] = self.nMarkedCell
        self.info["Subgrids"] = self.nSubgridCell
        self._dataOffSet = self._fp.tell()
        #print "info: ", self.info        
    
    def _readVer4Data(self, data = None, row = None, col = None):
        """
        Read version 4 cell data.
        If variables "row" and "col" are filled in, this function will return a single value. Otherwise it returns a chunk of data.
        """
        if (len(self._fileHeader) == 0 or data is None):
            return array()       
        dataOffset = self._dataOffSet          
        
        if (data == "Intensity" or data == "StdDev" or data == "Pixel"):
            self._fp.seek(dataOffset)                                                
            if (row != None and col != None): # Random access
                fmt = "<"+"ffh"                
                sizeToRead = 4+4+2                       
                #self._fp.seek(((self.info["Rows"]-1-col)*self.info["Cols"]+row)*sizeToRead, 1) #PF: changed so .getIntensity can be accessed as y,x
                self._fp.seek((int(row)*self.info["Cols"]+int(col))*sizeToRead, 1)
                dataArray = array(struct.unpack(fmt, self._fp.read(sizeToRead)))                
                if (data == "Intensity"):
                    dataArray = dataArray[0]
                elif (data == "StdDev"):
                    dataArray = dataArray[1]
                else:
                    dataArray = dataArray[2]
            else: # Load the entire chunk.
                fmt = "<"+"ffh" * self._fileHeader[4]
                sizeToRead = self._fileHeader[4]*(4+4+2)
                dataArray = array(struct.unpack(fmt, self._fp.read(sizeToRead)))    
                dataArray = reshape(dataArray, (self.size[0], self.size[1]*3))
        elif (data == "Outlier" or "Mask"):
            n = 0            
            if (data == "Outlier"):
                self._fp.seek(dataOffset+self._fileHeader[4]*(4+4+2))
                n = self.nOutlierCell  
            else:
                self._fp.seek(dataOffset+self._fileHeader[4]*(4+4+2+2+2))
                n = self.nMarkedCell
            fmt = "<" + "hh" * n           
            dataArray = array(struct.unpack(fmt, self._fp.read(n*(2+2))))            
            dataArray = reshape(dataArray, (n, 2))            
        elif (data == "SubGrid"):
            self._fp.seek(dataOffset+self._fileHeader[4]*(4+4+2+2+2+2+2))
            fmt = "<"+"iiffffffffiiii" * self.nSubgridCell
            dataArray = array(struct.unpack(fmt, self._fp.read(self.nSubgridCell*(4*14))))            
            dataArray = reshape(dataArray, (self.nSubgridCell, 14))       
        else:
            #print "[Warning] Undefined version 4 data query: ", data
            pass                   
        return dataArray
    
    def _readGenericDataHeader(self, parent = False):
        """
        Read Command Console Version 1 file header, algorithm parameters, and the DAT header.
        """
        if (not parent):
            self.info = dict() # Initialize the dictionary       
            self.info_exception = dict() # Initialize the exception dictionary 
            self._fp.seek(10)    
        # Parse Command Console version 1 file header ---------------------------------------------------  
        readStr = self._readString() # data type identifier            
        if (not parent):
             self.info["DataTypeId"] = readStr
        readStr = self._readString() # file identifier
        if (not parent):
            self.info["FileId"] = readStr
        readStr = self._readWString() # Date Time
        if (not parent):
            self.info["FileCreationDateTime"] = readStr        
        readStr = self._readWString() # Locale
        if (not parent):
            self.info["FileLocale"] = readStr        
        parameterCount = struct.unpack(">i", self._fp.read(4))[0]
        # Parse Command Console version 1 algorithm parameters ---------------------------------------               
        while (parameterCount > 0):
            pName = self._readWString()
            pValue = self._readString()
            pType = self._readWString()            
##            if (self.debug):               
##                print "v5 algo param: ", v5_to_global(pName), pValue 
            if (v5_to_global(pName) == "PartialDatHeader"):
                self._parseDatHeader(self._decodeMime(pValue, pType))            
            else:
                self._safeAssign(v5_to_global(pName), self._decodeMime(pValue, pType), 1)                
            parameterCount -= 1 
        self.size = (self.info["Cols"], self.info["Rows"])
        self.nCell = self.size[0] * self.size[1]   
        # Parse its parents' headers recursively -------------------------------------------------------------
        parentCount = struct.unpack(">i", self._fp.read(4))[0]                           
        while (parentCount > 0):
            self._readGenericDataHeader(True)
            parentCount -= 1
        

    def _readDataGroups(self):
        """
        Read Command Console Version 1 data group header.
        The real data section is skipped.
        """
        nextGroup = self._fileHeader[3]        
        while (nextGroup):
            self._fp.seek(nextGroup)
            groupInfo = struct.unpack(">IIi", self._fp.read(12))            
            groupName = self._readWString()  
            dataOffset = groupInfo[1]
            for i in xrange(groupInfo[2]):
                self._fp.seek(dataOffset)
                dataInfo = struct.unpack(">II", self._fp.read(8))
                dataName = self._readWString() 
                parameterCount = struct.unpack(">i", self._fp.read(4))[0]
                while (parameterCount > 0):
                    pName = self._readWString()                    
                    pValue = self._readString()
                    pType = self._readWString()
                    parameterCount -= 1
                    self._safeAssign(v5_to_global(pName), self._decodeMime(pValue, pType), 1)                                                     
                dataColumns = struct.unpack(">I", self._fp.read(4))[0]    
                rowSize = 0                
                for j in xrange(dataColumns):
                    colName = self._readWString()
                    colInfo = struct.unpack(">bi", self._fp.read(5))                    
                    rowSize += colInfo[1]
                dataRows = struct.unpack(">I", self._fp.read(4))[0]              
                dataOffset = dataInfo[1]
                
                 # identify the data name
                if (dataName == "Intensity"):
                    self._intensityOffset = dataInfo[0]                     
                    self._intensityCol = dataColumns
                elif (dataName == "StdDev"):
                    self._stdDevOffset = dataInfo[0]                    
                    self._stdDevCol = dataColumns
                elif (dataName == "Pixel"):
                    self._pixCntOffset = dataInfo[0]                    
                    self._pixCntCol = dataColumns
                elif (dataName == "Outlier"):
                    self._outlierOffset = dataInfo[0]
                    self.nOutlierCell = dataRows/dataColumns
                    self._outlierCol = dataColumns
                elif (dataName == "Mask"):
                    self._mkcellOffset = dataInfo[0]
                    self.nMarkedCell = dataRows/dataColumns
                    self._mkcellCol = dataColumns
                elif (dataName == "SubGrid"):
                    self._subgridOffset = dataInfo[0]
                    self.nSubgridCell = dataRows/dataColumns
                    self._subgridCol = dataColumns    
            nextGroup = groupInfo[0]
 
    def _readVer5Data(self, data = None, row = None, col = None):
        """
        Read Command Console version 1 cell data.
        If variables "row" and "col" are filled in, this function will return a single value. Otherwise it returns a chunk of data.
        """        
        randAcc = False       
        if (row != None and col != None):
            randAcc = True
                    
        # Compute the data offset from the beginning of the file.
        n = 0
        fmt = ""
        if (data == "Intensity"):
            dataOffset = self._intensityOffset              
            if (randAcc == True):
                n = 4
                dataSize = 1
            else:
                n = self._intensityCol*4          
                dataSize = self.nCell            
            fmt = ">%df" % dataSize
        elif (data == "StdDev"):
            dataOffset = self._stdDevOffset            
            if (randAcc == True):
                n = 4
                dataSize = 1
            else:
                n = self._stdDevCol*4          
                dataSize = self.nCell     
            fmt = ">%df" % dataSize
        elif (data == "Pixel"):
            dataOffset = self._pixCntOffset
            if (randAcc == True):
                n = 2
                dataSize = 1
            else:
                n = self._pixCntCol*2          
                dataSize = self.nCell
            fmt = ">%dh" % dataSize   
        elif (data == "Outlier"):
            dataOffset = self._outlierOffset
            dataSize = self.nOutlierCell * self._outlierCol
            n = self._outlierCol*1
            fmt = ">%dh" % dataSize
        elif (data == "Mask"):
            dataOffset = self._mkcellOffset
            dataSize = self.nMarkedCell * self._mkcellCol
            n = self._mkcellCol*1
            fmt = ">%dh" % dataSize   
        elif (data == "SubGrid"):
            dataOffset = self._subgridOffset
            dataSize = self.nSubgridCell * self._subgridCol            
            n = self._subgridCol*4/9            
            fmt = ">%di" % dataSize   
        else:
            #print "Undefined query: ", data
            return
        # Pull out the real data
        self._fp.seek(dataOffset)
        # Random access a single cell
        if (randAcc):
            #self._fp.seek(((self.info["Rows"]-1-int(col))*self.info["Cols"]+int(row))*n*dataSize, 1)            
            self._fp.seek((int(row)*self.info["Cols"]+int(col))*n*dataSize, 1)            
            dataArray = array(struct.unpack(fmt, self._fp.read(n*dataSize)))            
            return dataArray[0]        
        # Load the entire chunk.    
        dataArray = array(struct.unpack(fmt, self._fp.read(n*dataSize)))        
        if (data == "Intensity" or data == "StdDev" or data == "Pixel"):            
            dataArray = reshape(dataArray, (self.size[1], self.size[0]))            
        elif (data == "Outlier"):
            dataArray = reshape(dataArray, (self.nOutlierCell, self._outlierCol))
        elif (data == "Mask"):
            dataArray = reshape(dataArray, (self.nMarkedCell, self._mkcellCol))            
        elif (data == "SubGrid"):
            dataArray = reshape(dataArray, (self.nSubgridCell, self._subgridCol))        
        return dataArray
    
    def _readString(self):
        strLen = struct.unpack(">i", self._fp.read(4))[0]
        return self._fp.read(strLen)

    def _readWString(self):
        strLen = struct.unpack(">i", self._fp.read(4))[0]
        return self._fp.read(2*strLen).decode("utf-16-be")
    
    def _decodeMime(self, value, type):
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
    
    def loadIntensity(self):
        """ 
        Pre-fetch intensity data into the memory for fast access
        This code assumes the data won't change after the CEL object creation.
        """          
        if (self.version == 1):            
            self.intensity = self._readVer5Data("Intensity")                                                
            #self.intensity = fliplr(transpose(self.intensity))
        elif (self.version == 4):
            self.intensity = self._readVer4Data("Intensity")
            self.intensity = self.intensity[:, ::3]
            #self.intensity = fliplr(transpose(self.intensity))            
        else:
            return                    
        self._intensityInMem = 1      
        self.intensity=numpy.clip(self.intensity, 1, maxVal)   
        return   
    
    
    def loadStdDev(self):
        """ 
        Pre-fetch standard deviation into the memory for fast access
        This code assumes the data won't change after the CEL object creation.
        """        
        if (self.version == 1):            
            self.stdDev = self._readVer5Data("StdDev")
            #self.stdDev = fliplr(transpose(self.stdDev))           
        elif (self.version == 4):
            self.stdDev = self._readVer4Data("StdDev")
            self.stdDev = self.stdDev[:, 1::3]
            #self.stdDev = fliplr(transpose(self.stdDev))        
        else:
            #print "Wrong format!"
            return        
        self._stdDevInMem = 1
        return
    
    def loadPixelCount(self):
        """ 
        Pre-fetch pixel count information into the memory for fast access
        This code assumes the data won't change after the CEL object creation.
        """
        if (self.version == 1):            
            self.pixelCount = self._readVer5Data("Pixel")
            #self.pixelCount = fliplr(transpose(self.pixelCount))       
        elif (self.version == 4):
            self.pixelCount = self._readVer4Data("Pixel")
            self.pixelCount = self.pixelCount[:, 2::3]
            #self.pixelCount = fliplr(transpose(self.pixelCount))        
        else:
            #print "Wrong format!"
            return        
        self._pixCntInMem = 1
        return
    
    def loadOutlierCell(self):
        """ 
        Pre-fetch outlier coordinates into the memory for fast access
        This code assumes the data won't change after the CEL object creation.
        """
        if (self.version == 1):            
            self.outlierCell = self._readVer5Data("Outlier")
        elif (self.version == 4):
            self.outlierCell = self._readVer4Data("Outlier")                    
        else:
            #print "Wrong format!"
            return        
        self._outlierInMem = 1        
        return
    
    def loadMarkedCell(self):
        """ 
        Pre-fetch marked coordinates into the memory for fast access
        This code assumes the data won't change after the CEL object creation.
        """        
        if (self.version == 1):            
            self.markedCell = self._readVer5Data("Mask")        
        elif (self.version == 4):
            self.markedCell = self._readVer4Data("Mask")
            #self._readVer4Data("Mask")
            #return                     
        else:
            #print "Wrong format!"
            return
        self._mkcellInMem = 1
        return 
    
    def loadSubgridCell(self):
        """ 
        Pre-fetch subgrid information into the memory for fast access
        This code assumes the data won't change after the CEL object creation.
        """
        if (self.version == 1):            
            self.subgridCell = self._readVer5Data("SubGrid")          
        elif (self.version == 4):
            self.subgridCell = self._readVer4Data("SubGrid")                    
        else:
            #print "Wrong format!"
            return
        self._subgridInMem = 1
        return         
   
    def __getattr__( self, attr ):   
        """ 
        The overloading of the dot operator. It's called only when ".xxxx" doesn't exist. 
        """     
        if attr.startswith("__") and attr.endswith("__"):
            raise AttributeError, attr
        if (attr == "intensity"):            
            if (self._intensityInMem == 0):
##                if (self.debug):
##                    print "load intensity"
                self.loadIntensity()                
            return self.intensity            
        elif (attr == "stdDev"):            
            if (self._stdDevInMem == 0):
##                if (self.debug):
##                    print "load standard deviation"
                self.loadStdDev()                
            return self.stdDev
        elif (attr == "pixelCount"):            
            if (self._pixCntInMem == 0):
##                if (self.debug):
##                    print "load pixel count"
                self.loadPixelCount()                
            return self.pixelCount
        elif (attr == "outlierCell"):            
            if (self._outlierInMem == 0):
##                if (self.debug):
##                    print "load outlier cells"
                self.loadOutlierCell()                
            return self.outlierCell        
        elif (attr == "markedCell"):            
            if (self._mkcellInMem == 0):
##                if (self.debug):
##                    print "load marked cells"
                self.loadMarkedCell()                            
            return self.markedCell
        elif (attr == "subgridCell"):            
            if (self._subgridInMem == 0):
##                if (self.debug):
##                    print "load subgrid cells"
                self.loadSubgridCell()                 
            return self.subgridCell             
        
    def getIntensity(self, row, col):    
        """ 
        Random access the intensity value at cell(row, col). 
        """ 
        if (self.version == 1):
            if (self._intensityInMem == 0):
                return self._readVer5Data("Intensity", row, col)
            else:
                return self.intensity[row, col]
        elif (self.version == 4):     
            if (self._intensityInMem == 0):
                return self._readVer4Data("Intensity", row, col)
            else:
                return numpy.clip(self.intensity[row, col], 1, maxVal)
    
    def getStdDev(self, row, col):
        """ 
        Random access the standard deviation value at cell(row, col). 
        """     
        if (self.version == 1):
            if (self._stdDevInMem == 0):
                return self._readVer5Data("StdDev", row, col)
            else:
                return self.stdDev[row, col]
        elif (self.version == 4):     
            if (self._stdDevInMem == 0):
                return self._readVer4Data("StdDev", row, col)
            else:
                return self.stdDev[row, col]
            
    def getPixelCount(self, row, col):
        """ 
        Random access the pixel count value at cell(row, col). 
        """     
        if (self.version == 1):
            if (self._pixCntInMem == 0):
                return self._readVer5Data("Pixel", row, col)
            else:
                return self.pixelCount[row, col]
        elif (self.version == 4):     
            if (self._pixCntInMem == 0):
                return self._readVer4Data("Pixel", row, col)
            else:
                return self.pixelCount[row, col]
        
        
# API Begin ==========================================
def readCEL(celFileName):    
    """ 
    Read a CEL file. 
    """
    return CEL(celFileName)
# API End ===========================================