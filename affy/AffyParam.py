#!/usr/bin/python
import sys, string

def v3_to_global(param):
    """
    Define the mapping between version 3 and our definitions.
    """
    if (param == "GridCornerUL"):
        return "GridUL"
    elif (param == "GridCornerUR"):
        return "GridUR"
    elif (param == "GridCornerLR"):
        return "GridLR"
    elif (param == "GridCornerLL"):
        return "GridLL"    
    elif (param == "Axis-invertX"):
        return "AxisInvertX"
    elif (param == "swapXY"):
        return "SwapXY"
    elif (param == "Algorithm"):
        return "AlgorithmName"
    elif (param == "AlgorithmParameters"):
        return "V3AlgorithmParameters"
    else:
        return param        
    
def v5_to_global(param):
    """
    Define the mapping between version 3 and our definitions.
    """    
    ss = param.split('-')    
    param = ""
    for i in ss:
        if (len(ss) < 4 and (i == "affymetrix" or i == "param" or i == "cel")):
            continue
        elif (len(ss) >= 4 and (i == "affymetrix" or i == "param" or i == "cel" or i == "algorithm")):
            continue
        elif (i == "scanparameter"):
            param = "Scan"
            continue        
        elif (i == "ArrayParameter"):
            param = "Array"
            continue
        # Turn the unicode parameter name into a capitalized string
        if (i[0].isupper() == False):
            tmp_str = list(str(i))
            tmp_str[0] = str(tmp_str[0]).upper()
            i = "".join(tmp_str)
        else:             
            i = str(i)
        param += i
    return param

        
    
    
    