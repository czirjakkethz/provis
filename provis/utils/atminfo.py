import csv
import os
from numpy.core.fromnumeric import size
from numpy.core.numeric import outer
import provis

def import_atm_size_info(path, vw=0):
    """
    Funtion to load atomic radii from atmsize.csv, to be used in main
    
    :param name: dir - Directory path to base directory that contains data/data/atmsize.csv and data/data/atmvw.csv (downloadable from: https://github.com/czirjakkethz/provis/tree/main/data/data)
    :param type: str
    :param name: vw - Option to return vanderwaals radius.
    :param type: bool, optional
    
    :return: dict - dictionary of atomic radius by atom name
    :return: dict - dictionary of color by atom name
    :return: dict - !optional! return dictionary of vw rdius by atom name
    """
    
    # _path = os.path.dirname(provis.__file__)
    # _path points to provis/provis, we only want provis/
    # path = path[0:-6]
    
    
    filename_size = path + "data/data/atmsize.csv"
    filename_vw = path + "data/data/atmvw.csv"
    
    # opening the file using "with" 
    # statement

    firstline = []
    sizedict = {}
    with open(filename_size, newline='') as data:
        firstline = data.readline().rstrip().split(",")
        secondline = data.readline().rstrip().split(",")
        l = len(firstline)
        for i in range(0,l):
            sizedict[firstline[i]] = float(secondline[i])
    sizedict = {'H': 0.53, 
                'HE': 0.31, 
                'LI': 1.67, 
                'BE': 1.12, 
                'B': 0.87, 
                'C': 0.67, 
                'N': 0.56, 
                'O': 0.48, 
                'F': 0.42,
                'NE': 0.38,
                'NA': 1.9,
                'MG': 1.45,
                'AL': 1.18,
                'SI': 1.11,
                'P': 0.98,
                'S': 0.88,
                'CL': 0.79,
                'AR': 0.71,
                'K': 2.43,
                'CA': 1.94,
                'SC': 1.84,
                'TI': 1.76,
                'V': 1.71,
                'CR': 1.66,
                'MN': 1.61,
                'FE': 1.56,
                'CO': 1.52,
                'NI': 1.49,
                'CU': 1.45,
                'ZN': 1.42,
                'GA': 1.36,
                'GE': 1.25,
                'AS': 1.14,
                'SE': 1.03,
                'BR': 0.94,
                'KR': 0.88,
                'RB': 2.65,
                'SR': 2.19,
                'Y': 2.12,
                'ZR': 2.06,
                'NB': 1.98,
                'MO': 1.9,
                'TC': 1.83,
                'RU': 1.78,
                'RH': 1.73,
                'PD': 1.69,
                'AG': 1.65,
                'CD': 1.61,
                'IN': 1.56,
                'SN': 1.45,
                'SB': 1.33,
                'TE': 1.23,
                'I': 1.15,
                'XE': 1.08,
                'CS': 2.98,
                'BA': 2.53,
                'LA': 2.26,
                'CE': 2.1,
                'PR': 2.47,
                'ND': 2.06,
                'PM': 2.05,
                'SM': 2.38,
                'EU': 2.31,
                'GD': 2.33,
                'TB': 2.25,
                'DY': 2.28,
                'HO': 2.26,
                'ER': 2.26,
                'TM': 2.22,
                'YB': 2.22,
                'LU': 2.17,
                'HF': 2.08,
                'TA': 2.0,
                'W': 1.93,
                'RE': 1.88,
                'OS': 1.85,
                'IR': 1.8,
                'PT': 1.77,
                'AU': 1.74,
                'HG': 1.71,
                'TL': 1.56,
                'PB': 1.54,
                'BI': 1.43,
                'PO': 1.35,
                'AT': 1.27,
                'RN': 1.2,
                'FR': 0.0,
                'RA': 0.0,
                'AC': 0.0,
                'TH': 0.0,
                'PA': 0.0,
                'U': 0.0,
                'NP': 0.0,
                'PU': 0.0,
                'AM': 0.0,
                'CM': 0.0,
                'BK': 0.0,
                'CF': 0.0,
                'ES': 0.0,
                'FM': 0.0,
                'MD': 0.0,
                'NO': 0.0,
                'LR': 0.0,
                'RF': 0.0,
                'DB': 0.0,
                'SG': 0.0,
                'BH': 0.0,
                'HS': 0.0,
                'MT': 0.0,
                'DS': 0.0,
                'RG': 0.0,
                'CN': 0.0,
                'NH': 0.0,
                'FL': 0.0,
                'MC': 0.0,
                'LV': 0.0,
                'TS': 0.0,
                'OG': 0.0}
    
    
    color = {  'H': '#FFFFFF',
                        'C': '#333333',
                        'O': '#EA2128',
                        'N': '#2233FF',
                        'S': '#FFDC61',
                        'F': '#1FF01F',
                        'CL': '#1FF01F',
                        'BR': '#992200',
                        'I': '#6600BB',
                        'HE': '#00FFFF',
                        'NE': '#00FFFF',
                        'AR': '#00FFFF',
                        'KR': '#00FFFF',
                        'XE': '#00FFFF',
                        'P': '#FF9900',
                        'B': '#FFAA77',
                        'LI': '#7700FF',
                        'NA': '#7700FF',
                        'K': '#7700FF',
                        'RB': '#7700FF',
                        'CS': '#7700FF',
                        'FR': '#7700FF',
                        'BE': '#007700',
                        'MG': '#007700',
                        'CA': '#007700',
                        'SR': '#007700',
                        'BA': '#007700',
                        'RA': '#007700',
                        'BE': '#007700',
                        'TI': '#999999',
                        'FE': '#DD7700'
                    }
    
    if vw:
        vwdict = {'H': 1.2, 
                'HE': 1.4, 
                'LI': 1.82,
                'BE': 1.53, 
                'B': 1.92,
                'C': 1.7,
                'N': 1.55,
                'O': 1.52,
                'F': 1.47,
                'NE': 1.54,
                'NA': 2.27,
                'MG': 1.73,
                'AL': 1.84,
                'SI': 2.1,
                'P': 1.8,
                'S': 1.8,
                'CL': 1.75,
                'AR': 1.88,
                'K': 2.75,
                'CA': 2.31,
                'SC': 2.11,
                'TI': 0.0,
                'V': 0.0,
                'CR': 0.0,
                'MN': 0.0,
                'FE': 0.0,
                'CO': 0.0,
                'NI': 1.63,
                'CU': 1.4,
                'ZN': 1.39,
                'GA': 1.87,
                'GE': 2.11,
                'AS': 1.85,
                'SE': 1.9,
                'BR': 1.85,
                'KR': 2.02,
                'RB': 3.03,
                'SR': 2.49,
                'Y': 0.0,
                'ZR': 0.0,
                'NB': 0.0,
                'MO': 0.0,
                'TC': 0.0,
                'RU': 0.0,
                'RH': 0.0,
                'PD': 1.63,
                'AG': 1.72,
                'CD': 1.58,
                'IN': 1.93,
                'SN': 2.17,
                'SB': 2.06,
                'TE': 2.06,
                'I': 1.98,
                'XE': 2.16,
                'CS': 3.43,
                'BA': 2.68,
                'LA': 0.0,
                'CE': 0.0,
                'PR': 0.0,
                'ND': 0.0,
                'PM': 0.0,
                'SM': 0.0,
                'EU': 0.0,
                'GD': 0.0,
                'TB': 0.0,
                'DY': 0.0,
                'HO': 0.0,
                'ER': 0.0,
                'TM': 0.0,
                'YB': 0.0,
                'LU': 0.0,
                'HF': 0.0,
                'TA': 0.0,
                'W': 0.0,
                'RE': 0.0,
                'OS': 0.0,
                'IR': 0.0,
                'PT': 1.75,
                'AU': 1.66,
                'HG': 1.55,
                'TL': 1.96,
                'PB': 2.02,
                'BI': 2.07,
                'PO': 1.97,
                'AT': 2.02,
                'RN': 2.2,
                'FR': 3.48,
                'RA': 2.83,
                'AC': 0.0,
                'TH': 0.0,
                'PA': 0.0,
                'U': 1.86,
                'NP': 0.0,
                'PU': 0.0,
                'AM': 0.0,
                'CM': 0.0,
                'BK': 0.0,
                'CF': 0.0,
                'ES': 0.0,
                'FM': 0.0,
                'MD': 0.0,
                'NO': 0.0,
                'LR': 0.0,
                'RF': 0.0,
                'DB': 0.0,
                'SG': 0.0,
                'BH': 0.0,
                'HS': 0.0,
                'MT': 0.0,
                'DS': 0.0,
                'RG': 0.0,
                'CN': 0.0,
                'NH': 0.0,
                'FL': 0.0,
                'MC': 0.0,
                'LV': 0.0,
                'TS': 0.0,
                'OG': 0.0}

        return sizedict, color, vwdict
        
    return sizedict, color

  
def import_atm_mass_info():
    """
    Funtion to load atomic mass from atmmass.csv, to be used in main
    
    :return: dict - dictionary of atomic mass by atom name
    """
    path = os.path.dirname(provis.__file__)
    # path points to provis/provis, we only want provis/
    base_path = path[0:-6]
    
    filename_mass = base_path + "data/data/atmmass.csv"
    
    firstline = []
    massdict = {}
    with open(filename_mass, newline='') as data:
        firstline = data.readline().rstrip().split(",")
        secondline = data.readline().rstrip().split(",")
        l = len(firstline)
        for i in range(0,l):
            massdict[firstline[i].upper()] = float(secondline[i])

    return massdict
