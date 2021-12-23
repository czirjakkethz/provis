import csv
import os
import provis

def import_atm_size_info(vw=0):
    """
    Funtion to load atomic radii from atmsize.csv, to be used in main
    
    :param name: vw - Option to return vanderwaals radius.
    :param type: bool, optional
    
    :return: dict - dictionary of atomic radius by atom name
    :return: dict - dictionary of color by atom name
    :return: dict - !optional! return dictionary of vw rdius by atom name
    """
    
    path = os.path.dirname(provis.__file__)
    # path points to provis/provis, we only want provis/
    outer_path = path[0:-6]
    
    filename_size = outer_path + "data/data/atmsize.csv"
    filename_vw = outer_path + "data/data/atmvw.csv"
    
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
        firstline = []
        vwdict = {}
        with open(filename_vw, newline='') as data:
            firstline = data.readline().rstrip().split(",")
            secondline = data.readline().rstrip().split(",")
            l = len(firstline)
            for i in range(0,l):
                vwdict[firstline[i]] = float(secondline[i])
        return sizedict, color, vwdict
        
    return sizedict, color

# original dictionary used
# atoms_size_dict = {'C': pv.Sphere(radius=0.67, phi_resolution=phi_res, theta_resolution=theta_res),
#                         'O': pv.Sphere(radius=0.48, phi_resolution=phi_res, theta_resolution=theta_res),
#                         'N': pv.Sphere(radius=0.56, phi_resolution=phi_res, theta_resolution=theta_res),
#                         'S': pv.Sphere(radius=0.88, phi_resolution=phi_res, theta_resolution=theta_res),
#                         'CA': pv.Sphere(radius=1.94, phi_resolution=phi_res, theta_resolution=theta_res)
#                         }
# def main():
#     dictt = import_atm_size_info()
#     print(dictt)

# if __name__ == "__main__":
#     main()

  
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
