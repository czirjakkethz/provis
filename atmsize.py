import csv
  
def import_atm_info(vw):
    """
    Funtion to load atomic radii from atmsizeinfo.csv, to be used in main
    :param name: vw - option to return vanderwaals radius
    :param type: bool - returns if true
    :return: dict - dictionary of atomic radius by atom name
    :return: dict - dictionary of color by atom name
    """
    filename_size ="atmsizeinfo.csv"
    filename_vw ="atmvwinfo.csv"
    
    # opening the file using "with" 
    # statement

    firstline = []
    newdict = {}
    with open(filename_size, newline='') as data:
        firstline = data.readline().rstrip().split(",")
        secondline = data.readline().rstrip().split(",")
        l = len(firstline)
        for i in range(0,l):
            newdict[firstline[i]] = float(secondline[i])

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
        return newdict, color, vwdict


    return newdict, color

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