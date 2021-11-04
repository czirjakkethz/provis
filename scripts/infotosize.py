# python script to convert csv formatted atomic radii data from wikipedia to dict
# stores atomic names and radii in two lines of a csv file
import csv

# from https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_%28data_page%29
filename ="../data/atmsizeinfo.csv"
  
# opening the file using "with" 
# statement
atom_info_list = []
firstline = []
with open(filename, 'r') as data:
    firstline = data.readline().rstrip().split(",")
    data.seek(0, 0)
    for line in csv.DictReader(data):
        atom_info_list.append(line)

numatoms = len(atom_info_list)

options = [4, 5]

# for loop to run for atomic size and van der waals radius
for k in options:
    val_list = []
    key_list = []
    for i in range(numatoms):
        key_list.append(atom_info_list[i][firstline[1]].upper())
        if atom_info_list[i][firstline[k]].isdigit():
            val_list.append(float(atom_info_list[i][firstline[k]])/100.)
        else:
            val_list.append(0.0)


    # open file for writing, "w" is writing
    w = csv.writer(open("../data/atmsize.csv", "w", newline=''))
    if k == 4:
        w = csv.writer(open("../data/atmsize.csv", "w", newline=''))
    elif k == 5:
        w = csv.writer(open("../data/atmvw.csv", "w", newline=''))

    w.writerow(key_list)
    w.writerow(val_list)

