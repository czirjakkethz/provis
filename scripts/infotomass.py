# python script to convert csv formatted atomic mass data to dict
# stores atomic names and mass in two lines of a csv file
import csv

# from https://www.angelo.edu/faculty/kboudrea/periodic/structure_mass.htm
filename ="../data/atmmassinfo.csv"
  
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

val_list = []
key_list = []
for i in range(numatoms):
    key_list.append(atom_info_list[i][firstline[1]].upper())
    
    try:
        val_list.append(float(atom_info_list[i][firstline[3]]))
    except ValueError:
        val_list.append(0.0)

# open file for writing, "w" is writing
w = csv.writer(open("../data/atmmass.csv", "w", newline=''))
w.writerow(key_list)
w.writerow(val_list)

