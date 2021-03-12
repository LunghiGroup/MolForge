#import pymatgen as mg
import sys

# Count the arguments
arguments = len(sys.argv) - 1

# Output argument-wise
position = 1
while (arguments >= position):
    var = sys.argv[position]

    words = var.split(".")

    if words[len(words)-1] == "cif" :
        print ("Reading Cif file:",var) 
    else :
        print ("This is not a cif file, skipping file")
        exit

#    structure = mg.Structure.from_file(cifname".cif", primitive=False)
#    structure.to(filename=cifname".xyz")

    position = position + 1

