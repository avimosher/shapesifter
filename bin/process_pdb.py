#!/usr/bin/python

from Bio.PDB import *
from subprocess import Popen, PIPE, STDOUT

pdbl=PDBList()
pdb_filename=pdbl.retrieve_pdb_file('1FAT')

with open(pdb_filename,"r") as pdb_handle:
    pdb_content=pdb_handle.read()
    process=Popen(['./convert_structure'],stdout=PIPE,stdin=PIPE)
    stdout_data=process.communicate(input=pdb_content)
    print stdout_data[0]
