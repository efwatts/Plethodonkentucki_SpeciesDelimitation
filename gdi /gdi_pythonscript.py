# Script 3: Calculating GDI matrix directly from BPP A00-type analysis output
####
#python3

import math
import re
import sys
BPPresults = sys.argv[1]
outfile = sys.argv[2]
# Expand this to potentially include several input files and take GDI scores from all of them

# Required module for exponent e


# Regular expression pattern that finds the species name + number in the file
#pattern = re.compile(r"[0-9]+\s[A-Z]+")
pattern = re.compile(r"[0-9]+ [A-Z]\S*")

namedict = {}
thetameandict = {}
taumeandict = {}
thetamediandict = {}
taumediandict = {}


# Start by reading the population matrix from your output file and putting the species names in a dictionary
# For this code to work it's important you follow the naming convention of species in BPP we propose - e.g. when 26 species or less, simply name them A through Z.
# If you want to include more than 26 species, make sure each species name still starts with a capital letter, and add a number to each species' name to increase
# So A1 through Z1, and species number 27 could be 27, or number them A01, A02, A03, through A99 if you'd like
with open(BPPresults) as f1:
    for a in f1:
        a = a.strip()
        a = re.findall(pattern, a)
    if len(a) > 0:
        a = str(a)
        a = a[2:-2]
        a = a.split()
        namedict[int(a[0])] = a[1]

# the dictionary ##namedict## now contains an entry [0] that is the number used by BPP and entry [1] that is the species name given in the phylogeny and Imap file
# e.g.
# 1 = A
# 2 = B
# 18 = ABCFDEFGHIJKLMNOPQ

# Next we use the dictionary to output mean theta and tau values per species or species groups in excell file rows
# Here we calculate both the mean values of Theta and Tau, alternatively you can use the median values, for that you should use the thetamediandict and taumediandict later on in the code, for which a reminder will follow


# Thetas go in thetameandict
with open(BPPresults) as f2:
    for a in f2:
        if a[0:6] == "median":
            medianlist = a.split()
            limit = len(namedict) + 1
            medianlist = medianlist[:limit]
        for b in namedict:
            c = int(b)
        thetamediandict[c] = float(medianlist[c])
# Figure out Tau's, which we only want for composite species (more than 1 species in a group)
numbertau = 0
for a in namedict:
    a = namedict[a]
    if len(a) > 1:
        counter = 0
        indspecies = "ABCDE"
    for b in a:
        if b in indspecies:
            counter += 1
    if counter > 1:
        numbertau += 1

# So we want to fill our taumeandict with only those last columns that represent tau values
with open(BPPresults) as f2:
    for a in f2:
        if a[0:6] == "median":
            medianlist = a.split()
            medianlist = medianlist[-numbertau - 1:-1]
            taustart = len(namedict) - numbertau + 1
        for a in medianlist:
            taumediandict[taustart] = float(a)
            taustart += 1
# Veryifying dictionaries
print(thetamediandict)
print(taumediandict)

# Thetas go in thetameandict
with open(BPPresults) as f2:
    for a in f2:
        if a[0:4] == "mean":
            meanlist = a.split()
            limit = len(namedict) + 1
            meanlist = meanlist[:limit]
        for b in namedict:
            c = int(b)
            thetameandict[c] = float(meanlist[c])


# So we want to fill our taumeandict with only those last columns that represent tau values
with open(BPPresults) as f2:
    for a in f2:
        if a[0:4] == "mean":
            meanlist = a.split()
            meanlist = meanlist[-numbertau - 1:-1]
            taustart = len(namedict) - numbertau + 1
        for a in meanlist:
            taumeandict[taustart] = float(a)
            taustart += 1
 # Veryifying dictionaries
 # print(thetameandict)
 # print(taumeandict)

# Here we make a .csv file that has each species or their most recent common ancestor, their Theta, Tau and respective GDI values in a matrix
# REMEMBER to name your output file in the command line, e.g. 'GDI.csv'
# Write GDI values for each comparison to new file
with open(outfile, 'w') as f3:
    f3.write("Species,")
    counter = 1
    while counter <= len(namedict):
        entry = namedict[counter]
        f3.write(entry + ",")
        counter += 1
        f3.write("\n")

# This line writes the MeanThetas to the file, replace all occurences of thetameandict here by thetamediandict, or disable this section of code by placing # in front of each line and remove the # in front of each line in the next block
# with open(outfile, 'a') as f3:
#	f3.write("MeanTheta,")
#	counter = 1
#	while counter <= len(thetameandict):
#    	value = thetameandict[counter]
#    	f3.write(str(value) + ",")
#    	counter += 1
#	f3.write("\n")
#	f3.write("MeanTau,")
#	counter = 1
#	while counter <= len(thetameandict) - numbertau:
#    	f3.write(",")
#    	counter += 1
#	while counter <= len(thetameandict):
#    	value = taumeandict[counter]
#    	f3.write(str(value) + ",")
#    	counter += 1
#	f3.write("\n")


with open(outfile, 'a') as f3:
    f3.write("MedianTheta,")
    counter = 1
    while counter <= len(thetamediandict):
        value = thetamediandict[counter]
        f3.write(str(value) + ",")
        counter += 1
        f3.write("\n")
        f3.write("MedianTau,")
        counter = 1
    while counter <= len(thetamediandict) - numbertau:
        f3.write(",")
        counter += 1
    while counter <= len(thetamediandict):
        value = taumediandict[counter]
        f3.write(str(value) + ",")
        counter += 1
        f3.write("\n")


# Before we can calculate GDI, we need to make a dictionary that contains the most recent common ancestor for each speciesi pair:

name_key_list = list(namedict.keys())
name_value_list = list(namedict.values())


def RCA(a, b):
    fail = "N/A"
    templist = []
    for c in namedict:
        if a != b and a in namedict[c] and b in namedict[c] and len(namedict[c]) > 1:
            templist.append(namedict[c])
        sorted_list = sorted(templist, key=len)
    if len(sorted_list) > 0:
        return name_key_list[name_value_list.index(sorted_list[0])]
    else:
        return fail


with open(outfile, 'a') as f3:
    for a in namedict:
        f3.write(str(namedict[a]) + ',')
        counter = 1
    while counter <= len(namedict):
        #f3.write('GDI' + str(namedict[a]) + namedict[counter] + ',')
        ancestor = RCA(namedict[a], namedict[counter])
        if ancestor != "N/A":
            #GDI = 1 - (math.e**(-2 * taumeandict[ancestor] / thetameandict[a]))
            # ALTERNATIVE - to use the median values for tau and theta here, enable the next line by removing the # and add a # at the start of the previous line (containing taumeandict)
            GDI = 1 - (math.e**(-2 * taumediandict[ancestor] / thetamediandict[a]))
        #functionpart = taumeandict[namedict[counter]] / thetameandict[a]
        #GDI = (math.e**(-2))
            f3.write(str(GDI) + ",")
        #f3.write(RCA(namedict[a], namedict[counter]) + ',')
        else:
            f3.write(ancestor + ",")
            counter += 1
    f3.write("\n")


# 1-(math.e**(-2*meantaudict[counter]/meanthetadict[a]))
# 1-EXP(-2*Tau/Theta)
####
