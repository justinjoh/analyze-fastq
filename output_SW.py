import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import re

''' These are values that would be inconvenient to change at the command line '''
# The entire Cpf1(Cas12)-Repeat sequence. Terminates right before the SWSWSWS
cas12RepeatSequence = "ACAATTTCTACTGTTGTAGAT"
# This is the CRISPR-gate Promoter sequence, but terminated right before the SWSWSWS
CRISPRgatePromoter_first = "ttgacaacctcgtttg"


def readfastq(filename):
    ''' Returns all data from one fastq file'''
    f = open(filename)
    entirefastq = f.read()
    f.close()
    return entirefastq


def countlines(filename):
    '''just returns number of lines in data'''
    with open(filename) as f:
        num = len(f.readlines())
        return(num)


def createheatmap(fileR1, fileR2, maxlinenum):
    ''' Central function, return numpy "heatmaps"
     Ensures that files are of same length'''
    numlines1 = countlines(fileR1); numlines2 = countlines(fileR2)
    assert numlines1 == numlines2, 'These files are different lengths'
    vals1, numFaultyCpf1 = getvals(fileR1, maxlinenum, isR2=False)
    vals2, numFaultycrispr_gate = getvals(fileR2, maxlinenum, isR2=True)
    try:
        print("Number of faulty Cpf1/Cas12 repeat: " + str(numFaultyCpf1))
        print("Number of faulty crispr-gate promoter: " + str(numFaultycrispr_gate))
    except Exception as e:
        print(e)
        pass

    # Create histogram, heatmap for R1 and show histogram for R1
    heatmap1 = np.ones((1, 128), dtype=np.uint8)
    declist1 = []
    for i in range(len(vals1)):
        if vals1[i] != 'null':
            bin1 = vals1[i]
            dec1 = int(bin1, 2)
            declist1.append(dec1)
            heatmap1[0][dec1] += 1
        else:
            declist1.append('null')
    """ plt.hist(declist1, bins=128)
    plt.title(fileR1)
    plt.show()
"""
    # Create histogram, heatmap for R2 and show histogram for R2
    heatmap2 = np.zeros((1, 128), dtype=np.uint8)
    declist2 = []
    for i in range(len(vals2)):
        if (vals2[i] != 'null'):
            bin2 = vals2[i]
            dec2 = int(bin2, 2)
            declist2.append(dec2)
            heatmap2[0][dec2] += 1
        else:
            declist2.append('null')

    """plt.hist(declist2, bins=128)
    plt.title(fileR2)
    plt.show()
"""

    # Create 2d heatmap
    heatmap_2d = np.zeros((128, 128))
    try:
        for i in range(10000):
            for j in range(10000):
                r1val = declist1[i]
                r2val = declist2[j]
                if not(r1val == 'null' or r2val == 'null'):
                    heatmap_2d[r1val][r2val] += 1
        '''fig = plt.figure()
        fig.suptitle('heatmap')
        ax = fig.add_subplot(111)
        ax.set_title(fileR1)
        ax.set_ylabel(fileR2)
        ax.plot(heatmap_2d)'''
        heatmap_2d = np.log(heatmap_2d + .01)
        plt.imshow(heatmap_2d)
        plt.gca().invert_yaxis()
        plt.title('heatmap_2d')
        plt.xlabel(fileR2)
        plt.ylabel(fileR1)
        print('Displaying heatmap_2d; close it to exit program')
        plt.show()

        return heatmap1, heatmap2
    except Exception as e:
        print('*** ' + str(e))
        print('*** If got \"list index out of range\" error, is likely that no matches found')
        return

def reversecomplement(nucstring):
    '''Helper method for getvals: input is nucleotide sequence, returns its reverse complement'''
    bases_dictionary = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    revcomp = ''
    try:
        for pos in reversed(xrange(1, len(nucstring))):
            revcomp = revcomp + str((bases_dictionary[nucstring[pos-1].lower()]))
        return revcomp
    except Exception as e:
        # print (e)
        for pos in reversed(range(1, len(nucstring))):
            revcomp = revcomp + str((bases_dictionary[nucstring[pos-1].lower()]))
        return revcomp

def getvals(fastqfilename, maxlinenum, isR2):
    # TODO replace "constants" with variables that can be changed by user
    '''# Helper function for createHeatmap.
    Each possible SWSWSWS assigned binary location on heatmap
    Convention that can be tinkered with later: 'c' and 't' are binary 1
    Returns a list of 7-digit binary strings'''
    vals = []
    numFaultycrispr_gate = 0
    numFaultyCpf1 = 0
    with open(fastqfilename) as f:
        linenum = 0  # just to keep track of line
        if not isR2:  # this is the R1 case
            for line in f:
                if ((linenum-1) % 4) == 0 and linenum < maxlinenum:  # if is line index 1, 5, etc
                    try:  # Using this because str.find kept returning incorrect index
                        index1 = line.index(cas12RepeatSequence) + len(cas12RepeatSequence)  # If not found, go to except
                        index2 = index1 + 7

                        # create a boolean list, i.e. which positions match with ctctctc
                        heatIndexList = [line[i+index1].lower() == ['c', 't', 'c', 't', 'c', 't', 'c'][i] for i in range(index2-index1)]
                        heatIndex = np.reshape(heatIndexList, 7).astype(int).tolist()
                        # reform this boolean list into a string, effectively a binary identifier
                        heatIndex = ''.join([str(i) for i in heatIndex])
                        vals.append(heatIndex)

                    except Exception as e:  # If the Cpf1/Cas12-repeat does not show up in sequence
                        vals.append('null')
                        print('did not find matching Cpf1/Cas12-repeat: line ' + str(linenum))
                        numFaultyCpf1 += 1
                linenum += 1
        else:  # This is the R2 casekr
            linenum = 0
            for line in f:
                if ((linenum-1) % 4) == 0 and linenum < maxlinenum:  # if is line index 1, 5, etc
                    line = reversecomplement(line)
                    try:  # Using this because str.find kept returning incorrect index
                        index1 = line.index(CRISPRgatePromoter_first) + len(CRISPRgatePromoter_first)  # If not found, to except
                        index2 = index1 + 7
                        heatIndexList = [line[i + index1].lower() != ['c', 't', 'c', 't', 'c', 't', 'c'][i] for i in
                                         range(index2 - index1)]
                        heatIndex = np.reshape(heatIndexList, 7).astype(int).tolist()
                        heatIndex = ''.join([str(i) for i in heatIndex])
                        vals.append(heatIndex)
                    except Exception as e:  # If the crispr-gate promoter does not show up in sequence
                        # TODO need to append some sort of placeholder to vals
                        print('did not find matching crispr-gate promoter: line ' + str(linenum))
                        numFaultycrispr_gate += 1
                linenum += 1
    if numFaultyCpf1:
        return vals, numFaultyCpf1
    else:
        return vals, numFaultycrispr_gate

def getFastqFileList():
    '''Returns a list of all fastq filenames in current directory
    Known uses: in main'''
    fastqFileList = []
    for file in os.listdir("./"):
        if file.endswith(".fastq"):
            fastqFileList.append(file)
    return fastqFileList


if __name__ == "__main__":
    assert(len(sys.argv) == 1 or len(sys.argv) == 3), "must specify 0 or 2 fastq files"
    # If no filenames specified in cmdline call, run manual selection
    if len(sys.argv) == 1:
        print("Manual selection: ")
        fastqFileList = getFastqFileList()
        counter = 0
        for filename in fastqFileList:
            print("[" + str(counter) + "] " + filename + " [" + str(counter) + "]")
            counter += 1

        # Prompt user to pick the first file
        needToSelect = True
        while needToSelect:
            fileNum1 = int(input("Choose R1 file: "))
            if type(fileNum1) == int and 0 <= fileNum1 < len(fastqFileList):
                needToSelect = False
                fileR1 = fastqFileList[fileNum1]
                print("Using " + fileR1 + " as fileR1")
            else:
                print("Invalid input: choose number from list of filenames")

        # Prompt user to pick the second file
        needToSelect = True
        while needToSelect:
            fileNum2 = int(input("Choose R2 file: "))
            if type(fileNum2) == int and 0 <= fileNum2 < len(fastqFileList) and fileNum2 != fileNum1:
                needToSelect = False
                fileR2 = fastqFileList[fileNum2]
                print("Using " + fileR2 + " as fileR2")

        # Do the usual stuff
        print("Creating heatmaps")
        test_heatmap1, test_heatmap2 = createheatmap(fileR1, fileR2, maxlinenum=999999)
        print('Not printing single heatmaps')
        '''plt.imshow(test_heatmap1)
        plt.show()
        plt.imshow(test_heatmap2)
        plt.show()'''

    elif len(sys.argv) == 2:
        print("wrong # of cmdline args: specify either 0 or 2 fastq files")
    elif len(sys.argv) == 3:
        fileR1 = sys.argv[1]
        fileR2 = sys.argv[2]
        # do something
    # if there are 3 or more,
        # prompt for user input
    # if there are exactly two,
        # do the normal stuff
        # fileR1 = '8551_9805_58262_B87LH_SW-20p-Sucrose_TGACCA_R1.fastq'
        # fileR2 = '8551_9805_58262_B87LH_SW-20p-Sucrose_TGACCA_R2.fastq'
        test_heatmap1, test_heatmap2 = createheatmap(fileR1, fileR2, maxlinenum=99999)
        plt.imshow(test_heatmap1)
        plt.show()
        plt.imshow(test_heatmap2)
        plt.show()
