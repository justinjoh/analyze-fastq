import numpy as np
import matplotlib.pyplot as plt

def readfastq(filename):
    # returns all data from one fastq file
    f = open( filename )
    entirefastq = f.read()
    f.close()
    return entirefastq


def countlines(filename):
    #just returns number of lines in data
    with open(filename) as f:
        num = len(f.readlines())
        return(num)


def createHeatmap(fileR1, fileR2):
    # Main function: return numpy "heatmaps"
    # Ensure that files are of same length
    numlines1 = countlines(fileR1); numlines2 = countlines(fileR2)
    assert numlines1 == numlines2, "Files of different lengths"
    vals1 = getvals(fileR1, isR2 = False)
    vals2 = getvals(fileR2, isR2 = True)

    heatmap1 = np.ones((1, 128), dtype=np.uint8)

    declist1 = []
    for i in range(len(vals1)):
        bin1 = vals1[i]
        dec1 = int(bin1, 2)
        declist1.append(dec1)
        heatmap1[0][dec1] += 1
    plt.hist(declist1, bins=128)
    plt.show()

    heatmap2 = np.zeros((1, 128), dtype=np.uint8)
    declist2 = []
    for i in range(len(vals2)):
        bin2 = vals2[i]
        dec2 = int(bin2, 2)
        declist2.append(dec2)
        heatmap2[0][dec2] += 1
    plt.hist(declist2, bins=128)
    plt.show()
    return heatmap1, heatmap2

def reversecomplement(nucstring):
    dict = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    revcomp = ''
    for pos in reversed(xrange(1, len(nucstring))):
        revcomp = revcomp + str((dict[nucstring[pos-1].lower()]))
    return(revcomp)

def getvals(filename, isR2):
    '''# Helper function for createHeatmap.
    Each possible SWSWSWS assigned binary location on heatmap
    Convention that can be tinkered with later: 'c' and 't' are binary 1
    Returns a list of 7-digit binary strings'''
    maxlinenum = 999999999
    vals = []
    with open(filename) as f:
        linenum = 0  # just to keep track of line
        if not isR2:  # this is the R1 case
            for line in f:
                if ((linenum-1) % 4) == 0 and linenum < maxlinenum:  # if is line index 1, 5, etc
                    try:  # Using this because str.find kept returning incorrect index
                        index1 = line.index('ACAATTTCTACTGTTGTAGAT') + 21  # If not found, to except
                        index2 = index1 + 7

                        heatIndexList = [line[i+index1].lower() == ['c', 't', 'c', 't', 'c', 't', 'c'][i] for i in range(index2-index1)]
                        heatIndex = np.reshape(heatIndexList, 7).astype(int).tolist()
                        heatIndex = ''.join([str(i) for i in heatIndex])
                        vals.append(heatIndex)

                    except:  # If the Cpf1/Cas12-repeat does not show up in sequence
                        print('faulty Cpf1/Cas12 repeat: line ' + str(linenum))
                linenum += 1
        else:  # This is the R2 case
            linenum = 0
            for line in f:
                if ((linenum-1) % 4) == 0 and linenum < maxlinenum:  # if is line index 1, 5, etc
                    line = reversecomplement(line)
                    try:  # Using this because str.find kept returning incorrect index
                        index1 = line.index('ttgacaacctcgtttg') + 16  # If not found, to except
                        index2 = index1 + 7
                        heatIndexList = [line[i + index1].lower() != ['c', 't', 'c', 't', 'c', 't', 'c'][i] for i in
                                         range(index2 - index1)]
                        heatIndex = np.reshape(heatIndexList, 7).astype(int).tolist()
                        heatIndex = ''.join([str(i) for i in heatIndex])
                        vals.append(heatIndex)
                    except:  # If the crispr-gate promoter does not show up in sequence
                        print('faulty crispr-gate promoter ' + str(linenum))
                linenum += 1
    return vals


test_heatmap1, test_heatmap2 = createHeatmap('8551_9805_58262_B87LH_SW-20p-Sucrose_TGACCA_R1.fastq', '8551_9805_58262_B87LH_SW-20p-Sucrose_TGACCA_R2.fastq')

import matplotlib.pyplot as plt

plt.imshow(test_heatmap1)

plt.show()
plt.imshow(test_heatmap2)
plt.show()




