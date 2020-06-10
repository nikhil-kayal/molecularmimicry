#!/usr/bin/env python

####################################################
#title		get_lcs_parallel2.py
#description	Get the lcs between two fasta files
#authors	Praveen Anand, Nikhil 
#email		praveen@nference.net
#date		06/10/2020
#version	0.1
#pythonversion	3.6.0
####################################################

from __future__ import division
from __future__ import print_function

import sys
import os
import re
import json
import ast
import csv
import ast
import pickle
import string
import requests
from time import sleep
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import os
import functools
import operator as op
import numpy as np
from scipy.stats import binom


from itertools import islice

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


# Function to find Longest common substring of sequences X[0..m-1] and Y[0..n-1]
def LCS(X, Y, m, n):
    maxLength = 0 # stores the max length of LCS
    endingIndex = m # stores the ending index of LCS in X

    # lookup[i][j] stores the length of LCS of substring X[0..i-1], Y[0..j-1]
    lookup = [[0 for x in range(n + 1)] for y in range(m + 1)]

    # fill the lookup table in bottom-up manner
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            # if current character of X and Y matches
            if X[i - 1] == Y[j - 1]:
                lookup[i][j] = lookup[i - 1][j - 1] + 1

                # update the maximum length and ending index
                if lookup[i][j] > maxLength:
                    maxLength = lookup[i][j]
                    endingIndex = i

    # return Longest common substring having length maxLength
    return X[endingIndex - maxLength: endingIndex]

import configparser
import argparse
import re


if len(sys.argv) == 1:
   sys.argv[1:] = ["-h"]
   print("Wrong numberof arguments entered.")
   print("Usage:get_lcs_parallel2.py -f1 <coronavirusfasta> -f2 <humanfasta> -o <resultsoutfile>")

parser = argparse.ArgumentParser(prog='get_lcs_parallel2')
parser.add_argument('-f1',action='store',dest='coronavirusfasta',required=True,help='The fasta file containing all the coronavirus sequences')
parser.add_argument('-f2',action='store',dest='humanfasta',required=True, help='The fasta file containing all the human sequences')
parser.add_argument('-o',action='store',dest='resultsoutfile',required=True, help='The final output will be written to this file')


##Assigning arguments to variables
inputs=parser.parse_args()
coronavirusfasta=inputs.coronavirusfasta
humanfasta=inputs.humanfasta
outfile=inputs.resultsoutfile

coronarecords = list(SeqIO.parse(coronavirusfasta, "fasta"))
humanrecords = list(SeqIO.parse(humanfasta, "fasta"))


def singleseqhumanscanindex(index,coronarecords=coronarecords,humanrecords=humanrecords):
    s1name = str(coronarecords[index].name)
    s1 = str(coronarecords[index].seq)
    myresults = pd.DataFrame(columns = ['common','s1name','s1start','s2name','s2start'])
    for j in range(0,len(humanrecords)):
        s2 = str(humanrecords[int(j)].seq)
        x = LCS(s2,s1,len(s2),len(s1))
        while(len(x) >= 8):
            print(x)
            x = LCS(s2,s1,len(s2),len(s1))
            if len(x) >= 8:
                myresults = myresults.append({'common': x,'s1name':s1name,'s1start':s1.find(x),'s2name':str(humanrecords[j].name),'s2start':str(humanrecords[j].seq).find(x)}, ignore_index=True)
                s2 = s2.replace(x,' '*len(x))
    print("Done for " + str(coronarecords[index].name))
    return(myresults)

from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count()


results_chunk = Parallel(n_jobs=num_cores)(delayed(singleseqhumanscanindex)(i) for i in range(0,len(coronarecords)))

result = pd.concat(results_chunk)

result.to_csv(outfile, index=False)
