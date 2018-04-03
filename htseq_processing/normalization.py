# -*- coding: utf-8 -*-
from __future__ import division
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import sys
import argparse
import itertools
import warnings
import traceback
import os.path
import HTSeq
import math
import numpy as np
from collections import defaultdict
from time import gmtime, strftime
from common.files import FileProcessing
import os

class CoverageFileProcessor1:

    def __init__(self):

        self.genes_coverage = defaultdict(dict)

        self.line_type = 0

    def process_data(self, line):

        if (line != "\n"):
            list = line.split("\t")
            self.genes_coverage[list[0]] = int(list[1])

class CoverageFileProcessor2:

    def __init__(self):

        self.genes_coverage_in_points = defaultdict(dict)
        self.total_reads_in_sample = 0
        self.genes_exons = defaultdict(dict)

        self.line_type = 0

    def process_data(self, line):

       if (line != "\n"):
        list = line.split("\t")
        if(list[0]=="total_of_reads_in_sample"):
            self.total_reads_in_sample = int(list[1])
            self.line_type += 1

        elif(self.line_type==1):
            self.line_type += 1
            self.genes_exons[list[0]]["total_sum_of_exons"] = int(list[2])
            self.genes_exons[list[0]]["total_aligned_reads"] = int(list[1])

        elif(self.line_type==2):
            self.line_type = 1

            for k in range(0,100,10):
                i = int(k/10+1)
                self.genes_coverage_in_points[list[0]][k] = {"coverage": int(list[i])}


os.chdir("/home/kirill/bi/transcript/coverage/")
#положить все файлы в 1 папку



samples = list()

#get raw coverage in samples
for i in range(0,4,1):
    fp = CoverageFileProcessor2()
    file_path = str(i)+"_dict.txt"
    processor = FileProcessing(file_path, fp)
    processor.process_file()

    x = np.zeros(fp.genes_exons.__len__())

    samples.append(fp)

    i2 = 0
    max = 0
    for k, j in samples[i].genes_exons.iteritems():
        if (j["total_aligned_reads"]>0):
            x[i2] = j["total_aligned_reads"]
            i2 += 1
            if(j["total_aligned_reads"] > max):
                max=j["total_aligned_reads"]


    plt.hist(x, bins=100)

    plt.xlabel('Smarts')
    plt.ylabel('Probability')
    plt.title('Histogram of IQ')


    plt.grid(True)
    plt.show()



samples2 = list()
#we need genes lengths for normalization
for i in range(0,4,1):
    fp = CoverageFileProcessor1()
    file_path = str(i)+"_coverage.txt"
    processor = FileProcessing(file_path, fp)
    processor.process_file()

    samples2[i] = fp


genes_coverage_norm = []

for i in range(0,4,1):

    genes_coverage_norm[i] = defaultdict(dict)

    for gene_id, gene in samples2[i].genes_coverage.iteritems():
        genes_coverage_norm[i][gene_id] = (float(gene)/float(samples[i].genes_exons[gene_id]["total_sum_of_exons"]))/float(samples[i].total_reads_in_sample)



#TODO 10** ... convert?
#todo hist of normcoverage?, for each gene, need log normalization


#запись в файл нормализованных данных, для каждого образца отдельно!

#далее будем работать с ними для коэффициентов корреляции, и остального анализа!








