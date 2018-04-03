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



#TODO оформить в модуль!
class CoverageFileProcessor:

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


fp = CoverageFileProcessor()

sample = 0
colors = ["red", "blue", "green","yellow"]
handlers = []

#прочтет все файлы в папке
for i in range(0,3,1):
    file_path = str(i)+"_dict.txt"
    processor = FileProcessing(file_path, fp)
    processor.process_file()

    total_of_reads_in_sample = fp.total_reads_in_sample
    genes_coverage_in_points = fp.genes_coverage_in_points
    genes_exons = fp.genes_exons


    #можно делать график для одного образца
    #TODO можно будет засунуть в функцию и подключить из файла обработчика данных!
    #TODO как это все делать в numpy
    #TODO убрать из анализа гены с 0 покрытием, делить на количество генов с ненулевым общим покрытием, за размах будем брать 0 если он минимум и максимум
    #TODO среднеквадратичное отклонение делать надо а не max min
    y = np.zeros(10)
    x = np.arange(0, 100, 10)

    most_min = np.zeros(10)
    most_max = np.zeros(10)
    j = 0
    total_genes = genes_coverage_in_points.__len__()
    for gene_id, gene in genes_coverage_in_points.iteritems():

        i = 0
        if (genes_exons[gene_id]["total_aligned_reads"] == 0 or genes_exons[gene_id]["total_sum_of_exons"] == 0 or total_of_reads_in_sample == 0):
            total_genes-=1
            continue
        #10**11
        first_not_zero = np.zeros(10)
        for k, val in gene.iteritems():
            c = (((float(val["coverage"]) / float(genes_exons[gene_id]["total_aligned_reads"])) / float(genes_exons[gene_id]["total_sum_of_exons"])) / float(total_of_reads_in_sample))

            y[i] += c

            if (first_not_zero[i] == 0 and c != 0):  #first non zero value in given 10 percent interval  point
                most_min[i] = c
                first_not_zero[i] = 1

            if (most_max[i] < c):
                most_max[i] = c
            elif (most_min[i] > c and c != 0):
                most_min[i] = c

            i += 1
        j += 1

    y_means = np.zeros(10)
    i = 0
    for val in y:
        y_means[i] = val / total_genes
        i += 1

    #(каждый x в выборке на 10 % линии - среднее)^2 / количество генов, корень из этого
    ss = np.zeros(10)
    ss2 = np.zeros(10)

    for gene_id, gene in genes_coverage_in_points.iteritems():

        i = 0
        if (genes_exons[gene_id]["total_aligned_reads"] == 0 or genes_exons[gene_id]["total_sum_of_exons"] == 0 or total_of_reads_in_sample == 0):

            continue
        #10**11
        first_not_zero = np.zeros(10)
        for k, val in gene.iteritems():
            c = (((float(val["coverage"]) / float(genes_exons[gene_id]["total_aligned_reads"])) / float(genes_exons[gene_id]["total_sum_of_exons"])) / float(total_of_reads_in_sample))

            ss[i] += (c-y_means[i])**2


            i += 1

    for i in range(0, 10, 1):
        ss2[i] = math.sqrt(float(ss[i]/total_genes))


    """
    y_deviations = np.zeros(10)

    for i in range(0, 10, 1):
        y_deviations[i] = most_max[i] - most_min[i]
    """

    patch = mpatches.Patch(color=colors[sample])
    handlers.append(patch)

    # будет создан 1 график с четырьмя линиями!
    plt.errorbar(x, y_means, yerr=ss2, color=colors[sample], ls='--', marker='o', capsize=5, capthick=1,ecolor='black')

    sample += 1



plt.legend(handlers, ['Sample '+str(v) for v in range(0,sample,1)])
plt.title('Positions relative coverege')
plt.xlabel('5` -> 3` positions, %')
plt.ylabel('relative coverage')
plt.grid(True)

plt.savefig('/home/kirill/bi/transcript/covarage.png')
plt.show()
plt.close()

