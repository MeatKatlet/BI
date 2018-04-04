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
from calc_coverage import CalcCoverage



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




sample = 0
colors = ["red", "blue", "green","yellow"]
handlers = []

#прочтет все файлы в папке
for i in range(0,4,1):
    fp = CoverageFileProcessor()
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
    CalcCoverage.do_coverage(genes_coverage_in_points,genes_exons,total_of_reads_in_sample,colors,sample,handlers)

    sample += 1



plt.legend(handlers, ['Sample '+str(v) for v in range(0,sample,1)])
plt.title('Positions relative coverege')
plt.xlabel('5` -> 3` positions, %')
plt.ylabel('relative coverage')
plt.grid(True)

plt.savefig('/home/kirill/bi/transcript/covarage.png')
plt.show()
plt.close()

