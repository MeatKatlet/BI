# -*- coding: utf-8 -*-
from __future__ import division
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import math

class CalcCoverage():

    def do_coverage(self,genes_coverage_in_points,genes_exons,total_of_reads_in_sample,colors,sample,handlers):
        y = np.zeros(10)
        x = np.arange(0, 100, 10)

        most_min = np.zeros(10)
        most_max = np.zeros(10)
        j = 0
        total_genes = genes_coverage_in_points.__len__()
        for gene_id, gene in genes_coverage_in_points.iteritems():

            i = 0
            if (genes_exons[gene_id]["total_aligned_reads"] == 0 or genes_exons[gene_id][
                "total_sum_of_exons"] == 0 or total_of_reads_in_sample == 0):
                total_genes -= 1
                continue
            # 10**11
            first_not_zero = np.zeros(10)
            for k, val in gene.iteritems():
                c = (((float(val["coverage"]) / float(genes_exons[gene_id]["total_aligned_reads"])) / float(
                    genes_exons[gene_id]["total_sum_of_exons"])) / float(total_of_reads_in_sample))

                y[i] += c

                if (first_not_zero[i] == 0 and c != 0):  # first non zero value in given 10 percent interval  point
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
            y_means[i] = float(val) / total_genes
            i += 1

        # (каждый x в выборке на 10 % линии - среднее)^2 / количество генов, корень из этого
        ss = np.zeros(10)
        ss2 = np.zeros(10)

        for gene_id, gene in genes_coverage_in_points.iteritems():

            i = 0
            if (genes_exons[gene_id]["total_aligned_reads"] == 0 or genes_exons[gene_id][
                "total_sum_of_exons"] == 0 or total_of_reads_in_sample == 0):
                continue
            # 10**11
            first_not_zero = np.zeros(10)
            for k, val in gene.iteritems():
                c = (((float(val["coverage"]) / float(genes_exons[gene_id]["total_aligned_reads"])) / float(
                    genes_exons[gene_id]["total_sum_of_exons"])) / float(total_of_reads_in_sample))

                ss[i] += (c - y_means[i]) ** 2

                i += 1

        for i in range(0, 10, 1):
            ss2[i] = math.sqrt(float(ss[i]) / total_genes)

        """
        y_deviations = np.zeros(10)

        for i in range(0, 10, 1):
            y_deviations[i] = most_max[i] - most_min[i]
        """

        patch = mpatches.Patch(color=colors[sample])
        handlers.append(patch)

        y_means = y_means * (10 ** 12)
        ss2 = ss2 * (10 ** 12)

        # будет создан 1 график с четырьмя линиями!
        plt.errorbar(x, y_means, yerr=ss2, color=colors[sample], ls='--', marker='o', capsize=5, capthick=1,ecolor='black')