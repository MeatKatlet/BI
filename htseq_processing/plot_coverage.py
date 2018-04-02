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


#import HTSeq.scripts.count  as counts

#HTSeq.

#a = 1
#counts.count_reads_in_features()


# it determines for how many features(genes/exons ...) belongs certain intervsl

#exon - from - to - gene id
#exon - from - to - gene id

#ga = HTSeq.GenomicArray( [ "chr1", "chr2" ], stranded=False )

#ga[ HTSeq.GenomicInterval( "chr1", 250, 400, "+" ) ] = 20


#ga[ HTSeq.GenomicPosition( "chr1", 300, "+" ) ]
#>> 20

#gas = HTSeq.GenomicArrayOfSets( ["chr1", "chr2"], stranded=False )

#ivA = HTSeq.GenomicInterval( "chr1", 100, 300, "." )

#ivB = HTSeq.GenomicInterval( "chr1", 200, 500, "." )


#gas[ivA] += "gene A"
#gas[ivB] += "gene B"

#print([(st[0], sorted(st[1])) for st in gas[ HTSeq.GenomicInterval( "chr1", 0, 500, "." ) ].steps()])

#####
#fs = None
# take cigar intervals and query
"""
for iv3, fs2 in gas[ HTSeq.GenomicInterval( "chr1", 210, 290, "." ) ].steps():
    print(iv3)
    print(fs2)

    if fs is None:
        fs = fs2.copy()
    else:
        fs = fs.intersection(fs2)
"""


#my_set_1 = set(["gene A", "gene B", "gene C"])
#my_set_2 = set(["gene A", "gene D", "gene C"])

#print(my_set_1.intersection(my_set_2))





#print([(st[0], sorted(st[1])) for st in gas[ HTSeq.GenomicInterval( "chr1", 250, 290, "." ) ].steps()])




"""
1. sam file отсортирован по координате, там не обязательно могут идти риды в порядке возрастания координаты , иногда они могут идти и назад по координате

2. будем считать только для точек в 10%(интервал) - сколько ридов пересекает каждый интервал

нужна сумма всех экзонов для получения точек 10%
получим точки в абсолютных координатах исходя из 10% интервала, запишем их в массив
нам надо проверять пересекает ли рид эту точку и если да то +1, 
10% интервалы будем получать просто не учитывая в расчете интроны(не включая их в расстояние)
точки расчитаем 1 раз для каждого гена с учетом всех интронов, из координат ридов вычетать ничего не будем, а будем только искать пересечение(бинарным поиском, потом подумаем как можно ьулучшить) - хеш таблица?, фильтр блума?

3. нужно иметь список границ экзонов в каждом гене или уметь получать его(по границам гена из списка екзонов, значит надо запоминать только границы гена)
ген начало гена, конец гена, это сэкономит место в памяти
ссылки в питон?

4. приходит рид и по его координате надо найти экзон которому он принадлежит
получаем ген из уже существующего кода
получаем из массива с границами гена список интервалов с экзонами
собираем из кусков cigar левый и правый рид 
надо скоректировать левую границу начала левого рида (вычесть интроны)

для большинства ридов 
"""

class UnknownChrom(Exception):
    pass


def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2

##################################################################################################################
def count_reads_in_features(sam_filenames, gff_filename,
                            samtype,
                            order, max_buffer_size,
                            stranded, overlap_mode,
                            multimapped_mode,
                            secondary_alignment_mode,
                            supplementary_alignment_mode,
                            feature_type, id_attribute,
                            additional_attributes,
                            quiet, minaqual, samouts):

    def exists(obj, chain):
        _key = chain.pop(0)
        if _key in obj:
            return exists(obj[_key], chain) if chain else obj[_key]

    def check_overlapped_exons_and_calc_sum(gene):

        rightmost_value = gene["exons"][0][1]
        start = gene["exons"][0][0]
        new_exons = []
        total = 0
        for interval in gene["exons"]:

            if (interval[0] <= rightmost_value and interval[1] >= rightmost_value):

                rightmost_value = interval[1]

            elif (interval[0] > rightmost_value):
                new_exons.append([start, rightmost_value])
                total += (rightmost_value - start)
                start = interval[0]
                rightmost_value = interval[1]

        new_exons.append([start, rightmost_value])

        gene["exons"] = new_exons
        gene["total_sum_of_exons"] = total



    def check_and_count_points_coverage(gene_id, first_read, second_read):

        # определить какую из точек пересекает
        # вычесть из каждой координаты координату начала гена!
        if (first_read is None or second_read is None):
            return

        gene_begin = genes_exons[gene_id]["gene_begin"]

        fstart = first_read.iv.start - gene_begin
        fend = first_read.iv.end - gene_begin
        sstart = second_read.iv.start - gene_begin
        send = second_read.iv.end - gene_begin

        if (first_read.proper_pair == False or second_read.proper_pair == False):
            return

        #######test###############
        """
        if (gene_id == "ENSG00000000003.10" and first_read.iv.start > test_first_exon_start):

            if (first_read.iv.start > test_first_exon_start and second_read.iv.start > test_first_exon_start and first_read.iv.end < test_last_exon_end and second_read.iv.end < test_last_exon_end):
                cvg[first_read.iv] += 1
                cvg[second_read.iv] += 1
                test_n[0] += 1
        
        """
        ######################


        if (fend < sstart and fstart < fend and sstart < send):
            check2(gene_id, fstart, fend)
            check2(gene_id, sstart, send)

        elif (send < fstart and fstart < fend and sstart < send):
            check2(gene_id, fstart, fend)
            check2(gene_id, sstart, send)

        elif (fstart < fend and sstart < send and sstart >= fstart and send >= fend and sstart <= fend):
            check2(gene_id, fstart, send)

        elif (fstart < fend and sstart < send and sstart <= fstart and send >= fstart and send <= fend):
            check2(gene_id, sstart, fend)

        elif (fstart < sstart and send < fend):
            check2(gene_id, fstart, fend)
        elif (sstart < fstart  and fend < send):
            check2(gene_id, sstart, send)


    def check(gene_id, start, end):
        total = 100
        half = total / 2
        left_interval = right_interval = half

        while (left_interval >= 10):

            if (exists(genes_coverage_in_points, [gene_id, half]) == None):  # если точки нет то ищем ближаишую слева
                #half = math.ceil(half)
                half = int(math.floor(half / 10) * 10)
                point = genes_coverage_in_points[gene_id][half]["point"]
                right_interval += 5
                left_interval -= 5


            else:  # если точка есть,
                point = genes_coverage_in_points[gene_id][half]["point"]

            if (point < start):  # слева точка от рида, рид справой строны

                half = half + (right_interval / 2)
                left_interval = right_interval = right_interval / 2

            elif (point > end):  # точка справа от рида, рид слевой стороны

                half = half - (left_interval / 2)
                left_interval = right_interval = left_interval / 2


            elif (point > start and point < end):  # пересекает
                genes_coverage_in_points[gene_id][half]["coverage"] += 1
                return

    def check2(gene_id, start, end):
        #gene_begin = genes_exons[gene_id]["gene_begin"]
        for i in range(0,100,10):

            point = genes_coverage_in_points[gene_id][i]["point"]

            if (start < point and point < end):
                genes_coverage_in_points[gene_id][i]["coverage"] += 1
                return

    def clear_all_cov_points():
        for gene_id, gene in genes_coverage_in_points.iteritems():

            for k, val in gene.iteritems():
                val["coverage"] = 0


    def plot_gene_coverage():
        sys.stderr.write("ENSG00000000003.10 genes on: " + str(test_n[0]) + "\n")
        x = []
        y = []

        i = 0
        for k, val in enumerate(list(cvg[HTSeq.GenomicInterval("chrX", test_first_exon_start, test_last_exon_end)])):
            x.append(i)
            y.append(val)
            i += 1
        plt.plot(x, y)
        plt.show()
        """
         iv = HTSeq.GenomicInterval("chr3", 100, 200, "+")
        cvg[iv] += 1
        iv = HTSeq.GenomicInterval("chr3", 150, 250, "-")
        cvg[iv] += 1
        

        
        """



    if samouts != "":
        if len(samouts) != len(sam_filenames):
            raise ValueError('Select the same number of SAM input and output files')
        # Try to open samout files early in case any of them has issues
        for samout in samouts:
            with open(samout, 'w'):
                pass

    # Try to open samfiles to fail early in case any of them is not there
    if (len(sam_filenames) != 1) or (sam_filenames[0] != '-'):
        for sam_filename in sam_filenames:
            with open(sam_filename):
                pass

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ('M', '=', 'X')

    features = HTSeq.GenomicArrayOfSets("auto", stranded != "no")
    gff = HTSeq.GFF_Reader(gff_filename)

    #genes_coverage_in_points = {}
    genes_coverage_in_points = defaultdict(dict)
    #genes_exons = {}

    genes_exons = defaultdict(dict)
    #cvg = HTSeq.GenomicArray("auto", stranded != "no")



    test_n = [0]
    i= 0

    try:
        for f in gff:
            if f.type == feature_type:
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    raise ValueError("Feature %s does not contain a '%s' attribute" % (f.name, id_attribute))
                if stranded != "no" and f.iv.strand == ".":
                    raise ValueError("Feature %s at %s does not have strand information but you are ""running htseq-count in stranded mode. Use '--stranded=no'." % (f.name, f.iv))


                features[f.iv] += feature_id


                #counts[f.attr[id_attribute]] = 0

                #TODO экзоны не в порядке сортировки! координат
                #ген - граница экзона
                #здесь будут все интервалы и сумма всех интервалов
                gene_id = feature_id #f.attr[id_attribute]
                #if (count(genes_exons[gene_id])==0):
                if (exists(genes_exons, [gene_id]) == None):
                    #координата первого экзона

                    genes_exons[gene_id] = {"total_sum_of_exons": 0, "total_aligned_reads": 0, "gene_begin": 0, "exons": list([[f.iv.start, f.iv.end]])}

                    #genes_exons[gene_id].append({"coords": [f.iv.start, f.iv.end], "total_sum_of_exons": 0, "gene_begin": f.iv.start})

                else :

                    #genes_exons[gene_id]["total_sum_of_exons"] += f.iv.end - f.iv.start #last in list

                    genes_exons[gene_id]["exons"].append([f.iv.start, f.iv.end])

                    #genes_exons[gene_id].append({"coords": [f.iv.start, f.iv.end]})

                #10 точек для гена для которых будем считать покрытие(интроны вычтем)


            i += 1
            if i % 100000 == 0 and not quiet:
                sys.stderr.write("%d GFF lines processed.\n" % i)


    except:
        sys.stderr.write("Error occured when processing GFF file (%s):\n" % gff.get_line_number_string())
        raise


    if not quiet:
        sys.stderr.write("%d GFF lines processed.\n" % i)

    if len(genes_exons) == 0:
        sys.stderr.write("Warning: No features of type '%s' found.\n" % feature_type)



    #TODO подумать как можно улучшить алгоритм, чтобы не держать в памяти для сразу всех генов ряды по 10 точек


    #проход по всем генам и внутри каждого сортируем по первой координате экзона
    #в конце сортировки каждого гена назначаем крайнюю координату начала гена(первый экзон)
    #пересекающиеся экзоны надо склеивать и расширять границы
    #после склеивания будем получать сумму экзонов total_sum_of_exons, т.е. мы получим участки непокрытые ни на одном стренде


    for gene_id, gene in genes_exons.iteritems():

        gene["exons"].sort() #by first member
        gene["gene_begin"] = gene["exons"][0][0]
        #TODO слить все пересекающиеся экзоны и одновременно посчитать сумму длин без полученных промежутков

        check_overlapped_exons_and_calc_sum(gene)

        total = gene["total_sum_of_exons"] # длина всех экзонов

        for ten_interval in xrange(0,100,10):
            point = (total*ten_interval)/100 #точка в абсолютном исчислении % от длины экзона
            prev_exon_end = 0

            for exon_key, exon in enumerate(gene["exons"]):

                #prev_exon_length + exon.start +
                point += (exon[0] - prev_exon_end) #длина интрона

                if (point < exon[1]): #точка конца экзона
                    #пишем точку в конечный массив
                    genes_coverage_in_points[gene_id][ten_interval] = {"point": point-gene["gene_begin"], "coverage": 0}

                    break# переход на следующую точку 10%
                else:
                    #длину экзона не уложившегося записываем
                    #prev_exon_length += exon.end - exon.start
                    prev_exon_end = exon[1]


    """
    ##########################start test#################
    test_begin = genes_exons["ENSG00000000003.10"]["gene_begin"]
    test_first_exon_start = genes_exons["ENSG00000000003.10"]["exons"][0][0]
    last = len(genes_exons["ENSG00000000003.10"]["exons"]) - 1

    test_last_exon_end = genes_exons["ENSG00000000003.10"]["exons"][last][1]


    sys.stderr.write("ENSG00000000003.10 gene_begin: " + str(test_begin)+"\n")
    sys.stderr.write("ENSG00000000003.10 0 exon start: " + str(test_first_exon_start)+"\n")

    sys.stderr.write("ENSG00000000003.10 last exon end(end of gene): " + str(test_last_exon_end)+"\n")


    if(test_begin != test_first_exon_start):
        sys.stderr.write("not_equal!!!!!!\n")
    
    """

    ###########################end test######################




    if samtype == "sam":
        SAM_or_BAM_Reader = HTSeq.SAM_Reader
    elif samtype == "bam":
        SAM_or_BAM_Reader = HTSeq.BAM_Reader
    else:
        raise ValueError("Unknown input format %s specified." % samtype)

    #counts_all = []
    #empty_all = []
    #ambiguous_all = []
    #notaligned_all = []
    #lowqual_all = []
    #nonunique_all = []
    sample = 0


    colors = ["red", "blue", "green","yellow"]
    handlers = []
    sys.stderr.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + "\n")
    for isam, (sam_filename) in enumerate(sam_filenames):

        total_of_reads_in_sample = 0

        if samouts != '':
            samoutfile = open(samouts[isam], 'w')
        else:
            samoutfile = None

        try:
            if sam_filename != "-":
                read_seq_file = SAM_or_BAM_Reader(sam_filename)
                read_seq = read_seq_file
                first_read = next(iter(read_seq))
            else:
                read_seq_file = SAM_or_BAM_Reader(sys.stdin)
                read_seq_iter = iter(read_seq_file)
                first_read = next(read_seq_iter)
                read_seq = itertools.chain([first_read], read_seq_iter)
            pe_mode = first_read.paired_end
        except:
            sys.stderr.write(
                "Error occured when reading beginning of SAM/BAM file.\n")
            raise

        try:
            if pe_mode:
                if order == "name":
                    read_seq = HTSeq.pair_SAM_alignments(read_seq)
                elif order == "pos":
                    read_seq = HTSeq.pair_SAM_alignments_with_buffer(
                            read_seq,
                            max_buffer_size=max_buffer_size)
                else:
                    raise ValueError("Illegal order specified.")
            #empty = 0
            #ambiguous = 0
            notaligned = 0
            lowqual = 0
            #nonunique = 0
            i = 0
            for r in read_seq:
                #TODO 'NoneType' object has no attribute 'iv' raised in plot_coverage.py:169]
                total_of_reads_in_sample += 1
                if i > 0 and i % 100000 == 0 and not quiet:
                    sys.stderr.write(
                        "%d SAM alignment record%s processed.\n" %
                        (i, "s" if not pe_mode else " pairs"))
                    sys.stderr.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + "\n")
                i += 1
                if not pe_mode:
                    if not r.aligned:
                        #notaligned += 1
                        #write_to_samout(r, "__not_aligned", samoutfile)
                        continue
                    if ((secondary_alignment_mode == 'ignore') and
                       r.not_primary_alignment):
                        continue
                    if ((supplementary_alignment_mode == 'ignore') and
                       r.supplementary):
                        continue
                    try:
                        if r.optional_field("NH") > 1:
                            #nonunique += 1
                            #write_to_samout(r, "__alignment_not_unique", samoutfile)
                            if multimapped_mode == 'none':
                                continue
                    except KeyError:
                        pass
                    if r.aQual < minaqual:
                        lowqual += 1
                        #write_to_samout(r, "__too_low_aQual", samoutfile)
                        continue
                    if stranded != "reverse":
                        iv_seq = (co.ref_iv for co in r.cigar if co.type ==
                                  "M" and co.size > 0)
                    else:
                        iv_seq = (invert_strand(co.ref_iv)
                                  for co in r.cigar if (co.type in com and
                                                        co.size > 0))
                else:
                    if r[0] is not None and r[0].aligned:
                        if stranded != "reverse":
                            iv_seq = (co.ref_iv for co in r[0].cigar
                                      if co.type in com and co.size > 0)
                        else:
                            iv_seq = (invert_strand(co.ref_iv) for co in r[0].cigar
                                      if co.type in com and co.size > 0)
                    else:
                        iv_seq = tuple()
                    if r[1] is not None and r[1].aligned:
                        if stranded != "reverse":
                            iv_seq = itertools.chain(
                                    iv_seq,
                                    (invert_strand(co.ref_iv) for co in r[1].cigar
                                    if co.type in com and co.size > 0))
                        else:
                            iv_seq = itertools.chain(
                                    iv_seq,
                                    (co.ref_iv for co in r[1].cigar
                                     if co.type in com and co.size > 0))
                    else:
                        if (r[0] is None) or not (r[0].aligned):
                            #write_to_samout(r, "__not_aligned", samoutfile)
                            #notaligned += 1
                            continue
                    if secondary_alignment_mode == 'ignore':
                        if (r[0] is not None) and r[0].not_primary_alignment:
                            continue
                        elif (r[1] is not None) and r[1].not_primary_alignment:
                            continue
                    if supplementary_alignment_mode == 'ignore':
                        if (r[0] is not None) and r[0].supplementary:
                            continue
                        elif (r[1] is not None) and r[1].supplementary:
                            continue
                    try:
                        if ((r[0] is not None and r[0].optional_field("NH") > 1) or
                           (r[1] is not None and r[1].optional_field("NH") > 1)):
                            #nonunique += 1
                            #write_to_samout(r, "__alignment_not_unique", samoutfile)
                            if multimapped_mode == 'none':
                                continue
                    except KeyError:
                        pass
                    if ((r[0] and r[0].aQual < minaqual) or
                       (r[1] and r[1].aQual < minaqual)):
                        lowqual += 1
                        #write_to_samout(r, "__too_low_aQual", samoutfile)
                        continue

                try:
                    if overlap_mode == "union":
                        fs = set()
                        for iv in iv_seq:
                            if iv.chrom not in features.chrom_vectors:
                                raise UnknownChrom
                            for iv2, fs2 in features[iv].steps():
                                fs = fs.union(fs2)
                    elif overlap_mode in ("intersection-strict",
                                          "intersection-nonempty"):
                        fs = None
                        for iv in iv_seq:
                            if iv.chrom not in features.chrom_vectors:
                                continue
                                #raise UnknownChrom
                            for iv2, fs2 in features[iv].steps():
                                if ((len(fs2) > 0) or
                                   (overlap_mode == "intersection-strict")):
                                    if fs is None:
                                        fs = fs2.copy()
                                    else:
                                        fs = fs.intersection(fs2)
                    else:
                        sys.exit("Illegal overlap mode.")


                    if fs is not None and len(fs) > 0:
                        if multimapped_mode == 'none':
                            if len(fs) == 1:
                                #counts[list(fs)[0]] += 1
                                #read mapped only for one exon, (all cigar parts of both reads in pair mapped on one gene, but may be for several exons)
                                #we can take this read into account of analysis
                                #they must come in sorted order by coordinate!
                                #this is one unit of analysis. save it in memory and go throught it
                               
                                gene_name = list(fs)[0]# - имя гена

                                genes_exons[gene_name]["total_aligned_reads"] += 1


                                #if (total_of_reads_in_sample==100000):
                                #   break

                                check_and_count_points_coverage(gene_name, r[0], r[1])


                            """
                            elif multimapped_mode == 'all':
                                for fsi in list(fs):
                                    #counts[fsi] += 1 
                            """
                        else:
                            sys.exit("Illegal multimap mode.")


                except UnknownChrom:
                    #write_to_samout(r, "__no_feature", samoutfile)
                    #empty += 1
                    raise

        except:
            sys.stderr.write(
                "Error occured when processing SAM input (%s):\n" %
                read_seq_file.get_line_number_string())
            raise

        if not quiet:
            sys.stderr.write(
                "%d SAM %s processed.\n" %
                (i, "alignments " if not pe_mode else "alignment pairs"))

        if samoutfile is not None:
            samoutfile.close()

        #сохранить данные в таблицы чтобы работать с ними как угодно потом!




        outfile = open('/home/kirill/bi/transcript/'+str(sample)+'_dict.txt', 'w')
        outfile.write("total_of_reads_in_sample" + '\t' + str(total_of_reads_in_sample) + '\n')
        for gene_id, gene in genes_coverage_in_points.iteritems():

            outfile.write(str(gene_id) + '\t' + str(genes_exons[gene_id]["total_aligned_reads"]) + '\t' + str(genes_exons[gene_id]["total_sum_of_exons"]) + '\n')

            outfile.write(str(gene_id) + '\t')
            [outfile.write(str(val["coverage"]) + '\t') for k, val in gene.iteritems()]
            outfile.write('\n')

        outfile.close()

        #############test################

        #plot_gene_coverage()


        ################################

        #1. получить % от числа ридов картированных на ген в конкретной точке(сумма всех % на 10 точках = 100) - число ридов картированных на ген будем записывать в массив(это бывший массиыв count)
        #2 для каждой точки делим полученный процент на длину конкретного гена (total_sum of exons)
        #3. для каждой точки делим величину на общее число ридов в образце
        #4. deviance - min - max всех значений? точка на графике среднее между ними



        #TODO обнулять покрытие между файлами/образцами!!

        y = np.zeros(10)
        x = np.arange(0, 100, 10)


        most_min = np.zeros(10)
        most_max = np.zeros(10)
        j = 0
        for gene_id, gene in genes_coverage_in_points.iteritems():

            i = 0
            if (genes_exons[gene_id]["total_aligned_reads"]==0 or genes_exons[gene_id]["total_sum_of_exons"]==0 or total_of_reads_in_sample==0):
                continue

            for k, val in gene.iteritems():
                c = (((val["coverage"]/genes_exons[gene_id]["total_aligned_reads"])/genes_exons[gene_id]["total_sum_of_exons"])/total_of_reads_in_sample)

                y[i] += c

                if (j==0):# первый ген в выборке у каждого из 10
                    most_min[i] = c

                if (most_max[i] < c):
                    most_max[i] = c
                elif(most_min[i] > c):
                    most_min[i] = c

                i +=1
            j+=1


        y_means = np.zeros(10)
        i = 0
        for val in y:
            y_means[i] = val/genes_coverage_in_points.__len__()
            i += 1

        y_deviations = np.zeros(10)

        for i in range(0,10,1):
            y_deviations[i] = most_max[i] - most_min[i]


        patch = mpatches.Patch(color=colors[sample])
        handlers.append(patch)

        #будет создан 1 график с четырьмя линиями!
        plt.errorbar(x, y_means, yerr=y_deviations, color=colors[sample], ls='--', marker='o', capsize=5, capthick=1, ecolor='black')



        sample +=1

        #обнуление точек покрытия
        clear_all_cov_points()


    plt.legend(handlers, ['Sample '+str(v) for v in range(0,sample,1)])
    plt.title('Positions relative coverege')
    plt.xlabel('5` -> 3` positions, %')
    plt.ylabel('relative coverage')
    plt.grid(True)

    plt.savefig('/home/kirill/bi/transcript/covarage.png')
    plt.show()
    plt.close()



def my_showwarning(message, category, filename, lineno=None, file=None, line=None):
    sys.stderr.write("Warning: %s\n" % message)



def main():

    pa = argparse.ArgumentParser(
        usage="%(prog)s [options] alignment_file gff_file",
        description="This script takes one or more alignment files in SAM/BAM " +
        "format and a feature file in GFF format and calculates for each feature " +
        "the number of reads mapping to it. See " +
        "http://htseq.readthedocs.io/en/master/count.html for details.",
        epilog="Written by Simon Anders (sanders@fs.tum.de), " +
        "European Molecular Biology Laboratory (EMBL). (c) 2010. " +
        "Released under the terms of the GNU General Public License v3. " +
        "Part of the 'HTSeq' framework, version %s." % HTSeq.__version__)

    pa.add_argument(
            "samfilenames", nargs='+', type=str,
            help="Path to the SAM/BAM files containing the mapped reads. " +
            "If '-' is selected, read from standard input")

    pa.add_argument(
            "featuresfilename", type=str,
            help="Path to the file containing the features")

    pa.add_argument(
            "-f", "--format", dest="samtype",
            choices=("sam", "bam"), default="sam",
            help="type of <alignment_file> data, either 'sam' or 'bam' (default: sam)")

    pa.add_argument(
            "-r", "--order", dest="order",
            choices=("pos", "name"), default="name",
            help="'pos' or 'name'. Sorting order of <alignment_file> (default: name). Paired-end sequencing " +
            "data must be sorted either by position or by read name, and the sorting order " +
            "must be specified. Ignored for single-end data.")

    pa.add_argument(
            "--max-reads-in-buffer", dest="max_buffer_size", type=int,
            default=30000000,
            help="When <alignment_file> is paired end sorted by position, " +
            "allow only so many reads to stay in memory until the mates are " +
            "found (raising this number will use more memory). Has no effect " +
            "for single end or paired end sorted by name")

    pa.add_argument(
            "-s", "--stranded", dest="stranded",
            choices=("yes", "no", "reverse"), default="yes",
            help="whether the data is from a strand-specific assay. Specify 'yes', " +
            "'no', or 'reverse' (default: yes). " +
            "'reverse' means 'yes' with reversed strand interpretation")

    pa.add_argument(
            "-a", "--minaqual", type=int, dest="minaqual",
            default=10,
            help="skip all reads with alignment quality lower than the given " +
            "minimum value (default: 10)")

    pa.add_argument(
            "-t", "--type", type=str, dest="featuretype",
            default="exon", help="feature type (3rd column in GFF file) to be used, " +
            "all features of other type are ignored (default, suitable for Ensembl " +
            "GTF files: exon)")

    pa.add_argument(
            "-i", "--idattr", type=str, dest="idattr",
            default="gene_id", help="GFF attribute to be used as feature ID (default, " +
            "suitable for Ensembl GTF files: gene_id)")

    pa.add_argument(
            "--additional-attr", type=str, nargs='+',
            default=(), help="Additional feature attributes (default: none, " +
            "suitable for Ensembl GTF files: gene_name)")

    pa.add_argument(
            "-m", "--mode", dest="mode",
            choices=("union", "intersection-strict", "intersection-nonempty"),
            default="union", help="mode to handle reads overlapping more than one feature " +
            "(choices: union, intersection-strict, intersection-nonempty; default: union)")

    pa.add_argument(
            "--nonunique", dest="nonunique", type=str,
            choices=("none", "all"), default="none",
            help="Whether to score reads that are not uniquely aligned " +
            "or ambiguously assigned to features")

    pa.add_argument(
            "--secondary-alignments", dest="secondary_alignments", type=str,
            choices=("score", "ignore"), default="score",
            help="Whether to score secondary alignments (0x100 flag)")

    pa.add_argument(
            "--supplementary-alignments", dest="supplementary_alignments", type=str,
            choices=("score", "ignore"), default="score",
            help="Whether to score supplementary alignments (0x800 flag)")

    pa.add_argument(
            "-o", "--samout", type=str, dest="samouts", nargs='+',
            default="", help="write out all SAM alignment records into an output " +
            "SAM file called SAMOUT, annotating each line with its feature assignment " +
            "(as an optional field with tag 'XF')")

    pa.add_argument(
            "-q", "--quiet", action="store_true", dest="quiet",
            help="suppress progress report")  # and warnings" )

    args = pa.parse_args()

    warnings.showwarning = my_showwarning
    try:
        count_reads_in_features(
                args.samfilenames,
                args.featuresfilename,
                args.samtype,
                args.order,
                args.max_buffer_size,
                args.stranded,
                args.mode,
                args.nonunique,
                args.secondary_alignments,
                args.supplementary_alignments,
                args.featuretype,
                args.idattr,
                args.additional_attr,
                args.quiet,
                args.minaqual,
                args.samouts)
    except:
        sys.stderr.write("  %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write("  [Exception type: %s, raised in %s:%d]\n" %
                         (sys.exc_info()[1].__class__.__name__,
                          os.path.basename(traceback.extract_tb(
                              sys.exc_info()[2])[-1][0]),
                          traceback.extract_tb(sys.exc_info()[2])[-1][1]))
        sys.exit(1)


if __name__ == "__main__":
    main()









