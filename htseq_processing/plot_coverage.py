import sys
import argparse
import itertools
import warnings
import traceback
import os.path
import HTSeq
import HTSeq.scripts.count  as counts

#HTSeq.

a = 1
#counts.count_reads_in_features()


# it determines for how many features(genes/exons ...) belongs certain intervsl

#exon - from - to - gene id
#exon - from - to - gene id

#ga = HTSeq.GenomicArray( [ "chr1", "chr2" ], stranded=False )

#ga[ HTSeq.GenomicInterval( "chr1", 250, 400, "+" ) ] = 20


#ga[ HTSeq.GenomicPosition( "chr1", 300, "+" ) ]
#>> 20

gas = HTSeq.GenomicArrayOfSets( ["chr1", "chr2"], stranded=False )

ivA = HTSeq.GenomicInterval( "chr1", 100, 300, "." )

ivB = HTSeq.GenomicInterval( "chr1", 200, 500, "." )


gas[ivA] += "gene A"
gas[ivB] += "gene B"

#print([(st[0], sorted(st[1])) for st in gas[ HTSeq.GenomicInterval( "chr1", 0, 500, "." ) ].steps()])

#####
fs = None
# take cigar intervals and query

for iv3, fs2 in gas[ HTSeq.GenomicInterval( "chr1", 210, 290, "." ) ].steps():
    print(iv3)
    print(fs2)

    if fs is None:
        fs = fs2.copy()
    else:
        fs = fs.intersection(fs2)

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
    counts = {}
    genes_coverage_in_points = {}
    genes_exons = {}
    attributes = {}
    i = 0

    try:
        for f in gff:
            if f.type == feature_type:
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    raise ValueError("Feature %s does not contain a '%s' attribute" %
                                     (f.name, id_attribute))
                if stranded != "no" and f.iv.strand == ".":
                    raise ValueError("Feature %s at %s does not have strand information but you are "
                                     "running htseq-count in stranded mode. Use '--stranded=no'." %
                                     (f.name, f.iv))
                features[f.iv] += feature_id


                #counts[f.attr[id_attribute]] = 0

                #ген - граница экзона
                #здесь будут все интервалы и сумма всех интервалов
                genes_exons[f.attr[id_attribute]][] = [[f.interval.start, f.interval.end],сумма всех экзонов+=]

                #10 точек для гена для которых будем считать покрытие(интроны вычтем)
                #будем считать что экзоны приходят отсортированные, в этом надо будет убедиться!
                genes_coverage_in_points[f.attr[id_attribute]][] =




                attributes[f.attr[id_attribute]] = [
                        f.attr[attr] if attr in f.attr else ''
                        for attr in additional_attributes]

            i += 1
            if i % 100000 == 0 and not quiet:
                sys.stderr.write("%d GFF lines processed.\n" % i)
    except:
        sys.stderr.write(
            "Error occured when processing GFF file (%s):\n" %
            gff.get_line_number_string())
        raise


    if not quiet:
        sys.stderr.write("%d GFF lines processed.\n" % i)

    if len(genes_coverage_in_points) == 0:
        sys.stderr.write(
            "Warning: No features of type '%s' found.\n" % feature_type)


    for gene,key in genes_exons:

        for exon in gene:
            genes_coverage_in_points[key] =


    #TODO получаем здесь точки для измерения покрытия уже с учетом экзонов! 10 точек в каждом гене

    for gene, key in genes_exons:
        total = gene[total]# длина всех экзонов

        for ten_interval in range(0,100,10):
            point = (total*ten_interval)/100#точка в абсолютном исчислении
            prev_exon_end = 0
            for exon in gene:

                #prev_exon_length + exon.start +
                point += (exon.start - prev_exon_end) #длина интрона

                if (point < exon.end):
                    #пишем точку в конечный массив
                    break# переход на следующую точку 10%
                else:
                    #длину экзона не уложившегося записываем
                    #prev_exon_length += exon.end - exon.start
                    prev_exon_end = exon.end
    genes_coverage_in_points


    if samtype == "sam":
        SAM_or_BAM_Reader = HTSeq.SAM_Reader
    elif samtype == "bam":
        SAM_or_BAM_Reader = HTSeq.BAM_Reader
    else:
        raise ValueError("Unknown input format %s specified." % samtype)

    counts_all = []
    empty_all = []
    ambiguous_all = []
    notaligned_all = []
    lowqual_all = []
    nonunique_all = []
    for isam, (sam_filename) in enumerate(sam_filenames):
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
            empty = 0
            ambiguous = 0
            notaligned = 0
            lowqual = 0
            nonunique = 0
            i = 0
            for r in read_seq:
                if i > 0 and i % 100000 == 0 and not quiet:
                    sys.stderr.write(
                        "%d SAM alignment record%s processed.\n" %
                        (i, "s" if not pe_mode else " pairs"))

                i += 1
                if not pe_mode:
                    if not r.aligned:
                        notaligned += 1
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
                            nonunique += 1
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
                            notaligned += 1
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
                            nonunique += 1
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
                                raise UnknownChrom
                            for iv2, fs2 in features[iv].steps():
                                if ((len(fs2) > 0) or
                                   (overlap_mode == "intersection-strict")):
                                    if fs is None:
                                        fs = fs2.copy()
                                    else:
                                        fs = fs.intersection(fs2)
                    else:
                        sys.exit("Illegal overlap mode.")

                    if fs is None or len(fs) == 0:
                        #write_to_samout(r, "__no_feature", samoutfile)
                        empty += 1
                    elif len(fs) > 1:
                        #write_to_samout(r, "__ambiguous[" + '+'.join(fs) + "]",samoutfile)
                        ambiguous += 1
                    else:
                        #write_to_samout(r, list(fs)[0], samoutfile)

                    if fs is not None and len(fs) > 0:
                        if multimapped_mode == 'none':
                            if len(fs) == 1:
                                counts[list(fs)[0]] += 1
                                #read mapped only for one exon, (all cigar parts of both reads in pair mapped on one gene, but may be for several exons)
                                #we can take this read into account of analysis
                                #they must come in sorted order by coordinate!
                                #this is one unit of analysis. save it in memory and go throught it
                                """
                                храним  в памяти границы рида до тех пор пока его их не привысит левая граница очередного рида
                                
                                получаем массив с интервалами и покрытием каждого интервала , сначала в абсолютных координатах, потом в % с учетом вычитания интронов - это все храним для 1 гена, желательно делать для точек 10% разницей,
                                

                                сложение интервалов между двумя генами(будем складывать только точки в 10%, так проще и экономнее)
                                
                                """

                                #TODO определить пересекает ли рид(левый или правый или оба,тогда за 1) одну из 10 точек этого гена, если да то +1 на эту точку!
                                #вычесть из координвт рида координату начала гена(первого экзона?)
                                #gene_name = list(fs)[0] - имя гена
                                #genes_coverage_in_points[gene_name]
                                #todo функцию написать!



                        elif multimapped_mode == 'all':
                            for fsi in list(fs):
                                counts[fsi] += 1
                        else:
                            sys.exit("Illegal multimap mode.")


                except UnknownChrom:
                    #write_to_samout(r, "__no_feature", samoutfile)
                    empty += 1

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

        counts_all.append(counts.copy())
        for fn in counts:
            counts[fn] = 0
        empty_all.append(empty)
        ambiguous_all.append(ambiguous)
        lowqual_all.append(lowqual)
        notaligned_all.append(notaligned)
        nonunique_all.append(nonunique)

    pad = ['' for attr in additional_attributes]
    for fn in sorted(counts.keys()):
        print('\t'.join([fn] + attributes[fn] + [str(c[fn]) for c in counts_all]))
    print('\t'.join(["__no_feature"] + pad + [str(c) for c in empty_all]))
    print('\t'.join(["__ambiguous"] + pad + [str(c) for c in ambiguous_all]))
    print('\t'.join(["__too_low_aQual"] + pad + [str(c) for c in lowqual_all]))
    print('\t'.join(["__not_aligned"] + pad + [str(c) for c in notaligned_all]))
    print('\t'.join(["__alignment_not_unique"] + pad + [str(c) for c in nonunique_all]))







