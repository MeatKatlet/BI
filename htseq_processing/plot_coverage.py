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
ivB = HTSeq.GenomicInterval( "chr1", 400, 500, "." )

gas[ivA] += "gene A"
gas[ivB] += "gene B"

#print([(st[0], sorted(st[1])) for st in gas[ HTSeq.GenomicInterval( "chr1", 0, 500, "." ) ].steps()])

#####
fs = None
# take cigar intervals and query
for iv2, fs2 in gas[ HTSeq.GenomicInterval( "chr1", 150, 250, "." ) ].steps():
    print(iv2)
    print(fs2)

    if fs is None:
        fs = fs2.copy()
    else:
        fs = fs.intersection(fs2)

#my_set_1 = set(["gene A", "gene B", "gene C"])
#my_set_2 = set(["gene A", "gene D", "gene C"])

#print(my_set_1.intersection(my_set_2))





#print([(st[0], sorted(st[1])) for st in gas[ HTSeq.GenomicInterval( "chr1", 250, 290, "." ) ].steps()])


###########
#идем по всем экзонам, у каждого получаем id гена,
#array[id Гена] += длина экзона(end - start)
#получим массив со всеми генами и длинной каждого гена(сумма длин экзонов)
#длина гена - 100%
#позиция выравнивания рида это x %

#если все сигарные строки пересекаются с 1 геном(проверить на выход из границ гена, что будет?), только после того как мы убедимся в этом, то берем позицию первой сигарной строки за начало рида и от этого считаем % позиции в риде!
# в строке 279 оринтировочно проверяется окончательно картирование рида на (пары?) только на 1 ген, надо еще проверить не вылезает ли он за пределы гена(провести примерный эксперимент!)

#на графике по оси - x будет точки процентов, по оси y сколько ридов легло с таким %, надо посмотреть как строиться этот график,
#можно разделить на интервалы и посмотреть сколько попадает точек в интервал
#чтобы определить покрытие надо находить пересекающиеся риды по интервалам




