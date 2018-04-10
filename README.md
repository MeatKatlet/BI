# BI

\htseq_processing\plot_coverage.py - plot normalized coverage of genes, write raw absolute coverage 
data to file. It uses binary search to find intersection of reads on 10-percent interval control points. 
10 percent points of genes length. 

example of command run, arguments is the same as in HTseq-count, in this example we handle 4 bam files



python2.7 /home/kirill/bi/BI/htseq_processing/plot_coverage.py \
--stranded=no -f bam -m intersection-strict --order=pos \
/home/kirill/bi/transcript/Mock1-1/sorted_primary_proper.bam \
/home/kirill/bi/transcript/Mock2-1/sorted_primary_proper.bam \
/home/kirill/bi/transcript/ZIKV1-1/sorted_primary_proper.bam \
/home/kirill/bi/transcript/ZIKV2-1/sorted_primary_proper.bam \
/home/kirill/bi/mapping/HC2.gff3


\htseq_processing\merge_data_and_plot.py

merge txt files, produced by plot_coverage.py and plot all them on one plot


\String_Algorithms\Counting_DNA_Nucleotides\main.py - count nucliotides in different ways, process data by lines, 
by chars, and parallel reading of line
