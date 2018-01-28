class LogicProvider:

    def logic(self,line):
        """
        A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains.

An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."

Given: A DNA string ss of length at most 1000 nt.

Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in ss.

Sample Dataset
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
Sample Output
20 12 17 21

все сортировки работают не меньше чем за n*log n, поэтому быстрее чем за n никак!
        """

#TODO использовать итераторы для экономии памяти!
#TODO для того чтобы считывать символы поштучно, занимая в памяти только оодин символ за раз надо делать это через генератор, в генераторе будет считывание посимвольно из строки(побайтова)


    def method1(self, line):
        #simple
        result = {"A" : 0, "C": 0, "G": 0, "T" : 0}

        try:
            for char in line:
                result[char] += 1

        except:
            print("Символ не из алфавита ДНК!")

        return " ".join(result)

    def method2(self,line):

        #распаралеливание строки на части!
        #

        pass





