

#TODO тесты на свой код, использовать фреймворк для этого!
#профилировщик для кода найти, чтобы определять что работает быстро а что нет
#инструменты для памяти какие есть?

#надо все задачи делать с расчетом на большие данные!,(большие файлы) не помещающиеся в оперативную память

#распаралеливание кода?

#будет внедряться нужный класс конфигурации и внем будет конкретная логика по работе со строкой

from common.files import FileProcessing, FileProcessing_by_chars

from String_Algorithms.Counting_DNA_Nucleotides import DNAFileProcessor
from String_Algorithms.Counting_DNA_Nucleotides import DataProvider
from String_Algorithms.Counting_DNA_Nucleotides.LogicProvider import LinesLogicProvider, CharsLogicProvider

dp = DataProvider
lp = LinesLogicProvider
fp = DNAFileProcessor(dp,lp)
processor = FileProcessing("text.txt",fp)# построчно читает большой файл, без загрузки всего файла в память todo проверить это на большом файле

processor.process_file()

lp.print_result()

cl = CharsLogicProvider
fp.logic_provider = cl
processor = FileProcessing_by_chars("text.txt",fp)

processor.process_file()

cl.print_result()













