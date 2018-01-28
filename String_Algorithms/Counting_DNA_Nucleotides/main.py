

#TODO тесты на свой код, использовать фреймворк для этого!
#профилировщик для кода найти, чтобы определять что работает быстро а что нет

#надо все задачи делать с расчетом на большие данные!,(большие файлы) не помещающиеся в оперативную память

#распаралеливание кода?

#TODO написать универсальный класс для чтения больших файлов , чтобы использовать его в любых задачах(сделать это через dependency injection)
#будет внедряться нужный класс конфигурации и внем будет конкретная логика по работе со строкой

from common.files import FileProcessing
from String_Algorithms.Counting_DNA_Nucleotides import DNAFileProcessor
from String_Algorithms.Counting_DNA_Nucleotides import DataProvider
from String_Algorithms.Counting_DNA_Nucleotides import LogicProvider

dp = DataProvider
lp = LogicProvider
fp = DNAFileProcessor(dp,lp)

FileProcessing("text.txt",fp)# построчно читает большой файл, без загрузки всего файла в память todo проверить это на большом файле











