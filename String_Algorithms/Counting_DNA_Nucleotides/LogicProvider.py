from math import floor
from multiprocessing import Pool

class LinesLogicProvider:

    def __init__(self):

        self.char_total_counter = {"A" : 0, "C": 0, "G": 0, "T" : 0}


    def print_result(self):
        print(" ".join(self.char_total_counter))

    def logic(self,line):
        """
            все сортировки работают не меньше чем за n*log n, поэтому быстрее чем за n никак!
        """

        result = self.method1(line)

        #для всего файла целиком считаем здесь
        self.char_total_counter["A"] += result["A"]
        self.char_total_counter["C"] += result["C"]
        self.char_total_counter["G"] += result["G"]
        self.char_total_counter["T"] += result["T"]

        return " ".join(result)#только для строки

    def method1(self, line):
        #simple
        result = {"A" : 0, "C": 0, "G": 0, "T" : 0}

        try:
            for char in line:
                result[char] += 1

            return result
        except Exception as e:
            print("Символ не из алфавита ДНК!")
            raise e

        #return " ".join(result)



    def method2(self,line):

        #TODO распаралеливание строки на части! в памяти
        #т.е. каждую строку в памяти паралельно!
        #middle = floor(len(line)/2);

        #простой вариант - делим строку пополам, на каждую часть будет worker
        #или делим на несколько частей
        #

        #task1 = line[0:middle]
        #task2 = line[middle:]

        #Это увеличивает память, надо наверное посимвольно отправлять буквы в worker ы


        if __name__ == '__main__':
            p = Pool(2)
            p.map(LinesLogicProvider.line_worker(), line) #здесь мы пока что только получаем статистику для 1 воркера

        def line_worker(line):
            self.method1(line)

            #как передать очередь?
            #как получить результат вычислений наружу
            #TODO как получить результат вычислений наружу и передать его обратно для следующих вычислений?




class CharsLogicProvider(LinesLogicProvider):

    def __init__(self):

        self.char_total_counter = {"A" : 0, "C": 0, "G": 0, "T" : 0}

    def logic(self,char):

        try:
            self.char_total_counter[char] += 1

        except:
            print("Not correct symbol!")

