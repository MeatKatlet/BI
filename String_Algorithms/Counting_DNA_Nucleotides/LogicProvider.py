from math import floor
from multiprocessing import Pool, Lock, Array
import multiprocessing as mp

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


        #согластно ответу на SO память расти не должна

        if __name__ == '__main__':
            a = mp.Array('i', 4)
            self.counter_lock = Lock()
            p = Pool(processes=2,initializer=self.init, initargs=(a))
            p.map(self.line_worker, line,2) #здесь мы пока что только получаем статистику для 1 воркера

            p.close()
            p.join()

            dict = {"A" : a[0], "C": a[1], "G": a[2], "T" : a[3]}

            return dict


    def init(self, aa):
        global a
        a = aa


    def line_worker(self, line):
        result = self.method1(line)

        with self.counter_lock:
            a[0] += result["A"]
            a[1] += result["C"]
            a[2] += result["G"]
            a[3] += result["T"]

        # {"A" : 0, "C": 0, "G": 0, "T" : 0}

            # пока вижу 2 пути для распаралеливания
            # 1. каждый worker независимо производит свой результат, потом результаты суммиируются
            # 2. каждый worker разделяет shared_memory
            #   - manager(list, dict, и другие)
            #   - Value/Array + lock т.к. операции например += не атомарны, но как говорят на SO быстрее работает

            #написать все варианты для этого!
            # как передавать shared_memory object в worker



class CharsLogicProvider(LinesLogicProvider):

    def __init__(self):

        self.char_total_counter = {"A" : 0, "C": 0, "G": 0, "T" : 0}

    def logic(self,char):

        try:
            self.char_total_counter[char] += 1

        except:
            print("Not correct symbol!")

