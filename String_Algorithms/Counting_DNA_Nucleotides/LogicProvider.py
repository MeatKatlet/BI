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

        except:
            print("Символ не из алфавита ДНК!")

        #return " ".join(result)
        return result

    def method2(self,line):

        #TODO распаралеливание строки на части! в памяти
        #т.е. каждую строку в памяти паралельно!

        pass


class CharsLogicProvider(LinesLogicProvider):

    def __init__(self):

        self.char_total_counter = {"A" : 0, "C": 0, "G": 0, "T" : 0}

    def logic(self,char):

        try:
            self.char_total_counter[char] += 1

        except:
            print("Not correct symbol!")

