import random

class DataProvider:
    def __init__(self):
        pass
    def write_line(self,f):
        dict = {0: "A", 1: "C", 2: "G", 3: "T"}

        for i in range(1000000):
            for i in range(1000):
                f.write(dict[random.randint(0, 4)])
            f.write('\n')