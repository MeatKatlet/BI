#класс - конкретная реализация логики чтения/записи в файл
class DNAFileProcessor:

    def __init__(self,data_provider,logic_provider):

        self.data_provider = data_provider#класс - провайдер данных для записи

        self.logic_provider = logic_provider

    def process_data(self, data):

        #здесь должна быть логика работы со строкой/символом
        if (data != "\n"):
            self.logic_provider.logic(data)


    def write_to_file(self,f):

        self.data_provider.write_line(f)



