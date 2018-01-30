#класс - конкретная реализация логики чтения/записи в файл
class DNAFileProcessor:

    def __init__(self,data_provider,logic_provider):

        self.data_provider = data_provider#класс - провайдер данных для записи

        self.logic_provider = logic_provider

    def process_data(self, data):

        #здесь должна быть логика работы со строкой/символом, с учетом того что в памяти они одновременно не встречаются, а частями
        if (data != "\n"):
            self.logic_provider.logic(data)

        #подумать как делать посимвольное считывание файла - не будет ли это медленно работать или все таки лучше построчно, по сути генератор и не нужен если мы посимвольно собираемся читать
    def write_to_file(self,f):
        #генерить данные на лету
        #либо заранее подготовленный список данных использовать
        #Это зависит от класса провайдера

        data = self.data_provider.get_line()

        f.write(data+'\n')

