class FileProcessing:

    def __init__(self, path,fp):

        self.file_path = path

        self.concrete_fileprocessor = fp

    def process_file(self):

        try:
            with open(self.file_path, 'r') as f:
                for line in f:
                    self.concrete_fileprocessor.process_data(line)
        except OSError:
            print('No such file.')



    #TODO распаралеливание чтения файла построчно?

    def write_to_file(self):
        try:
            with open(self.file_path, 'w') as f:
                self.concrete_fileprocessor.write_to_file(f)

        except OSError:
            print('No such file.')


class FileProcessing_by_chars(FileProcessing):

    def process_file(self):
        с = ""
        try:
            with open(self.file_path, 'r') as f:
                while True:
                    c = f.read(1)# 1 символ считываем
                    if not c:
                        print("End of file")
                        break

                    self.concrete_fileprocessor.process_data(c)

        except OSError:
                print('No such file.')


    #TODO распаралеливание чтения файла посимвольно?
