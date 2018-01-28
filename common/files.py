class FileProcessing:

    def __init__(self, path,fp):

        self.file_path = path

        self.concrete_fileprocessor = fp

    def process_file(self):
        with open(self.file_path) as f:
            for line in f:
                self.concrete_fileprocessor.process_line(line)


