class FileProcessing:

    def __init__(self, path,fp):

        self.file_path = path

        self.concrete_fileprocessor = fp

    def process_file(self):
        with open(self.file_path, 'a') as f:
            for line in f:
                self.concrete_fileprocessor.process_line(line)

    def write_to_file(self):
        with open(self.file_path, 'w') as f:
            self.concrete_fileprocessor.write_to_file(f)



