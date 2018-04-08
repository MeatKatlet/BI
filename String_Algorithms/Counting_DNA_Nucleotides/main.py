
from common.files import FileProcessing, FileProcessing_by_chars

from String_Algorithms.Counting_DNA_Nucleotides.DNAFileProcessor import DNAFileProcessor
from String_Algorithms.Counting_DNA_Nucleotides.DataProvider import DataProvider
from String_Algorithms.Counting_DNA_Nucleotides.LogicProvider import LinesLogicProvider, CharsLogicProvider , ParalellLogicProvider



if __name__ == '__main__':

    file_path = "../files/text.txt"

    #*******settings**************************************************
    dp = DataProvider
    lines_logic_provider = LinesLogicProvider
    fp = DNAFileProcessor(dp,lines_logic_provider)
    processor = FileProcessing(file_path, fp)# построчно читает большой файл, без загрузки всего файла в память
    #*********write to file*******************************************

    if os.path.exists(file_path)==False:
        processor.write_to_file()
    # ****************************************************************
    processor.process_file()

    lines_logic_provider.print_result()
    #*********settings************************************************
    char_logic_provider = CharsLogicProvider
    fp.logic_provider = char_logic_provider
    processor = FileProcessing_by_chars(file_path, fp)# читает файл посимвольно

    #******************************************************************
    processor.process_file()

    char_logic_provider.print_result()
    #****settings for paralell ****************************************

    fp.logic_provider = ParalellLogicProvider

    processor = FileProcessing(file_path, fp)
    #******************************************************************
    processor.process_file()

    fp.logic_provider.print_result()












