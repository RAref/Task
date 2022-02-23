import xlrd
from Bio import SeqIO

class Parser():
    fastafile=""
    excel_file=""
    colnumber=1



    def __init__(self,fastafile,excel_file,colnumber):
        self.fastafile=fastafile
        self.excel_file=excel_file
        self.colnumber=colnumber


    def parser(self):
        list=[0,self.colnumber]
        dic={}
        file=xlrd.open_workbook(self.excel_file)
        sheet=file.sheet_by_index(0)
        for i in range(sheet.nrows):
            dic[sheet.cell_value(i,0)]=sheet.cell_value(i,self.colnumber)

        return dic


    def read_data(self):
        sequnce = SeqIO.parse(self.fastafile, "fasta")
        list = []
        for record in sequnce:
            list.append(record.seq)

        return list




    def Genarate_value(self,sequnce):
        file = SeqIO.parse(self.fastafile, "fasta")
        list = []
        value_list=[]
        for record in file:
            if sequnce==str(record.seq):
                list.append(record.description)


        for i in list:

            a = i.split("reference=")
            b = a[1].split(" ")
            value_list.append(b[0])

        return value_list


    def search_value(self,sequnce):
        genarte_value=self.Genarate_value(sequnce)
        dic=self.parser()
        for i in genarte_value:
            if i in dic:
               # print(i,dic[i])
                return dic[i]






