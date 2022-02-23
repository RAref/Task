import math
import numpy as np
import pylab
from matplotlib import cm
import matplotlib.pyplot as plt
from  Bio import SeqIO
from collections import defaultdict
import sys

class CGR():
    k = 0
    Data = ""
    i = 0
    fasta_file=""

    def empty_dict(self):
        """
    	None type return vessel for defaultdict
    	:return:
    	"""
        return None

    def __init__(self, a, fasta_file,k):
        self.fasta_file=fasta_file
        self.i = a
        self.k=k


    def read_data(self):
        sequnce = SeqIO.parse(self.fasta_file, "fasta")
        list = []
        for record in sequnce:
            list.append(record.seq)

        return list

    # base on splite for our fcgr
    def fcgr_splite(self ,sequnce):

        j = 0
        new_list=[]
        for i in range(math.ceil(len(sequnce) / self.k)):
            new_list.append(sequnce[j:j+self.k])
            j=j+self.k
        return new_list

    # base on slide for our fcgr
    def fcgr_slide(self ,sequence):
        new_list=[]
        for i in range(len(sequence)-(self.k-1)):
                new_list.append(sequence[i:i+self.k])

        return new_list



    #base on slide for classic fcgr
    def count_kmers(self, sequence):
        d = {}
        for i in range(len(sequence) - (self.k - 1)):
            if str(sequence)[i:i + self.k] not in d:
                d[str(sequence)[i:i + self.k]]=0
            d[str(sequence)[i:i + self.k]] += 1
        d.pop("N", None)
        print(d)
        input("Press Enter to continue...")
        print()
        return d

    #base on split for classic fcgr
    def count_kmers2(self, sequence):
        d = {}
        j=0
        for i in range(math.ceil(len(sequence)/self.k)):
            if str(sequence)[j:j + self.k] not in d:
                d[str(sequence)[j:j + self.k]]=0

            d[str(sequence)[j:j + self.k]] += 1
            j=j+self.k

        d.pop("N", None)
        print(d)
        input("Press Enter to continue...")
        print()
        return d


    def probabilities(self, kmer_count):
        probabilities = {}
        N = len(self.Data)
        for key, value in kmer_count.items():
            probabilities[key] = float(value) / (N - self.k + 1)
        return probabilities

    def classic_frequency_chaos_game_representation(self, probabilities):
        array_size = int(math.sqrt(4 ** self.k))
        chaos = []
        for i in range(array_size):
            chaos.append([0] * array_size)
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
        for key, value in probabilities.items():
            for char in key:
                if char == "T":
                    posx += maxx / 2
                elif char == "C":
                    posy += maxy / 2
                elif char == "G":
                    posx += maxx / 2
                    posy += maxy / 2
                maxx /= 2
                maxy /= 2

            chaos[int(posy - 1)][int(posx - 1)] = value
            maxx = array_size
            maxy = array_size
            posx = 1
            posy = 1
            print("this is the kmers:")
            print()
            print(chaos)
            input("Press Enter to continue...")
            print()
        m = float(np.amax(chaos))

        c = np.array(chaos) / m
        print("Ending.................................")
        print("result of given sequnce : ")
        return c

    def our_frequency_chaos_game_representation(self,kmer):
        CGR_X_MAX = 1 * (2 ** (self.k - 1))
        CGR_Y_MAX = 1 * (2 ** (self.k - 1))
        CGR_X_MIN = -1 * (2 ** (self.k - 1))
        CGR_Y_MIN = -1 * (2 ** (self.k - 1))
        CGR_A = (CGR_X_MIN, CGR_Y_MAX)
        CGR_T = (CGR_X_MIN, CGR_Y_MIN)
        CGR_G = (CGR_X_MAX, CGR_Y_MIN)
        CGR_C = (CGR_X_MAX, CGR_Y_MAX)
        CGR_CENTER = (0, 0)
        CGR_DICT = defaultdict(
            self.empty_dict,
            [
                ('A', CGR_A),  # Adenine
                ('T', CGR_T),  # Thymine
                ('G', CGR_G),  # Guanine
                ('C', CGR_C),  # Cytosine
                ('U', CGR_T),  # Uracil demethylated form of thymine
                ('a', CGR_A),  # Adenine
                ('t', CGR_T),  # Thymine
                ('g', CGR_G),  # Guanine
                ('c', CGR_C),  # Cytosine
                ('u', CGR_T)  # Uracil/Thymine
            ]
        )

        cgr = []
        cgr_marker = CGR_CENTER[:]
        for s in kmer:
            cgr_corner = CGR_DICT[s]
            if cgr_corner:
                cgr_marker = (
                    (cgr_corner[0] + cgr_marker[0]) / 2,
                    (cgr_corner[1] + cgr_marker[1]) / 2
                )
                cgr.append([s, cgr_marker])
            else:
                sys.stderr.write("Bad Nucleotide: " + s + " \n")

        fcgr = np.zeros(shape=(2 ** self.k, 2 ** self.k))
        size = 2 ** self.k
        for list in cgr:
            x = list[1][0]
            y = list[1][1]
            i, j = math.ceil((size / 2 - y)), math.ceil((size / 2 + x))
            fcgr[i - 1][j - 1] = 1
        return fcgr



    #base on slide
    def our_fcgr_load_fasta(self):
        data=self.read_data()
        for i in data:
            print("your sequnce is :")
            print(i)
            sequnce_slides=self.fcgr_slide(i)
            for kmer in sequnce_slides:
                result=self.our_frequency_chaos_game_representation(kmer)
                yield {kmer:result}



    #base on splite

    def our_fcgr_load_fasta2(self):
        data = self.read_data()
        for i in data:
            print("your sequnce is :")
            print(i)
            sequnce_splites = self.fcgr_splite(i)
            for kmer in sequnce_splites:
                print(kmer)
                result = self.our_frequency_chaos_game_representation(kmer)
                yield {kmer: result}




    #base on slide
    def classic_fcgr_load_fasta(self):

        data = self.read_data()
        for sequnce in data:
            self.Data = sequnce
            f4 = self.count_kmers(sequnce)
        

            f4_prob = self.probabilities(f4)
            chaos_k4 = self.classic_frequency_chaos_game_representation(f4_prob)


            yield chaos_k4



    #base on split

    def classic_fcgr_load_fasta2(self):
        data=self.read_data()
        for sequnce in data:
            self.Data = sequnce
            f4 = self.count_kmers2(sequnce)

            f4_prob = self.probabilities(f4)
            chaos_k4 = self.classic_frequency_chaos_game_representation(f4_prob)

            yield chaos_k4


    def fcgr_show(self):
        pylab.figure(figsize=(8, 8))
        pylab.title('CGR of ' + str(self.K) )
        pylab.imshow(self.c, cmap=cm.gray_r)  # ,interpolation = "spline36")
        pylab.savefig(str(self.i) + ".PNG")
        pylab.show()

    def liner_show(self):
        y=np.array(self.c)
        plt.plot(y,marker='o')
        plt.show()












