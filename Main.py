import Parser
from FCGR import CGR
from RNN import RNN
import numpy as np
from keras.preprocessing.text import Tokenizer
from tensorflow.keras.utils import to_categorical
# You Can  Change the value of  k to 3 , 5 , 7 ................................................
k = 3
fcgr = CGR(1,"16S-reads.fa",k)
X_train = []
Y_train = []
Data = {}
fastafile = "16S-reads.fa"
excelfile = "type.xlsx"
# You Can Change the colnumber to 1,2,3,4
colnumber = 1
parser = Parser.Parser(fastafile,excelfile,colnumber)
list = parser.read_data()

for i in list:
    Data[str(i)]=parser.search_value(i)

f = open('output_col_{}'.format(colnumber)+'_k_{}.txt'.format(k),'w')
f.write(str(Data))
f.write('\n')

print(Data)

for i in Data:
    # You Can Change the method to split , here we use sliding method
    X_train.append(fcgr.fcgr_slide(i))
    Y_train.append(Data[i])
X_train = np.array(X_train,dtype=object)
Y_train = np.array(Y_train)
f.write(str(X_train))
f.write('\n')
f.write(str(Y_train))
f.write('\n')
print(X_train)
print(Y_train)


L = []
for list in X_train:

    L2 = []
    for kmer in list:
        sum = 0
        z = fcgr.our_frequency_chaos_game_representation(kmer)
        for i in range(len(z)):
            for j in range(len(z[i])):
                sum = sum+z[i][j]*(2**(i+j)/(i+k))
        L2.append(int(sum))
    L.append(L2)

if colnumber == 1:
    for i in range(len(Y_train)):
        if Y_train[i] == 'Alphaproteobacteria':
            Y_train[i] = 0
        elif Y_train[i] == 'Gammaproteobacteria':
            Y_train[i] = 2
        else:
            Y_train[i] = 1


# TODO: baraye sotoon haye 2,3,4 dg ham shartharo ezafe kon


print(len(X_train))
print(len(L))
print(X_train[0])
print(L[0])
print(Y_train)
f.write(str(len(X_train)))
f.write('\n')
f.write(str(len(L)))
f.write('\n')
f.write(str(X_train[0]))
f.write('\n')
f.write(str(L[0]))
f.write('\n')
f.write(str(Y_train))
f.write('\n')

# You Can Change the max_word size
max_word = 10000
tokenizer = Tokenizer(num_words=max_word)
X_train = tokenizer.sequences_to_matrix(L,mode='binary')
Y_train = to_categorical(Y_train,3)
X_test = X_train[19000:]
Y_test = Y_train[19000:]
X_train = X_train[0:19000]
Y_train = Y_train[0:19000]
f.write("this is X_train:")
print("this is X_train:")
f.write('\n')
print()
f.write(str(X_train))
print(X_train)
f.write('\n')
f.write("this is Y_train:")
print("this is Y_train:")
f.write('\n')
print()
f.write(str(Y_train))
f.write('\n')
print(Y_train)
f.write("this is x test :")
print("this is x test :")
f.write('\n')
f.write(str(X_test))
print(X_test)
print()
f.write('\n')
f.write("this is y test :")
print("this is y test :")
print()
f.write('\n')
f.write(str(Y_test))
print(Y_test)
f.write('\n')
f.write(str(X_train.shape))
f.write('\n')
f.write(str(Y_train.shape))
f.write('\n')
f.write(str(X_test.shape))
f.write('\n')
f.write(str(Y_test.shape))
f.write('\n')
print(X_train.shape)
print(Y_train.shape)
print(X_test.shape)
print(Y_test.shape)
rnn = RNN(X_train,Y_train,X_test,Y_test,max_word)
run = rnn.Make_RNN(f)
f.close()


