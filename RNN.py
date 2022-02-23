from keras.models import Sequential
from keras.layers import Dense , LSTM , Dropout , Activation
from keras.losses import categorical_crossentropy


class RNN():

    def __init__(self, X_train,Y_train,X_test,Y_test,max_word):
        self.X_train=X_train
        self.Y_train=Y_train
        self.X_test=X_test
        self.Y_test=Y_test
        self.max_word=max_word

    def Make_RNN(self,file):
        model=Sequential()
        model.add(Dense(512,input_shape=(self.max_word,)))
        model.add(Activation('relu'))
        model.add(Dropout(0.5))
        model.add(Dense(3))
        model.add(Activation('softmax'))
        model.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])
        batch_size=32
        epoch=19000
        train_network = model.fit(self.X_train,self.Y_train,batch_size=batch_size,epochs=epoch,verbose=1,validation_split=0.1)
        score=model.evaluate(self.X_test,self.Y_test,batch_size=batch_size,verbose=1)
        print('Test loss:{}'.format(score[0]))
        print('Test accuracy: {}'.format(score[1]))
        predict_test_data = model.predict(self.X_test)
        print(self.Y_test)
        print(predict_test_data)

        file.write('loss:'+str(train_network.history['loss']))
        file.write('\n')
        file.write('accuracy:'+str(train_network.history['accuracy']))
        file.write('\n')
        file.write('val_loss:'+str(train_network.history['val_loss']))
        file.write('\n')
        file.write('val_accuracy:'+str(train_network.history['val_accuracy']))
        file.write('\n')
        file.write('Test loss:{}'.format(score[0]))
        file.write('\n')
        file.write('Test accuracy: {}'.format(score[1]))
        file.write('\n')
        file.write(str(self.Y_test))
        file.write('\n')
        file.write(str(predict_test_data))






