import numpy as np
from tensorflow import keras
from tensorflow.keras import layers
class nn_db:
    def __init__(self,qrange,minvalu,maxvalu):
        self.all_x = []
        self.all_y = []
        self.model = keras.Sequential([])
        self.minvalu = np.array(minvalu)
        self.maxvalu = np.array(maxvalu)

        self.qrange = qrange

    def build_model(self,nnode,input_size):
        model = keras.Sequential([])
        model.add(layers.Dense(nnode[0],activation='relu',input_shape =
                               [input_size]))       
        for i in range(1,len(nnode)):
            model.add(layers.Dense(nnode[i],activation='relu'))
            model.add(layers.Dense(1))
            model.compile(optimizer='adam',loss='mse',metrics=["mae"])
        self.model = model

    def add_data(self,new_params,new_iq,data_qrange):
        for pi,to_push_x in enumerate(new_params):
            to_push_y = new_iq[pi]
            for qi,q in enumerate(data_qrange):
                self.all_x.append(np.append(to_push_x,q))
                self.all_y.append(to_push_y[qi])


    def preprocess(self,xs):

        minvalu_wq = np.append(self.minvalu,self.qrange[0])
        maxvalu_wq = np.append(self.maxvalu,self.qrange[-1])
        return np.array([(x-minvalu_wq)/(maxvalu_wq-minvalu_wq) for x in xs])
    
    def train_model(self,epoch,validation_split,minvalu,maxvalu,preshuffle = True):
        all_x_array = np.array(self.all_x)
        all_y_array = np.array(self.all_y)
        for i in range(len(all_y_array)):
            if all_y_array[i] < 1e-9:
                all_y_array[i] = 1e-9
        
        all_x_array = self.preprocess(all_x_array)
        all_y_array = np.log10(all_y_array)

        if preshuffle:
            rng_state = np.random.get_state()
            np.random.shuffle(all_x_array)
            np.random.set_state(rng_state)
            np.random.shuffle(all_y_array)

        self.model.fit(all_x_array,all_y_array,epochs = epoch, 
                       validation_split = validation_split,
                       shuffle = True,verbose = 0)



