import json
import kink_finder
from keras.models import Model
from keras.layers import Input, Dense, Flatten, Dropout
from keras.optimizers import Adam
from keras.initializers import RandomUniform
import numpy as np
from random import shuffle

def get_model(shape):

    # Model parameters
    input_shape = [shape]

    output_init = RandomUniform(minval=0,maxval = 0.2)

    hidden_size = 10

    inp = Input(shape=input_shape)
    hidden_1 = Dense(hidden_size, activation='relu')(inp)
    dropout_1 = Dropout(0.3)(hidden_1)
    out = Dense(2, activation='softmax',kernel_initializer = output_init)(dropout_1)

    model = Model(inputs=inp, outputs=out)

    print(model.summary())

    return model



data = json.load(open('data_big_succinate.json'))

y_plot = np.array(list(map(lambda i : i['obj1'],data)))
x_plot = np.array(list(map(lambda i : i['obj2'],data)))

x_lim = kink_finder.get_kink_point(x_plot,y_plot)[0]

unshuffled = []
for item in data:
	unshuffled.append(item)

shuffle(data)


gene_size = len(data[0]['gene_set'])

model = get_model(gene_size)

X = np.array(list(map(lambda i : i['gene_set'],data)))
Y = np.array(list(map(lambda i : [0,1] if i['obj2']>x_lim else [1,0], data)))


tot = len(data)
test_size = int(tot*0.1)

X_test = X[:test_size]
y_test = Y[:test_size]
X_train = X[test_size:]
y_train = Y[test_size:]




batch_size = 50
nb_epoch = 200



model.compile(loss='categorical_crossentropy', optimizer=Adam(),metrics=['accuracy'])

model.fit(X_train, y_train, batch_size=batch_size, epochs=nb_epoch,verbose=1, validation_data=(X_test, y_test))

score = model.evaluate(X_test, y_test, verbose=1)

model.save('NN_succinate_2.h5')

print("Accuracy:", score[1])
