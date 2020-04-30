#os.environ['THEANO_FLAGS'] = 'optimizer=None'
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten
from keras.utils import np_utils
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping
import numpy as np
import pandas as pd
import sys
import glob
import matplotlib.pyplot as plt
import time
import datetime

# fix random seed for reproducibility
seed = 7
np.random.seed(seed)

## load root files
import uproot
import h5py

filename = {}
upfile = {}
treename = {}

filename['sig'] = '../../../output_GluGluToHHTo2B2G_numEvent1000.root'
filename['bkg'] = '../../../output_diphojets_1bjet_numEvent5000.root'

upfile['sig'] = uproot.open(filename['sig'])
upfile['bkg'] = uproot.open(filename['bkg'])

treename['sig'] = 'tagsDumper/trees/GluGluToHHTo2B2G_13TeV_DoubleHTag_0'
treename['bkg'] = 'tagsDumper/trees/diphojets_1bjet_13TeV_DoubleHTag_0'

#print(upfile['sig'][treename['sig']].show())
#print(upfile['bkg'][treename['bkg']].show())


## Convert tree to pandas DataFrames
import pandas as pd

VARS = ['absCosTheta_gg','customLeadingPhotonIDMVA','weight'] # choose which vars to use 

df = {}
df['bkg'] = upfile['bkg'][treename['bkg']].pandas.df(branches=VARS)
df['sig'] = upfile['sig'][treename['sig']].pandas.df(branches=VARS)

# print first entry
print 'Printing the first entry of the signal dataframe'
print(df['sig'].iloc[:1])

# add isSignal variable
#df['sig']['isSignal'] = np.ones(len(df['sig'])) 
#df['bkg']['isSignal'] = np.zeros(len(df['bkg']))

## Define the model
# baseline keras model
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten
from keras.utils import np_utils

NDIM = len(VARS[:-1])
inputs = Input(shape=(NDIM,), name = 'input')  
outputs = Dense(1, name = 'output', kernel_initializer='normal', activation='sigmoid')(inputs)

# creae the model
model = Model(inputs=inputs, outputs=outputs)
# compile the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
# print the model summary
model.summary()


## Divide in training and test dataset
X_sig = df['sig'][VARS[:-1]].values
X_bkg = df['bkg'][VARS[:-1]].values
Y_sig = np.ones(len(X_sig))
Y_bkg = np.zeros(len(X_bkg))
W_sig = df['sig']['weight'].values
W_bkg = df['bkg']['weight'].values

X = np.vstack([X_sig, X_bkg])
Y = np.hstack([Y_sig, Y_bkg])
W = np.hstack([W_sig, W_bkg])

from sklearn.model_selection import train_test_split
# X = features for training, Y = expected outcome. In this case, Y = 0 (for signal) or 1 (for background)
X_train_val, X_test, Y_train_val, Y_test, W_train_val, W_test = train_test_split(X, Y, W, test_size=0.2, random_state=7)
print X_train_val.shape
print Y_train_val.shape


# preprocessing: standard scalar
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler().fit(X_train_val)
X_train_val = scaler.transform(X_train_val)
X_test = scaler.transform(X_test)


# early stopping callback
from keras.callbacks import EarlyStopping
early_stopping = EarlyStopping(monitor='val_loss', patience=10)

# model checkpoint callback
# this saves our model architecture + parameters into dense_model.h5
from keras.callbacks import ModelCheckpoint
model_checkpoint = ModelCheckpoint('dense_model.h5', monitor='val_loss', 
                                   verbose=0, save_best_only=True, 
                                   save_weights_only=False, mode='auto', 
                                   period=1)


## Train classifier
print 'Start training the classifier...'
start_time = time.time()

history = model.fit(X_train_val, 
                    Y_train_val, 
                    epochs=100, 
                    batch_size=1024, 
                    verbose=0, # switch to 1 for more verbosity 
                    callbacks=[early_stopping, model_checkpoint], 
                    validation_split=0.20)

print 'Training done. It took', time.time()-start_time, 'seconds.'


# Evaluation
print(">>> Computing AUC...")
from sklearn.metrics import roc_auc_score, roc_curve
pred = model.predict(X_test)
roc_auc = roc_auc_score(Y_test,pred)
fp, tp, th = roc_curve(Y_test, pred)
print("AUC score: " + str(roc_auc))

## Plot classifier output
fig, ax1 = plt.subplots(figsize=(7,6), dpi=100)
plt.hist(pred[Y_test==0], bins=50,density=True, label="false", histtype="step")
plt.hist(pred[Y_test==1], bins=50, density=True, label="true", histtype="step")
plt.legend()
plt.savefig('./DNN_output.pdf')

## plot loss vs epoch
plt.clf()
fig, ax1 = plt.subplots(figsize=(7,6), dpi=100)
plt.plot(history.epoch, history.history["val_loss"], label="validation loss")
plt.plot(history.epoch, history.history["loss"], label="training loss")
plt.legend(loc="upper right")
ax1.set_xlabel('epoch')
ax1.set_ylabel('loss')
plt.savefig('./DNN_loss_vs_epoch.pdf')

## plot accuracy vs epoch
plt.clf()
fig,ax1 = plt.subplots(figsize=(7,6))
ax1.plot(history.history['acc'], label='acc')
ax1.plot(history.history['val_acc'], label='val_acc')
ax1.legend(loc="upper left")
ax1.set_xlabel('epoch')
ax1.set_ylabel('accuracy')
plt.savefig('./DNN_accuracy_vs_epoch.pdf')

## plot roc
plt.clf()
fig,ax1 = plt.subplots(figsize=(7,6))
ax1.plot(fp, tp, lw=2, color='cyan', label='auc = %.3f' % (roc_auc))
ax1.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k', label='random chance')
ax1.set_xlim([0, 1.0])
ax1.set_ylim([0, 1.0])
ax1.set_xlabel('false positive rate')
ax1.set_ylabel('true positive rate')
ax1.set_title('receiver operating curve')
ax1.legend(loc="lower right")
plt.savefig('./DNN_roc.pdf')



