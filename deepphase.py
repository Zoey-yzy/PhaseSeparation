
import os
import numpy as np
import pandas as pd
import keras as K
from keras.preprocessing.image import ImageDataGenerator
from keras.preprocessing import image
from keras.layers import Activation,Dropout,Flatten,Dense,Convolution2D,MaxPooling2D,GlobalAveragePooling2D
from keras.applications import ResNet50
from keras.models import Model,Sequential
import requests


def greenonly(x):
    return x*[0.,1.,0.]

#Setting training and validation data

sz = 512
batch_size = 64
train_data_dir = 'data_train/'
validation_data_dir = 'data_test/'
train_datagen = ImageDataGenerator(rescale=1./255,
                rotation_range=180,
                horizontal_flip = True)
                                   # ,preprocessing_function=greenonly)
test_datagen = ImageDataGenerator(rescale=1./255)
                                  #,preprocessing_function=greenonly)
train_generator = train_datagen.flow_from_directory(train_data_dir,shuffle = True,
                                                    target_size = (sz,sz),
                                                   batch_size = batch_size,class_mode = 'binary')
validation_generator = test_datagen.flow_from_directory(validation_data_dir,
                                                       shuffle = False,
                                                       target_size = (sz,sz),
                                                       batch_size = batch_size,
                                                       class_mode = 'binary')

# if K.backend.image_dim_ordering() == 'th':
if K.backend.image_data_format() == 'channels_first':
    input_shape = (3, sz, sz)
else:
    input_shape = (sz, sz, 3)

#
#Building CNN

model = Sequential()
model.add(Convolution2D(32, (3, 3), input_shape=input_shape))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))

model.add(Convolution2D(32, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))

model.add(Convolution2D(32, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))

model.add(Convolution2D(32, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))

model.add(Convolution2D(64, (3, 3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))

model.add(GlobalAveragePooling2D())
model.add(Dense(64))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(1))
model.add(Activation('sigmoid'))

model.compile(loss='binary_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])

# Or using pretrained weight
model.load_weights('pretrained.h5')

import requests
import pandas as pd 
def download_image(url, save_path):
    try:
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(save_path, 'wb') as file:
                for chunk in response.iter_content(1024):
                    file.write(chunk)
            print(f"Downloaded {save_path}")
        else:
            print(f"Failed to download {url}")
    except Exception as e:
        print(f"Error downloading {url}: {e}")

"""
import pandas as pd
import numpy as np
urls0 = pd.read_csv("/home/dinglab/Desktop/yzy/IF/image_url.csv")
for ii in range(urls0.shape[0]):
    for jj in range(urls0.shape[1]):
        if type(urls0.iloc[ii,jj])!=str:
            if np.isnan(urls0.iloc[ii,jj]):
                urls0.iloc[ii,jj] = 0
import operator as op
urls = []
for ii in range(len(urls0)):
    url1 = list(urls0.iloc[ii,:])[2:]
    name = urls0.iloc[ii,1]
    site0 = [x for x in url1 if type(x)==str]
    site = [x for x in site0 if op.contains(x,'blue_red_green.jpg')]
    if site:
        site.insert(0,name)
        urls.append(site)
df = pd.DataFrame(urls)
df.to_csv("/home/dinglab/Desktop/yzy/IF/IFimage_url.csv")
"""
urls = pd.read_csv("/Users/zoey/Desktop/PS/IF-deepphase/IFimage_url.csv")
dir = "/Users/zoey/Desktop/PS/IF-deepphase/picture"
for kk in range(len(urls)):
    urlp = urls.iloc[kk,2:].dropna()
    for ii in range(len(urlp)):
        url = urlp[ii]
        download_image(url, dir + "/" + urls.iloc[kk,1] + "_" + str(ii) + ".jpg")

    # Calculating scores for all images
    file=os.listdir('data_real/').copy()
    img=[]
    for i in file:
        try:
            x = image.load_img('data_real/'+i, target_size=(512,512))
            x = image.img_to_array(x)
            x = greenonly(x)
            x = x.reshape((1,) + x.shape)
            x = x/255.
            img.append([i,model.predict(x)[0][0]])
        except:pass
    pd.DataFrame(img).to_csv(urls.iloc[kk,1] + '.csv')
    for filename in os.listdir(dir):
        file_path = os.path.join(dir, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))