#!/usr/bin/env python

import numpy as np
import scipy
import matplotlib.pyplot as plt
import imageio
import plotly
import plotly.express as px

from functions import find_labels, filter_by_size

###############
#Exercise 1
###############

week10_path = "/Users/cmdb/qbb2024-answers/week10/"
genes = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
fields = ["field0", "field1"]
stains = ["DAPI", "nascentRNA", "PCNA"]

#print(imageio.v3.imread('/Users/cmdb/qbb2024-answers/week10/APEX1_field0_DAPI.tif').shape)

rgbimg_dict = {}

for gene in genes:
    for field in fields:
        rgbimg = np.zeros((520, 616, 3), np.uint16)
        for i, stain in enumerate(stains):
            img_path = week10_path + gene + "_" + field + "_" + stain + ".tif"
            rgbimg[:, :, i] = imageio.v3.imread(img_path)
        condition = gene + '_' + field
        rgbimg_dict[condition] = rgbimg

###############
#Exercise 2
###############

#Step 2.1
dapi_mask_dict = {} #Dictionary storing masks of DAPI channels
for key, value in rgbimg_dict.items():
    dapi_mean = np.mean(value[:,:,0])
    dapi_mask = value[:,:,0] >= dapi_mean
    dapi_mask_dict[key] = dapi_mask

#Step 2.2
dapi_labels_dict = {}
for key, value in dapi_mask_dict.items():
    dapi_labels_dict[key] = find_labels(value)

#Step 2.3
#First filtering
dapi_labels_filtered_dict = {}
for key, value in dapi_labels_dict.items():
    minsize = 100
    maxsize = 1000000
    dapi_labels_filtered_dict[key] = filter_by_size(value, minsize, maxsize)

#Find label size mean and sd
dapi_labels_filtered_mean_std_dict = {}
for key, value in dapi_labels_filtered_dict.items():
    label_sizes = np.bincount(value.ravel())[1:]
    label_sizes_mean = np.mean(label_sizes)
    label_sizes_std = np.std(label_sizes)
    dapi_labels_filtered_mean_std_dict[key] = (label_sizes_mean, label_sizes_std)

#Second filtering
dapi_labels_filtered_dict2 = {}
for key, value in dapi_labels_filtered_dict.items(): #should this be after first filter or before?
    mean = dapi_labels_filtered_mean_std_dict[key][0]
    std = dapi_labels_filtered_mean_std_dict[key][1]
    minsize = mean - std
    maxsize = mean + std
    dapi_labels_filtered_dict2[key] = filter_by_size(value, minsize, maxsize)

###############
#Exercise 3
###############

condition_list = []
for key, value in rgbimg_dict.items():
    nucleus_list = []
    label_array = dapi_labels_filtered_dict2[key]
    
    for i in range(np.amax(label_array)+1):
        #skip label = 0 (background)
        if i == 0:
            continue
        #append gene name to the list
        list = []
        list.append(key.split('_')[0]) #only include gene name, not field
        
        where = np.where(label_array == i)
        #append nascent RNA mean to the list
        RNA_values = value[:,:,1]
        nucleus_i_RNA_mean = np.mean(RNA_values[where])
        list.append(nucleus_i_RNA_mean)
        #append PCNA mean to the list
        PCNA_values = value[:,:,2]
        nucleus_i_PCNA_mean = np.mean(PCNA_values[where])
        list.append(nucleus_i_PCNA_mean)
        #append ratio to the list
        log_ratio = np.log2(nucleus_i_RNA_mean/nucleus_i_PCNA_mean)
        list.append(log_ratio)
        #append the list to the master list
        nucleus_list.append(list)

    condition_list.append(nucleus_list)

output_path = '/Users/cmdb/qbb2024-answers/week10/output.txt'

with open(output_path, 'w') as file:
    file.write('Gene,nascentRNA,PCNA,ratio\n')
    for condition in condition_list:
        for nucleus in condition:
            condition = nucleus[0]
            RNA_mean = nucleus[1]
            PCNA_mean = nucleus[2]
            log2_ratio = nucleus[3]
            line = f"{condition},{RNA_mean},{PCNA_mean},{log2_ratio}\n"
            file.write(line)