#!/usr/bin/env python

import numpy as np
import scipy
import matplotlib.pyplot as plt
import imageio
import plotly
import plotly.express as px

###############
#Exercise 1
###############

week10_path = "/Users/cmdb/qbb2024-answers/week10/"
genes = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
fields = ["field0", "field1"]
stains = ["DAPI", "nascentRNA", "PCNA"]

rgbimg_dict = {}
j = 0
for gene in genes:

    for field in fields:
        rgbimg = np.zeros((520, 616, 3), np.uint16)

        for i, stain in enumerate(stains):
            img_path = week10_path + gene + "_" + field + "_" + stain + ".tif"
            rgbimg[:, :, i] = imageio.v3.imread(img_path)

        j += 1
        rgbimg_dict[f"rgbimg_{j}"] = rgbimg
        
plt.imshow(rgbimg_dict["rgbimg_8"][:, :, 0])
plt.show()
plt.imshow(rgbimg_dict["rgbimg_8"][:, :, 1])
plt.show()
plt.imshow(rgbimg_dict["rgbimg_8"][:, :, 2])
plt.show()