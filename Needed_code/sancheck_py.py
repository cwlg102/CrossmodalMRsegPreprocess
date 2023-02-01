import os
import shutil
import numpy as np
from PIL import Image

bp = r"D:\!Crossmodality_2023\Contour_polyvol"
dp = r"D:\!Crossmodality_2023\Contour_polyvol_re"
plist = os.listdir(bp)

for i in range(len(plist)):
    dlist = os.listdir(os.path.join(bp, plist[i]))
    path2 = os.path.join(bp, plist[i])

    try: os.mkdir(r"D:\!Crossmodality_2023\Contour_polyvol_re\P%03d" %(i+1))
    except: pass

    for j, dest in enumerate(dlist):
        arr = np.load(os.path.join(path2, dest))
        np.save(os.path.join(dp + "/P%03d" %(i+1), "P%03d_%03dslice_ctr.npy" %(i+1, j+1)), arr)
        
        #arr = arr.astype("float64")
        #arr += 1024
        #arr /= 4095
        #arr *= 255
        #arr = arr.astype("uint8")
        #img = Image.fromarray(arr)
        #img.save(r"D:\!Crossmodality_2023\sancheck_CT\%d_%d.jpg" %(i+1, j+1))