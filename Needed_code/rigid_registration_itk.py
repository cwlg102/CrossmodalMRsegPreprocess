import os
import shutil
import numpy as np
import pydicom
import SimpleITK as sitk
from skimage import draw
from PIL import Image

def convert_to_sitk(img_volume, spacing, origin, direction=None):
    """convert numpy volume to sitk image"""
    # numpy [z, y, x] -> sitk [x, y, z]
    sitk_volume = sitk.GetImageFromArray(img_volume.astype(np.float64))
    sitk_volume.SetOrigin(origin)
    sitk_volume.SetSpacing(spacing)
    if direction:
        sitk_volume.SetDirection(direction)
    return sitk_volume

def convert_to_numpy(sitk_volume):
    """convert sitk image to numpy volume"""
    img_volume = sitk.GetArrayFromImage(sitk_volume)
    
    return img_volume

ct_path = r"D:\!Crossmodality_2023\CT_npy"
mr_path = r"D:\!Crossmodality_2023\new_MR_npy"
ct_path_list = os.listdir(ct_path)
mr_path_list = os.listdir(mr_path)

for idx in range(len(ct_path_list)):
    ct_path2 = os.path.join(ct_path, ct_path_list[idx])
    ct_dcm_list = os.listdir(ct_path2)
    mr_path2 = os.path.join(mr_path, mr_path_list[idx])
    mr_dcm_list = os.listdir(mr_path2)
    
    ct_vol = np.zeros((len(ct_dcm_list), 512, 512))
    mr_vol = np.zeros((len(mr_dcm_list), 512, 512))

    for i in range(len(ct_dcm_list)):
        ct_vol[i] = np.load(os.path.join(ct_path2, ct_dcm_list[i]))
    for i in range(len(mr_dcm_list)):
        mr_vol[i] = np.load(os.path.join(mr_path2, mr_dcm_list[i]))
    ct_vol = ct_vol.astype("int16")
    mr_vol = mr_vol.astype("int32")
    
    ct_itk = convert_to_sitk(ct_vol, [0.9765625, 0.9765625, 3], [0, 0, 0])
    mr_itk = convert_to_sitk(mr_vol, [0.9765625, 0.9765625, 3], [0, 0, 0])


    

    """1) Set ElastixImageFilter"""
    elastixImageFilter = sitk.ElastixImageFilter()

    """2) Set Parameters"""
    elastixImageFilter.SetFixedImage(ct_itk)
    elastixImageFilter.SetMovingImage(mr_itk)
    parametermap = sitk.GetDefaultParameterMap('translation')
    parametermap['MaximumNumberOfIterations'] = ['7000']
    parametermap['MaximumNumberOfSamplingAttempts'] = ['16']
    elastixImageFilter.SetParameterMap(parametermap) #method

    """3) Execute"""
    elastixImageFilter.Execute()
    resultimage = elastixImageFilter.GetResultImage()

    #save parameter map
    sitk.WriteParameterFile(elastixImageFilter.GetTransformParameterMap()[0], 
    r"D:\!Crossmodality_2023\parameter_map\%d.txt" %(idx+1))

    mr_regi_vol = convert_to_numpy(resultimage)
    try:
        os.mkdir(r"D:\!Crossmodality_2023\new_MR_npy_regi\P%03d" %(idx+1))
    except:
        pass
    for i, arr in enumerate(mr_regi_vol):
        np.save(os.path.join(r"D:\!Crossmodality_2023\new_MR_npy_regi\P%03d" %(idx+1), r"P%03d_%03d_MRslice.npy" %(idx+1, i+1)), arr)
    
    