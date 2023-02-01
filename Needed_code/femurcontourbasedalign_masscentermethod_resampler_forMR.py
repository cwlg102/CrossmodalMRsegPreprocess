import os
import numpy as np
import pydicom
import matplotlib.pyplot as plt
import SimpleITK as sitk
from skimage import draw
from PIL import Image

def getRTST(dcm_fol_path):
    dcm_path_list = os.listdir(dcm_fol_path)

    rtst_list = []
    for rtst_dcm in dcm_path_list:
        if rtst_dcm[-3:] == "dcm" or rtst_dcm[-3:] == 'DCM':
            pass
        else:continue
        temp = pydicom.dcmread(os.path.join(dcm_fol_path, rtst_dcm), force=True)
        if temp.Modality == 'RTSTRUCT':
            rtst_list.append(temp)
    # assert len(rtst_list) == 1, "number of rtst should be one"

    return rtst_list[0]

def getsortedMRdcm(dcm_fol_path): #reverse: False
    dcm_path_list = os.listdir(dcm_fol_path)
    mr_dcms = []
    for mr_dcm in dcm_path_list:
        #print(mr_dcm)
        if mr_dcm[-3:] == "dcm" or mr_dcm[-3:] == 'DCM':
            pass
        else:continue
        temp = pydicom.dcmread(os.path.join(dcm_fol_path, mr_dcm), force=True)
        
        
        if temp.Modality == 'RTPLAN' or temp.Modality == 'RTSTRUCT':
            continue
        #print(temp.ImagePositionPatient)
        mr_dcms.append((temp.InstanceNumber, temp))
        
    mr_dcms.sort(key=lambda x: x[0]) #오름차순 정렬
    sorted_mr_dcms = []
    for i in range(len(mr_dcms)):
        sorted_mr_dcms.append(mr_dcms[i][1])
    #last_idx = max(sorted_mr_dcms[0].InstanceNumber, sorted_mr_dcms[-1].InstanceNumber)
    return sorted_mr_dcms, sorted_mr_dcms[0], len(sorted_mr_dcms)   

def getContourNumbers(rtst):
        ###총 Item이 몇개 있는지 가져옴###
        item_num = 0
        while True:
            try:
                rtst.ROIContourSequence[item_num]
                item_num += 1
            except:
                break

        roi_slice_num = []

        ### Item마다 slice 몇개 있는지 가져옴 ###
        for i in range(item_num):
            for j in range(10000):
                try: rtst.ROIContourSequence[i].ContourSequence[j]
                except:break
            roi_slice_num.append(j)

        return item_num, roi_slice_num

def setZeroandScale(sorted_mr_dcms):
    inf_ct = sorted_mr_dcms[0]
    sup_ct = sorted_mr_dcms[-1]
    x1, y1, z1 = inf_ct.ImagePositionPatient
    x2, y2, z2 = sup_ct.ImagePositionPatient
    real_world_mm = abs(z1-z2)

    xzp, yzp, zzp = inf_ct.ImagePositionPatient; zero_set = (float(xzp), float(yzp), float(zzp)) 
    xs, ys = inf_ct.PixelSpacing; zs = real_world_mm/len(sorted_mr_dcms); scale_set = (float(xs), float(ys), float(zs)) 
    
    return zero_set, scale_set

def ContourtoPixel(ctr_data, zero_set, scale_set):
        ### Contour 좌표 inferior CT position(zero_set) 기준으로 픽셀 공간으로 변환 ###
        x_zr, y_zr , z_zr = zero_set
        x_sca, y_sca, z_sca = scale_set
        ctr_int = []
        for i in range(0, len(ctr_data), 3):
            raw_x, raw_y, raw_z = ctr_data[i], ctr_data[i+1], ctr_data[i+2]
            x = round((raw_x - x_zr)/x_sca)
            y = round((raw_y - y_zr)/y_sca)
            z = round((raw_z - z_zr)/z_sca)
            ctr_int.append((x, y, z))

        return ctr_int

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

def resample(sitk_volume, new_spacing, new_size, default_value=0):
    """1) Create resampler"""
    resample = sitk.ResampleImageFilter() 
    
    """2) Set parameters"""
    #set interpolation method, output direction, default pixel value
    resample.SetInterpolator(sitk.sitkLinear)
    resample.SetOutputDirection(sitk_volume.GetDirection())
    resample.SetDefaultPixelValue(default_value)
    
    #set output spacing
    new_spacing = np.array(new_spacing)
    resample.SetOutputSpacing(new_spacing)
    
    #set output size and origin
    old_size = np.array(sitk_volume.GetSize())
    old_spacing = np.array(sitk_volume.GetSpacing())
    new_size_no_shift = np.int16(np.ceil(old_size*old_spacing/new_spacing))
    old_origin = np.array(sitk_volume.GetOrigin())
    
    shift_amount = np.int16(np.floor((new_size_no_shift - new_size)/2))*new_spacing
    new_origin = old_origin + shift_amount
    
    new_size = [int(s) for s in new_size]
    resample.SetSize(new_size)
    resample.SetOutputOrigin(new_origin)
    
    """3) execute"""
    new_volume = resample.Execute(sitk_volume)
    return new_volume

def main_func(fol_path,ind):
    # rtst = getRTST(fol_path)
    MR_list, inf_ct, last_idx = getsortedMRdcm(fol_path)
    # item_num, roi_slice_num = getContourNumbers(rtst)
    zero_set, scale_set = setZeroandScale(MR_list)
    
    # roi = rtst.ROIContourSequence
    # structureset = rtst.StructureSetROISequence

    # ctr_vol = np.zeros((last_idx, inf_ct.Rows, inf_ct.Columns, 3), dtype="float64")
    # ctr_vol_1c = np.zeros((last_idx, inf_ct.Rows, inf_ct.Columns), dtype= "uint8")
    # z_ctr_data = {"Femur_Head_L": [], "Femur_Head_R": []}
    # for idx in range(item_num):
        # max_slices = roi_slice_num[idx]
        
        # if str(structureset[idx].ROIName)[0] == 'Z' or 'TV' in str(structureset[idx].ROIName) or 'nbl' in str(structureset[idx].ROIName):
            
            # continue
        # # if structureset[idx].ROIName == "Cortex_R" or structureset[idx].ROIName == "Cortex_L":
            # # pass
        # # else:
            # # continue
        # #print(structureset[idx].ROIName)
        # for jdx in range(max_slices):
            # ctr_data = roi[idx].ContourSequence[jdx].ContourData
            # color_chan = roi[idx].ROIDisplayColor
            # ctr_int_data = ContourtoPixel(ctr_data, zero_set, scale_set)
            # if structureset[idx].ROIName == "Femur_Head_L":
                # z_ctr_data["Femur_Head_L"].append(ctr_int_data)
            # elif structureset[idx].ROIName == "Femur_Head_R":
                # z_ctr_data["Femur_Head_R"].append(ctr_int_data)

            # poly_r = []
            # poly_c = []

            # for i in range((len(ctr_int_data)-1)):

                # x1, y1, z1 = ctr_int_data[i]
                # x2, y2, z2 = ctr_int_data[i+1]
                

                # draw_y, draw_x = draw.line(y1, x1, y2, x2)
                
                # ctr_vol[z1][draw_y, draw_x] = color_chan
                # if structureset[idx].ROIName == "Femur_Head_L":
                    # ctr_vol_1c[z1][draw_y, draw_x] = 1
                # elif structureset[idx].ROIName == "Femur_Head_R":
                    # ctr_vol_1c[z1][draw_y, draw_x] = 2

                     

                # poly_r.append(y1)
                # poly_c.append(x1)
            # poly_r.append(y2)
            # poly_c.append(x2)
            # r_arr = np.array(poly_r)
            # c_arr = np.array(poly_c)
            
            #draw polygon

            # rr, cc = draw.polygon(r_arr, c_arr)
            # ctr_vol[z1][rr, cc] = color_chan 
            # if structureset[idx].ROIName == "Femur_Head_L":
                # ctr_vol_1c[z1][rr, cc] = 1
            # elif structureset[idx].ROIName == "Femur_Head_R":
                # ctr_vol_1c[z1][rr, cc] = 2
            # draw_y, draw_x = draw.line(y2, x2, ctr_int_data[0][1], ctr_int_data[0][0])
            
            
            # ctr_vol[z1][draw_y, draw_x] = color_chan
            # if structureset[idx].ROIName == "Femur_Head_L":
                # ctr_vol_1c[z1][draw_y, draw_x] = 1
            # elif structureset[idx].ROIName == "Femur_Head_R":
                # ctr_vol_1c[z1][draw_y, draw_x] = 2
            
        # pivot = []
        # if structureset[idx].ROIName == 'Femur_Head_L' or structureset[idx].ROIName == "Femur_Head_R":
            # sorted_ctr_int_data = sorted(ctr_int_data, key=lambda x: x[1])
            # pivot.append(sorted_ctr_int_data[0])
        # print(pivot)

            # print(ctr_int_data[0]) #여기서 호출하면 제일 머리 쪽임
    #################################################################################
    # Femur Head 기반 contour
    # volume 하고 질량중심잡으면됨... 이 생각을 왜못했지 ㅋㅋ
    
    # 없는 MR 도있음...
    # assert z_ctr_data["Femur_Head_L"] != [] and z_ctr_data["Femur_Head_R"] != [], "이름 check"
    # femur_l = np.transpose(np.where(ctr_vol_1c == 1), (1, 0))
    # femur_r = np.transpose(np.where(ctr_vol_1c == 2), (1, 0))
    # m_center_l = np.squeeze(np.sum(femur_l, axis=0, keepdims=True)//(len(femur_l)))
    # m_center_r = np.squeeze(np.sum(femur_r, axis=0, keepdims=True)//(len(femur_r)))
    # print(scale_set)
    # print("z, y, x order:", m_center_l, m_center_r)
    # #################################################################################
    
  
    # #이동 좌표 구하기
    # cur_im_center = (m_center_l[2]+m_center_r[2])//2, (m_center_l[1]+m_center_r[1])//2 # x y z
    # dest_center = (239, 234)
    # d_x, d_y = dest_center[0]-cur_im_center[0], dest_center[1]-cur_im_center[1]
    # pad_param = 100

    # #contour 이동 부분
    # ctr_pad_vol = np.pad(ctr_vol, ((0, 0), (pad_param, pad_param), (pad_param, pad_param), (0, 0)), "constant")
    
    # ctr_moved_vol = ctr_pad_vol[:, pad_param-d_y:pad_param-d_y+inf_ct.Rows, pad_param-d_x:pad_param-d_x+inf_ct.Columns, :]
    # ctr_moved_vol = ctr_moved_vol.astype("uint8")

    #ct 이동 부분
    mr_vol = np.zeros((last_idx, inf_ct.Rows, inf_ct.Columns), dtype="float32")
    mr_jpg_vol = np.zeros((last_idx, inf_ct.Rows, inf_ct.Columns, 3), dtype="float64")
    for idx in range(len(MR_list)):
        MR_list[idx].file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
        # pix_arr = MR_list[idx].pixel_array
        # mr_pad_arr = np.pad(pix_arr, ((pad_param, pad_param), (pad_param, pad_param)), "constant" )
        
        # mr_moved_arr = mr_pad_arr[pad_param-d_y:pad_param-d_y+inf_ct.Rows, pad_param-d_x:pad_param-d_x+inf_ct.Columns]
        mr_moved_arr = MR_list[idx].pixel_array
        rescale = 1
        intersep = 0
        mr_moved_arr = mr_moved_arr * rescale + intersep
 
        rgb_arr = np.stack((mr_moved_arr,) *3, axis=-1)
        mr_vol[idx] = mr_moved_arr
        mr_jpg_vol[idx] = rgb_arr
        #print(np.max(ct_arr), np.min(ct_arr))
    
    mr_one_dim = np.ravel(mr_vol)

    #####normalization#####
    max_pixel_value = 2000
    mrnorm_hyperparam = 1000
    counts_list, bin_locations, patches = plt.hist(mr_one_dim, max_pixel_value, (0, max_pixel_value))

    
    for idx_val in range(max_pixel_value-1, -1, -1):
        if counts_list[idx_val] > mrnorm_hyperparam:
            val_norm = idx_val+1
            
            # norm_val_list.append(val_norm)
            break
    
    mr_vol = np.where(mr_vol > val_norm, val_norm, mr_vol)
    mr_jpg_vol = np.where(mr_jpg_vol > val_norm, val_norm, mr_jpg_vol)
    
    mr_jpg_vol = 255*(mr_jpg_vol - np.min(mr_jpg_vol))/(np.max(mr_jpg_vol)-np.min(mr_jpg_vol))
    mr_jpg_vol = mr_jpg_vol.astype("uint8")

    total_vol = np.zeros_like(mr_jpg_vol)
    
    #total_vol = np.where(ctr_moved_vol == [0, 0, 0], mr_jpg_vol, ctr_moved_vol)
    #total_vol = total_vol.astype("uint8")

    rev_mr_vol = np.flip(mr_vol, axis=0)
    #rev_t_vol = np.flip(total_vol, axis=0)
    
    #inf_ct로 slicethickne
    sitk_vol = convert_to_sitk(rev_mr_vol, spacing=
    [
        inf_ct.PixelSpacing[0], 
        inf_ct.PixelSpacing[1], 
        scale_set[2]

        ], 
        origin=[0, 0, 0])
    # sitk_color_vol = convert_to_sitk(rev_t_vol, spacing=
    # [
        # inf_ct.PixelSpacing[0], 
        # inf_ct.PixelSpacing[1], 
        # scale_set[2]

        # ], 
        # origin=[0, 0, 0])

    rescaled_sitk_vol = resample(sitk_vol, [0.9765625, 0.9765625, 3], [512, 512, len(rev_mr_vol)])
    #rescaled_color_sitk_vol = resample(sitk_color_vol, [0.9765625, 0.9765625, 3], [512, 512, len(rev_mr_vol)])
    
    
    rev_mr_vol = convert_to_numpy(rescaled_sitk_vol)
    rev_mr_vol = rev_mr_vol.astype("int16")
    #rev_t_vol = convert_to_numpy(rescaled_color_sitk_vol)
    #rev_t_vol = rev_t_vol.astype("uint8")

    try:
        os.mkdir(r"D:\!Crossmodality_2023\new_MR_npy\P%03d" %(ind))
    except:
        pass

    for i in range(len(rev_mr_vol)):
        if np.all(rev_mr_vol[i] == 0):
            continue

        #im = Image.fromarray(rev_mr_vol[i])
        np.save(r"D:\!Crossmodality_2023\new_MR_npy/P%03d/P%03d_%03dslice.npy" %(ind, ind, i+1) ,rev_mr_vol[i])
        gray_arr = 255*(rev_mr_vol[i]-np.min(rev_mr_vol))/(np.max(rev_mr_vol)-np.min(rev_mr_vol))
        gray_arr = gray_arr.astype("uint8")
        gray_im = Image.fromarray(gray_arr)
        #im.save(os.path.join(r"D:\!Crossmodality_2023\mrjpgtest_scaling","%d_%d.jpg" %(ind, i+1)))
        gray_im.save(os.path.join(r"D:\!Crossmodality_2023\mrjpg_gray_scaling", "%d_%d.jpg" %(ind, i+1)))
        

basepath = r"D:\!Crossmodality_2023\!Pelvic_Unity_Data"
file_list = os.listdir(basepath)

# for i in range(len(file_list)):
    # print(file_list[i])
    # try:main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\%s\mr1" %(file_list[i]), i+1); print(1)
    # except:pass
    # try:main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\%s\mr1 ats" %(file_list[i]), i+1);print(1)
    # except:pass
    # try:main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\%s\mr1 rtst" %(file_list[i]), i+1);print(1)
    # except:pass
    # try:main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\%s\mr10 ats" %(file_list[i]), i+1);print(1)
    # except:pass
    # print(0)
# main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\20210916 kyc\mr3", 14)
# main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\20210923 lka\mr11", 18)