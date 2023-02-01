import os
import numpy as np
import pydicom
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
    assert len(rtst_list) == 1, "number of rtst should be one"

    return rtst_list[0]

def getsortedCTdcm(dcm_fol_path):
    dcm_path_list = os.listdir(dcm_fol_path)
    ct_dcms = []
    for ct_dcm in dcm_path_list:
        #print(ct_dcm)
        if ct_dcm[-3:] == "dcm" or ct_dcm[-3:] == 'DCM':
            pass
        else:continue
        temp = pydicom.dcmread(os.path.join(dcm_fol_path, ct_dcm), force=True)
        
        
        if temp.Modality == 'RTPLAN' or temp.Modality == 'RTSTRUCT':
            continue
        #print(temp.ImagePositionPatient)
        ct_dcms.append((temp.InstanceNumber, temp))
        
    ct_dcms.sort(key=lambda x: x[0], reverse=True) #내림차순 정렬
    sorted_ct_dcms = []
    for i in range(len(ct_dcms)):
        sorted_ct_dcms.append(ct_dcms[i][1])
    
    return sorted_ct_dcms, sorted_ct_dcms[0]    

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

def setZeroandScale(inf_ct):
    xzp, yzp, zzp = inf_ct.ImagePositionPatient; zero_set = (float(xzp), float(yzp), float(zzp)) 
    xs, ys = inf_ct.PixelSpacing; zs = inf_ct.SliceThickness; scale_set = (float(xs), float(ys), float(zs)) 
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


def main_func(fol_path,ind):
    rtst = getRTST(fol_path)
    CT_list, inf_ct = getsortedCTdcm(fol_path)
    item_num, roi_slice_num = getContourNumbers(rtst)
    zero_set, scale_set = setZeroandScale(inf_ct)
    roi = rtst.ROIContourSequence
    structureset = rtst.StructureSetROISequence

    ctr_vol = np.zeros((inf_ct.InstanceNumber, inf_ct.Rows, inf_ct.Columns, 3), dtype="float64")
    ctr_vol_1c = np.zeros((inf_ct.InstanceNumber, inf_ct.Rows, inf_ct.Columns), dtype= "uint8")
    z_ctr_data = {"Femur_Head_L": [], "Femur_Head_R": []}
    for idx in range(item_num):
        max_slices = roi_slice_num[idx]
        
        if str(structureset[idx].ROIName)[0] == 'Z' or 'TV' in str(structureset[idx].ROIName) or 'nbl' in str(structureset[idx].ROIName):
            
            continue
        # if structureset[idx].ROIName == "Cortex_R" or structureset[idx].ROIName == "Cortex_L":
            # pass
        # else:
            # continue
        #print(structureset[idx].ROIName)
        for jdx in range(max_slices):
            ctr_data = roi[idx].ContourSequence[jdx].ContourData
            color_chan = roi[idx].ROIDisplayColor
            ctr_int_data = ContourtoPixel(ctr_data, zero_set, scale_set)
            if structureset[idx].ROIName == "Femur_Head_L":
                z_ctr_data["Femur_Head_L"].append(ctr_int_data)
            elif structureset[idx].ROIName == "Femur_Head_R":
                z_ctr_data["Femur_Head_R"].append(ctr_int_data)

            poly_r = []
            poly_c = []

            for i in range((len(ctr_int_data)-1)):

                x1, y1, z1 = ctr_int_data[i]
                x2, y2, z2 = ctr_int_data[i+1]
                
                draw_y, draw_x = draw.line(y1, x1, y2, x2)
                try:
                    ctr_vol[z1][draw_y, draw_x] = color_chan
                    if structureset[idx].ROIName == "Femur_Head_L":
                        ctr_vol_1c[z1][draw_y, draw_x] = 1
                    elif structureset[idx].ROIName == "Femur_Head_R":
                        ctr_vol_1c[z1][draw_y, draw_x] = 2

                except:
                    pass
                     

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
            draw_y, draw_x = draw.line(y2, x2, ctr_int_data[0][1], ctr_int_data[0][0])
            
            try:
                ctr_vol[z1][draw_y, draw_x] = color_chan
                if structureset[idx].ROIName == "Femur_Head_L":
                    ctr_vol_1c[z1][draw_y, draw_x] = 1
                elif structureset[idx].ROIName == "Femur_Head_R":
                    ctr_vol_1c[z1][draw_y, draw_x] = 2
            except:
                pass
        # pivot = []
        # if structureset[idx].ROIName == 'Femur_Head_L' or structureset[idx].ROIName == "Femur_Head_R":
            # sorted_ctr_int_data = sorted(ctr_int_data, key=lambda x: x[1])
            # pivot.append(sorted_ctr_int_data[0])
        # print(pivot)

            # print(ctr_int_data[0]) #여기서 호출하면 제일 머리 쪽임
    #################################################################################
    # Femur Head 기반 contour
    # volume 하고 질량중심잡으면됨... 이 생각을 왜못했지 ㅋㅋ

    assert z_ctr_data["Femur_Head_L"] != [] and z_ctr_data["Femur_Head_R"] != [], "이름 check"
    femur_l = np.transpose(np.where(ctr_vol_1c == 1), (1, 0))
    femur_r = np.transpose(np.where(ctr_vol_1c == 2), (1, 0))
    m_center_l = np.squeeze(np.sum(femur_l, axis=0, keepdims=True)//(len(femur_l)))
    m_center_r = np.squeeze(np.sum(femur_r, axis=0, keepdims=True)//(len(femur_r)))
    print(inf_ct.PixelSpacing, inf_ct.SliceThickness)
    print("z, y, x order:", m_center_l, m_center_r)
    #################################################################################
    
  
    #이동 좌표 구하기
    cur_im_center = (m_center_l[2]+m_center_r[2])//2, (m_center_l[1]+m_center_r[1])//2 # x y z
    dest_center = (255, 250)
    d_x, d_y = dest_center[0]-cur_im_center[0], dest_center[1]-cur_im_center[1]
    pad_param = 100

    #contour 이동 부분
    ctr_pad_vol = np.pad(ctr_vol, ((0, 0), (pad_param, pad_param), (pad_param, pad_param), (0, 0)), "constant")
    
    ctr_moved_vol = ctr_pad_vol[:, pad_param-d_y:pad_param-d_y+512, pad_param-d_x:pad_param-d_x+512, :]
    ctr_moved_vol = ctr_moved_vol.astype("uint8")

    #ct 이동 부분
    ct_vol = np.zeros((inf_ct.InstanceNumber, inf_ct.Rows, inf_ct.Columns), dtype="float32")
    ct_jpg_vol = np.zeros((inf_ct.InstanceNumber, inf_ct.Rows, inf_ct.Columns, 3), dtype="float64")
    for idx in range(len(CT_list)):
        CT_list[idx].file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
        pix_arr = CT_list[idx].pixel_array
        ct_pad_arr = np.pad(pix_arr, ((pad_param, pad_param), (pad_param, pad_param)), "constant" )
        
        ct_moved_arr = ct_pad_arr[pad_param-d_y:pad_param-d_y+512, pad_param-d_x:pad_param-d_x+512]

        rescale = CT_list[idx].RescaleSlope
        intersep = CT_list[idx].RescaleIntercept
        ct_moved_arr = ct_moved_arr * rescale + intersep
 
        rgb_arr = np.stack((ct_moved_arr,) *3, axis=-1)
        ct_vol[idx] = ct_moved_arr
        ct_jpg_vol[idx] = rgb_arr
        #print(np.max(ct_arr), np.min(ct_arr))
    img_max = 3071
    img_min = -1024
    w_max = 3071; w_min = -1024
    ct_vol[ct_vol>img_max] = img_max
    ct_vol[ct_vol<img_min] = img_min
    ct_jpg_vol[ct_jpg_vol>w_max] = w_max
    ct_jpg_vol[ct_jpg_vol<w_min] = w_min
    
    ct_jpg_vol = 255*(ct_jpg_vol - w_min)/(w_max-w_min)
    ct_jpg_vol = ct_jpg_vol.astype("uint8")

    total_vol = np.zeros_like(ct_jpg_vol)
    
    total_vol = np.where(ctr_moved_vol == [0, 0, 0], ct_jpg_vol, ctr_moved_vol)
    total_vol = total_vol.astype("uint8")

    rev_ct_vol = np.flip(ct_vol, axis=0)
    rev_t_vol = np.flip(total_vol, axis=0)
    
    #inf_ct로 slicethickne

    for i in range(len(rev_t_vol)):
        im = Image.fromarray(rev_t_vol[i])
        gray_arr = 255*(rev_ct_vol[i]+1024)/4095
        gray_arr = gray_arr.astype("uint8")
        gray_im = Image.fromarray(gray_arr)
        im.save(os.path.join(r"D:\!Crossmodality_2023\jpgtest","%d_%d.jpg" %(ind, i+1)))
        gray_im.save(os.path.join(r"D:\!Crossmodality_2023\test_jpg_gray", "%d_%d.jpg" %(ind, i+1)))

basepath = r"D:\!Crossmodality_2023\!Pelvic_Unity_Data"
file_list = os.listdir(basepath)

for i in range(len(file_list)):
    print(file_list[i])
    try:main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\%s\ct" %(file_list[i]), i+1)
    except:pass
    try:main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\%s\ct1" %(file_list[i]), i+1)
    except:pass
    try:main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\%s\ct2" %(file_list[i]), i+1)
    except:pass
    try:main_func(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\%s\ct3" %(file_list[i]), i+1)
    except:pass