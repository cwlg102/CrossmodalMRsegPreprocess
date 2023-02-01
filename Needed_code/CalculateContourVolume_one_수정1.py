import os
from skimage import draw
import numpy as np
import pydicom
from PIL import Image
import cv2
#폴더에 한 환자의 ct.dcm 과 해당 dcm에 맞는 rtst 한꺼번에 넣어놓으면 됨.
def CalContourVolSize(fol_path):
    dcm_path_list = os.listdir(fol_path)
    ct_list = []; rt_list = []
    
    for i, dcm_path in enumerate(dcm_path_list):
        
        dcm_data = pydicom.dcmread(os.path.join(fol_path, dcm_path), force=True)
        dcm_data.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
        
        if dcm_data.Modality == "RTSTRUCT":
            rt_list.append(dcm_data)
        elif dcm_data.Modality == 'RTPLAN':
            pass
        #Modality가 CT 혹은 MR일 경우
        else:
            ct_list.append((int(dcm_data.InstanceNumber), dcm_data))

    #CT 혹은 MR InstanceNumber에 따른 정렬
    ct_list.sort(key=lambda x: x[0])

    #가장 Inferior 쪽에 있는 CT
    xzp, yzp, zzp = ct_list[-1][1].ImagePositionPatient; zero_set = (float(xzp), float(yzp), float(zzp)) 
    #픽셀 공간으로 보내기 위한 scale value
    xs, ys = ct_list[-1][1].PixelSpacing; zs = ct_list[-1][1].SliceThickness; scale_set = (float(xs), float(ys), float(zs)) 

    #RTst에 있는 roi item 갯수와 roi item마다 존재하는 slice 갯수 가져옴
    item_num, roi_slice_num = getContourNumbers(rt_list[0])
    
    roi = rt_list[0].ROIContourSequence
    #ie이름 파악하기 위해서
    structureset = rt_list[0].StructureSetROISequence

    total_vol_val = {} #각 roi 마다 volume 몇인지 dictionary 형태로 저장
    I_NEED_IT = ["Cortex_R", "Cortex_L"] #여기에 있는 것들만 Volume 계산

    for idx in range(item_num):
        ctr_vol = np.zeros((ct_list[-1][1].InstanceNumber, ct_list[-1][1].Rows, ct_list[-1][1].Columns), dtype="float64")
        NAME = str(structureset[idx].ROIName)
        
        #필요한 roi면 반복문 진행
        if NAME in I_NEED_IT:
            pass
        else:continue
        numvol = 0
        max_slices = roi_slice_num[idx]
        for jdx in range(max_slices):
            ctr_data = roi[idx].ContourSequence[jdx].ContourData
            #color_chan = roi[idx].ROIDisplayColor
            
            ctr_int_data = ContourtoPixel(ctr_data, zero_set, scale_set)
            poly_r = []
            poly_c = []
            for i in range(len(ctr_int_data)-1):
                x1, y1, z1 = ctr_int_data[i]
                x2, y2, z2 = ctr_int_data[i+1]
                #draw_y, draw_x = draw.line(y1, x1, y2, x2)
                #ctr_vol[z1][draw_y, draw_x] = 1
                poly_r.append(y1)
                poly_c.append(x1)
            poly_r.append(y2)
            poly_c.append(x2)
            r_arr = np.array(poly_r)
            c_arr = np.array(poly_c)
            
            #skimage-draw.polygon : 좋은 라이브러리긴 하나 느림

            rr, cc = draw.polygon(r_arr, c_arr)
            ctr_vol[z1][rr, cc] = 255
            
            arr = ctr_vol[z1]
            arr = arr.astype("uint8")
            #arr_im = Image.fromarray(arr)
            #ret,thresh = cv2.threshold(arr_im,127,255,0)
            contours,hierarchy = cv2.findContours(arr, cv2.RETR_TREE, \
                                            cv2.CHAIN_APPROX_SIMPLE)
            
            for i, contour in enumerate(contours):
                M = cv2.moments(contour)
                h = np.squeeze(hierarchy)
                if hierarchy[0][i][3] == -1:
                    numvol+= M["m00"]
                else:
                    numvol -= M["m00"]
            
            #draw_y, draw_x = draw.line(y2, x2, ctr_int_data[0][1], ctr_int_data[0][0])
            #ctr_vol[z1][draw_y, draw_x] = 1
        num_vox = len(np.where(ctr_vol==255)[0])
        #mm^3 -> cc로 변환
        vol_size = num_vox*scale_set[0]*scale_set[1]*scale_set[2]/1000 
        total_vol_val[str(structureset[idx].ROIName)] = vol_size
        print(numvol*scale_set[0]*scale_set[1]*scale_set[2]/1000)
    print(total_vol_val) 

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
            for j in range(1000):
                try: rtst.ROIContourSequence[i].ContourSequence[j]
                except:break
            roi_slice_num.append(j)

        return item_num, roi_slice_num

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
            

CalContourVolSize(r"E:\test_이식외과")
