import os
from skimage import draw
import numpy as np
import pydicom
from PIL import Image

#폴더에 한 환자의 ct.dcm 과 해당 dcm에 맞는 rtst 한꺼번에 넣어놓으면 됨.
def getContourVolume(fol_path):
    dcm_path_list = os.listdir(fol_path)
    ct_list = []; rt_list = []
    
    for i, dcm_path in enumerate(dcm_path_list):
        
        dcm_data = pydicom.dcmread(os.path.join(fol_path, dcm_path), force=True)
        dcm_data.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
        try:
            if dcm_data.Modality == "RTSTRUCT":
                rt_list.append(dcm_data)
            elif dcm_data.Modality == 'RTPLAN':
                pass
            #Modality가 CT 혹은 MR일 경우
            else:
                ct_list.append((int(dcm_data.InstanceNumber), dcm_data))
        except:
            pass

    #CT 혹은 MR InstanceNumber에 따른 정렬
    ct_list.sort(key=lambda x: x[0])

    #가장 Inferior 쪽에 있는 CT
    xzp, yzp, zzp = ct_list[-1][1].ImagePositionPatient; zero_set = (float(xzp), float(yzp), float(zzp)) 
    #픽셀 공간으로 보내기 위한 scale value
    xs, ys = ct_list[-1][1].PixelSpacing; zs = ct_list[-1][1].SliceThickness; scale_set = (float(xs), float(ys), float(zs)) 
    img_orient = np.array(ct_list[-1][1].ImageOrientationPatient)
    
    axis_x = img_orient[:3] 
    axis_y = img_orient[3:] 
    axis_z = np.cross(axis_x, axis_y)
    #axis_z = axis_z / np.linalg.norm(axis_z)
    transform = np.vstack([axis_x, axis_y, axis_z]).T
    sca_mat = np.array([[xs, xs, xs], [ys, ys, ys], [zs, zs, zs]])
    # transform *= sca_mat
    transform_inv = np.linalg.inv(transform)
    #RTst에 있는 roi item 갯수와 roi item마다 존재하는 slice 갯수 가져옴
    item_num, roi_slice_num = getContourNumbers(rt_list[0])
    
    roi = rt_list[0].ROIContourSequence
    #ie이름 파악하기 위해서
    structureset = rt_list[0].StructureSetROISequence
    
    total_vol_val = {} #각 roi 마다 volume 몇인지 dictionary 형태로 저장
    tvv2 = {}
    I_NEED_IT = ["Cortex_R", "Cortex_L"] #여기에 있는 것들만 Volume 계산
    ctr_vol = np.zeros((200, ct_list[-1][1].Rows, ct_list[-1][1].Columns, 3), dtype="float32")
    for idx in range(item_num):
        
        NAME = str(structureset[idx].ROIName)
        
        #필요한 roi면 반복문 진행
        # if NAME in I_NEED_IT:
            # pass
        # else:continue
        numvol = 0
        max_slices = roi_slice_num[idx]
        # if NAME == "Cortex_R":
            # color_chan = [255, 0, 0]
        # elif NAME == "Cortex_L":
            # color_chan = [0, 0, 255]
        for jdx in range(max_slices):
            ctr_data = roi[idx].ContourSequence[jdx].ContourData
            color_chan = roi[idx].ROIDisplayColor
            
            ctr_int_data = ContourtoPixel(ctr_data, zero_set, scale_set, transform_inv)

    
            poly_r = []
            poly_c = []
            for i in range(len(ctr_int_data)-1):
                x1, y1, z1 = ctr_int_data[i]
                x2, y2, z2 = ctr_int_data[i+1]
                ctr_vol[z1][y1][x1] = color_chan
                draw_y, draw_x = draw.line(y1, x1, y2, x2)
                ctr_vol[z1][draw_y, draw_x] = color_chan
                # poly_r.append(y1)
                # poly_c.append(x1)

            # poly_r.append(y2)
            # poly_c.append(x2)
            # r_arr = np.array(poly_r)
            # c_arr = np.array(poly_c)
            

            #skimage-draw.polygon : 좋은 라이브러리긴 하나 느림

            # rr, cc = draw.polygon(r_arr, c_arr)
            # ctr_vol[z1][rr, cc] = color_chan
            
            

            
            draw_y, draw_x = draw.line(y2, x2, ctr_int_data[0][1], ctr_int_data[0][0])
            ctr_vol[z1][draw_y, draw_x] = color_chan
    ctr_vol = ctr_vol.astype("uint8")
    return ctr_vol

def getCTVolume(fol_path, RGB=False):
    dcm_path_list = os.listdir(fol_path)
    ct_list = []; rt_list = []
    
    for i, dcm_path in enumerate(dcm_path_list):
        
        dcm_data = pydicom.dcmread(os.path.join(fol_path, dcm_path), force=True)
        dcm_data.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
        try:
            if dcm_data.Modality == "RTSTRUCT":
                rt_list.append(dcm_data)
            elif dcm_data.Modality == 'RTPLAN':
                pass
            #Modality가 CT 혹은 MR일 경우
            else:
                ct_list.append((int(dcm_data.InstanceNumber), dcm_data))
        except:
            pass

    #CT 혹은 MR InstanceNumber에 따른 정렬
    temptemp = 200
    ct_list.sort(key=lambda x: -x[0])
    ct_vol = np.zeros((temptemp, 512, 512), dtype="float32")
    for i in range(len(ct_list)):
        temp_arr = \
        ct_list[i][1].pixel_array *\
        ct_list[i][1].RescaleSlope +\
        ct_list[i][1].RescaleIntercept
        
        temp_arr = temp_arr.astype("float32")
        ct_vol[i] = temp_arr
    win_min = -160
    win_max = 350
    ct_vol[ct_vol < -160] = -160
    ct_vol[ct_vol > 350] = 350
    ct_vol = 255*(ct_vol-win_min)/(win_max-win_min)
    ct_vol = ct_vol.astype("uint8")
    
    if RGB:
        ct_vol = np.stack((ct_vol, )*3, axis=-1)
    return ct_vol

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

def ContourtoPixel(ctr_data, zero_set, scale_set, transform_inv):
        ### Contour 좌표 inferior CT position(zero_set) 기준으로 픽셀 공간으로 변환 ###
        x_zr, y_zr , z_zr = zero_set
        x_sca, y_sca, z_sca = scale_set
        zr_arr = np.array(zero_set)
        sca_arr = np.array(scale_set)
        ctr_int = []
        
        for i in range(0, len(ctr_data), 3):
            raw_x, raw_y, raw_z = ctr_data[i], ctr_data[i+1], ctr_data[i+2]
            coord_inf = np.array([raw_x-x_zr, raw_y-y_zr, raw_z-z_zr])
            coord_inf /= sca_arr
            coord_ct = np.matmul(transform_inv, coord_inf)
            x = round(coord_ct[0])
            y = round(coord_ct[1])
            z = round(coord_inf[2])
            ctr_int.append((x, y, z))

        return ctr_int
bp = r"C:\Users\cwl98\Downloads\CT"
ctr_vol = getContourVolume(os.path.join(bp, "KIM SOON BONG_10399539_CT_2022-02-03_101434_Plann_(rotated).PRISM.ABD.Body.3.0....PRISM.ABD_n160__00000"))

ct_vol = getCTVolume(os.path.join(bp, "KIM SOON BONG_10399539_CT_2022-02-03_101434_Plann_(rotated).PRISM.ABD.Body.3.0....PRISM.ABD_n160__00000"), RGB=True)
total_vol = np.where(ctr_vol==[0, 0, 0], ct_vol, ctr_vol)
for i in range(len(total_vol)):
    img = Image.fromarray(total_vol[i])
    img.save(r"C:\Users\cwl98\Downloads\CT/%d.jpg" %i)
quit()

fol_paths = os.listdir(bp)

for i in range(113, len(fol_paths)):
    ctr_vol = getContourVolume(os.path.join(bp, fol_paths[i]))
    ct_vol = getCTVolume(os.path.join(bp, fol_paths[i]), RGB=True)
    total_vol = np.where(ctr_vol==[0, 0, 0], ct_vol, ctr_vol)
    for idx, arr in enumerate(total_vol):
        img = Image.fromarray(arr)
        img.save(os.path.join(r"D:\!Oncosoft\이식외과_done\ck_img_polygon", fol_paths[i] + "_slice%03d.jpg" %(idx + 1)))
                 
