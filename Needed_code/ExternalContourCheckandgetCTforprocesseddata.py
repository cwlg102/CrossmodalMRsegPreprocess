import os
import numpy as np
import pydicom
from skimage import draw
from PIL import Image
import nibabel as nib
from skimage.morphology import convex_hull_image
slash = '/'
#SFD_path = r"D:\brain_segmentation\!Sorted Full Dataset"
#SFD_list = os.listdir(SFD_path)
#label_table = {"brain_stem": 1, "spinal_cord": 2, "spinal cord": 2, "R_cochlea": 3, "L_cochlea": 4, "R_TMJ": 5, "L_TMJ": 6, "R_parotidG": 7, "L_parotidG": 8,
#"oral_cavity": 9, "mandible":10, "R_SMG":11, "L_SMG": 12, "pharynx": 13, "phaynx": 13, "larynx": 14, "thyroid": 15, "esophagus": 16, "R_eye": 17, "L_eye": 18,
#"R_lens": 19, "L_lens": 20, "R_optic_nerve": 21, "R_optic nerve": 21, "L_optic_nerve": 22, "L_optic nerve": 22, "optic_chiasm": 23, "optic_pathway": 24, "Line 1" : 0, "Line 2" : 0}
def run(ct_fol_path, rtst_fol_path,run_order, pnum):

    # ct_fol_path = r"D:\MRSIM_abdomen_data\Patient_data\001_5978364\CT_aligned\CHANG MYUNG HWAN_5978364_CT_2022-02-07_133956_Planning.CT.Abd.+.Pelvis.(contrast).carbon-311R_Align_n50__00000"
    # rtst_fol_path = r"D:\MRSIM_abdomen_data\Patient_data\001_5978364\MR_RTst\CHANG MYUNG HWAN_5978364_RTst_2022-02-07_140653_MR.limited.Abdomen.(noncontrast)_._n1__00000"

    ct_dcm_path = os.listdir(ct_fol_path)
    rtst_dcm_path = os.listdir(rtst_fol_path)

    def getsortedCTdcm(ct_fol_path ,ct_dcm_path):
        ct_dcms = []
        for ct_dcm in ct_dcm_path:
            if 'dcm' not in ct_dcm:
                continue
            temp = pydicom.dcmread(os.path.join(ct_fol_path, ct_dcm))
            #print(temp.ImagePositionPatient)
            ct_dcms.append((temp.InstanceNumber, temp))

        ct_dcms.sort(key=lambda x: x[0], reverse=True)
        sorted_ct_dcms = []
        for i in range(len(ct_dcms)):

            sorted_ct_dcms.append(ct_dcms[i][1])
        
        return sorted_ct_dcms, sorted_ct_dcms[0]

    ct_dcms, inf_ct = getsortedCTdcm(ct_fol_path, ct_dcm_path)

    x_zr, y_zr , z_zr = inf_ct.ImagePositionPatient[0], inf_ct.ImagePositionPatient[1], inf_ct.ImagePositionPatient[2]
    x_sca, y_sca, z_sca = inf_ct.PixelSpacing[0], inf_ct.PixelSpacing[1], inf_ct.SliceThickness
    
    zero_set = (x_zr, y_zr, z_zr)
    scale_set = (x_sca, y_sca, z_sca)
    print(x_zr, y_zr, z_zr)
    print(x_sca, y_sca, z_sca)
    #################################################################
    for file in rtst_dcm_path:
        if 'dcm' not in file:
            continue
        rtst_path = file
    rtst = pydicom.dcmread(os.path.join(rtst_fol_path, rtst_path))

    def getContourNumbers(rtst):
        item_num = 0
        while True:
            try:
                rtst.ROIContourSequence[item_num]
                item_num += 1
            except:
                break

        roi_slice_num = []

        for i in range(item_num):
            for j in range(1000):
                try: rtst.ROIContourSequence[i].ContourSequence[j]
                except:break
            roi_slice_num.append(j)

        return item_num, roi_slice_num

    item_num, roi_slice_num = getContourNumbers(rtst)

    #################################################################

    roi = rtst.ROIContourSequence
    structureset = rtst.StructureSetROISequence

    def ContourtoPixel(ctr_data, zero_set, scale_set):
        x_zr, y_zr , z_zr = zero_set
        x_sca, y_sca, z_sca = scale_set

        ctr_int = []
        for i in range(0, len(ctr_data), 3):
            raw_x, raw_y, raw_z = ctr_data[i], ctr_data[i+1], ctr_data[i+2]
            #print("raw:", raw_x, raw_y, raw_z)
            #print("zero_point", x_zr, y_zr, z_zr)
            #print("scale", x_sca, y_sca, z_sca)
            x = round((raw_x - x_zr)/x_sca)
            y = round((raw_y - y_zr)/y_sca)
            z = round((raw_z - z_zr)/z_sca)
            #print(x, y, z)
            ctr_int.append((x, y, z))
    
        return ctr_int
    def ContourtoPixelraw(ctr_data):
        ctr_raw = []
        for i in range(0, len(ctr_data), 3):
            ctr_raw.append((ctr_data[i], ctr_data[i+1], ctr_data[i+2]))
        return ctr_raw

    ctr_vol = np.zeros((inf_ct.InstanceNumber, inf_ct.Rows, inf_ct.Columns, 3), dtype="float64")
    ctr_gvol = np.zeros((inf_ct.InstanceNumber, inf_ct.Rows, inf_ct.Columns), dtype="float64")
    for idx in range(item_num):
        
        if structureset[idx].ROIName == "External" or structureset[idx].ROIName == 'external':
            pass
        else:
            continue
        max_slices = roi_slice_num[idx]
        for jdx in range(max_slices):
            ctr_data = roi[idx].ContourSequence[jdx].ContourData
            color_chan = roi[idx].ROIDisplayColor
            #color_chan = label_table[str(structureset[idx].ROIName)]
            
            ctr_int_data = ContourtoPixel(ctr_data, zero_set, scale_set)
            # ctr_raw_data = ContourtoPixelraw(ctr_data)
            poly_r = []
            poly_c = []
            sx, sy, sz = ctr_int_data[0]
            if sz > inf_ct.InstanceNumber:
                break

            for i in range((len(ctr_int_data)-1)):

                x1, y1, z1 = ctr_int_data[i]
                x2, y2, z2 = ctr_int_data[i+1]
                # Px, Py, Pz = ctr_raw_data[i]
                draw_y, draw_x = draw.line(y1, x1, y2, x2)
                ctr_vol[z1][draw_y, draw_x] = color_chan
                ctr_gvol[z1][draw_y, draw_x] = 1
                #Xx, Yx, Xy, Yy, Xz, Yz = ct_dcms[z1].ImageOrientationPatient
                #Sx, Sy, Sz = ct_dcms[z1].ImagePositionPatient
                #
                #del_i = x_sca
                #del_j = y_sca
                #cvt_mat = np.array([[Xx*del_i, Yx*del_j, 0, Sx], 
                #                    [Xy*del_i, Yy*del_j, 0, Sy],
                #                    [Xz*del_i, Yz*del_j, 0, Sz],
                #                    [0    ,    0,        0, 1]])
                #inv_cvt_mat = np.linalg.pinv(cvt_mat)
                #target = np.dot(inv_cvt_mat, np.array([Px, Py, Pz, 1]))
                #print(target)                

                poly_r.append(y1)
                poly_c.append(x1)
            poly_r.append(y2)
            poly_c.append(x2)
            r_arr = np.array(poly_r)
            c_arr = np.array(poly_c)
            
            
            rr, cc = draw.polygon(r_arr, c_arr)
            ctr_vol[z1][rr, cc] = color_chan
            ctr_gvol[z1][rr, cc] = 1
            draw_y, draw_x = draw.line(y2, x2, ctr_int_data[0][1], ctr_int_data[0][0])
            ctr_vol[z1][draw_y, draw_x] = color_chan
            ctr_gvol[z1][rr, cc] = 1
            

    ctr_vol[ctr_vol<0] = 0
    ctr_vol = ctr_vol.astype("uint8")
    ctr_gvol = ctr_gvol.astype("uint8")
    #######################draw contour on CT########################################
    ct_vol = np.zeros((inf_ct.InstanceNumber, inf_ct.Rows, inf_ct.Columns), dtype="float32")
    ct_jpg_vol = np.zeros((inf_ct.InstanceNumber, inf_ct.Rows, inf_ct.Columns, 3), dtype="float64")
    
    print(inf_ct.InstanceNumber)
    print(len(ct_dcms))
    for idx in range(len(ct_dcms)):
        ct_arr = ct_dcms[idx].pixel_array
        rescale = ct_dcms[idx].RescaleSlope
        intersep = ct_dcms[idx].RescaleIntercept
        ct_arr = ct_arr * rescale + intersep
        circle_r, circle_c = draw.disk((255, 255), 252)
        mask = np.zeros_like(ct_arr)
        mask[circle_r, circle_c] = 1
        
        ct_arr = np.where(mask==1, ct_arr, -1000)

        #img_max = 1000
        #img_min = -1000
        #ct_arr[ct_arr>img_max] = img_max
        #ct_arr[ct_arr<img_min] = img_min
        #ct_arr += 2047 - abs(np.min(ct_arr))
        #max_ = np.max(ct_arr)
        #min_ = np.min(ct_arr)
        #ct_arr = 255*(ct_arr - min_)/(max_-min_)
        #ct_arr = ct_arr.astype("uint8")
        #ct_jpg_vol[idx] = ct_arr
        rgb_arr = np.stack((ct_arr,) *3, axis=-1)
        ct_vol[idx] = ct_arr
        ct_jpg_vol[idx] = rgb_arr

    img_max = 1200
    img_min = -1100
    ct_jpg_vol[ct_jpg_vol>img_max] = img_max
    ct_jpg_vol[ct_jpg_vol<img_min] = img_min
    ct_jpg_vol += 2047 - abs(np.min(ct_jpg_vol))
    max_ = np.max(ct_jpg_vol)
    min_ = np.min(ct_jpg_vol)
    ct_jpg_vol = 255*(ct_jpg_vol - min_)/(max_-min_)
    ct_jpg_vol = ct_jpg_vol.astype("uint8")

    

    total_vol = np.zeros_like(ct_jpg_vol)
    #total_vol = np.where(ctr_vol != [0, 0, 0], ctr_vol, ct_jpg_vol)
    total_vol = np.where(ctr_vol != [0, 0, 0], ct_jpg_vol, 0)
    ct_vol = np.where(ctr_gvol != 0, ct_vol, -1000)
    ###################################################################################

    ct_256_vol = np.zeros((inf_ct.InstanceNumber, 256, 256))
    ctr_256_vol = np.zeros((inf_ct.InstanceNumber, 256, 256))
    
    for i in range(len(total_vol)):
        if np.all(total_vol[i] == 0) :
            continue
        #np.save(r"D:\DLnetwork\traindata_dicomct_each\each_patients\test_slices\P%03d_CT_slice%03d.npy" %(pnum, len(total_vol)-i), ct_vol[i])
        #img = Image.fromarray(total_vol[i])
        #img.save(r"D:\DLnetwork\traindata_dicomct_each\each_patients\tpic\P%03d_CT_slice%03d.jpg" %(pnum, len(total_vol)-i))
        
        np.save(r"D:\DLnetwork\traindata_dicomct_each\ctr_vol\npy\P%03d_CT_ctr%03d.npy" %(pnum, len(total_vol)-i), ctr_gvol[i])
        ctr_gvol[i] *= 255
        ctr_gimg = Image.fromarray(ctr_gvol[i])
        ctr_gimg.save(r"D:\DLnetwork\traindata_dicomct_each\ctr_vol\jpg\P%03d_CT_ctr%03d.jpg" %(pnum, len(total_vol)-i))
    # ct_vol = np.transpose(ct_vol, (2, 1, 0)) 
    # ct_vol = np.rot90(ct_vol, 2, (1, 2))

    # ctr_vol = np.transpose(ctr_vol, (2, 1, 0))
    # ctr_vol = np.rot90(ctr_vol, 2, (1, 2))

    

    #label = nib.Nifti1Image(ctr_vol, np.eye(4))
    #img = nib.Nifti1Image(ct_vol, np.eye(4))
    #img.header.get_xyzt_units()
    #label.header.get_xyzt_units()
    #img.to_filename(r"D:\brain_segmentation\ImT/la_%03d_0000.nii.gz" %(index+1))
    #label.to_filename(r"D:\brain_segmentation\laT/la_%03d.nii.gz" %(index+1))

#for i in range(len(ctr_vol)):
#    a= Image.fromarray(ctr_vol[i])
#    a.save(r"D:\brain_segmentation\lab\%d.jpg" %i)

#import shutil
#basepath = r"D:\DLnetwork\OncostudioContour_exceptP032\OncostudioContour_exceptP032"
#
#file_paths = os.listdir(basepath)
#for i in range(len(file_paths)):
#    if file_paths[i] == "test":
#        target_path = os.path.join(basepath, file_paths[i])
#        break
#
#rt_paths = os.listdir(target_path)
#for i in range(len(rt_paths)):
#    num = int(rt_paths[i][1:4])
#    try:os.mkdir(os.path.join(basepath, "test_pr") + "/P%03d" %(num))
#    except:pass
#    shutil.copy(os.path.join(target_path ,rt_paths[i]), os.path.join(basepath, "test_pr") + "/P%03d" %(num))

ct_path = r"D:\DLnetwork\traindata_dicomct_each\each_patients\test"
rt_path = r"D:\DLnetwork\traindata_dicomct_each\OncostudioContour_exceptP032\OncostudioContour_exceptP032\test_pr"
ct_paths = os.listdir(ct_path)
rt_paths = os.listdir(rt_path)

for i in range(len(ct_paths)):
    ct_fol_path = os.path.join(ct_path, ct_paths[i])
    rt_fol_path = os.path.join(rt_path, rt_paths[i])
    
    pnum = int(ct_paths[i][1:4])
    print(pnum)
    run(ct_fol_path, rt_fol_path, i+1, pnum)



        
