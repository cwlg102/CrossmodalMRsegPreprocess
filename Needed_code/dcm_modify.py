import os
import pydicom
dcm_path = r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\20211014 lbo\ct\9643504_StrctrSets.dcm"
rtst = pydicom.dcmread(dcm_path, force=True)

rtst.StructureSetROISequence[5].ROIName = "Femur_Head_L"
rtst.StructureSetROISequence[4].ROIName = "Femur_Head_R"
print(rtst.StructureSetROISequence[4])
rtst.save_as(r"D:\!Crossmodality_2023\!Pelvic_Unity_Data\20211014 lbo\ct\9643504_StrctrSets.dcm")