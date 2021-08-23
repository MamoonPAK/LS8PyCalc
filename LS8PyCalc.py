# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a LST calculation script file.
"""

import os
from glob import glob
#import numpy as np
import tarfile
import re
import fiona
import rasterio
import rasterio.mask
import numpy as np
import rasterio
import rasterio.plot
from rio_toa import radiance
import tkinter
from tkinter import filedialog
from osgeo import gdal,ogr
import matplotlib.pyplot as plt
import shutil

# import shutil
# import sh

srcDataPath = ""
TempWorkingFolder = ""
fileDate = ""
ExtractedFolder = ""
MaskedFolder = ""
NDVIFolder = ""
NDVIRatioFolder = ""
StackedFolder = ""
AOIFile = ""
DataPath = ""
OutputPath = ""
count = -1


def GUIDeveloper():
    global mystringAOI
    global mystringInput
    global mystringOutput
    global mystringMessage
    global myStringPMLocation
    
    window = tkinter.Tk()
    window.title("LST Calculator")
    window.geometry("680x250")
    # window.resizable(0, 0)
    
    mystringAOI = tkinter.StringVar()       
    
    tkinter.Label(window, text = "AOI File (.shp)").grid(row = 0)  #label
    tkinter.Entry(window, textvariable = mystringAOI,width=70).grid(row = 0, column = 1,padx=5,pady=5,ipady=3)
    tkinter.Button(window, text="Browse", command=getAOIFile).grid(row=0, column=2, sticky=tkinter.W) #button
    
    mystringInput = tkinter.StringVar()       
    
    tkinter.Label(window, text = "Input Folder Path * ").grid(row = 2)  #label
    tkinter.Entry(window, textvariable = mystringInput,width=70).grid(row = 2, column = 1,padx=5,pady=5,ipady=3)
    tkinter.Button(window, text="Browse", command=getInputFolderPath).grid(row=2, column=2, sticky=tkinter.W) #button
        
    mystringOutput = tkinter.StringVar()       
    
    tkinter.Label(window, text = "Output Folder Path *").grid(row = 7)  #label
    tkinter.Entry(window, textvariable = mystringOutput,width=70).grid(row = 7, column = 1,padx=5,pady=5,ipady=3)
    tkinter.Button(window, text="Browse", command=getOutputFolderPath).grid(row=7, column=2, sticky=tkinter.W) #button
    
    myStringPMLocation = tkinter.StringVar()       
    
    tkinter.Label(window, text = "Hand Held Thermometer File (.shp)").grid(row = 8)  #label
    tkinter.Entry(window, textvariable = myStringPMLocation,width=70).grid(row = 8, column = 1,padx=5,pady=5,ipady=3)
    tkinter.Button(window, text="Browse", command=getPMLocationFile).grid(row=8, column=2, sticky=tkinter.W) #button
    
    tkinter.Button(window, text="Process", command=processMethod).grid(row=10, column=1,padx=135, sticky=tkinter.W,columnspan=3) #button    
    
    mystringMessage = tkinter.StringVar() 
    mystringMessage.set("")
    
    tkinter.Label(window, textvariable = mystringMessage,fg = "red",font = "Times").grid(row = 11,column = 0)  #label    
        
    window.mainloop()

def getAOIFile():
    global AOIFile
    input = filedialog.askopenfile(initialdir="/")
    AOIFile = input.name
    mystringAOI.set(input.name)
    
def getInputFolderPath():
    global DataPath
    input = filedialog.askdirectory()    
    DataPath = input
    mystringInput.set(input)    
    
def getOutputFolderPath():
    global OutputPath
    input = filedialog.askdirectory()
    OutputPath = input
    mystringOutput.set(input)
    
def getPMLocationFile():
    global PMLocationFile
    input = filedialog.askopenfile(initialdir="/")
    PMLocationFile = input.name
    myStringPMLocation.set(input.name)
    
def initilizeGlobalVariables():
    
    global ExtractedFolder
    global MaskedFolder
    global LSTFolder
    global LSEFolder
    global NDVIFolder
        
    ExtractedFolder = "/1-Extracted"
    MaskedFolder = "/2-Masked"
    LSTFolder = "/3-LST"
    LSEFolder = "/4-LSE"
    NDVIFolder = "/5-NDVI"
    
def createDirectories():
    
    # Make new folder Extracted    
    if not os.path.exists(OutputPath + ExtractedFolder):         
        os. mkdir(OutputPath + ExtractedFolder)
        
    # Make new folder Masked
    if not os.path.exists(OutputPath + MaskedFolder):         
        os. mkdir(OutputPath + MaskedFolder)
    
    # Make new folder LST
    if not os.path.exists(OutputPath + LSTFolder):         
        os. mkdir(OutputPath + LSTFolder)
        
    # Make new folder LSE
    if not os.path.exists(OutputPath + LSEFolder):         
        os. mkdir(OutputPath + LSEFolder)
    
    # Make new folder NDVI
    if not os.path.exists(OutputPath + NDVIFolder):         
        os. mkdir(OutputPath + NDVIFolder)
        

def extractAndMasking():
    Count = 0
    global fileDate
    global fileOriginalName
    
    # DataPath = 'D:/Mamoon(Personal)/PhD/1-WorkingFolder/9-Practical/0-LSTCalculator_WithMasking/Data'
    
    # bands_grep = re.compile(".*_(B4|B5)\.TIF")
    bands_grep = re.compile(".*_(B4|B5|B10|MTL)\.TIF|.txt")

    # List files with the tar.gz extension in the data/ directory
    # tar_list = glob(DataPath + '*.tar.gz')
    tar_list = glob(DataPath + "/" + '*.tar')

    for tarFile in tar_list:
        # tarFile = tar_list[0]
    
        fileDate = tarFile.split("\\")[1].split(".")[0].split("_")[3]
        fileOriginalName = tarFile.split("\\")[1].split(".")[0]
        # Open a connection with the first archive (in fact there's only one)
        tar = tarfile.open(tarFile)
    
    
        # Retrieve the file names (Bands Name) from the archive
        file_list = tar.getnames()
        
        # Filter and keep only useful files (bands 4,5,10 and MTL to compute LST)
        bands = filter(lambda x: bands_grep.search(x), file_list)
        for item in bands:
            tar.extract(item, path=OutputPath + ExtractedFolder)
            # output = Basic_Path + ExtractedFolder 
            if  (AOIFile != ""):
                if item.split(".")[1] == "TIF":
                    imageMasking(item,ExtractedFolder)        
            
            
        # Calculate LST
        LSTCal(fileDate,fileOriginalName)
        
        # Close the connection
        tar.close()  
        
        deleteFolderContents(OutputPath+ExtractedFolder)
        deleteFolderContents(OutputPath+MaskedFolder)
        Count = Count + 1
        # time.sleep(15)
        # mystringMessage.set("<<<<<-----  " + str(Count) + "  LST Generated ----->>>>>>")
        # print("<<<<<-----  " + str(Count) + "  LST Generated ----->>>>>>" )       
    
    # mystringMessage.set("Completed")
    mystringMessage.set(str(Count) + "  LST Generated")

def imageMasking(IBandName,ExtractedFolderName):
            
    with fiona.open(AOIFile, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
        
    # print(Basic_Path+ExtractedFolderName+"/"+IBandName)
    pp = OutputPath + ExtractedFolderName + "/" + IBandName
    # with rasterio.open(Basic_Path + ExtractedFolderName + "/" + IBandName) as src:
    with rasterio.open(pp) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta
        
    out_meta.update({"driver": "GTiff",
                      "height": out_image.shape[1],
                      "width": out_image.shape[2],
                      "transform": out_transform})
    
    maskedFileName = IBandName.split("_")[3] + "_" + IBandName.split("_")[-1]
    
    with rasterio.open(OutputPath + MaskedFolder + "/" + maskedFileName, "w", **out_meta) as dest:
        dest.write(out_image)  

def movesDownloadedFiles(strDestFolderName,strSrcPath,strFileName):
    # src_folder = r"C:\Users\Lenovo\.spyder-py3\Downloaded_Files\\"    
    # dst_folder = r"C:\Users\Lenovo\.spyder-py3\CLSM\\"
    
    src_folder = strSrcPath + "\Downloaded_Files\\" + strFileName
    dst_folder = strSrcPath + "/" + strDestFolderName + "\\"
    
            
    if not os.path.exists(strSrcPath + "/" + strDestFolderName):         
        os. mkdir(strSrcPath + "/" + strDestFolderName)
        
    # move file whose name starts with string 'emp'
    # pattern = src_folder + "\GLDAS_CLSM10*"
    pattern = src_folder
    for file in glob.iglob(pattern, recursive=True):
        # extract file name form file path
        file_name = os.path.basename(file)
        shutil.move(file, dst_folder + file_name)
        print('Moved:', file)  
        
def LSTCal(fileDate,OriginalName):
    strFileName = ''
    #Get RADIANCE_MULT_BAND_10 and RADIANCE_ADD_BAND_10 value from text file of Band 10
    filePath = OutputPath + ExtractedFolder + "/" + OriginalName + "_MTL.txt"
    # print ("File Name:   ",filePath)
    f = open(filePath, "r")
    for line in f:
        if 'RADIANCE_MULT_BAND_10' == str(line.split("=")[0].strip()):
    #        print ("RADIANCE_MULT_BAND_10: ",line.split("=")[1].strip())            
    #        print ("RADIANCE_MULT_BAND_10: ",float(line.split("=")[1].strip()))            
            Radiance_Mult_Band_10 = float(line.split("=")[1].strip())
            
        if 'RADIANCE_ADD_BAND_10' == str(line.split("=")[0].strip()):
    #        print ("RADIANCE_ADD_BAND_10: ",line.split("=")[1].strip())            
            Radiance_Add_Band_10 = float(line.split("=")[1].strip())
            
        if 'K1_CONSTANT_BAND_10' == str(line.split("=")[0].strip()):
    #        print ("K1_CONSTANT_BAND_10: ",line.split("=")[1].strip())            
            K1_Const_Band_10 = float(line.split("=")[1].strip())
            
        if 'K2_CONSTANT_BAND_10' == str(line.split("=")[0].strip()):
    #        print ("K2_CONSTANT_BAND_10: ",line.split("=")[1].strip())            
            K2_Const_Band_10 = float(line.split("=")[1].strip())
            
    
    if  (AOIFile != ""):
        operationalPath = OutputPath + MaskedFolder
        strFileName = fileDate
    else:
        operationalPath = OutputPath + ExtractedFolder
        strFileName = OriginalName
        
    #Read Band 10    
    dsBand10 = rasterio.open(operationalPath + '/' + strFileName + "_B10.TIF")
    Band10 = dsBand10.read(1).astype('float64')
    
    # print("Display Band10   ",dsBand10)
    
    #set Output image dimensions
    IHeight = dsBand10.height
    IWidth = dsBand10.width
    Icrs =  dsBand10.crs        
    
    #1.- Calculation of TOA (Top of Atmospheric) spectral radiance.
    #    TOA (L) = M L  * Q cal  + A L
    #    M L  = RADIANCE_MULT_BAND_x, where x is the band number
    #    Q cal  = corresponds to band 10.
    #    A L  = RADIANCE_ADD_BAND_x, where x is the band number
    
    #toa = radiance.radiance(tif, ML, AL, src_nodata=0)
    L = radiance.radiance(Band10, Radiance_Mult_Band_10, Radiance_Add_Band_10, src_nodata=0)
    
    # L = np.where(tmpL==-272.339017,0, tmpL)
    # print("Display TOA     ",L)
    
    with rasterio.open(
        operationalPath + '/' + 'TOA.tif', 'w',
        driver='GTiff', width=IWidth, height=IHeight, count=1,
        # dtype=L.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstTOA:
        dtype=L.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstTOA:
        dstTOA.write(L, 1)
    # print (dstTOA.nodatavals)    
        
    # print ("Step 1 Completed")

    ##2.- TOA to Brightness Temperature conversion
    #    BT = (K 2  / ln ((K 1  / L) + 1)) − 273.15
    #    K 1  = K1_CONSTANT_BAND_x, where x is the thermal band number
    #    K 2  = K2_CONSTANT_BAND_x, where x is the thermal band number
    #    L = TOA        
    
    # tmpVal = 
    Bt = (K2_Const_Band_10 / np.log ((K1_Const_Band_10 / L) + 1)) - 273.15
    
    # print ("BT.............." , Bt)
    
    with rasterio.open(
        operationalPath + '/'+'BT.tif', 'w',
        driver='GTiff', width=IWidth, height=IHeight, count=1,
        # dtype=Bt.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstBT:
        dtype=Bt.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=None) as dstBT:
        dstBT.write(Bt, 1)
    # print ("Step 2 COmpleted")    
    ##3.- Calculate the NDVI
    #    #NDVI = (Band 5 – Band 4) / (Band 5 + Band 4)
    #    #NDVI = Float(Band 5 – Band 4) / Float(Band 5 + Band 4)
    #
    #Open Band4 as red band
    dsBand4 = rasterio.open(operationalPath + '/'+ strFileName +"_B4.TIF")
    red = dsBand4.read(1).astype('float64')
    # print("Display Band4   ",red)
    
    #Open Band5 as nir band
    dsBand5 = rasterio.open(operationalPath + '/'+ strFileName +"_B5.TIF")
    nir = dsBand5.read(1).astype('float64')
    # print("Display Band5   ",nir)
    
    #Calculate NDVI
    ndvi=np.where((nir+red)==0.,0, (nir-red)/(nir+red))
    # print("NDVI", ndvi)
    
    with rasterio.open(
        OutputPath + NDVIFolder + "/" + strFileName +'_NDVI.tif', 'w',
        driver='GTiff', width=IWidth, height=IHeight, count=1,
        # dtype=ndvi.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstNDVI:
        dtype=ndvi.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstNDVI:
        dstNDVI.write(ndvi, 1)
    
    # plt.imsave(operationalPath + '/'+'NDVIColor.png', ndvi, cmap=plt.cm.summer)        
    
    # print ("Step 3 COmpleted")    
    #4.- Calculate the proportion of vegetation P v
    #P v  = Square ((NDVI – NDVI min ) / (NDVI max  – NDVI min ))
    #P v  = Square(“NDVI_Image” – 0.216901) / (0.632267 – 0.216901)        
    
    p = np.square(ndvi-0.216901)/(0.632267-0.216901)
    # print("proportion of vegetation P",p)
    
    with rasterio.open(
        operationalPath + '/'+'P.tif', 'w',
        driver='GTiff', width=IWidth, height=IHeight, count=1,
        # dtype=ndvi.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstP:
        dtype=ndvi.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstP:
        dstP.write(ndvi, 1)
    
    # print ("Step 4 COmpleted")    
    #5.- Calculate Emissivity ε
    #   ε = 0.004 * P v  + 0.986
    Ee = 0.004*p+0.986
    # print("Emissivity ε",Ee)
    
    with rasterio.open(
        OutputPath + LSEFolder + "/" + strFileName +'_Emissivity.tif', 'w',    
        # operationalPath + '/4-LSE/'+'Emissivity.tif', 'w',
        driver='GTiff', width=IWidth, height=IHeight, count=1,
        # dtype=ndvi.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstE:
        dtype=ndvi.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstE:
        dstE.write(ndvi, 1)
        
    # print ("Step 5 COmpleted")
    #6.- Calculate the Land Surface Temperature
    #LST = (BT / (1 + (0.00115 * BT / 1.4388) * Ln(ε)))
    #lst = (Bt/(1+(0.00115*Bt/1.4388)*np.log(ee)))
    tmpLST = (Bt/(1+(0.00115*Bt/1.4388)*np.log(Ee)))
    lst = np.where(tmpLST==-272.339,0, tmpLST)
    
    #Save LST image
    with rasterio.open(
        OutputPath+LSTFolder+'/'+ fileDate +'_'+'LST.tif', 'w',
        driver='GTiff', width=IWidth, height=IHeight, count=1,
        # dtype=lst.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstLST:
        dtype=lst.dtype,crs=Icrs,resolution=dsBand10.res,transform=dsBand10.transform,nodatavals=0) as dstLST:
        dstLST.write(lst, 1)
def training_points(raster, shp):
    lsx=[]
    lsy=[]
    
    # Open Raster Image
    src_ds=gdal.Open(raster)
    gt=src_ds.GetGeoTransform()
    rb=src_ds.GetRasterBand(1)
    
    # Open shape file
    ds=ogr.Open(shp)
    lyr=ds.GetLayer()    
    count = 0
    for feat in lyr:
        geom = feat.GetGeometryRef()
        mx,my=geom.GetX(), geom.GetY()  #coord in map units
        #Convert from map to pixel coordinates.
        #Only works for geotransforms with no rotation.
        px = int((mx - gt[0]) / gt[1]) #x pixel
        py = int((my - gt[3]) / gt[5]) #y pixel
        
        intval=rb.ReadAsArray(px,py,1,1)
        
        value = float(intval[0]) #### this is the value of the pixel, forcing it to a float
        
        lsx.append(str(count))
        lsy.append(value)
        count +=1
    # return lsx, lsy, lsz
    return lsx, lsy
    
def getLayerFeatureCount(strShapeFile):
    ds=ogr.Open(strShapeFile)
    lyr=ds.GetLayer()
    
    return lyr.GetFeatureCount()

def getReflectance():
    
    shp = PMLocationFile        
    source = OutputPath + LSTFolder    
    
    fig,ax = plt.subplots(1,1)
    ax.set_xlabel('Locations') ; ax.set_ylabel('Reflectance')
    # ax.set_title('pcolormesh')
        
    ax.set_xlim(0,getLayerFeatureCount(shp)) ; 
    legendContents=[]
    ## loop through my source directory
    for root, dirs, filenames in os.walk(source):
        for file in filenames:
            ## looking for tiles ending with TIF
            if file.endswith('tif'):
                raster = os.path.join(source + '/', file)
                
                lsx, lsy = training_points(raster, shp) 
                
                legendVal = str(raster.split('/')[len(raster.split('/'))-1]).split('_')[0]
                
                legendContents.append(legendVal)                
                ax.plot(lsx, lsy)
                fig.canvas.draw()                           
    
    plt.legend(legendContents)
    # fig = plt.figure(figsize=(18,12))
    plt.savefig(source + '/myfig.png')
    ## show my plot
    plt.show()

def deleteFolderContents(dirpath):    
    
    for i in os.listdir(dirpath):
        os.remove(dirpath+'/'+i)

def processMethod():
    
    initilizeGlobalVariables()
    createDirectories()
    extractAndMasking()
    getReflectance()
    
def main():    
    GUIDeveloper()
main()

