#encoding=utf-8
import wave
#from pymedia import *
import pymedia.muxer as muxer
import pymedia.audio.acodec as acodec
#import pymedia.audio.sound as sound

import time
import os
import gc
import math
import sys
import random as rnd

import pylab as pl
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import *
from array import array
from scipy import fftpack
import pywt

from ctypes import *
import ctypes

# global control
flag = ""
outputWave = False
showDetails = False

TOTAL_TIME = 0

dll = cdll.LoadLibrary("AudioAnalysis.dll")

#PROJECT_DIR = r"H:\音乐\项目目录"
#TOTAL_MUSIC = 478 # 8 songs for musical instrument testing
PROJECT_DIR = r"H:\音乐\项目目录\音乐基准数据集\606201"
PROJECT_DIR = unicode(PROJECT_DIR, "utf8")
TOTAL_MUSIC = 1892 # 8 songs for musical instrument testing
#PROJECT_DIR = r"H:\音乐\项目目录\音乐基准数据集\606201\instruments"
#TOTAL_MUSIC = 149 # 8 songs for musical instrument testing
PROJECT_TEST_DIR = PROJECT_DIR + r"\test"
PROJECT_DOC_DIR = PROJECT_DIR + r"\doc"
PROJECT_INFO_DIR = PROJECT_DIR + r"\info"
PROJECT_PCASPACE_DIR = r"\pca_space.txt"
PROJECT_RAWFEAT_DIR = r"\data_freqIt5_db4_melOnce_SVDfeatureWithMel_2596.txt"
PROJECT_COSDIS_DIR = r"\cosDis.txt"
EXCEL_DATA_NAME = r"\data.xls"
INFO_TXT_NAME = r"\info.txt"
AudioFileCount = 0
currentAudioFileCount = 0
dataToWrite = []
MFCC_TOP_FEAT = 320
WLT_TOP_FEAT = 200
BEAT_TOP_FEAT = 4
TOP_FEAT = 100
PROJECT_PCADATA_DIR = u"lowDData_" + str(MFCC_TOP_FEAT) + "_" + str(WLT_TOP_FEAT)\
                      + "_" + str(BEAT_TOP_FEAT) + "_pca.txt"
#PROJECT_PCADATA_DIR = u"lowDData_" + str(TOP_FEAT) + "_pca.txt"

#data = np.zeros(0, dtype = np.short)

nchannels = 2
sampwidth = 2
framerate = 44100
nframes = 10000
comptype = "NONE"
compname = "compressed"

params = [nchannels, sampwidth, framerate, nframes, comptype, compname]

# calculate the Mel scale
max_frequency = framerate / 2
max_melFrequency = 2595 * (math.log10(1 + float(max_frequency) / 700))
#MUSIC_CHUNK_NUM = 4
#TIME_SLIP = 0.5
CHANNELS = 2
beats_dimension = 8 # define the dimension of Beats
melWindow = 690
MFCCs_dimension = 25 # define the dimension of MFCCs
dwtIterations_time = 13 # not used
dwtIterations_freq = 5 # this value should bigger than 0
CLength = 827
dwtParts = dwtIterations_time + 1 # not used
dwtModel = "db4"
chunkNum = 1

getPCATarget = False
targetEigVals = np.zeros(0)
targetEigVects_L = np.zeros(0)
targetEigVects_R = np.zeros(0)

melParaArr = np.arange(0, dtype = ctypes.c_int) # define with 0

# melParaArr is used for weighting wlt parts in mel scale
def initMelParaArr(dwtParts):
    global melParaArr
    
    melParaArr = np.arange(dwtParts, dtype = ctypes.c_int)
    for i in xrange(dwtParts):
        melParaArr[i] = (max_melFrequency / dwtParts) * (i + 1)
#    print melParaArr
    for i in xrange(dwtParts):
        melParaArr[i] = (math.pow(10, float(melParaArr[i]) / 2595) - 1) * 700
    maxDiff = melParaArr[-1] - melParaArr[-2]
#    print melParaArr
#    print maxDiff
    for i in xrange(dwtParts):
        if i == 0:
            temp = melParaArr[i]
            melParaArr[i] = maxDiff / melParaArr[i]
        else:
            res = maxDiff / (melParaArr[i] - temp)
            temp = melParaArr[i]
            melParaArr[i] = res
#    print melParaArr

wltAnalysisLength = 60 * framerate * 2 # in seconds
MFCCPara_dimension = chunkNum * MFCCs_dimension * pow(2, dwtIterations_freq) * \
                     CHANNELS  * 2
#wltPara_dimension = CLength * 12
wltPara_dimension = 2048
frequencyPara_dimension = MFCCPara_dimension + wltPara_dimension
print "MFCCPara_dimension: ", frequencyPara_dimension - wltPara_dimension
print "wltPara_dimension: ", wltPara_dimension
print "beats_dimension: ", beats_dimension


melScale = np.arange(melWindow, dtype = ctypes.c_int)

for i in xrange(melWindow):
    melScale[i] = (max_melFrequency / melWindow) * (i + 1)
#print melScale
    
# find the corresponding frequency scale of Mel scale
for i in xrange(melWindow):
    melScale[i] = (math.pow(10, float(melScale[i]) / 2595) - 1) * 700
#print melScale


melScale = np.arange(melWindow, dtype = ctypes.c_int)
for i in xrange(melWindow):
    melScale[i] = (max_frequency / melWindow) * (i + 1)
#print melScale


MEL_INTENSITY = np.arange(melWindow, dtype = ctypes.c_float)
MFCCs = np.arange(MFCCs_dimension, dtype = ctypes.c_float)
FIRST_DIFF = np.zeros(MFCCs_dimension)
#MEL_INTENSITY_VAR = np.zeros(melWindow)

initMelParaArr(pow(2, dwtIterations_freq))
#initMelParaArr(2596)

musicAttributes = np.arange(20, dtype = ctypes.c_float)
musicAttributes[0] = 0.0# AVE_BIG_WINDOW_LENGTH
musicAttributes[1] = 0.0# AVERAGE_INTENSITY_L
musicAttributes[2] = 0.0# AVERAGE_INTENSITY_R
musicAttributes[3] = 0.0# AVERAGE_INTENSITY
musicAttributes[4] = -1000000.0# MAX_VALUE_L
musicAttributes[5] = -1000000.0# MAX_VALUE_R
musicAttributes[6] = -1000000.0# MAX_VALUE
musicAttributes[7] = 1000000.0# MIN_VALUE_L
musicAttributes[8] = 1000000.0# MIN_VALUE_R
musicAttributes[9] = 1000000.0# MIN_VALUE
musicAttributes[10] = 0.0# MAX_L / AVERAGE_L
musicAttributes[11] = 0.0# MAX_R / AVERAGE_R
musicAttributes[12] = 0.0# MAX / AVERAGE
musicAttributes[13] = 0.0# NONZERO_EXCEPTSILENCE_LENGTH
musicAttributes[14] = 0.0# AVE_BEAT_INTENSITY_L
musicAttributes[15] = 0.0# AVE_BEAT_INTENSITY_H
musicAttributes[16] = 0.0# BEATS_L
musicAttributes[17] = 0.0# BEATS_H
musicAttributes[18] = 0.0# LOW_PASS_AVE_REDUCE_RATIO
musicAttributes[19] = 0.0# HIGH_PASS_AVE_REDUCE_RATIO

# AVE_BIG_WINDOW_LENGTH = musicAttributes[0]
# AVERAGE_INTENSITY_L = musicAttributes[1]
# AVERAGE_INTENSITY_R = musicAttributes[2]
# AVERAGE_INTENSITY = musicAttributes[3]
# MAX_VALUE_L = musicAttributes[4]
# MAX_VALUE_R = musicAttributes[5]
# MAX_VALUE = musicAttributes[6]
# MIN_VALUE_L = musicAttributes[7]
# MIN_VALUE_R = musicAttributes[8]
# MIN_VALUE = musicAttributes[9]
# MAX_L / AVERAGE_L = musicAttributes[10]
# MAX_R / AVERAGE_R = musicAttributes[11]
# MAX / AVERAGE = musicAttributes[12]
# NONZERO_EXCEPTSILENCE_LENGTH = musicAttributes[13]// NOT USED
# AVE_BEAT_INTENSITY_L = musicAttributes[14]
# AVE_BEAT_INTENSITY_H = musicAttributes[15]
# BEATS_L = musicAttributes[16]
# BEATS_H = musicAttributes[17]
# LOW_PASS_AVE_REDUCE_RATIO = musicAttributes[18]
# HIGH_PASS_AVE_REDUCE_RATIO = musicAttributes[19]

def errorCheck(error, num = -1):
    if error is 0:
        if (showDetails):
            if num is 0:
                print "-------------------------------------------------------------------"
                print "make silence finish"
                print "-------------------------------------------------------------------"
            if num is 1:
                print "-------------------------------------------------------------------"
                print "difference detection finish"
                print "-------------------------------------------------------------------"
            if num is 2:
                print "-------------------------------------------------------------------"
                print "low pass filter finish"
                print "-------------------------------------------------------------------"
            if num is 3:
                print "-------------------------------------------------------------------"
                print "high pass filter finish"
                print "-------------------------------------------------------------------"
            if num is 4:
                print "-------------------------------------------------------------------"
                print "intensity conversion finish"
                print "-------------------------------------------------------------------"
            if num is 5:
                print "-------------------------------------------------------------------"
                print "get attributes finish"
                print "-------------------------------------------------------------------"
            if num is 6:
                print "-------------------------------------------------------------------"
                print "beatFilter finish"
                print "-------------------------------------------------------------------"
            if num is 7:
                print "-------------------------------------------------------------------"
                print "getBeats finish"
                print "-------------------------------------------------------------------"
            if num is 8:
                print "-------------------------------------------------------------------"
                print "countBeats finish"
                print "-------------------------------------------------------------------"
            if num is 9:
                print "-------------------------------------------------------------------"
                print "fft finish"
                print "-------------------------------------------------------------------"
            if num is 10:
                print "-------------------------------------------------------------------"
                print "getMFCCs finish"
                print "-------------------------------------------------------------------"
            if num is 11:
                print "-------------------------------------------------------------------"
                print "getAMPlitude finish"
                print "-------------------------------------------------------------------"
            if num is 12:
                print "-------------------------------------------------------------------"
                print "getMelIntensity finish"
                print "-------------------------------------------------------------------"
        
    if error is 1:
        print "memery alloction ERROR"
    if error is 2:
        print "NULL input ERROR"
    if error is 3:
        print "open file ERROR"

def getPath(fileName):
    global PROJECT_DIR
    # gb18030 to unicode
    fileName = fileName.decode("gb18030")
    # unicode to utf-8
    fileName = fileName.encode("UTF-8")
    PROJECT_DIR = PROJECT_DIR.encode("UTF-8")
    path = PROJECT_DIR + fileName
    path = unicode(path, "utf8")
    return path, fileName

def getData(fileName):
    global PROJECT_DIR
    global TOTAL_TIME
#    global data
    
    unicodeFlag = isinstance(fileName, unicode)
    # the path is typed in command line
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        # the path is return in unicode
        path = fileName
        fileName = fileName.strip().split("\\")[-1]

    #data = np.zeros(0, dtype = np.short)
    if ((fileName.strip().split(".")[-1]) == "wav") or \
       ((fileName.strip().split(".")[-1]) == "WAV"):
        if (showDetails):
            print "get data from .wav"
        T1 = time.time()
        data = getWAVData(path)
        T2 = time.time()
        TOTAL_TIME = TOTAL_TIME + (T2- T1)
        print "get data in .wav takes %fs" % (T2- T1)
        print fileName
    if ((fileName.strip().split(".")[-1]) == "mp3") or \
       ((fileName.strip().split(".")[-1]) == "MP3"):
        if (showDetails):
            print "get data from .mp3"
        T1 = time.time()
        data = getMP3Data(path)
        T2 = time.time()
        TOTAL_TIME = TOTAL_TIME + (T2- T1)
        print "get data from .mp3 takes %fs" % (T2- T1)
        print fileName

    del unicodeFlag
    if path:
        del path
    del fileName
    
    return data

def getWAVData(path):
#    global data
    
    f = wave.open(path, "rb")
    params = f.getparams()
    nframes = params[3]
    data = f.readframes(nframes)
    data = np.fromstring(data, dtype = np.short)
    f.close()

    del f
    del params
    del nframes
    del path
    
    return data

def getMP3Data(path):
#    global data
    
    fr = open(path, "rb")

    #data = np.zeros(0, dtype = np.short)
    
    rawData = fr.read(3) # the first 3 bytes is "ID3", if it has
#    print rawData
    if rawData == "ID3":
#        print "get ID3"
        del rawData
        rawData = fr.read(3) # contains the revision and flag
        del rawData
        rawData = fr.read(4) # contains the size of the ID3V2
                             # each byte has only 7 bits is available x0000000
        ID3Size = np.zeros(5, dtype = np.int)
        for i in xrange(4):
            ID3Size[i] =  np.fromstring(rawData[i] + chr(0) * 3, \
                                        dtype = np.int)
            ID3Size[i] = ID3Size[i] & 0x7F # calculate each byte
        ID3Size[0] = ID3Size[0] * 0x200000
        ID3Size[1] = ID3Size[1] * 0x4000
        ID3Size[2] = ID3Size[2] * 0x80
        ID3Size[3] = ID3Size[3]
        ID3Size[4] = ID3Size.sum()
        fr.seek(ID3Size[4]) # skip the ID3V2 field
#        print "skip ID3V2 in %dbytes" % ID3Size[4]
        del ID3Size
    else:
#        print "seek to 0"
        fr.seek(0)

    dm = None
    dec = None
    index = 0
    readFrame = 1
    startDec = False
    readLength = 0
    while rawData:
        if startDec == False:
            #readFrame = 50
            del readLength
            readLength = 4096 * readFrame
            del rawData
            rawData = fr.read(readLength)
            index += readLength
        if startDec:
            del rawData
            rawData = fr.read()
        if dm == None:
            dm = muxer.Demuxer("mp3")
        frames = dm.parse(rawData)
        numOfFrames = len(frames)
        if frames:
            if dec == None:
                del dec
                dec = acodec.Decoder(dm.streams[0])
            if startDec:
                frameData = []
                for i in xrange(numOfFrames):
                    frameData.append(frames[i][1])
                del frames
                frameDataCon = "".join(frameData)
                del frameData
                r = dec.decode(frameDataCon)
                del frameDataCon
                data = np.fromstring(r.data, dtype = np.short)
                del r

            if startDec == False:
                if ((np.fromstring(frames[0][1][0] + chr(0) * 3, dtype = np.int) \
                     == 255) and \
                     np.fromstring(frames[0][1][1] + chr(0) * 3, dtype = np.int) \
                     == 251):
                    startDec = True
                    index = index - readLength
                    fr.seek(index)
                del frames
                
    del rawData                
    gc.collect()
    if startDec == False:
#        print "cannot find header"
        fr.seek(0)
        rawData = fr.read()
        del dm
        dm = muxer.Demuxer("mp3")
        frames = dm.parse(rawData)
        del rawData
        numOfFrames = len(frames)
        frameData = []
        for i in xrange(numOfFrames):
            frameData.append(frames[i][1])
        del frames
        frameDataCon = "".join(frameData)
        del frameData
        r = dec.decode(frameDataCon)
        del frameDataCon
        data = np.fromstring(r.data, dtype = np.short)
        del r
    fr.close()
    
    del dm
    del dec
    del fr
    
    del path
    del index
    del readFrame
    del startDec
    del readLength

    return data

#delta is the factor of average of the data
def makeSilence(data, fileName, delta_L = 1, delta_R = 1):

#        Type code   C Type             Minimum size in bytes
#        'c'         character         　　　　　　　  1
#        'b'         signed integer     　　　　　    1
#        'B'         unsigned integer   　　　　　    1
#        'u'         Unicode character  　　　　　    2
#        'h'         signed integer     　　　　　　　 2
#        'H'         unsigned integer   　　　　　　   2
#        'i'         signed integer     　　　　　　　 2
#        'I'         unsigned integer  　　　　　　    2
#        'l'         signed integer     　　　　　　　 4
#        'L'         unsigned integer 　　　　　　     4
#        'f'         floating point    　　　　　　　　 4
#        'd'         floating point    　　　　　　　　 8

    global musicAttributes
    global PROJECT_DIR
    global PROJECT_TEST_DIR
    
    unicodeFlag = isinstance(fileName, unicode)
    #人为输入文件名
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        #读取程序返回的unicode编码的路径
        path = fileName
        fileName = fileName.strip().split("\\")[-1]

    tempDelta_L = delta_L
    tempDelta_R = delta_R
    getMusicAttributes(data)
    delta_L = delta_L * int(musicAttributes[1])# AVERAGE_INTENSITY_L
    delta_R = delta_R * int(musicAttributes[2])# AVERAGE_INTENSITY_R
    
    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(makeSilence) music loaded!"
        print "-------------------------------------------------------------------"

        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]

    if (outputWave):
        fileName = fileName[0 : len(fileName) - \
                        len(fileName.strip().split(".")[-1]) - 1] \
                        + "_makeSilence(" + \
                        str(tempDelta_L) + "-" + str(tempDelta_R) + ").wav"
        path = PROJECT_TEST_DIR + fileName
        if os.path.exists(path):
            os.remove(path)
            print "same file detected, delete"
        fw = wave.open(path, "w")
        fw.setparams(params)
        
    if (showDetails):    
        print "-------------------------------------------------------------------"
        print "start make silence"
        print "-------------------------------------------------------------------"

    dataLength = len(data)
    error = -1
    T1 = time.time()
    error = dll.makeSilence(data.ctypes.data_as(POINTER(ctypes.c_int16)), \
                            ctypes.c_int(dataLength), \
                            ctypes.c_float(delta_L), \
                            ctypes.c_float(delta_R), \
                            musicAttributes.ctypes.data_as(\
                                POINTER(ctypes.c_float)), \
                            ctypes.c_bool(True))
    T2 = time.time()

    if (outputWave):
        fw.writeframes(data[0 : ].tostring())
        fw.close()
    
    
    errorCheck(error, 0)
    if (showDetails):
        print "%fs" % float(T2 - T1)
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]
        print "\n"
    del data
    return path

# diff = data[i] - data[i - 1]
def differenceDetection(data, fileName, ratio = 1.0):
    global musicAttributes
    global PROJECT_DIR
    global PROJECT_TEST_DIR
    
    unicodeFlag = isinstance(fileName, unicode)
    #人为输入文件名
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        #读取程序返回的unicode编码的路径
        path = fileName
        fileName = fileName.strip().split("\\")[-1]

    getMusicAttributes(data)
    
    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(differenceDetection) music loaded!"
        print "-------------------------------------------------------------------"

        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]

    if (outputWave):
        fileName = fileName[0 : len(fileName) - \
                        len(fileName.strip().split(".")[-1]) - 1] \
                        + "_diffDetec(" + \
                       str(ratio) + ").wav"
        path = PROJECT_TEST_DIR + fileName
        if os.path.exists(path):
            os.remove(path)
            print "same file detected, delete"
        fw = wave.open(path, "w")
        fw.setparams(params)

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start difference detection"
        print "-------------------------------------------------------------------"

    dataLength = len(data)
    error = -1
    T1 = time.time()
    error = dll.differenceDetection(data.ctypes.data_as(POINTER(ctypes.c_int16)), \
                            ctypes.c_int(dataLength), \
                            ctypes.c_float(ratio), \
                            musicAttributes.ctypes.data_as(\
                                POINTER(ctypes.c_float)), \
                            ctypes.c_bool(True))
    T2 = time.time()

    if (outputWave):
        fw.writeframes(data[0 : ].tostring())
        fw.close()
    
    errorCheck(error, 1)
    if (showDetails):
        print "%fs" % float(T2 - T1)
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]
        print "\n"
    del data
    return path

def lowPassFilter(data, fileName, alpha = 0.005, sigma = 1):
    global musicAttributes
    global PROJECT_DIR
    global PROJECT_TEST_DIR

    unicodeFlag = isinstance(fileName, unicode)
    #人为输入文件名
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        #读取程序返回的unicode编码的路径
        path = fileName
        fileName = fileName.strip().split("\\")[-1]

    getMusicAttributes(data)
        
    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(lowPassFilter) music loaded!"
        print "-------------------------------------------------------------------"

        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]

    if (outputWave):
        fileName = fileName[0 : len(fileName) - \
                        len(fileName.strip().split(".")[-1]) - 1] \
                        + "_lowPass(" + \
                       str(alpha) + "-" + str(sigma) + ").wav"
        path = PROJECT_TEST_DIR + fileName
        if os.path.exists(path):
            os.remove(path)
            print "same file detected, delete"
        fw = wave.open(path, "w")
        fw.setparams(params)
        
    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start lowPassFilter"
        print "-------------------------------------------------------------------"

    dataLength = len(data)
    error = -1
    T1 = time.time()
    error = dll.lowPassFilter(data.ctypes.data_as(POINTER(ctypes.c_int16)), \
                            ctypes.c_int(dataLength), \
                            ctypes.c_float(alpha), \
                            ctypes.c_float(sigma), \
                            musicAttributes.ctypes.data_as(\
                                POINTER(ctypes.c_float)), \
                            ctypes.c_bool(True))
    T2 = time.time()

    if (outputWave):
        fw.writeframes(data[0 : ].tostring())
        fw.close()
    
    errorCheck(error, 2)
    if (showDetails):
        print "%fs" % float(T2 - T1)
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]
        print "\n"
    del data
    return path

def highPassFilter(data, fileName, alpha = 0.005):
    global musicAttributes
    global PROJECT_DIR
    global PROJECT_TEST_DIR

    unicodeFlag = isinstance(fileName, unicode)
    #人为输入文件名
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        #读取程序返回的unicode编码的路径
        path = fileName
        fileName = fileName.strip().split("\\")[-1]

    getMusicAttributes(data)
    
    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(highPassFilter) music loaded!"
        print "-------------------------------------------------------------------"

        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]

    if (outputWave):
        fileName = fileName[0 : len(fileName) - \
                        len(fileName.strip().split(".")[-1]) - 1] \
                        + "_highPass(" + \
                        str(alpha) + ").wav"
        path = PROJECT_TEST_DIR + fileName
        if os.path.exists(path):
            os.remove(path)
            print "same file detected, delete"
        fw = wave.open(path, "w")
        fw.setparams(params)

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start highPassFilter"
        print "-------------------------------------------------------------------"

    dataLength = len(data)
    error = -1
    T1 = time.time()
    error = dll.highPassFilter(data.ctypes.data_as(POINTER(ctypes.c_int16)), \
                               ctypes.c_int(dataLength), \
                               ctypes.c_float(alpha), \
                               musicAttributes.ctypes.data_as(\
                                   POINTER(ctypes.c_float)), \
                               ctypes.c_bool(True))
    T2 = time.time()

    if (outputWave):
        fw.writeframes(data[0 : ].tostring())
        fw.close()
    
    errorCheck(error, 3)
    if (showDetails):
        print "%fs" % float(T2 - T1)
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]
        print "\n"
    del data
    return path

def convertToIntensity(data, fileName):
    global musicAttributes
    global PROJECT_DIR
    global PROJECT_TEST_DIR
    
    unicodeFlag = isinstance(fileName, unicode)
    #人为输入文件名
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        #读取程序返回的unicode编码的路径
        path = fileName
        fileName = fileName.strip().split("\\")[-1]

    getMusicAttributes(fileName)

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(convertToIntensity) music loaded!"
        print "-------------------------------------------------------------------"
    
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]

    if (outputWave):
        fileName = fileName[0 : len(fileName) - \
                        len(fileName.strip().split(".")[-1]) - 1] \
                        + "_intensityConvert.wav"
        path = PROJECT_TEST_DIR + fileName
        if os.path.exists(path):
            os.remove(path)
            print "same file detected, delete"
        fw = wave.open(path, "w")
        fw.setparams(params)
        
    valuesInEachFrame = len(data) / nframes

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start intensity conversion"
        print "-------------------------------------------------------------------"

    dataLength = len(data)
    error = -1
    T1 = time.time()
    error = dll.convertToIntensity(data.ctypes.data_as(POINTER(ctypes.c_int16)), \
                            ctypes.c_int(dataLength), \
                            musicAttributes.ctypes.data_as(\
                                POINTER(ctypes.c_float)), \
                            ctypes.c_bool(True))
    T2 = time.time()

    if (outputWave):
        fw.writeframes(data[0 : ].tostring())
        fw.close()

    errorCheck(error, 4)
    if (showDetails):
        print "%fs" % float(T2 - T1)
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]
        print "\n"
    del data
    return path

def getMusicAttributes(data):
    global musicAttributes

    dataLength = len(data)
    error = -1
    T1 = time.time()
    error = dll.getMusicAttributes(data.ctypes.data_as(POINTER(ctypes.c_int16)), \
                                   ctypes.c_int(dataLength), \
                                   musicAttributes.ctypes.data_as(\
                                       POINTER(ctypes.c_float)))
    T2 = time.time()
    
    errorCheck(error, 5)

def beatFilter(data, fileName, beatFilterStyle):
    global musicAttributes
    global PROJECT_DIR
    global PROJECT_TEST_DIR
    
    unicodeFlag = isinstance(fileName, unicode)
    #人为输入文件名
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        #读取程序返回的unicode编码的路径
        path = fileName
        fileName = fileName.strip().split("\\")[-1]

    getMusicAttributes(data)

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(beatFilter) music loaded!"
        print "-------------------------------------------------------------------"
    
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]

    if (outputWave):
        fileName = fileName[0 : len(fileName) - \
                        len(fileName.strip().split(".")[-1]) - 1] \
                        + "_beatFilter(" + str(beatFilterStyle) + ").wav"
        path = PROJECT_TEST_DIR + fileName
        if os.path.exists(path):
            os.remove(path)
            print "same file detected, delete"
        fw = wave.open(path, "w")
        fw.setparams(params)

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start beatFilter"
        print "-------------------------------------------------------------------"

    dataLength = len(data)
    error = -1
    T1 = time.time()
    error = dll.beatFilter(data.ctypes.data_as(POINTER(ctypes.c_int16)), \
                           ctypes.c_int(dataLength), \
                           ctypes.c_bool(beatFilterStyle), \
                           musicAttributes.ctypes.data_as(\
                               POINTER(ctypes.c_float)), \
                               ctypes.c_bool(True))
    T2 = time.time()

    if (outputWave):
        fw.writeframes(data[0 : ].tostring())
        fw.close()
    
    errorCheck(error, 6)
    if (showDetails):
        print "%fs" % float(T2 - T1)
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]
        print "\n"
    del data
    return path

def getBeats(data, fileName, alpha_h = 0.5, alpha_l = 0.0001):
    global musicAttributes
    global PROJECT_DIR
    global TOTAL_TIME
    global PROJECT_TEST_DIR
    global outputWave
    global showDetails
    
    unicodeFlag = isinstance(fileName, unicode)
    #人为输入文件名
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        #读取程序返回的unicode编码的路径
        path = fileName
        fileName = fileName.strip().split("\\")[-1]

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(getBeats) music loaded!"
        print "-------------------------------------------------------------------"
    
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]

    if (outputWave):
        fileName = fileName[0 : len(fileName) - \
                        len(fileName.strip().split(".")[-1]) - 1] \
                        + "_getBeats(" + str(alpha_h) + \
                        "-" + str(alpha_l) + ").wav"
        path = PROJECT_TEST_DIR + fileName
        if os.path.exists(path):
            os.remove(path)
            print "same getBeats.wav file detected, delete"
        fw = wave.open(path, "w")
        fw.setparams(params)

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start getBeats"
        print "-------------------------------------------------------------------"

    dataLength = len(data)
    error = -1
    T1 = time.time()
    error = dll.getBeats(data.ctypes.data_as(POINTER(ctypes.c_int16)), \
                         ctypes.c_int(dataLength), \
                         ctypes.c_float(alpha_h), \
                         ctypes.c_float(alpha_l), \
                         musicAttributes.ctypes.data_as(\
                             POINTER(ctypes.c_float)), \
                         ctypes.c_bool(True))
    T2 = time.time()

    if (outputWave):
        fw.writeframes(data[0 : ].tostring())
        fw.close()
    
    errorCheck(error, 7)
    if (showDetails):
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]
        print "\n"
    del data
    return path

def countBeats(data, fileName, windowSize, style): # true: lowPass false: highPass
    global musicAttributes
    global PROJECT_DIR
    global PROJECT_TEST_DIR
    
    unicodeFlag = isinstance(fileName, unicode)
    #人为输入文件名
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        #读取程序返回的unicode编码的路径
        path = fileName
        fileName = fileName.strip().split("\\")[-1]

    getMusicAttributes(data)
    
    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(countBeats) music loaded!"
        print "-------------------------------------------------------------------"

        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]

    if (outputWave):
        fileName = fileName[0 : len(fileName) - \
                        len(fileName.strip().split(".")[-1]) - 1] \
                        + "_countBeats(" + str(windowSize) + \
                        "-" + str(style) + ").wav"
        path = PROJECT_TEST_DIR + fileName
        if os.path.exists(path):
            os.remove(path)
            print "same file detected, delete"
        fw = wave.open(path, "w")
        fw.setparams(params)

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start countBeats"
        print "-------------------------------------------------------------------"

    dataLength = len(data)
    error = -1
    T1 = time.time()
    error = dll.countBeats(data.ctypes.data_as(POINTER(ctypes.c_int16)), \
                           ctypes.c_int(dataLength), \
                           ctypes.c_int(windowSize), \
                           ctypes.c_bool(style), \
                           musicAttributes.ctypes.data_as(\
                               POINTER(ctypes.c_float)), \
                           ctypes.c_bool(True))
    T2 = time.time()

    if (outputWave):
        fw.writeframes(data[0 : ].tostring())
        fw.close()
    
    errorCheck(error, 8)
    if (showDetails):
        print "%fs" % float(T2 - T1)
        print "FILE_DIR: ", PROJECT_DIR + fileName
        print "AVE_BIG_WINDOW_LENGTH: ", musicAttributes[0]
        print "AVERAGE_INTENSITY_L: ", musicAttributes[1]
        print "AVERAGE_INTENSITY_R: ", musicAttributes[2]
        print "AVERAGE_INTENSITY: ", musicAttributes[3]
        print "MAX_VALUE_L: ", musicAttributes[4]
        print "MAX_VALUE_R: ", musicAttributes[5]
        print "MAX_VALUE: ", musicAttributes[6]
        print "MIN_VALUE_L: ", musicAttributes[7]
        print "MIN_VALUE_R: ", musicAttributes[8]
        print "MIN_VALUE: ", musicAttributes[9]
        print "MAX_L / AVERAGE_L: ", musicAttributes[10]
        print "MAX_R / AVERAGE_R: ", musicAttributes[11]
        print "MAX / AVERAGE: ", musicAttributes[12]
        print "NONZERO_EXCEPTSILENCE_LENGTH: ", musicAttributes[13]
        print "AVE_BEAT_INTENSITY_L: ", musicAttributes[14]
        print "AVE_BEAT_INTENSITY_H: ", musicAttributes[15]
        print "BEATS_L: ", musicAttributes[16]
        print "BEATS_H: ", musicAttributes[17]
        print "LOW_PASS_AVE_REDUCE_RATIO: ", musicAttributes[18]
        print "HIGH_PASS_AVE_REDUCE_RATIO: ", musicAttributes[19]
        print "\n"
    del data
    return path

def dwt(data, model):
    CA, CD = pywt.dwt(data, model)
    return CA, CD

def fft(data): 
    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(fft) music loaded!"
        print "-------------------------------------------------------------------"

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start fft"
        print "-------------------------------------------------------------------"

    fftSize = 1
    expIndex = 0
    while (fftSize <= len(data) and expIndex <= 24):
        fftSize = fftSize * 2
        expIndex += 1
    fftSize = fftSize / 2
    expIndex -= 1

    error = -1

    dataLeft = len(data) - fftSize
    partData = data[dataLeft / 2 : fftSize + dataLeft / 2]
    del data

    smoothLen = int(1.0 / 20 * float(fftSize))

    stableLen = fftSize - (smoothLen * 2)
    
    factor = np.linspace(0, 1, smoothLen)
    factor = np.append(factor, np.ones(stableLen))
    factor = np.append(factor, np.linspace(1, 0, smoothLen))

    partData = np.short(partData * factor)
    del factor

    fftData = np.fft.rfft(partData) / fftSize
    del partData
    
    dataLength = len(fftData)

    amplitude = np.zeros(dataLength, dtype = ctypes.c_float)

    real = fftData.real
    imag = fftData.imag
    del fftData
    
    error = dll.getAmplitude(real.ctypes.data_as(POINTER(ctypes.c_double)), \
                             imag.ctypes.data_as(POINTER(ctypes.c_double)), \
                             ctypes.c_int(dataLength), \
                             amplitude.ctypes.data_as(POINTER(ctypes.c_float)))

    errorCheck(error, 11)

    '''
    freqs = np.linspace(0, 44100 / 2, dataLength)
    pl.figure(figsize=(8,4))
    pl.subplot(111)
    pl.plot(freqs, amplitude)
    pl.xlabel(u"频率(Hz)")
    pl.show()
    '''

    del real
    del imag
    del fftSize
    del expIndex
    del error
    del dataLeft
    del smoothLen
    del stableLen
    
    return amplitude

def getMelIntensity(data, framerate): 
    global MEL_INTENSITY
    global melWindow
    global melScale
    
    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(getMelIntensity) music loaded!"
        print "-------------------------------------------------------------------"

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start getMelIntensity"
        print "-------------------------------------------------------------------"
    
    error = -1
    dataLength = len(data)
    overlap = 0.6
    error = dll.getMelIntensity(data.ctypes.data_as(POINTER(ctypes.c_float)), \
                         ctypes.c_int(dataLength), \
                         melScale.ctypes.data_as(\
                            POINTER(ctypes.c_int)), \
                         ctypes.c_int(melWindow), \
                         MEL_INTENSITY.ctypes.data_as(\
                            POINTER(ctypes.c_float)), \
                         ctypes.c_int(framerate), \
                         ctypes.c_float(overlap))
    
    errorCheck(error, 12)

    if (showDetails):
        print "getMFCCs takes %fs" % float(T2 - T1)
        print "\n"
        
    del data
    del error
    del dataLength

def getMFCCs(data): 
    global MFCCs
    global MFCCs_dimension
    global melWindow
    
    if (showDetails):
        print "-------------------------------------------------------------------"
        print "(getMFCCs) music loaded!"
        print "-------------------------------------------------------------------"

    if (showDetails):
        print "-------------------------------------------------------------------"
        print "start getMFCCs"
        print "-------------------------------------------------------------------"
    
    error = -1
    error = dll.getMFCCs(data.ctypes.data_as(POINTER(ctypes.c_float)), \
                         ctypes.c_int(MFCCs_dimension), \
                         ctypes.c_int(melWindow), \
                         MFCCs.ctypes.data_as(\
                            POINTER(ctypes.c_float)))
    
    errorCheck(error, 10)

    if (showDetails):
        print "getMFCCs takes %fs" % float(T2 - T1)
        print "\n"
        
    del data
    del error

def getWltFeatures(incomingData, showFigure):
    global CHANNEL
    global dataToWrite
    global musicAttributes

    smoothPara = 9 # odd number
    featureLen = -1

    getMusicAttributes(incomingData)

    incomingData.shape = -1, CHANNELS
    incomingData = incomingData.T
    
    for channel in xrange(CHANNELS):
        # 80hz
        ori = incomingData[channel]        
        CA, CD = dwt(ori, dwtModel)#CA
        CA, CD = dwt(CA, dwtModel)#CA
        CA, CD = dwt(CA, dwtModel)#CA
        CA, CD = dwt(CA, dwtModel)#CA
        CA, CD = dwt(CA, dwtModel)#CA
        CA, CD = dwt(CA, dwtModel)#CA
#        CA, CD = dwt(CA, dwtModel)#CA
#        CA, CD = dwt(CA, dwtModel)#CA
#        CA, CD = dwt(CD, dwtModel)#CD
#        CA, CD = dwt(CD, dwtModel)#CD
#        CA, CD = dwt(CA, dwtModel)#CA
#        CA, CD = dwt(CD, dwtModel)#CA
#        CA, CD = dwt(CA, dwtModel)#CA
#        print len(CA)
#        CA = CA / musicAttributes[channel + 1]
        CA = CA / 20
        temp = np.sqrt(CA ** 2)
        
        feature = temp
#        feature = feature / (musicAttributes[channel + 1] / 35)
        feature.shape = -1, 10
        var = np.zeros(0)
        for line in feature:
            var = np.append(var, line.var())
        var.shape = -1, 5
        varVar = np.zeros(0)
        for line in var:
            varVar = np.append(varVar, line.var())
        temp = varVar
        if featureLen == -1: # feature has same length
            featureLen = len(varVar)
        for i in xrange(featureLen):
            if i < int(smoothPara / 2) or i >= featureLen - int(smoothPara / 2):
                varVar[i] = varVar[i]
            else:
                for smoothNum in xrange((smoothPara - 1) / 2):
                    varVar[i] = varVar[i] + temp[i - smoothNum - 1]
                    varVar[i] = varVar[i] + temp[i + smoothNum + 1]
                varVar[i] = varVar[i] / smoothPara
            dataToWrite.append(str(varVar[i]) + "\t")
        

        
        print len(dataToWrite)
        
        if showFigure:
            freqs = np.linspace(0, 60, len(varVar))
            pl.figure(figsize=(8,4))
            pl.subplot(111)
            pl.plot(freqs, varVar)
            pl.xlabel(u"频率(Hz)")
            pl.show()
        
        del CA
        del CD
        del feature
    
    for channel in xrange(CHANNELS):
        # 690hz
        ori = incomingData[channel]        
        CA, CD = dwt(ori, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        print len(CD)
#        CD = CD / musicAttributes[channel + 1]
        CD = CD / 8
        temp = np.sqrt(CD ** 2)
        
        feature = temp
#        feature = feature / (musicAttributes[channel + 1] / 35)
        feature.shape = -1, 10
        var = np.zeros(0)
        for line in feature:
            var = np.append(var, line.var())
        var.shape = -1, 5
        varVar = np.zeros(0)
        for line in var:
            varVar = np.append(varVar, line.var())
        temp = varVar
        if featureLen == -1: # feature has same length
            featureLen = len(varVar)
        for i in xrange(featureLen):
            if i < int(smoothPara / 2) or i >= featureLen - int(smoothPara / 2):
                varVar[i] = varVar[i]
            else:
                for smoothNum in xrange((smoothPara - 1) / 2):
                    varVar[i] = varVar[i] + temp[i - smoothNum - 1]
                    varVar[i] = varVar[i] + temp[i + smoothNum + 1]
                varVar[i] = varVar[i] / smoothPara
            dataToWrite.append(str(varVar[i]) + "\t")
        
        if showFigure:
            freqs = np.linspace(0, 60, len(CA))
            pl.figure(figsize=(8,4))
            pl.subplot(111)
            pl.plot(freqs, CD)
            pl.xlabel(u"频率(Hz)")
            pl.show()     
        del CA
        del CD
        del feature

    for channel in xrange(CHANNELS):
        # 1000hz
        ori = incomingData[channel]
        CA, CD = dwt(ori, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        print len(CD)
#        CD = CD / musicAttributes[channel + 1]
        CD = CD / 4
        temp = np.sqrt(CD ** 2)
        
        feature = temp
#        feature = feature / (musicAttributes[channel + 1] / 35)
        feature.shape = -1, 10
        var = np.zeros(0)
        for line in feature:
            var = np.append(var, line.var())
        var.shape = -1, 5
        varVar = np.zeros(0)
        for line in var:
            varVar = np.append(varVar, line.var())
        temp = varVar
        if featureLen == -1: # feature has same length
            featureLen = len(varVar)
        for i in xrange(featureLen):
            if i < int(smoothPara / 2) or i >= featureLen - int(smoothPara / 2):
                varVar[i] = varVar[i]
            else:
                for smoothNum in xrange((smoothPara - 1) / 2):
                    varVar[i] = varVar[i] + temp[i - smoothNum - 1]
                    varVar[i] = varVar[i] + temp[i + smoothNum + 1]
                varVar[i] = varVar[i] / smoothPara
            dataToWrite.append(str(varVar[i]) + "\t")
        
        if showFigure:
            freqs = np.linspace(0, 60, len(CD))
            pl.figure(figsize=(8,4))
            pl.subplot(111)
            pl.plot(freqs, CD)
            pl.xlabel(u"频率(Hz)")
            pl.show()     
        del CA
        del CD
        del feature

    for channel in xrange(CHANNELS):
        # 4500hz
        ori = incomingData[channel]
        CA, CD = dwt(ori, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CD, dwtModel)
        CA, CD = dwt(CD, dwtModel)
        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        print len(CA)
#        CA = CA/ musicAttributes[channel + 1]
        CA = CA / 1
        temp = np.sqrt(CA ** 2)
        
        feature = temp
#        feature = feature / (musicAttributes[channel + 1] / 35)
        feature.shape = -1, 10
        var = np.zeros(0)
        for line in feature:
            var = np.append(var, line.var())
        var.shape = -1, 5
        varVar = np.zeros(0)
        for line in var:
            varVar = np.append(varVar, line.var())
        temp = varVar
        if featureLen == -1: # feature has same length
            featureLen = len(varVar)
        for i in xrange(featureLen):
            if i < int(smoothPara / 2) or i >= featureLen - int(smoothPara / 2):
                varVar[i] = varVar[i]
            else:
                for smoothNum in xrange((smoothPara - 1) / 2):
                    varVar[i] = varVar[i] + temp[i - smoothNum - 1]
                    varVar[i] = varVar[i] + temp[i + smoothNum + 1]
                varVar[i] = varVar[i] / smoothPara
            dataToWrite.append(str(varVar[i]) + "\t")
        
        if showFigure:
            freqs = np.linspace(0, 60, len(CA))
            pl.figure(figsize=(8,4))
            pl.subplot(111)
            pl.plot(freqs, CA)
            pl.xlabel(u"频率(Hz)")
            pl.show()     
        del CA
        del CD
        del feature

    for channel in xrange(CHANNELS):
        # 16000hz
        ori = incomingData[channel]
        CA, CD = dwt(ori, dwtModel)
        CA, CD = dwt(CD, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CD, dwtModel)
        CA, CD = dwt(CD, dwtModel)
        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        print len(CA)
#        CA = CA / musicAttributes[channel + 1]
        CA = CA / 0.2
        temp = np.sqrt(CA ** 2)
        
        feature = temp
#        feature = feature / (musicAttributes[channel + 1] / 35)
        feature.shape = -1, 10
        var = np.zeros(0)
        for line in feature:
            var = np.append(var, line.var())
        var.shape = -1, 5
        varVar = np.zeros(0)
        for line in var:
            varVar = np.append(varVar, line.var())
        temp = varVar
        if featureLen == -1: # feature has same length
            featureLen = len(varVar)
        for i in xrange(featureLen):
            if i < int(smoothPara / 2) or i >= featureLen - int(smoothPara / 2):
                varVar[i] = varVar[i]
            else:
                for smoothNum in xrange((smoothPara - 1) / 2):
                    varVar[i] = varVar[i] + temp[i - smoothNum - 1]
                    varVar[i] = varVar[i] + temp[i + smoothNum + 1]
                varVar[i] = varVar[i] / smoothPara
            dataToWrite.append(str(varVar[i]) + "\t")
        
        if showFigure:
            freqs = np.linspace(0, 60, len(CA))
            pl.figure(figsize=(8,4))
            pl.subplot(111)
            pl.plot(freqs, CA)
            pl.xlabel(u"频率(Hz)")
            pl.show()     
        del CA
        del CD
        del feature

    for channel in xrange(CHANNELS):
        # 18000hz
        ori = incomingData[channel]
        CA, CD = dwt(ori, dwtModel)
        CA, CD = dwt(CD, dwtModel)
        CA, CD = dwt(CD, dwtModel)
        CA, CD = dwt(CA, dwtModel)
        CA, CD = dwt(CD, dwtModel)
        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        CA, CD = dwt(CA, dwtModel)
#        CA, CD = dwt(CD, dwtModel)
#        print len(CA)
#        CA = CA / musicAttributes[channel + 1]
        CA = CA / 0.5
        temp = np.sqrt(CA ** 2)
        
        feature = temp
#        feature = feature / (musicAttributes[channel + 1] / 35)
        feature.shape = -1, 10
        var = np.zeros(0)
        for line in feature:
            var = np.append(var, line.var())
        var.shape = -1, 5
        varVar = np.zeros(0)
        for line in var:
            varVar = np.append(varVar, line.var())
        temp = varVar
        if featureLen == -1: # feature has same length
            featureLen = len(varVar)
        for i in xrange(featureLen):
            if i < int(smoothPara / 2) or i >= featureLen - int(smoothPara / 2):
                varVar[i] = varVar[i]
            else:
                for smoothNum in xrange((smoothPara - 1) / 2):
                    varVar[i] = varVar[i] + temp[i - smoothNum - 1]
                    varVar[i] = varVar[i] + temp[i + smoothNum + 1]
                varVar[i] = varVar[i] / smoothPara
            dataToWrite.append(str(varVar[i]) + "\t")
        
        if showFigure:
            freqs = np.linspace(0, 60, len(CA))
            pl.figure(figsize=(8,4))
            pl.subplot(111)
            pl.plot(freqs, CA)
            pl.xlabel(u"频率(Hz)")
            pl.show()     
        del CA
        del CD
        del feature
    

def frequencyAnalysis(incomingData, fileName):
    global PROJECT_DIR
    global DATA_NAME
    global AudioFileCount
    global currentAudioFileCount
    global EXCEL_LINE
    global PROJECT_TEST_DIR
    global PROJECT_PCASPACE_DIR
    global showDetails
    global dataToWrite
    global TOTAL_TIME
    global CHANNELS
    global MFCCs
    global FIRST_DIFF
    global max_frequency
    global max_melFrequency
    global melScale
    global framerate
    global dwtParts
    global dwtIterations_time
    global dwtIterations_freq
    global dwtModel
    global MEL_INTENSITY
    global wltAnalysisLength
    global chunkNum
    global getPCATarget
    global targetEigVals
    global targetEigVects_L
    global targetEigVects_R

    dataToWrite = []
    fileNameToTXT = fileName.decode("gb18030")
    fileNameToTXT = fileNameToTXT.encode("UTF-8")# for Korean
    dataToWrite.append(fileNameToTXT + "\t")
    resLen = 0
    CListLen = 0

    txtPath = PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR

    dataLength = len(incomingData)
#    print dataLength
    chunkSize = int(dataLength / chunkNum)
#    incomingData.shape = -1, CHANNELS
#    incomingData = incomingData.T
    if chunkSize %2 != 0:
        chunkSize -= 1
#    print chunkSize

    for chunk in xrange(chunkNum):
        chunkData = incomingData[chunk * chunkSize : (chunk + 1) * chunkSize]
#        print chunkData
#        print len(chunkData)
        chunkData.shape = -1, CHANNELS
        chunkData = chunkData.T
        for channel in xrange(CHANNELS):
            wltParts = 1
            res = np.zeros(0)
            T1 = time.time()
#            print len(chunkData[channel])
            for iteration in xrange(dwtIterations_freq):
    #            print "iteration: ", iteration
                CAList = np.zeros(0)
                CDList = np.zeros(0)
#                print "iterations", iteration
#                print "wltParts", wltParts
                for wltPart in xrange(wltParts):
                    if iteration == 0:
                        dataToDwt = chunkData[channel]
                    else:
                        if wltPart != wltParts - 1:
                            dataToDwt = \
                                res[int(float (wltPart) / float (wltParts) * resLen) : \
                                int(float (wltPart + 1) / float (wltParts) * resLen)]
                        if wltPart == wltParts - 1:
                            dataToDwt = \
                                res[int(float (wltPart) / float (wltParts) * resLen) : ]
                    CA, CD = dwt(dataToDwt, dwtModel)
                    CAList = np.append(CAList, CA)
                    CDList = np.append(CDList, CD)
                    del CA
                    del CD
                    
                del dataToDwt   
                del res # release res
                res = np.zeros(0)
                CListLen = len(CAList)
#                print "CListLen", CListLen
                for wltPart in xrange(wltParts):
                    if iteration == 0:
                        res = np.append(res, CAList)
                        res = np.append(res, CDList)
                        
                    else:
                        if wltPart != wltParts - 1:
                            res = np.append(res, CAList[ \
                                int(float (wltPart) / float (wltParts) * CListLen) : \
                                int(float(wltPart + 1) / float (wltParts) * CListLen)])
                            res = np.append(res, CDList[ \
                                int(float (wltPart) / float (wltParts) * CListLen) : \
                                int(float (wltPart + 1) / float (wltParts) * CListLen)])
                        if wltPart == wltParts - 1:
                            res = np.append(res, CAList[ \
                                int(float (wltPart) / float (wltParts) * CListLen) : \
                                int(float(wltPart + 1) / float (wltParts) * CListLen)])
                            res = np.append(res, CDList[ \
                                int(float (wltPart) / float (wltParts) * CListLen) : \
                                int(float (wltPart + 1) / float (wltParts) * CListLen)])
                        
                
                            
                resLen = len(res)
                wltParts = wltParts * 2
                del CAList
                del CDList
                CAList = np.zeros(0)
                CDList = np.zeros(0)
                           
#            print resLen   
#            print res
#            print wltParts / 2

            T2 = time.time()
            print "Wlt takes %fs for L/R channel" % float(T2 - T1)

            T1 = time.time()
#            wltParts = wltParts / 2
#            print wltParts
            for wltPart in xrange(wltParts):
                dataToFFT = res[int(float (wltPart) / float (wltParts) * resLen) : \
                                int(float(wltPart + 1) / float (wltParts) * resLen)]
                fftDataInDB = fft(dataToFFT)
                getMelIntensity(fftDataInDB, framerate)
                getMFCCs(MEL_INTENSITY)
                MFCCs = MFCCs * melParaArr[wltPart]
                for i in xrange(MFCCs_dimension):
                    if i == 0:
                        FIRST_DIFF[i] = MFCCs[i]
                    if i >= 1:
                        FIRST_DIFF[i] = MFCCs[i] - MFCCs[i - 1]
                for i in xrange(MFCCs_dimension):
                    dataToWrite.append(str(MFCCs[i]) + "\t")
                for i in xrange(MFCCs_dimension):
                    dataToWrite.append(str(FIRST_DIFF[i]) + "\t")
                    
            T2 = time.time()
            print "FFT takes %fs for L/R channel" % float(T2 - T1)
            
            del res
            del dataToFFT
            del fftDataInDB
            del CAList
            del CDList

    T1 = time.time()    
    lengthDiff = dataLength - wltAnalysisLength
#    print str(lengthDiff) + " = " + str(dataLength)\
#          + " + " + str(wltAnalysisLength)
    while lengthDiff < 0:
        if abs(lengthDiff) <= dataLength:
            incomingData = np.append(incomingData, \
                                incomingData[dataLength - abs(lengthDiff) : ])
            dataLength = len(incomingData)
            lengthDiff = dataLength - wltAnalysisLength
        else:
            incomingData = np.append(incomingData, incomingData)
            dataLength = len(incomingData)
            lengthDiff = dataLength - wltAnalysisLength

    if lengthDiff >= 0:
        lengthDiff = int(lengthDiff / 2)
        incomingData = incomingData[lengthDiff : wltAnalysisLength + lengthDiff]
#    print "len(incomingData): ", len(incomingData)

#    getWltFeatures(incomingData, False)

#    getPCAFeature(incomingData)

    
    #get SVD features
    incomingData.shape = -1, CHANNELS
    incomingData = incomingData.T
    SVDFeature = zeros(0)
    scale = initMelScale(pow(2, 10))
    for channel in xrange(CHANNELS):
        dwtData = dwtMusic(incomingData[channel], 10)
        for i in xrange(len(dwtData)):
            dwtData[i] = dwtData[i] * scale[i]
#        dwtData = incomingData[channel]
#        dwtData.shape = 1000, -1
        print "dwtData.shape: ", dwtData.shape
        dwtMat = mat(dwtData)
        U, Sigma, VT = np.linalg.svd(np.mat(dwtMat))
        print "Sigma.shape: ", Sigma.shape
        for i in xrange(len(Sigma)):
            SVDFeature = np.append(SVDFeature, Sigma[i])
    print "SVDFeature.shape: ", SVDFeature.shape
    print len(SVDFeature)
    for i in xrange(len(SVDFeature)):
        SVDFeature[i] = SVDFeature[i] / 2000.0
        dataToWrite.append(str(SVDFeature[i]) + "\t")
    
    '''
    #get PCA features
    incomingData.shape = -1, CHANNELS
    incomingData = incomingData.T
    PCAFeature = zeros(0)
#    scale = initMelScale(pow(2, 10))
    for channel in xrange(CHANNELS):
#        dwtData = dwtMusic(incomingData[channel], 10)
#        for i in xrange(len(dwtData)):
#            dwtData[i] = dwtData[i] * scale[i]
        dwtData = incomingData[channel]
        dwtData.shape = 1000, -1
        print "dwtData.shape: ", dwtData.shape
        dwtMat = mat(dwtData)
        covMat = np.cov(np.mat(dwtData), rowvar = 0)
        eigVals, eigVects = np.linalg.eig(np.mat(covMat))
        print "eigVals.shape: ", eigVals.shape
        for i in xrange(len(eigVals)):
            PCAFeature = np.append(PCAFeature, eigVals[i])
    print "PCAFeature.shape: ", PCAFeature.shape
    print len(PCAFeature)
    for i in xrange(len(PCAFeature)):
        PCAFeature[i] = PCAFeature[i]
        dataToWrite.append(str(PCAFeature[i]) + "\t")
    '''
    T2 = time.time()
    print "get Wlt features takes %fs for L/R channel" % float(T2 - T1)

    del dataLength
    del incomingData
    del txtPath
    del fileNameToTXT

def initMelScale(dwtParts):
    global max_melFrequency
    scale = np.arange(dwtParts, dtype = ctypes.c_float)
    for i in xrange(dwtParts):
        scale[i] = (float(max_melFrequency) / float(dwtParts)) * (i + 1)
#    print melParaArr
    for i in xrange(dwtParts):
        scale[i] = (math.pow(10, float(scale[i]) / 2595) - 1) * 700
    maxDiff = scale[-1] - scale[-2]
#    print melParaArr
#    print maxDiff
    for i in xrange(dwtParts):
        if i == 0:
            temp = scale[i]
            scale[i] = maxDiff / scale[i]
        else:
            res = maxDiff / (scale[i] - temp)
            temp = scale[i]
            scale[i] = res
    scale = scale / scale[0]
#    print scale
    return scale

def getPCAFeature(incomingData):
    incomingData.shape = -1, CHANNELS
    incomingData = incomingData.T

    testMusic = 59
    data = incomingData[0]
    CA, CD = dwt(data, dwtModel)
#    CA, CD = dwt(CA, dwtModel)
#    CA, CD = dwt(CA, dwtModel)
#    CA, CD = dwt(CD, dwtModel)
#    CA, CD = dwt(CD, dwtModel)
    data = CD
    fr = open(PROJECT_DOC_DIR + u"meanVals.txt", "r")
    for line in fr:
        line = line.strip().split("\t")
    fr.close()
    meanVals = line
    del line
    meanVals = np.array(meanVals, np.float)
    data = np.absolute(data)
    data = np.mat(data)
    meanVals = np.mat(meanVals)
#    print "data.shape: ", data.shape
#    print "meanVals.shape: ", meanVals.shape
#    meanRemoved = data - meanVals
    meanRemoved = data
    meanRemoved = np.array(meanRemoved)
    meanRemoved = meanRemoved[0]
    del data
    fr = open(PROJECT_DOC_DIR + PROJECT_PCASPACE_DIR, "r")
    PCASpace = np.zeros(0)
    for line in fr:
        line = line.strip().split("\t")
    PCASpace = np.append(PCASpace, np.array(line, np.float))
    PCASpace.shape = -1, testMusic
    del line
#    print "PCASpace.shape:", PCASpace.shape
    PCASpace = PCASpace.T
    Omega = []
    for i in xrange(len(PCASpace)):
        Omega = np.inner(PCASpace[i], meanRemoved)
        dataToWrite.append(str(Omega) + "\t")
    print len(dataToWrite)

def dwtMusic(data, level):
    dwtData = np.zeros(0)
    tmpData = data
    tmpData.shape = -1, len(tmpData)
    Clen = 0
    for dwtLevel in xrange(level):
        for i in xrange(len(tmpData)):
            CA, CD = dwt(tmpData[i], "db4")
            CA = np.short(CA)
            CD = np.short(CD)
            Clen = len(CA)
            dwtData = np.append(dwtData, CA)
            dwtData = np.append(dwtData, CD)
        dwtData.shape = -1, Clen
        tmpData = dwtData
        dwtData = np.zeros(0)
    dwtData = tmpData
#    print dwtData.shape
    return dwtData

def beatsDetection(incomingData, fileName, flag = ""):
    global PROJECT_DIR
    global DATA_NAME
    global AudioFileCount
    global currentAudioFileCount
    global EXCEL_LINE
    global PROJECT_TEST_DIR
    global dataToWrite
    global MFCCs_dimension
    global CHANNELS
    global dwtParts

    txtPath = PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR
    fw = open(txtPath, "a")

#    data = incomingData
#    del incomingData
    
    getMusicAttributes(incomingData)
    # AVERAGE_INTENSITY
    AVE = musicAttributes[3]
    dataToWrite.append(str(AVE) + "\t")

    print "beatsDetection (%d/%d)" % (currentAudioFileCount, AudioFileCount)

    if flag is "":
        path = getBeats(incomingData, fileName, 0.5, 0.005)

    if flag is "low": # for test
        path = lowPassFilter(incomingData, fileName, 0.005, 2)
        path = beatFilter(incomingData, path, True) #True: lowPass False: highPass
        path = countBeats(incomingData, path, int(44100 * 0.08), True)#True: lowPass False: highPass

    if flag is "high": # for test
        path = highPassFilter(incomingData, fileName, 0.9)
        path = beatFilter(incomingData, path, False) #True: lowPass False: highPass
        path = countBeats(incomingData, path, int(44100 * 0.08), False)#True: lowPass False: highPass
    
    #AVE_BIG_WINDOW_LENGTH
    dataToWrite.append(str(musicAttributes[0] * 1000) + "\t")
    
    # LOW_PASS_AVE_REDUCE_RATIO
    dataToWrite.append(str(musicAttributes[18] * 1000) + "\t")

    # HIGH_PASS_AVE_REDUCE_RATIO
    dataToWrite.append(str(musicAttributes[19] * 1000) + "\t")

    # BEATS_L
    dataToWrite.append(str(musicAttributes[16] * 1000) + "\t")

    # BEATS_H
    dataToWrite.append(str(musicAttributes[17] * 1000) + "\t")

    # AVE_BEAT_INTENSITY_L
    dataToWrite.append(str(musicAttributes[14] * 0.2) + "\t")
    
    # AVE_BEAT_INTENSITY_H
    dataToWrite.append(str(musicAttributes[15] * 0.2) + "\n")

    print len(dataToWrite)
    fw.writelines(dataToWrite)
    fw.close()
    del dataToWrite

    del fw
    del AVE
    del txtPath

def testFrequencyDomain(fileName, OPWave = False, SDetails = False):
    global outputWave
    global showDetails

    outputWave = OPWave
    showDetails = SDetails

    tempFileName = unicode(fileName, "UTF-8")
    tempFileName = tempFileName.encode("gb18030")

    frequencyAnalysis(getData(fileName), tempFileName)

def testTimeDomain(fileName, OPWave = True, SDetails = False, word = ""):
    global outputWave
    global showDetails
    global flag

    flag = word # flag is the parameter shows the testing type, whether low pass
                # or high pass

    outputWave = OPWave
    showDetails = SDetails

    fileName = unicode(fileName, "utf-8")
    fileName = fileName.encode("gb18030")

    beatsDetection(getData(fileName), fileName, flag)

def testAll(OPWave = False, SDetails = False):
    global PROJECT_DIR
    global PROJECT_RAWFEAT_DIR
    global AudioFileCount
    global currentAudioFileCount
    global TOTAL_TIME
    global outputWave
    global showDetails
    global MFCCs_dimension

    outputWave = OPWave
    showDetails = SDetails

    path = PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR
    lineInTXT = 0
    if os.path.exists(path):
        fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
        completedMusicName = []
        for line in fr:
            lineInTXT += 1
            line = line.strip().split("\t")
            completedMusicName.append(line[0])
        fr.close()

    fileList = os.listdir(PROJECT_DIR[:len(PROJECT_DIR) - 1])
    '''
    fw = open("H://fileList.txt", "w")
    for i in xrange(len(fileList)):
        fw.write(str(fileList[i]) + "\n")
    fw.close()
    '''

    for audioFile in fileList:
        if (audioFile.strip().split(".")[-1] == "wav") or \
           (audioFile.strip().split(".")[-1] == "WAV") or \
           (audioFile.strip().split(".")[-1] == "mp3") or \
           (audioFile.strip().split(".")[-1] == "MP3"):
            AudioFileCount += 1
            
    if os.path.exists(path):
        currentAudioFileCount += len(completedMusicName)
    
    fileCount = 0
    finish = False
    for audioFile in fileList:
        if (audioFile.strip().split(".")[-1] == "wav") or \
           (audioFile.strip().split(".")[-1] == "WAV") or \
           (audioFile.strip().split(".")[-1] == "mp3") or \
           (audioFile.strip().split(".")[-1] == "MP3"):
            fileCount += 1
            if(fileCount > currentAudioFileCount):
                currentAudioFileCount += 1
                data = getData(audioFile.encode("gb18030"))
                tempData = np.copy(data)
                T1 = time.time()
                frequencyAnalysis(data, audioFile.encode("gb18030"))
                #data = None
                T2 = time.time()
                print "frequency Analysis takes %fs" % (T2 - T1)
                TOTAL_TIME = TOTAL_TIME + (T2- T1)
                del T1
                del T2
                
                T1 = time.time()
                beatsDetection(tempData, audioFile.encode("gb18030"))
                #tempData = None
                T2 = time.time()
                print "beats analysis takes %fs\n"  % float(T2 - T1)
                TOTAL_TIME = TOTAL_TIME + (T2- T1)
                del T1
                del T2
                #tempData = None
                del tempData
                #data = None
                del data
                #gc.collect()
            else:
                # last audio file
                # -1 ~ -6 are folders
                for i in xrange(AudioFileCount):
                    if fileList[-1 * i - 1].strip().split(".")[-1] == "wav" or \
                       fileList[-1 * i - 1].strip().split(".")[-1] == "WAV" or \
                       fileList[-1 * i - 1].strip().split(".")[-1] == "mp3" or \
                       fileList[-1 * i - 1].strip().split(".")[-1] == "MP3":
                        if audioFile == fileList[-1 * i - 1]: 
                            print "All songs are analysized"
                            finish = True
                            break
                    if finish:
                        break
    if(AudioFileCount - lineInTXT) != 0:
        AVE_time = TOTAL_TIME / (AudioFileCount - lineInTXT)
        print "each music takes %fs." % AVE_time
    #fileList = None
    del fileList
    #gc.collect()
    
def similarity(fileName, word = ""):
#    global MUSIC_CHUNK_NUM
    global CHANNELS
    global dwtParts
    global TOTAL_MUSIC
    global MFCCPara_dimension
    global wltPara_dimension
    global beats_dimension
    global frequencyPara_dimension
    
    path, fileName = getPath(fileName)
    fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
    dataSet = []
    lineCount = 0
    for line in fr:
        line = line.strip().split("\t")
        dataSet.append(line)
        lineCount += 1
    fr.close()

    dataDimention = 1 + frequencyPara_dimension + beats_dimension
    res = []
    target = []
    for eachLine in xrange(lineCount):
        if (str)(dataSet[eachLine][0]) == (str)(fileName):
            print fileName
            for i in xrange(dataDimention):
                target.append(dataSet[eachLine][i])
            del dataSet[eachLine]
            lineCount -= 1
            break
#    print target

    if word == "BeatsMFCCs":
        AVE_para = 1500
        lowReduce_para = 8
        highReduce_para = 4
        lowBeat_para = 1
        highBeat_para = 1
        lowBeatsInt_para = 1000
        highBeatsInt_para = 2000
        for eachLine in xrange(lineCount):
            simi = 0
            simi = simi + (float(dataSet[eachLine][frequencyPara_dimension + 1]) / AVE_para - float(target[frequencyPara_dimension+ 1]) / AVE_para) * \
                          (float(dataSet[eachLine][frequencyPara_dimension + 1]) / AVE_para - float(target[frequencyPara_dimension + 1]) / AVE_para)
            simi = simi + (float(dataSet[eachLine][frequencyPara_dimension + 3]) * lowReduce_para - float(target[frequencyPara_dimension + 3]) * lowReduce_para) * \
                          (float(dataSet[eachLine][frequencyPara_dimension + 3]) * lowReduce_para - float(target[frequencyPara_dimension + 3]) * lowReduce_para)
            simi = simi + (float(dataSet[eachLine][frequencyPara_dimension + 4]) * highReduce_para - float(target[frequencyPara_dimension + 4]) * highReduce_para) * \
                          (float(dataSet[eachLine][frequencyPara_dimension + 4]) * highReduce_para - float(target[frequencyPara_dimension + 4]) * highReduce_para) 
            simi = simi + (float(dataSet[eachLine][frequencyPara_dimension + 5]) * lowBeat_para - float(target[frequencyPara_dimension + 5])) * lowBeat_para * \
                          (float(dataSet[eachLine][frequencyPara_dimension + 5]) * lowBeat_para - float(target[frequencyPara_dimension + 5])) * lowBeat_para
            simi = simi + (float(dataSet[eachLine][frequencyPara_dimension + 6]) * highBeat_para - float(target[frequencyPara_dimension + 6])) * highBeat_para * \
                          (float(dataSet[eachLine][frequencyPara_dimension + 6]) * highBeat_para - float(target[frequencyPara_dimension + 6])) * highBeat_para
            simi = simi + (float(dataSet[eachLine][frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target[frequencyPara_dimension + 7]) / lowBeatsInt_para) * \
                          (float(dataSet[eachLine][frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target[frequencyPara_dimension + 7]) / lowBeatsInt_para)
            simi = simi + (float(dataSet[eachLine][frequencyPara_dimension + 8]) / highBeatsInt_para - float(target[frequencyPara_dimension + 8]) / highBeatsInt_para) * \
                          (float(dataSet[eachLine][frequencyPara_dimension + 8]) / highBeatsInt_para - float(target[frequencyPara_dimension + 8]) / highBeatsInt_para)
            simi = math.sqrt(simi)
            dataSet[eachLine].extend([simi])
            res.append(dataSet[eachLine])

        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 1])
        dataSet = None
        del dataSet

        dataSetForMFCCs = []
        top = 10
        for i in xrange(top):
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension \
                              + beats_dimension + 1])).encode("UTF-8")
            dataSetForMFCCs.append(res[i])
        res = None
        del res

        res = []
        for eachLine in xrange(len(dataSetForMFCCs)):
            for i in xrange(frequencyPara_dimension - CLength * 12):
                diff = math.pow((float(dataSetForMFCCs[eachLine][i + 1]) \
                                 - float(target[i + 1])), 2)
                simi += diff
            simi = math.sqrt(simi) * 1
            dataSetForMFCCs[eachLine].extend([simi])
            dataSetForMFCCs[eachLine].extend([(simi + \
                            dataSetForMFCCs[eachLine]\
                            [frequencyPara_dimension \
                             + beats_dimension + 1]) / 2])
            res.append(dataSetForMFCCs[eachLine])
            
        res.sort(key = lambda x:x[frequencyPara_dimension \
                                  + beats_dimension + 3])
        dataSetForMFCCs = None
        del dataSetForMFCCs

        top = 3
        for i in xrange(top):
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension \
                              + beats_dimension + 3])).encode("UTF-8") + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 1])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 2])).encode("UTF-8") + ")"

    if word == "MFCCsBeats":
        for eachLine in xrange(lineCount):
            simi = 0
            for i in xrange(frequencyPara_dimension - CLength * 12):
                diff = math.pow((float(dataSet[eachLine][i + 1]) \
                                 - float(target[i + 1])), 2)
                simi += diff
            simi = math.sqrt(simi) * 0.1
            dataSet[eachLine].extend([simi])
            res.append(dataSet[eachLine])

        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 1])
        dataSet = None
        del dataSet

        dataSetForBeats = []
        top = 20
        for i in xrange(top):
            
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension
                              + beats_dimension + 1])).encode("UTF-8")
            
            dataSetForBeats.append(res[i])
        res = None
        del res

        res = []
        AVE_para = 1500
        lowReduce_para = 8
        highReduce_para = 4
        lowBeat_para = 1
        highBeat_para = 1
        lowBeatsInt_para = 1000
        highBeatsInt_para = 1000
        for eachLine in xrange(len(dataSetForBeats)):
            simi = 0
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 1]) / AVE_para - float(target[frequencyPara_dimension + 1]) / AVE_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 1]) / AVE_para - float(target[frequencyPara_dimension + 1]) / AVE_para)
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 3]) * lowReduce_para - float(target[frequencyPara_dimension + 3]) * lowReduce_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 3]) * lowReduce_para - float(target[frequencyPara_dimension + 3]) * lowReduce_para)
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 4]) * highReduce_para - float(target[frequencyPara_dimension + 4]) * highReduce_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 4]) * highReduce_para - float(target[frequencyPara_dimension + 4]) * highReduce_para)
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 5]) * lowBeat_para - float(target[frequencyPara_dimension + 5])) * lowBeat_para * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 5]) * lowBeat_para - float(target[frequencyPara_dimension + 5])) * lowBeat_para
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 6]) * highBeat_para - float(target[frequencyPara_dimension + 6])) * highBeat_para * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 6]) * highBeat_para - float(target[frequencyPara_dimension + 6])) * highBeat_para
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target[frequencyPara_dimension + 7]) / lowBeatsInt_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target[frequencyPara_dimension + 7]) / lowBeatsInt_para)
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 8]) / highBeatsInt_para - float(target[frequencyPara_dimension + 8]) / highBeatsInt_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 8]) / highBeatsInt_para - float(target[frequencyPara_dimension + 8]) / highBeatsInt_para)

            simi = math.sqrt(simi)
            dataSetForBeats[eachLine].extend([simi])
            dataSetForBeats[eachLine].extend([(simi + \
                            dataSetForBeats[eachLine]\
                            [frequencyPara_dimension \
                             + beats_dimension + 1]) / 2])
            res.append(dataSetForBeats[eachLine])
            
        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 3])
        dataSetForBeats = None
        del dataSetForBeats

        top = 5
        for i in xrange(top):
            
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension \
                   + beats_dimension + 3])).encode("UTF-8") + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 1])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 2])).encode("UTF-8") + ")"
            
    if word == "MFCCsWlt":
        for eachLine in xrange(lineCount):
            CosSimi = getCosineDis(dataSet[eachLine][1 : MFCCPara_dimension + 1], \
                                    target[1 : MFCCPara_dimension + 1])
            eucSimi = 0
            for i in xrange(frequencyPara_dimension - wltPara_dimension):
                diff = math.pow((float(dataSet[eachLine][i + 1]) \
                                 - float(target[i + 1])), 2)
                eucSimi += diff
            eucSimi = math.sqrt(eucSimi) * 8
            simi = eucSimi * (1 - CosSimi)
            dataSet[eachLine].extend([simi])
            res.append(dataSet[eachLine])

        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 1])
        dataSet = None
        del dataSet

        dataSetForWlt = []
        top = int(TOTAL_MUSIC * 0.6)
        for i in xrange(top):

            if i < 5:
                print "----------  NO.%d  ----------" % (i + 1)
                print (str(res[i][0])).encode("gb18030")
                print (str(res[i][frequencyPara_dimension
                              + beats_dimension + 1])).encode("UTF-8")
            
            dataSetForWlt.append(res[i])
        del res

        res = []
        for eachLine in xrange(len(dataSetForWlt)):
            CosSimi = getCosineDis(dataSetForWlt[eachLine][MFCCPara_dimension + 1 : frequencyPara_dimension + 1], \
                                target[MFCCPara_dimension + 1 : frequencyPara_dimension + 1])
            
            eucSimi = 0
            for i in xrange(wltPara_dimension):
                diff = math.pow((float(dataSetForWlt[eachLine][frequencyPara_dimension - wltPara_dimension + i + 1]) / 10 \
                                 - float(target[frequencyPara_dimension - wltPara_dimension + i + 1]) / 10), 2)
                eucSimi += diff
            eucSimi = math.sqrt(eucSimi) * 7000
            simi = eucSimi * (1 - CosSimi)
            dataSetForWlt[eachLine].extend([simi])
            res.append(dataSetForWlt[eachLine])
            
        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 2])
        del dataSetForWlt
        dataSetForBeats = []
        top = int(TOTAL_MUSIC * 0.6)
        for i in xrange(top):
            '''
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension \
                   + beats_dimension + 3])).encode("UTF-8") + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 1])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 2])).encode("UTF-8") + ")"
            '''
            dataSetForBeats.append(res[i])
        del res
        
        res = []
        AVE_para = 500
        lowReduce_para = 10
        highReduce_para = 10
        lowBeat_para = 10
        highBeat_para = 10
        lowBeatsInt_para = 500
        highBeatsInt_para = 500
        for eachLine in xrange(len(dataSetForBeats)):
            CosSimi = getCosineDis(dataSetForBeats[eachLine][frequencyPara_dimension + 1 : frequencyPara_dimension + 9], \
                                target[frequencyPara_dimension + 1 : frequencyPara_dimension + 9])
            eucSimi = 0            
            eucSimi = eucSimi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 1]) / AVE_para - float(target[frequencyPara_dimension + 1]) / AVE_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 1]) / AVE_para - float(target[frequencyPara_dimension + 1]) / AVE_para)
            eucSimi = eucSimi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 2]) - float(target[frequencyPara_dimension + 2])) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 2]) - float(target[frequencyPara_dimension + 2]))
            eucSimi = eucSimi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 3]) * lowReduce_para - float(target[frequencyPara_dimension + 3]) * lowReduce_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 3]) * lowReduce_para - float(target[frequencyPara_dimension + 3]) * lowReduce_para)
            eucSimi = eucSimi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 4]) * highReduce_para - float(target[frequencyPara_dimension + 4]) * highReduce_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 4]) * highReduce_para - float(target[frequencyPara_dimension + 4]) * highReduce_para)            
            eucSimi = eucSimi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 5]) * lowBeat_para - float(target[frequencyPara_dimension + 5])) * lowBeat_para * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 5]) * lowBeat_para - float(target[frequencyPara_dimension + 5])) * lowBeat_para
            eucSimi = eucSimi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 6]) * highBeat_para - float(target[frequencyPara_dimension + 6])) * highBeat_para * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 6]) * highBeat_para - float(target[frequencyPara_dimension + 6])) * highBeat_para
            eucSimi = eucSimi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target[frequencyPara_dimension + 7]) / lowBeatsInt_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target[frequencyPara_dimension + 7]) / lowBeatsInt_para)
            eucSimi = eucSimi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 8]) / highBeatsInt_para - float(target[frequencyPara_dimension + 8]) / highBeatsInt_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 8]) / highBeatsInt_para - float(target[frequencyPara_dimension + 8]) / highBeatsInt_para)

            eucSimi = math.sqrt(eucSimi) * 0.006
            simi = eucSimi * (1 - CosSimi)
            
            dataSetForBeats[eachLine].extend([simi])
            dataSetForBeats[eachLine].extend([(simi + \
                            dataSetForBeats[eachLine]\
                            [frequencyPara_dimension \
                             + beats_dimension + 1] + \
                            dataSetForBeats[eachLine]\
                            [frequencyPara_dimension \
                             + beats_dimension + 2]) / 3])
            res.append(dataSetForBeats[eachLine])
            
        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 4])
        dataSetForBeats = None
        del dataSetForBeats

        top = 5
        for i in xrange(top):

            if(i < top):
                print "----------  NO.%d  ----------" % (i + 1)
                print (str(res[i][0])).encode("gb18030")
            
            if(i < top):
                # MFCCs
                # Wlt
                # Beats
                # in total
                print (str(res[i][frequencyPara_dimension \
                   + beats_dimension + 4])).encode("UTF-8") + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 1])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 2])).encode("UTF-8") + ")"+ \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 3])).encode("UTF-8") + ")"+ \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 4])).encode("UTF-8") + ")"
                '''
                freqs = np.linspace(0, 120, CLength * 12)
                pl.figure(figsize=(8,4))
                pl.subplot(111)
                pl.plot(freqs, target\
                    [frequencyPara_dimension - CLength * 12:frequencyPara_dimension],
                    alpha = 0.7)
                pl.plot(freqs, res[i]\
                    [frequencyPara_dimension - CLength * 12:frequencyPara_dimension],
                    alpha = 0.7)
                pl.xlabel(u"相似度")
                pl.show()
                '''
        '''
        errorCount = 0
        for i in xrange(top):
            if "@" in res[i][0]:
                if("acoustic" not in res[i][0].strip().split("@")[1].strip().split("#")):
                    errorCount += 1
        accuracy = float(1) - float(errorCount) / float(top)
        print "error: %d" % errorCount
        print "accuracy: %f" % accuracy
        return accuracy
        '''
            
    if word == "WltMFCCs":
        for eachLine in xrange(lineCount):
            simi = 0
            for i in xrange(CLength * 2):
                diff = math.pow((float(dataSet[eachLine][frequencyPara_dimension - CLength * 6 + i + 1]) \
                                 - float(target[frequencyPara_dimension - CLength * 6 + i + 1])), 2)
            simi += diff
            simi = math.sqrt(simi) * 1
            dataSet[eachLine].extend([simi])
            res.append(dataSet[eachLine])

        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 1])
        del dataSet

        dataSetForMid = []
        top = int(TOTAL_MUSIC * 0.3)
        for i in xrange(top):
            '''
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension
                              + beats_dimension + 1])).encode("UTF-8")
            '''
            dataSetForMid.append(res[i])
        res = None
        del res

        res = []
        for eachLine in xrange(len(dataSetForMid)):
            simi = 0
            for i in xrange(CLength * 2):
                diff = math.pow((float(dataSetForMid[eachLine][frequencyPara_dimension - CLength * 4 + i + 1]) \
                                 - float(target[frequencyPara_dimension - CLength * 4 + i + 1])), 2)
                simi += diff
            simi = math.sqrt(simi)
            dataSetForMid[eachLine].extend([simi])
            res.append(dataSetForMid[eachLine])
            
        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 2])
        del dataSetForMid

        dataSetForHigh = []
        top = int(TOTAL_MUSIC * 0.2)
        for i in xrange(top):
            '''
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension
                              + beats_dimension + 1])).encode("UTF-8")
            '''
            dataSetForHigh.append(res[i])
        res = None
        del res

        res = []
        for eachLine in xrange(len(dataSetForHigh)):
            simi = 0
            for i in xrange(CLength * 2):
                diff = math.pow((float(dataSetForHigh[eachLine][frequencyPara_dimension - CLength * 2 + i + 1]) \
                                 - float(target[frequencyPara_dimension - CLength * 2 + i + 1])), 2)
                simi += diff
            simi = math.sqrt(simi)
            dataSetForHigh[eachLine].extend([simi])
            res.append(dataSetForHigh[eachLine])
            
        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 3])
        del dataSetForHigh

        dataSetForMFCCs = []
        top = int(TOTAL_MUSIC * 0.1)
        for i in xrange(top):
            '''
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension
                              + beats_dimension + 1])).encode("UTF-8")
            '''
            dataSetForMFCCs.append(res[i])
        res = None
        del res

        res = []
        for eachLine in xrange(len(dataSetForMFCCs)):
            simi = 0
            for i in xrange(frequencyPara_dimension - CLength * 12):
                diff = math.pow((float(dataSetForMFCCs[eachLine][i + 1]) \
                                 - float(target[i + 1])), 2)
                simi += diff
            simi = math.sqrt(simi)
            dataSetForMFCCs[eachLine].extend([simi])
            res.append(dataSetForMFCCs[eachLine])
            
        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 4])
        del dataSetForMFCCs

        dataSetForBeats = []
        top = int(TOTAL_MUSIC * 0.01)
        for i in xrange(top):
            '''
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension \
                   + beats_dimension + 3])).encode("UTF-8") + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 1])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 2])).encode("UTF-8") + ")"
            '''
            dataSetForBeats.append(res[i])
            
        res = []
        AVE_para = 150
        lowReduce_para = 8
        highReduce_para = 4
        lowBeat_para = 1
        highBeat_para = 1
        lowBeatsInt_para = 50
        highBeatsInt_para = 50
        for eachLine in xrange(len(dataSetForBeats)):
            simi = 0
            '''
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 1]) / AVE_para - float(target[frequencyPara_dimension + 1]) / AVE_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 1]) / AVE_para - float(target[frequencyPara_dimension + 1]) / AVE_para)
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 3]) * lowReduce_para - float(target[frequencyPara_dimension + 3]) * lowReduce_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 3]) * lowReduce_para - float(target[frequencyPara_dimension + 3]) * lowReduce_para)
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 4]) * highReduce_para - float(target[frequencyPara_dimension + 4]) * highReduce_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 4]) * highReduce_para - float(target[frequencyPara_dimension + 4]) * highReduce_para)
            '''
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 5]) * lowBeat_para - float(target[frequencyPara_dimension + 5])) * lowBeat_para * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 5]) * lowBeat_para - float(target[frequencyPara_dimension + 5])) * lowBeat_para
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 6]) * highBeat_para - float(target[frequencyPara_dimension + 6])) * highBeat_para * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 6]) * highBeat_para - float(target[frequencyPara_dimension + 6])) * highBeat_para
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target[frequencyPara_dimension + 7]) / lowBeatsInt_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target[frequencyPara_dimension + 7]) / lowBeatsInt_para)
            simi = simi + (float(dataSetForBeats[eachLine][frequencyPara_dimension + 8]) / highBeatsInt_para - float(target[frequencyPara_dimension + 8]) / highBeatsInt_para) * \
                          (float(dataSetForBeats[eachLine][frequencyPara_dimension + 8]) / highBeatsInt_para - float(target[frequencyPara_dimension + 8]) / highBeatsInt_para)

            simi = math.sqrt(simi)
            dataSetForBeats[eachLine].extend([simi])
            dataSetForBeats[eachLine].extend([(simi + \
                            dataSetForBeats[eachLine]\
                            [frequencyPara_dimension \
                             + beats_dimension + 1] + \
                            dataSetForBeats[eachLine]\
                            [frequencyPara_dimension \
                             + beats_dimension + 2] + \
                            dataSetForBeats[eachLine]\
                            [frequencyPara_dimension \
                             + beats_dimension + 3] + \
                            dataSetForBeats[eachLine]\
                            [frequencyPara_dimension \
                             + beats_dimension + 4]) / 5])
            res.append(dataSetForBeats[eachLine])
            
        res.sort(key = lambda x:x[frequencyPara_dimension\
                                  + beats_dimension + 6])
        dataSetForBeats = None
        del dataSetForBeats

        top = 5
        for i in xrange(top):
            print "----------  NO.%d  ----------" % (i + 1)
            print (str(res[i][0])).encode("gb18030")
            print (str(res[i][frequencyPara_dimension \
                   + beats_dimension + 6])).encode("UTF-8") + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 1])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 2])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 3])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 4])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 5])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][frequencyPara_dimension \
                                    + beats_dimension + 6])).encode("UTF-8") + ")"


            freqs = np.linspace(0, 60, CLength * 6)
            pl.figure(figsize=(8,4))
            pl.subplot(111)
            pl.plot(freqs, res[eachLine]\
                    [frequencyPara_dimension - CLength * 6:frequencyPara_dimension],
                    alpha = 1)
            pl.plot(freqs, target\
                    [frequencyPara_dimension - CLength * 6:frequencyPara_dimension],
                    alpha = 0.7)
            pl.xlabel(u"相似度")
            pl.show()
    
def compare(fileName1, fileName2, field = ""):
    global MFCCs_dimension
    
    path1, fileName1 = getPath(fileName1)
    path2, fileName2 = getPath(fileName2)
    dataDimention =  1 + frequencyPara_dimension \
                    + beats_dimension
    
    fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
    dataSet = []
    lineCount = -1
    for line in fr:
        line = line.strip().split("\t")
        dataSet.append(line)
        lineCount += 1

    target1 = []
    target2 = []
    for eachLine in xrange(lineCount):
        if (str)(dataSet[eachLine][0]) == (str)(fileName1):
            for i in xrange(dataDimention):
                target1.append(dataSet[eachLine][i])
            break
    for eachLine in xrange(lineCount):
        if (str)(dataSet[eachLine][0]) == (str)(fileName2):
            for i in xrange(dataDimention):
                target2.append(dataSet[eachLine][i])
            break

    #print target1
    #print target2

    if field == "MFCCs":
        simi = 0
        for i in xrange(frequencyPara_dimension):
            diff = math.pow((float(target1[i + 1]) \
                            - float(target2[i + 1])), 2)
            print str(i + 1) + "  " + str(diff)
            simi += diff
        
    if field == "Beats":
        AVE_para = 1500
        lowReduce_para = 8
        highReduce_para = 4
        lowBeat_para = 1
        highBeat_para = 1
        lowBeatsInt_para = 1000
        highBeatsInt_para = 2000
        simi = 0
        
        simi = simi + (float(target1[frequencyPara_dimension + 1]) / AVE_para - float(target2[frequencyPara_dimension + 1]) / AVE_para) * \
                      (float(target1[frequencyPara_dimension + 1]) / AVE_para - float(target2[frequencyPara_dimension + 1]) / AVE_para)
        
        simi = simi + (float(target1[frequencyPara_dimension + 3]) * lowReduce_para - float(target2[frequencyPara_dimension + 3]) * lowReduce_para) * \
                      (float(target1[frequencyPara_dimension + 3]) * lowReduce_para - float(target2[frequencyPara_dimension + 3]) * lowReduce_para)
        simi = simi + (float(target1[frequencyPara_dimension + 4]) * highReduce_para - float(target2[frequencyPara_dimension + 4]) * highReduce_para) * \
                      (float(target1[frequencyPara_dimension + 4]) * highReduce_para - float(target2[frequencyPara_dimension + 4]) * highReduce_para)
        
        simi = simi + (float(target1[frequencyPara_dimension + 5]) * lowBeat_para - float(target2[frequencyPara_dimension + 5]) * lowBeat_para) * \
                      (float(target1[frequencyPara_dimension + 5]) * lowBeat_para - float(target2[frequencyPara_dimension + 5]) * lowBeat_para)
        simi = simi + (float(target1[frequencyPara_dimension + 6]) * highBeat_para - float(target2[frequencyPara_dimension + 6]) * highBeat_para) * \
                      (float(target1[frequencyPara_dimension + 6]) * highBeat_para - float(target2[frequencyPara_dimension + 6]) * highBeat_para)
        
        simi = simi + (float(target1[frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target2[frequencyPara_dimension + 7]) / lowBeatsInt_para) * \
                      (float(target1[frequencyPara_dimension + 7]) / lowBeatsInt_para - float(target2[frequencyPara_dimension + 7]) / lowBeatsInt_para)
        simi = simi + (float(target1[frequencyPara_dimension + 8]) / highBeatsInt_para - float(target2[frequencyPara_dimension + 8]) / highBeatsInt_para) * \
                      (float(target1[frequencyPara_dimension + 8]) / highBeatsInt_para - float(target2[frequencyPara_dimension + 8]) / highBeatsInt_para)

    simi = math.sqrt(simi)

    print simi

def serializeFileName():
    global PROJECT_DIR
    global PROJECT_INFO_DIR
    global INFO_TXT_NAME

    path = PROJECT_INFO_DIR + INFO_TXT_NAME

    infoData = []

    fr = None
    fw = None

    if os.path.exists(path):
        fileList = os.listdir(PROJECT_DIR[:len(PROJECT_DIR) - 1])
        fr = open(path, "r")
        infoData = fr.read()
        fr.close()
        infoData = infoData.strip().split("\n")
        fileNum = infoData[0]
        fileName = infoData[1:]
        AudioFileCount = 0
        dataToWrite = []
        for audioFile in fileList:
            if (audioFile.strip().split(".")[-1] == "wav") or \
               (audioFile.strip().split(".")[-1] == "WAV") or \
               (audioFile.strip().split(".")[-1] == "mp3") or \
               (audioFile.strip().split(".")[-1] == "MP3"):
                AudioFileCount += 1
        dataToWrite.extend(fileName)
        for audioFile in fileList:
            if (audioFile.strip().split(".")[-1] == "wav") or \
               (audioFile.strip().split(".")[-1] == "WAV") or \
               (audioFile.strip().split(".")[-1] == "mp3") or \
               (audioFile.strip().split(".")[-1] == "MP3"):
                if ("@" not in audioFile):
                    srcFile = os.path.join(PROJECT_DIR, audioFile)
                    destFile = os.path.join(PROJECT_DIR, \
                                            str(100000 + int(fileNum)) + "@" + audioFile)
                    os.rename(srcFile, destFile)
                    dataToWrite.append(destFile.strip().split("\\")[-1])           
                    fileNum = int(fileNum) + 1
        dataToWrite.append(fileNum)
        os.remove(path)
        print "same inof.txt file detected, delete"
        fw = open(path, "w")
        fw.write(str(dataToWrite[-1]) + "\n")
        for i in xrange(len(dataToWrite) - 1):
            fw.write(dataToWrite[i] + "\n")
    else:
        fileList = os.listdir(PROJECT_DIR[:len(PROJECT_DIR) - 1])
        AudioFileCount = 0
        for audioFile in fileList:
            if (audioFile.strip().split(".")[-1] == "wav") or \
               (audioFile.strip().split(".")[-1] == "WAV") or \
               (audioFile.strip().split(".")[-1] == "mp3") or \
               (audioFile.strip().split(".")[-1] == "MP3"):
                AudioFileCount += 1
        fw = open(path, "w")
        fw.write(str(AudioFileCount) + "\n")
        count = 0
        for audioFile in fileList:
            if (audioFile.strip().split(".")[-1] == "wav") or \
               (audioFile.strip().split(".")[-1] == "WAV") or \
               (audioFile.strip().split(".")[-1] == "mp3") or \
               (audioFile.strip().split(".")[-1] == "MP3"):
                if (audioFile.count("@") == 1):
                    srcFile = os.path.join(PROJECT_DIR, audioFile)
                    destFile = os.path.join(PROJECT_DIR, \
                                        str(100000 + count) + "@" + audioFile)
                    os.rename(srcFile, destFile)
                count += 1
                fw.write(destFile.strip().split("\\")[-1] + "\n")
        
    if fw != None:
        fw.close()

def resetFileName():
    global PROJECT_DIR

    path = PROJECT_INFO_DIR + INFO_TXT_NAME
    fileList = os.listdir(PROJECT_DIR[:len(PROJECT_DIR) - 1])
    if os.path.exists(path):
        os.remove(path)
        print "inof.txt file delete"
    for fileName in fileList:
        if "@" in fileName:
            srcFile = os.path.join(PROJECT_DIR, fileName)
            #destFile = os.path.join(PROJECT_DIR, fileName.strip().split("@")[-1])
            destFile = os.path.join(PROJECT_DIR, \
                fileName.strip().split("@")[-2] + "@" + fileName.strip().split("@")[-1])
            os.rename(srcFile, destFile)

def test():
    global PROJECT_DIR
    fileList = os.listdir(PROJECT_DIR[:len(PROJECT_DIR) - 1])
    print fileList[:5]

def readXml(tag):
    import xml.dom.minidom
    global PROJECT_DIR

    if(tag == "genre"):
        xmlList = []
        fileName = []
        genre = []
    
        fileList = os.listdir(PROJECT_DIR[ : len(PROJECT_DIR) - 1])
        print len(fileList)

        # get names of .xml file
        for f in fileList:
            if (f.strip().split(".")[-1] == "xml" and f.strip().split(".")[-2] != \
                "user_taxonomies"):
                xmlList.append(f)

        for i in xrange(len(xmlList)):
            print xmlList[i]
            dom = xml.dom.minidom.parse(PROJECT_DIR + xmlList[i])
            root = dom.documentElement
            itemList = root.getElementsByTagName("song")

            for item in itemList:
                srcName = item.getAttribute("path")
                desName = item.getAttribute("genre") + "@" + srcName
                if(srcName in fileList):
                    os.rename(os.path.join(PROJECT_DIR, srcName), \
                          os.path.join(PROJECT_DIR, desName))

    if(tag == "instruments"):
        import shutil
        xmlList = []
        fileName = []
        nameList = []
        instruments = "Rap"
    
        fileList = os.listdir(PROJECT_DIR[ : len(PROJECT_DIR) - 1])
        print len(fileList)

        # get names of .xml file
        for f in fileList:
            if (f.strip().split(".")[-1] == "xml" and f.strip().split(".")[-2] == \
                "user_taxonomies"):
                xmlList.append(f)

        print xmlList[0]
        dom = xml.dom.minidom.parse(PROJECT_DIR + xmlList[0])
        root = dom.documentElement
        #clusterList = root.getElementsByTagName("cluster")
        propertyTag = root.getElementsByTagName("property")
        
        for pTag in propertyTag:
            if pTag.getAttribute("title") == "Alternative / Rock":
                songs = pTag.getElementsByTagName("object")
                for eachSong in songs:
                        if eachSong.getAttribute("id") not in nameList:
                            nameList.append(eachSong.getAttribute("id"))    

        '''
        for cluster in clusterList:
            if cluster.getAttribute("title") == instruments:
                if(False):
                    detaiList = cluster.getElementsByTagName("cluster")
                    for detail in detaiList:
                        if detail.getAttribute("title") == "acoustic" or \
                           detail.getAttribute("title") != "acoustic guitar":
                            songs = detail.getElementsByTagName("object")
                            for eachSong in songs:
                                if eachSong.getAttribute("id") not in nameList:
                                    nameList.append(eachSong.getAttribute("id"))
                else:
                    songs = cluster.getElementsByTagName("object")
                    for eachSong in songs:
                        if eachSong.getAttribute("id") not in nameList:
                            nameList.append(eachSong.getAttribute("id"))
        '''
        
        for i in xrange(len(nameList)):
            nameList[i] = nameList[i][:-4] + ".mp3"
            print nameList[i]
            
            for eachFile in fileList:
                if (eachFile.replace("-", "_") == nameList[i]):
                    srcDir = PROJECT_DIR + eachFile
                    desDir = PROJECT_DIR + "instruments\\" + "rock" + "@" + eachFile
            if(os.path.exists(desDir) == False):
                shutil.copy(srcDir, desDir)
            
def perfermance(start, end):
    maxAcc = -1
    name = []
    aveAcc = -1
    sumAcc = 0
    
    fileList = os.listdir(PROJECT_DIR[ : len(PROJECT_DIR) - 1])
    tempRemove = []
    for i in xrange(len(fileList)):
        if(fileList[i].strip().split(".")[-1] != "mp3"):
            tempRemove.append(fileList[i])
    print len(tempRemove)
    for i in xrange(len(tempRemove)):
        fileList.remove(tempRemove[i])
    print len(fileList)
    for i in xrange(end - start):
        acc = similarity(fileList[start + i], "MFCCsWlt")
        sumAcc = sumAcc + acc
        if(acc > maxAcc):
            maxAcc = acc
            name = fileList[start + i]
    aveAcc = sumAcc / (end - start)
    return maxAcc, name, aveAcc

def drawSample(tag):
    fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
    if(tag == "instruments"):
        global MFCCPara_dimension
        global wltPara_dimension
        global beats_dimension
        global frequencyPara_dimension
        maxD = 0
        acoustic = np.zeros(0)
        piano = np.zeros(0)
        hiphop = np.zeros(0)
        relaxing = np.zeros(0)
        vocal = np.zeros(0)
        rock = np.zeros(0)
        tmp = np.zeros(0)
        getDim = False
        dim = []
        freqVec = np.zeros(0)
        wltVec = np.zeros(0)
        beatVec = np.zeros(0)
        writeVec = []
        fw = open(PROJECT_DOC_DIR + u"Vec.txt", "a")
        for line in fr:
            line = line.strip().split("\t")
            if getDim == False:
                #fVec = open(PROJECT_DOC_DIR + u"Vec.txt", "r")
                #vec = []
                #for vecLine in fVec:
                #    vec.append(vecLine)
                dim.append(len(line[1:frequencyPara_dimension - wltPara_dimension]))
#                dim.append(len(line[frequencyPara_dimension - wltPara_dimension : frequencyPara_dimension]))
                dim.append(40)
                freqVec = np.append(freqVec, np.random.randint(1, 100, dim[0]))
                wltVec = np.append(wltVec, np.random.randint(1, 100, dim[1]))
                beatVec = np.append(beatVec, np.random.randint(1, 100, 8))

                '''
                vec[0] = vec[0].strip().split("\t")
                for i in xrange(len(vec[0])):
                    freqVec = np.append(freqVec, np.array(vec[0][i]))
                print len(freqVec)
                vec[1] = vec[1].strip().split("\t")
                for i in xrange(len(vec[1])):
                    wltVec = np.append(wltVec, np.array(vec[1][i]))
                '''
                print "get dimension of frequency and wlt"
                getDim = True
            if "acoustic" in line[0]:
                tmpLine = []
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:MFCCPara_dimension + 1]))
                tmpLine.extend(line[MFCCPara_dimension + 1:MFCCPara_dimension + 21])
                tmpLine.extend(line[MFCCPara_dimension + 442:MFCCPara_dimension + 462])
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, np.array(tmpLine, np.float)))
                tmp = np.append(tmp, getCosineDis( \
                    beatVec, line[frequencyPara_dimension + 1 : ]))
                acoustic = np.append(acoustic, tmp)
                tmp = np.zeros(0)
            if "piano" in line[0]:
                tmpLine = []
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:MFCCPara_dimension + 1]))
                tmpLine.extend(line[MFCCPara_dimension + 1:MFCCPara_dimension + 21])
                tmpLine.extend(line[MFCCPara_dimension + 442:MFCCPara_dimension + 462])
                tmp = np.append(tmp, getCosineDis( \
                   wltVec, np.array(tmpLine, np.float)))
                tmp = np.append(tmp, getCosineDis( \
                    beatVec, line[frequencyPara_dimension + 1 : ]))
                piano = np.append(piano, tmp)
                tmp = np.zeros(0)
            if "hiphop" in line[0]:
                tmpLine = []
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:MFCCPara_dimension + 1]))
                tmpLine.extend(line[MFCCPara_dimension + 1:MFCCPara_dimension + 21])
                tmpLine.extend(line[MFCCPara_dimension + 442:MFCCPara_dimension + 462])
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, np.array(tmpLine, np.float)))
                tmp = np.append(tmp, getCosineDis( \
                    beatVec, line[frequencyPara_dimension + 1 : ]))
                hiphop = np.append(hiphop, tmp)
                tmp = np.zeros(0)
            if "relaxing" in line[0]:
                tmpLine = []
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:MFCCPara_dimension + 1]))
                tmpLine.extend(line[MFCCPara_dimension + 1:MFCCPara_dimension + 21])
                tmpLine.extend(line[MFCCPara_dimension + 442:MFCCPara_dimension + 462])
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, np.array(tmpLine, np.float)))
                tmp = np.append(tmp, getCosineDis( \
                    beatVec, line[frequencyPara_dimension + 1 : ]))
                relaxing = np.append(relaxing, tmp)
                tmp = np.zeros(0)
            if "vocal" in line[0]:
                tmpLine = []
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:MFCCPara_dimension + 1]))
                tmpLine.extend(line[MFCCPara_dimension + 1:MFCCPara_dimension + 21])
                tmpLine.extend(line[MFCCPara_dimension + 442:MFCCPara_dimension + 462])
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, np.array(tmpLine, np.float)))
                tmp = np.append(tmp, getCosineDis( \
                    beatVec, line[frequencyPara_dimension + 1 : ]))
                vocal = np.append(vocal, tmp)
                tmp = np.zeros(0)
            if "rock" in line[0]:
                tmpLine = []
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:MFCCPara_dimension + 1]))
                tmpLine.extend(line[MFCCPara_dimension + 1:MFCCPara_dimension + 21])
                tmpLine.extend(line[MFCCPara_dimension + 442:MFCCPara_dimension + 462])
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, np.array(tmpLine, np.float)))
                tmp = np.append(tmp, getCosineDis( \
                    beatVec, line[frequencyPara_dimension + 1 : ]))
                rock = np.append(rock, tmp)
                tmp = np.zeros(0)

        acoustic.shape = -1, 3
        acoustic = acoustic.T
        piano.shape = -1, 3
        piano = piano.T
        hiphop.shape = -1, 3
        hiphop = hiphop.T
        relaxing.shape = -1, 3
        relaxing = relaxing.T
        vocal = vocal.T
        vocal.shape = -1, 3
        rock.shape = -1, 3
        rock = rock.T

        draw3D = False
        mark = ['s','o','^','v','>','<','d','p','h','8','+','*']

        if draw3D == False:
            plt.scatter(acoustic[0], acoustic[1], marker = mark[2], color = np.array([[0,0,1,0.5]]), label = str("acoustic"))
            plt.scatter(acoustic[0].mean(), acoustic[1].mean(), s = 100, marker = mark[2], color = np.array([[0,0,1,1]]))            
            plt.scatter(piano[0], piano[1],  marker = mark[1], color = np.array([[0,1,0, 0.5]]), label = str("piano"))
            plt.scatter(piano[0].mean(), piano[1].mean(), s = 100, marker = mark[1], color = np.array([[0,1,0, 1]]))
            plt.scatter(hiphop[0], hiphop[1],  marker = mark[6], color = np.array([[1,0,0, 0.5]]), label = str("hiphop"))
            plt.scatter(hiphop[0].mean(), hiphop[1].mean(), s = 100, marker = mark[6], color = np.array([[1,0,0, 1]]))
            plt.scatter(relaxing[0], relaxing[1],  marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.5]]), label = str("relaxing"))
            plt.scatter(relaxing[0].mean(), relaxing[1].mean(), s = 100, marker = mark[3], color = np.array([[0.7,0.3,0.2, 1]]))
#            plt.scatter(vocal[0], vocal[1], marker = mark[10], color = np.array([[0.27,0.37865,0.889, 1]]), label = str("vocal"))
#            plt.scatter(vocal[0].mean(), vocal[1].mean(), s = 100, marker = mark[10], color = np.array([[0.27,0.37865,0.889, 1]]))
            plt.scatter(rock[0], rock[1],  marker = mark[11], color = np.array([[0,0,0, 0.5]]), label = str("rock"))
            plt.scatter(rock[0].mean(), rock[1].mean(), s = 100, marker = mark[11], color = np.array([[0,0,0, 1]]))           

            pl.xlabel("Frequency Feature", fontsize = 15)
#            pl.xlabel("Auditory Perceptual Feature", fontsize = 15)
            pl.ylabel("Auditory Perceptual Feature", fontsize = 15)
#            pl.ylabel("Statistical Characteristic of Beat", fontsize = 15)
            plt.legend()
            plt.show()

        if draw3D == True:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = "3d")
            ax.scatter(acoustic[0], acoustic[1], acoustic[2],  marker = mark[2], color = np.array([[0,0,1,0.5]]), label = str("acoustic"))
            ax.scatter(piano[0], piano[1], piano[2],  marker = mark[1], color = np.array([[0,1,0, 0.5]]), label = str("piano"))
            ax.scatter(hiphop[0], hiphop[1], hiphop[2],  marker = mark[6], color = np.array([[1,0,0, 0.5]]), label = str("hiphop"))
            ax.scatter(relaxing[0], relaxing[1], relaxing[2],  marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.5]]), label = str("relaxing"))
            ax.scatter(vocal[0], vocal[1], vocal[2],  marker = mark[10], color = np.array([[0.27,0.37865,0.889, 1]]), label = str("vocal"))
            ax.scatter(rock[0], rock[1], rock[2],  marker = mark[11], color = np.array([[0,0,0, 0.5]]), label = str("rock"))

            ax.set_xlabel("Frequency Feature", fontsize = 15)
            ax.set_ylabel("Auditory Perceptual Feature", fontsize = 15)
            ax.set_zlabel("Statistical Characteristic of Beat", fontsize = 15)
            plt.show()

    if(tag == "all"):
        fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
        alternative = np.zeros(0)
        blues = np.zeros(0)
        electronic = np.zeros(0)
        folkcountry = np.zeros(0)
        funksoulrnb = np.zeros(0)
        jazz = np.zeros(0)
        pop = np.zeros(0)
        raphiphop = np.zeros(0)
        rock = np.zeros(0)
        tmp = np.zeros(0)

        getDim = False
        dim = []
        freqVec = np.zeros(0)
        wltVec = np.zeros(0)
        
        for line in fr:
            line = line.strip().split("\t")
            if getDim == False:
                dim.append(len(line[1:frequencyPara_dimension - wltPara_dimension]))
                dim.append(len(line[frequencyPara_dimension - wltPara_dimension : frequencyPara_dimension]))
                freqVec = np.append(freqVec, np.random.randint(1, 100, dim[0]))
                wltVec = np.append(wltVec, np.random.randint(1, 100, dim[1]))
                print "get dimension of frequency and wlt"
                getDim = True
            if "alternative" in line[0]:
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:frequencyPara_dimension - CLength * 12]))
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, line[frequencyPara_dimension - CLength * 12:frequencyPara_dimension]))
                alternative = np.append(alternative, tmp)
                tmp = np.zeros(0)
            if "blues" in line[0]:
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:frequencyPara_dimension - CLength * 12]))
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, line[frequencyPara_dimension - CLength * 12:frequencyPara_dimension]))
                blues = np.append(blues, tmp)
                tmp = np.zeros(0)
            if "electronic" in line[0]:
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:frequencyPara_dimension - CLength * 12]))
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, line[frequencyPara_dimension - CLength * 12:frequencyPara_dimension]))
                electronic = np.append(electronic, tmp)
                tmp = np.zeros(0)
            if "folkcountry" in line[0]:
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:frequencyPara_dimension - CLength * 12]))
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, line[frequencyPara_dimension - CLength * 12:frequencyPara_dimension]))
                folkcountry = np.append(folkcountry, tmp)
                tmp = np.zeros(0)
            if "funksoulrnb" in line[0]:
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:frequencyPara_dimension - CLength * 12]))
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, line[frequencyPara_dimension - CLength * 12:frequencyPara_dimension]))
                funksoulrnb = np.append(funksoulrnb, tmp)
                tmp = np.zeros(0)
            if "jazz" in line[0]:
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:frequencyPara_dimension - CLength * 12]))
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, line[frequencyPara_dimension - CLength * 12:frequencyPara_dimension]))
                jazz = np.append(jazz, tmp)
                tmp = np.zeros(0)
            if "pop" in line[0]:
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:frequencyPara_dimension - CLength * 12]))
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, line[frequencyPara_dimension - CLength * 12:frequencyPara_dimension]))
                pop = np.append(pop, tmp)
                tmp = np.zeros(0)
            if "raphiphop" in line[0]:
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:frequencyPara_dimension - CLength * 12]))
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, line[frequencyPara_dimension - CLength * 12:frequencyPara_dimension]))
                raphiphop = np.append(raphiphop, tmp)
                tmp = np.zeros(0)
            if "rock" in line[0]:
                tmp = np.append(tmp, getCosineDis( \
                    freqVec, line[1:frequencyPara_dimension - CLength * 12]))
                tmp = np.append(tmp, getCosineDis( \
                    wltVec, line[frequencyPara_dimension - CLength * 12:frequencyPara_dimension]))
                rock = np.append(rock, tmp)
                tmp = np.zeros(0)


        alternative.shape = -1, 2
        alternative = alternative.T
        blues.shape = -1, 2
        blues = blues.T
        electronic.shape = -1, 2
        electronic = electronic.T
        folkcountry.shape = -1, 2
        folkcountry = folkcountry.T
        funksoulrnb.shape = -1, 2
        funksoulrnb = funksoulrnb.T
        jazz.shape = -1, 2
        jazz = jazz.T
        pop.shape = -1, 2
        pop = pop.T
        raphiphop.shape = -1, 2
        raphiphop = raphiphop.T
        rock.shape = -1, 2
        rock = rock.T

        mark = ['s','o','^','v','>','<','d','p','h','8','+','*']
        
        plt.scatter(alternative[0], alternative[1], marker = mark[0], color = np.array([[0,0,1,0.5]]), label = str("alternative"))
        plt.scatter(blues[0], blues[1], marker = mark[1], color = np.array([[0,1,0, 0.5]]), label = str("blues"))
        plt.scatter(electronic[0], electronic[1], marker = mark[6], color = np.array([[1,0,0, 0.5]]), label = str("electronic"))
        plt.scatter(folkcountry[0], folkcountry[1], marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.5]]), label = str("folkcountry"))
        plt.scatter(funksoulrnb[0], funksoulrnb[1], marker = mark[11], color = np.array([[0.3,0.1,0.9, 0.5]]), label = str("funksoulrnb"))
        plt.scatter(jazz[0], jazz[1], marker = mark[11], color = np.array([[0.65,0.145,0.8, 0.5]]), label = str("jazz"))
        plt.scatter(pop[0], pop[1], marker = mark[11], color = np.array([[0.2,0.5789,0.67, 0.5]]), label = str("pop"))
        plt.scatter(raphiphop[0], raphiphop[1], marker = mark[11], color = np.array([[0.57,0.2,0.3, 0.5]]), label = str("raphiphop"))
        plt.scatter(rock[0], rock[1], marker = mark[11], color = np.array([[0,0,0, 0.5]]), label = str("rock"))
        

#        plt.scatter(linspace(0,1,len(alternative[1])), alternative[0], marker = mark[0], color = np.array([[0,0,1,0.5]]), label = str("alternative"))
#        plt.scatter(linspace(0,1,len(blues[1])), blues[0], marker = mark[1], color = np.array([[0,1,0, 0.5]]), label = str("blues"))
#        plt.scatter(linspace(0,1,len(electronic[1])), electronic[0], marker = mark[6], color = np.array([[1,0,0, 0.5]]), label = str("electronic"))
#        plt.scatter(linspace(0,1,len(folkcountry[1])), folkcountry[0], marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.5]]), label = str("folkcountry"))
#        plt.scatter(linspace(0,1,len(funksoulrnb[1])), funksoulrnb[0], marker = mark[11], color = np.array([[0.3,0.1,0.9, 0.5]]), label = str("funksoulrnb"))
#        plt.scatter(linspace(0,1,len(jazz[1])), jazz[0], marker = mark[11], color = np.array([[0.65,0.145,0.8, 0.5]]), label = str("jazz"))
#        plt.scatter(linspace(0,1,len(pop[1])), pop[0], marker = mark[11], color = np.array([[0.2,0.5789,0.67, 0.5]]), label = str("pop"))
#        plt.scatter(linspace(0,1,len(raphiphop[1])), raphiphop[0], marker = mark[11], color = np.array([[0.57,0.2,0.3, 0.5]]), label = str("raphiphop"))
#        plt.scatter(linspace(0,1,len(rock[1])), rock[0], marker = mark[11], color = np.array([[0,0,0, 0.5]]), label = str("rock"))
        plt.legend()
        plt.show()
    fr.close()

def getCosineDis(srcVec, disVec):
    cosAlpha = 2
    up = 0
    down = 0
    sqSumSrc = 0
    sqSumDis = 0
    for i in xrange(len(srcVec)):
        up = up + float(srcVec[i]) * float(disVec[i])
        sqSumSrc = float(sqSumSrc) + float(srcVec[i]) ** 2
        sqSumDis = float(sqSumDis) + float(disVec[i]) ** 2
    down = sqrt(sqSumSrc) * sqrt(sqSumDis)
    cosAlpha = up / down
    return cosAlpha

def drawASample():
    fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
    mark = ['s','o','^','v','>','<','d','p','h','8','+','*']
    vector = np.ones(1)
    vector = np.append(vector, np.zeros(12798))

    for i in xrange(0):
        fr.readline()
    line = fr.readline().strip().split("\t")[1:frequencyPara_dimension - CLength * 12]
    plt.scatter(1, getCosineDis(vector, line), marker = mark[0], color = np.array([[0,0,1,0.5]]), label = str("1"))
    for i in xrange(145):
        fr.readline()
    line = fr.readline().strip().split("\t")[1:]
    plt.scatter(1, getCosineDis(vector, line), marker = mark[1], color = np.array([[0,1,0,0.5]]), label = str("2"))
    for i in xrange(122):
        fr.readline()
    line = fr.readline().strip().split("\t")[1:]
    plt.scatter(1, getCosineDis(vector, line), marker = mark[2], color = np.array([[1,0,0,0.5]]), label = str("3"))
    for i in xrange(120):
        fr.readline()
    line = fr.readline().strip().split("\t")[1:]
    plt.scatter(1, getCosineDis(vector, line), marker = mark[3], color = np.array([[0.46,0.267,0.67,0.5]]), label = str("4"))
    for i in xrange(590):
        fr.readline()
    line = fr.readline().strip().split("\t")[1:]
    plt.scatter(1, getCosineDis(vector, line), marker = mark[4], color = np.array([[0.7,0.5,0.1,0.5]]), label = str("5"))
    
    plt.legend()
    plt.show()

def PCA(dataMat, topNFeat):
    meanVals = mean(dataMat, axis = 0) #axis = 0 mean(list) axis = 1 mean(line)
    meanRemoved = dataMat - meanVals
    covMat = np.cov(meanRemoved, rowvar = 0)
    print "covMat.shape: ", covMat.shape
    eigVals, eigVects = np.linalg.eig(mat(covMat))
    del covMat
    eigValIndex = argsort(eigVals)
    eigValIndex = eigValIndex[:-(topNFeat + 1):-1]
    topEigVects = eigVects[:,eigValIndex]
    resData = meanRemoved * topEigVects
    sumEig = 0
    for i in xrange(len(eigVals)):
        sumEig += eigVals[i]
    partition = 0
    printForOnce = False
    for i in xrange(len(eigVals)):
        partition += eigVals[i]
        ratio = float(partition) / float(sumEig)
        if ratio < 0.95:
            print "%d, %f" % (i, ratio)
        else:
            if printForOnce == False:
                print "%d, %f" % (i, ratio)
                printForOnce = True
    return resData, eigVals

def writePCAData(tag):
    global MFCC_TOP_FEAT
    global WLT_TOP_FEAT
    global BEAT_TOP_FEAT
    global MFCCPara_dimension
    global wltPara_dimension
    global beats_dimension
    global frequencyPara_dimension
    
    MFCCDataMat = np.zeros(0)
    wltDataMat = np.zeros(0)
    beatDataMat = np.zeros(0)
    fileName = []

    if tag == "separately": # PCA for each dimension
        print "MFCC_TOP_FEAT: ", MFCC_TOP_FEAT
        print "WLT_TOP_FEAT: ", WLT_TOP_FEAT
        print "BEAT_TOP_FEAT", BEAT_TOP_FEAT
        fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
        for line in fr:
            line = line.strip().split("\t")
            fileName.append(line[0])
            MFCCTmp = []
            wltTmp = []
            beatTmp = []
            MFCCTmp.append(line[1 : MFCCPara_dimension + 1])
            wltTmp.append(line[MFCCPara_dimension + 1 : MFCCPara_dimension + \
                               wltPara_dimension + 1])
            for i in xrange(beats_dimension):
                if line[i + 1 + frequencyPara_dimension] == "nan":
                    line[i + 1 + frequencyPara_dimension] = "0"
                    print "fix a nan"
                beatTmp.append(line[i + 1 + frequencyPara_dimension])
            MFCCTmp = np.array(MFCCTmp, np.float)
            wltTmp = np.array(wltTmp, np.float)
            beatTmp = np.array(beatTmp, np.float)
            MFCCDataMat = np.append(MFCCDataMat, MFCCTmp)
            wltDataMat = np.append(wltDataMat, wltTmp)
            beatDataMat = np.append(beatDataMat, beatTmp)
            MFCCTmp = np.zeros(0)
            wltTmp = np.zeros(0)
            beatTmp = np.zeros(0)
        MFCCDataMat.shape = -1, MFCCPara_dimension
        wltDataMat.shape = -1, wltPara_dimension
        beatDataMat.shape = -1, beats_dimension
        fr.close()
                
        #   print len(dataMat)
        #   print len(dataMat[0])
                
        MFCCDataMat = mat(MFCCDataMat)
        wltDataMat = mat(wltDataMat)
        beatDataMat = mat(beatDataMat)
        print MFCCDataMat.shape
        print wltDataMat.shape
        print beatDataMat.shape
               
        MFCCLowDData, MFCCVals = PCA(MFCCDataMat, MFCC_TOP_FEAT)
        print "MFCC PCA finish"
        wltLowDData, wltVals = PCA(wltDataMat, WLT_TOP_FEAT)
        print "wlt PCA finish"
        beatLowDData, beatVals = PCA(beatDataMat, BEAT_TOP_FEAT)
        print "beat PCA finish"

        MFCCLowDData = np.array(MFCCLowDData)
        wltLowDData = np.array(wltLowDData)
        beatLowDData = np.array(beatLowDData)

        fw = open(PROJECT_DOC_DIR + PROJECT_PCADATA_DIR, "w")
        dataToDisk = []
        for i in xrange(len(MFCCLowDData)):
            dataToDisk.append(fileName[i] + "\t")
            MFCCLowDDataTmp = MFCCLowDData[i].real 
            for j in xrange(len(MFCCLowDDataTmp)):
                dataToDisk.append(str(MFCCLowDDataTmp[j]) + "\t")
            wltLowDDataTmp = wltLowDData[i].real
            for j in xrange(len(wltLowDDataTmp)):
                dataToDisk.append(str(wltLowDDataTmp[j]) + "\t")
            beatLowDDataTmp = beatLowDData[i].real
            for j in xrange(len(beatLowDDataTmp)):
                dataToDisk.append(str(beatLowDDataTmp[j]) + "\t")
            dataToDisk.append("\n")
            
        fw.writelines(dataToDisk)
        dataToDisk = []
        fw.close()
        
    if tag == "all": # PCA for all elements in feature
        print "TOP_FEAT: ", TOP_FEAT
        fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
        lineCount = 0
        data = []
        for line in fr:
            line = line.strip().split("\t")
            fileName.append(line[0])
            for i in xrange(beats_dimension):
                if line[i + 1 + frequencyPara_dimension] == "nan":
                    line[i + 1 + frequencyPara_dimension] = "0"
                    print "fix a nan"
            for i in xrange(len(line) - 1):
                data.append(float(line[1 + i]))
            lineCount += 1
        data = np.array(data)
        print data.shape
        data.shape = lineCount, -1
        fr.close()
                
        #   print len(dataMat)
        #   print len(dataMat[0])
                
        data = mat(data)
        print data.shape
               
        lowDData, vals = PCA(data, TOP_FEAT)
        print "data PCA finish"
        print lowDData.shape

        lowDData = np.array(lowDData)

        fw = open(PROJECT_DOC_DIR + PROJECT_PCADATA_DIR, "w")
        dataToDisk = []
        for i in xrange(len(lowDData)):
            dataToDisk.append(fileName[i] + "\t")
            lowDDataTmp = lowDData[i].real 
            for j in xrange(len(lowDDataTmp)):
                dataToDisk.append(str(lowDDataTmp[j]) + "\t")
            dataToDisk.append("\n")
            
        fw.writelines(dataToDisk)
        dataToDisk = []
        fw.close()
    
def drawPCASample(tag):
    global MFCC_TOP_FEAT
    global WLT_TOP_FEAT
    global BEAT_TOP_FEAT
    global MFCCPara_dimension
    global wltPara_dimension
    global beats_dimension
    global frequencyPara_dimension
    global PROJECT_COSDIS_DIR
    if tag == "all":
        alternative = np.zeros(0)
        blues = np.zeros(0)
        electronic = np.zeros(0)
        folkcountry = np.zeros(0)
        funksoulrnb = np.zeros(0)
        jazz = np.zeros(0)
        pop = np.zeros(0)
        raphiphop = np.zeros(0)
        rock = np.zeros(0)
        tmp = np.zeros(0)
        MFCCDataMat = np.zeros(0)
        wltDataMat = np.zeros(0)
        beatDataMat = np.zeros(0)
        MFCCTmp = np.zeros(0)
        wltTmp = np.zeros(0)
        beatTmp = np.zeros(0)
        fileName = []
                
        MFCCLowDData = np.zeros(0)
        wltLowDData = np.zeros(0)
        beatLowDData = np.zeros(0)
            
        frLowD = open(PROJECT_DOC_DIR + PROJECT_PCADATA_DIR, "r")
        for line in frLowD:
            line = line.strip().split("\t")
            fileName.append(line[0])
            MFCCLowDData = np.append(MFCCLowDData, np.array(line[1 : MFCC_TOP_FEAT + 1]))
            wltLowDData = np.append(wltLowDData, np.array(line[MFCC_TOP_FEAT + 1 : MFCC_TOP_FEAT + WLT_TOP_FEAT + 1]))
            beatLowDData = np.append(beatLowDData, np.array(line[MFCC_TOP_FEAT + WLT_TOP_FEAT + 1 : \
                            MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 1]))
            
        MFCCLowDData.shape = -1, MFCC_TOP_FEAT
        wltLowDData.shape = -1, WLT_TOP_FEAT
        beatLowDData.shape = -1, BEAT_TOP_FEAT

        randIndex = random.randint(TOTAL_MUSIC)
        MFCCVec = MFCCLowDData[randIndex]
        wltVec = wltLowDData[randIndex]
        beatVec = beatLowDData[randIndex]
        print fileName[randIndex]

        fw = open(PROJECT_DOC_DIR + PROJECT_COSDIS_DIR, "w")
        cosData =[]
        tmp = np.zeros(0)
        for i in xrange(len(MFCCLowDData)):
            if "alternative" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                                - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                                - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                alternative = np.append(alternative, tmp)
                tmp = np.zeros(0)
            if "blues" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                                - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                                 - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                                - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                blues = np.append(blues, tmp)
                tmp = np.zeros(0)
            if "electronic" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                electronic = np.append(electronic, tmp)
                tmp = np.zeros(0)
            if "folkcountry" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                folkcountry = np.append(folkcountry, tmp)
                tmp = np.zeros(0)
            if "funksoulrnb" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                funksoulrnb = np.append(funksoulrnb, tmp)
                tmp = np.zeros(0)
            if "jazz" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                jazz = np.append(jazz, tmp)
                tmp = np.zeros(0)
            if "pop" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                pop = np.append(pop, tmp)
                tmp = np.zeros(0)
            if "raphiphop" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                raphiphop = np.append(raphiphop, tmp)
                tmp = np.zeros(0)
            if "rock" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                rock = np.append(rock, tmp)
                tmp = np.zeros(0)

        fw.writelines(cosData)
        fw.close()
        alternative.shape = -1, 3
        alternative = alternative.T
        blues.shape = -1, 3
        blues = blues.T
        electronic.shape = -1, 3
        electronic = electronic.T
        folkcountry.shape = -1, 3
        folkcountry = folkcountry.T
        funksoulrnb.shape = -1, 3
        funksoulrnb = funksoulrnb.T
        jazz.shape = -1, 3
        jazz = jazz.T
        pop.shape = -1, 3
        pop = pop.T
        raphiphop.shape = -1, 3
        raphiphop = raphiphop.T
        rock.shape = -1, 3
        rock = rock.T

        draw3D = True
        mark = ['s','o','^','v','>','<','d','p','h','8','+','*']

        if draw3D == True:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = "3d")
#           ax.scatter(alternative[0], alternative[1], alternative[2], marker = mark[0], color = np.array([[0,0,1,0.3]]), label = str("alternative"))
            ax.scatter(blues[0], blues[1], blues[2], marker = mark[1], color = np.array([[0,1,0, 0.3]]), label = str("blues"))
#           ax.scatter(electronic[0], electronic[1], electronic[2], marker = mark[6], color = np.array([[1,0,0, 0.3]]), label = str("electronic"))
            ax.scatter(folkcountry[0], folkcountry[1], folkcountry[2], marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.3]]), label = str("folkcountry"))
            ax.scatter(funksoulrnb[0], funksoulrnb[1], funksoulrnb[2], marker = mark[11], color = np.array([[0.3,0.1,0.9, 0.3]]), label = str("funksoulrnb"))
            ax.scatter(jazz[0], jazz[1], jazz[2], marker = mark[11], color = np.array([[0.65,0.145,0.8, 0.3]]), label = str("jazz"))
            ax.scatter(pop[0], pop[1], pop[2], marker = mark[11], color = np.array([[0.2,0.5789,0.67, 0.3]]), label = str("pop"))
#           ax.scatter(raphiphop[0], raphiphop[1], raphiphop[2], marker = mark[11], color = np.array([[0.57,0.2,0.3, 0.3]]), label = str("raphiphop"))
            ax.scatter(rock[0], rock[1], rock[2], marker = mark[11], color = np.array([[0,0,0, 0.3]]), label = str("rock"))

            ax.set_xlabel("MFCC")
            ax.set_ylabel("wlt")
            ax.set_zlabel("beat")
            plt.show()
        
    if tag == "instruments":
        acoustic = np.zeros(0)
        piano = np.zeros(0)
        hiphop = np.zeros(0)
        relaxing = np.zeros(0)
        rock = np.zeros(0)
        MFCCDataMat = np.zeros(0)
        wltDataMat = np.zeros(0)
        beatDataMat = np.zeros(0)
        MFCCTmp = np.zeros(0)
        wltTmp = np.zeros(0)
        beatTmp = np.zeros(0)
        fileName = []
            
        fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
        for line in fr:
            line = line.strip().split("\t")
            fileName.append(line[0])
                
        MFCCLowDData = np.zeros(0)
        wltLowDData = np.zeros(0)
        beatLowDData = np.zeros(0)
            
        frLowD = open(PROJECT_DOC_DIR + PROJECT_PCADATA_DIR, "r")
        for line in frLowD:
            line = line.strip().split("\t")
            fileName.append(line[0])
            MFCCLowDData = np.append(MFCCLowDData, np.array(line[1 : MFCC_TOP_FEAT + 1]))
            wltLowDData = np.append(wltLowDData, np.array(line[MFCC_TOP_FEAT + 1 : MFCC_TOP_FEAT + WLT_TOP_FEAT + 1]))
            beatLowDData = np.append(beatLowDData, np.array(line[MFCC_TOP_FEAT + WLT_TOP_FEAT + 1 : \
                            MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 1]))
            
        MFCCLowDData.shape = -1, MFCC_TOP_FEAT
        wltLowDData.shape = -1, WLT_TOP_FEAT
        beatLowDData.shape = -1, BEAT_TOP_FEAT
        '''
        randIndex = random.randint(TOTAL_MUSIC)
        MFCCVec = MFCCLowDData[randIndex]
        wltVec = wltLowDData[randIndex]
        beatVec = beatLowDData[randIndex]
        print fileName[randIndex]
        '''
        MFCCVec = np.ones(MFCC_TOP_FEAT)
        wltVec = np.ones(WLT_TOP_FEAT)
        beatVec = np.ones(BEAT_TOP_FEAT)
        fw = open(PROJECT_DOC_DIR + PROJECT_COSDIS_DIR, "w")
        cosData = []
        tmp = np.zeros(0)
        for i in xrange(TOTAL_MUSIC):
            if "acoustic" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                acoustic = np.append(acoustic, tmp)
                tmp = np.zeros(0)
            if "piano" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                piano = np.append(piano, tmp)
                tmp = np.zeros(0)
            if "hiphop" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                hiphop = np.append(hiphop, tmp)
                tmp = np.zeros(0)
            if "relaxing" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                cosData.append(str(simi) + "\t")
                tmp = np.append(tmp, simi)
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                relaxing = np.append(relaxing, tmp)
                tmp = np.zeros(0)
            if "rock" in fileName[i]:
                cosData.append(fileName[i] + "\t")
                cosSimi = getCosineDis(MFCCVec, MFCCLowDData[i])
                eucSimi = 0
                for j in xrange(MFCC_TOP_FEAT):
                    diff = math.pow((float(MFCCLowDData[i][j]) \
                            - float(MFCCVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(wltVec, wltLowDData[i])
                eucSimi = 0
                for j in xrange(WLT_TOP_FEAT):
                    diff = math.pow((float(wltLowDData[i][j]) \
                            - float(wltVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\t")
                cosSimi = getCosineDis(beatVec, beatLowDData[i])
                eucSimi = 0
                for j in xrange(BEAT_TOP_FEAT):
                    diff = math.pow((float(beatLowDData[i][j]) \
                            - float(beatVec[j])), 2)
                    eucSimi += diff
                eucSimi = math.sqrt(eucSimi)
                simi = eucSimi * (1 - cosSimi)
                tmp = np.append(tmp, simi)
                cosData.append(str(simi) + "\n")
                rock = np.append(rock, tmp)
                tmp = np.zeros(0)

        fw.writelines(cosData)
        fw.close()
        acoustic.shape = -1, 3
        acoustic = acoustic.T
        piano.shape = -1, 3
        piano = piano.T
        hiphop.shape = -1, 3
        hiphop = hiphop.T
        relaxing.shape = -1, 3
        relaxing = relaxing.T
        rock.shape = -1, 3
        rock = rock.T

        draw3D = True
        mark = ['s','o','^','v','>','<','d','p','h','8','+','*']
             
        if draw3D == False:
        #   plt.scatter(acoustic[0], acoustic[2], marker = mark[0], color = np.array([[0,0,1,0.5]]), label = str("acoustic"))
            plt.scatter(acoustic[0].mean(), acoustic[2].mean(), s = 100, marker = mark[0], color = np.array([[0,0,1,1]]))            
        #   plt.scatter(piano[0], piano[2],  marker = mark[1], color = np.array([[0,1,0, 0.5]]), label = str("piano"))
            plt.scatter(piano[0].mean(), piano[2].mean(), s = 100, marker = mark[1], color = np.array([[0,1,0, 1]]))
            plt.scatter(hiphop[0], hiphop[2],  marker = mark[6], color = np.array([[1,0,0, 0.5]]), label = str("hiphop"))
            plt.scatter(hiphop[0].mean(), hiphop[2].mean(), s = 100, marker = mark[6], color = np.array([[1,0,0, 1]]))
        #   plt.scatter(relaxing[0], relaxing[2],  marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.5]]), label = str("relaxing"))
            plt.scatter(relaxing[0].mean(), relaxing[2].mean(), s = 100, marker = mark[3], color = np.array([[0.7,0.3,0.2, 1]]))
            plt.scatter(rock[0], rock[2],  marker = mark[11], color = np.array([[0,0,0, 0.5]]), label = str("rock"))
            plt.scatter(rock[0].mean(), rock[2].mean(), s = 100, marker = mark[11], color = np.array([[0,0,0, 1]]))           
                
            '''
        #        plt.scatter(acoustic[0], acoustic[1], marker = mark[0], color = np.array([[0,0,1,0.5]]), label = str("acoustic"))
                plt.scatter(acoustic[0].mean(), acoustic[1].mean(), s = 100, marker = mark[0], color = np.array([[0,0,1,1]]))            
        #        plt.scatter(piano[0], piano[1],  marker = mark[1], color = np.array([[0,1,0, 0.5]]), label = str("piano"))
                plt.scatter(piano[0].mean(), piano[1].mean(), s = 100, marker = mark[1], color = np.array([[0,1,0, 1]]))
                plt.scatter(hiphop[0], hiphop[1],  marker = mark[6], color = np.array([[1,0,0, 0.5]]), label = str("hiphop"))
                plt.scatter(hiphop[0].mean(), hiphop[1].mean(), s = 100, marker = mark[6], color = np.array([[1,0,0, 1]]))
        #        plt.scatter(relaxing[0], relaxing[1],  marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.5]]), label = str("relaxing"))
                plt.scatter(relaxing[0].mean(), relaxing[1].mean(), s = 100, marker = mark[3], color = np.array([[0.7,0.3,0.2, 1]]))
                plt.scatter(rock[0], rock[1],  marker = mark[11], color = np.array([[0,0,0, 0.5]]), label = str("rock"))
                plt.scatter(rock[0].mean(), rock[1].mean(), s = 100, marker = mark[11], color = np.array([[0,0,0, 1]]))
            '''
            '''
        #        plt.scatter(acoustic[1], acoustic[2], marker = mark[0], color = np.array([[0,0,1,0.5]]), label = str("acoustic"))
                plt.scatter(acoustic[1].mean(), acoustic[2].mean(), s = 100, marker = mark[0], color = np.array([[0,0,1,1]]))            
        #        plt.scatter(piano[1], piano[2],  marker = mark[1], color = np.array([[0,1,0, 0.5]]), label = str("piano"))
                plt.scatter(piano[1].mean(), piano[2].mean(), s = 100, marker = mark[1], color = np.array([[0,1,0, 1]]))
                plt.scatter(hiphop[1], hiphop[2],  marker = mark[6], color = np.array([[1,0,0, 0.5]]), label = str("hiphop"))
                plt.scatter(hiphop[1].mean(), hiphop[2].mean(), s = 100, marker = mark[6], color = np.array([[1,0,0, 1]]))
        #        plt.scatter(relaxing[1], relaxing[2],  marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.5]]), label = str("relaxing"))
                plt.scatter(relaxing[1].mean(), relaxing[2].mean(), s = 100, marker = mark[3], color = np.array([[0.7,0.3,0.2, 1]]))
                plt.scatter(rock[1], rock[2],  marker = mark[11], color = np.array([[0,0,0, 0.5]]), label = str("rock"))
                plt.scatter(rock[1].mean(), rock[2].mean(), s = 100, marker = mark[11], color = np.array([[0,0,0, 1]]))
            '''
            plt.legend()
            plt.show()

        if draw3D == True:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = "3d")
            ax.scatter(acoustic[0], acoustic[1], acoustic[2],  marker = mark[0], color = np.array([[0,0,1,0.5]]), label = str("acoustic"))
            ax.scatter(piano[0], piano[1], piano[2],  marker = mark[1], color = np.array([[0,1,0, 0.5]]), label = str("piano"))
            ax.scatter(hiphop[0], hiphop[1], hiphop[2],  marker = mark[6], color = np.array([[1,0,0, 0.5]]), label = str("hiphop"))
            ax.scatter(relaxing[0], relaxing[1], relaxing[2],  marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.5]]), label = str("relaxing"))
            ax.scatter(rock[0], rock[1], rock[2],  marker = mark[11], color = np.array([[0,0,0, 0.5]]), label = str("rock"))

            ax.set_xlabel("MFCC")
            ax.set_ylabel("wlt")
            ax.set_zlabel("beat")
            plt.show()

        

def PCASimilarity(fileName, word):
    global TOTAL_MUSIC
    global PROJECT_DOC_DIR
    global PROJECT_PCADATA_DIR
    global PROJECT_RAWFEAT_DIR
    global MFCC_TOP_FEAT
    global WLT_TOP_FEAT
    global BEAT_TOP_FEAT
    fileName = (fileName.decode("gb18030")).encode("utf8")
    nameList = []

    #read PCA data and insert the music title
    MFCCLowDData = np.zeros(0)
    wltLowDData = np.zeros(0)
    beatLowDData = np.zeros(0)
    PCAData = []
    frLowD = open(PROJECT_DOC_DIR + PROJECT_PCADATA_DIR, "r")
    for line in frLowD:
        line = line.strip().split("\t")
        nameList.append(line[0])
        MFCCLowDData = np.append(MFCCLowDData, np.array(line[1 : MFCC_TOP_FEAT + 1]))
        wltLowDData = np.append(wltLowDData, np.array(line[MFCC_TOP_FEAT + 1 : MFCC_TOP_FEAT + WLT_TOP_FEAT + 1]))
        beatLowDData = np.append(beatLowDData, np.array(line[MFCC_TOP_FEAT + WLT_TOP_FEAT + 1 : \
                            MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 1]))

    MFCCLowDData.shape = -1, MFCC_TOP_FEAT
    wltLowDData.shape = -1, WLT_TOP_FEAT
    beatLowDData.shape = -1, BEAT_TOP_FEAT

    PCATmpInLine = []
    for i in xrange(TOTAL_MUSIC):
        PCATmpInLine.append(nameList[i])
        PCATmpInLine.extend(MFCCLowDData[i].tolist())
        PCATmpInLine.extend(wltLowDData[i].tolist())
        PCATmpInLine.extend(beatLowDData[i].tolist())
        PCAData.append(PCATmpInLine)
        PCATmpInLine = []
   
    #find target
    target = []
    count = 0
    for each in PCAData:
        if fileName == each[0]:
            for ele in each:
                target.append(ele)
            del PCAData[count]
            break
        count += 1
    if word == "MFCCsWlt":
        res = []
        for lineCount in xrange(len(PCAData)):
            cosSimi = getCosineDis(PCAData[lineCount][1 : MFCC_TOP_FEAT + 1], \
                                target[1 : MFCC_TOP_FEAT + 1])
            eucSimi = 0
            for i in xrange(MFCC_TOP_FEAT):
                diff = math.pow((float(PCAData[lineCount][i + 1]) \
                                 - float(target[i + 1])), 2)
                eucSimi += diff
            eucSimi = math.sqrt(eucSimi)
            simi = eucSimi * (1 - cosSimi)
            simi = simi * 0.003
            PCAData[lineCount].extend([simi])
            res.append(PCAData[lineCount])
        res.sort(key = lambda x : x[MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 1])
        del PCAData

        dataForWlt = []
        top = int(TOTAL_MUSIC * 0.6)
        for i in xrange(top):
            if i < 5:
                print "----------  NO.%d  ----------" % (i + 1)
                print (str(res[i][0])).encode("gb18030")
                print (str(res[i][MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 1])).encode("UTF-8")
            dataForWlt.append(res[i])
        del res

        res = []
        for lineCount in xrange(len(dataForWlt)):
            cosSimi = getCosineDis(dataForWlt[lineCount][MFCC_TOP_FEAT + 1 : MFCC_TOP_FEAT + WLT_TOP_FEAT + 1], \
                                target[MFCC_TOP_FEAT + 1 : MFCC_TOP_FEAT + WLT_TOP_FEAT + 1])
            eucSimi = 0
            for i in xrange(WLT_TOP_FEAT):
                diff = math.pow((float(dataForWlt[lineCount][MFCC_TOP_FEAT + i + 1]) \
                                 - float(target[MFCC_TOP_FEAT + i + 1])), 2)
                eucSimi += diff
            eucSimi = math.sqrt(eucSimi)
            simi = eucSimi * (1 - cosSimi)
            simi = simi * 0.1
            dataForWlt[lineCount].extend([simi])
            res.append(dataForWlt[lineCount])
        res.sort(key = lambda x : x[MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 2])
        del dataForWlt

        dataForBeats = []
        top = int(TOTAL_MUSIC * 0.6)
        for i in xrange(top):
            dataForBeats.append(res[i])
        del res

        res = []
        for lineCount in xrange(len(dataForBeats)):
            cosSimi = getCosineDis(dataForBeats[lineCount][MFCC_TOP_FEAT + WLT_TOP_FEAT + 1 : MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 1], \
                                target[MFCC_TOP_FEAT + WLT_TOP_FEAT + 1 : MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 1])
            eucSimi = 0
            for i in xrange(BEAT_TOP_FEAT):
                eucSimi = eucSimi + (float(dataForBeats[lineCount][MFCC_TOP_FEAT + WLT_TOP_FEAT + i + 1]) - float(target[MFCC_TOP_FEAT + WLT_TOP_FEAT + i + 1])) * \
                              (float(dataForBeats[lineCount][MFCC_TOP_FEAT + WLT_TOP_FEAT + i + 1]) - float(target[MFCC_TOP_FEAT + WLT_TOP_FEAT + i + 1]))

            eucSimi = math.sqrt(eucSimi)
            simi = eucSimi * (1 - cosSimi)
            simi = simi * 30
            dataForBeats[lineCount].extend([simi])
            dataForBeats[lineCount].extend([(simi + \
                        dataForBeats[lineCount]\
                        [MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 1] + \
                        dataForBeats[lineCount]\
                        [MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 2]) / 3])
            res.append(dataForBeats[lineCount])
        res.sort(key = lambda x : x[MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 4])
        del dataForBeats

        top = 5
        for i in xrange(top):
            if(i < top):
                print "----------  NO.%d  ----------" % (i + 1)
                print (str(res[i][0])).encode("gb18030")
            
            if(i < top):
                # MFCCs
                # Wlt
                # Beats
                # in total
                print (str(res[i][MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 4])).encode("UTF-8") + \
                  "(" + (str(res[i][MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 1])).encode("UTF-8") + ")" + \
                  "(" + (str(res[i][MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 2])).encode("UTF-8") + ")"+ \
                  "(" + (str(res[i][MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 3])).encode("UTF-8") + ")"+ \
                  "(" + (str(res[i][MFCC_TOP_FEAT + WLT_TOP_FEAT + BEAT_TOP_FEAT + 4])).encode("UTF-8") + ")"
                '''
                freqs = np.linspace(0, 120, CLength * 12)
                pl.figure(figsize=(8,4))
                pl.subplot(111)
                pl.plot(freqs, target\
                    [frequencyPara_dimension - CLength * 12:frequencyPara_dimension],
                    alpha = 0.7)
                pl.plot(freqs, res[i]\
                    [frequencyPara_dimension - CLength * 12:frequencyPara_dimension],
                    alpha = 0.7)
                pl.xlabel(u"相似度")
                pl.show()
                '''
def getEighenMusicSpace():
    testMusic = 59
    dataMat = np.zeros(0, dtype = ctypes.c_int)
    fileList = os.listdir(PROJECT_TEST_DIR[:len(PROJECT_TEST_DIR) - 1])
    for audioFile in fileList:
        if (audioFile.strip().split(".")[-1] == "wav") or \
           (audioFile.strip().split(".")[-1] == "WAV") or \
           (audioFile.strip().split(".")[-1] == "mp3") or \
           (audioFile.strip().split(".")[-1] == "MP3"):
            #form the music data to 60s length
            incomingData = getData(audioFile.encode("gb18030"))
            dataLength = len(incomingData)
            lengthDiff = dataLength - wltAnalysisLength
            while lengthDiff < 0:
                if abs(lengthDiff) <= dataLength:
                    incomingData = np.append(incomingData, \
                                incomingData[dataLength - abs(lengthDiff) : ])
                    dataLength = len(incomingData)
                    lengthDiff = dataLength - wltAnalysisLength
                else:
                    incomingData = np.append(incomingData, incomingData)
                    dataLength = len(incomingData)
                    lengthDiff = dataLength - wltAnalysisLength
            if lengthDiff >= 0:
                lengthDiff = int(lengthDiff / 2)
                incomingData = incomingData[lengthDiff : \
                                wltAnalysisLength + lengthDiff]
            incomingData.shape = -1, CHANNELS
            incomingData = incomingData.T

            incomingData = incomingData[0]
            CA, CD = dwt(incomingData, dwtModel)
#            CA, CD = dwt(CA, dwtModel)
#            CA, CD = dwt(CA, dwtModel)
#            CA, CD = dwt(CD, dwtModel)
#            CA, CD = dwt(CD, dwtModel)

            #use time to save memory
            dataMat = np.append(dataMat, CD)

    #start PCA
    topNFeat = testMusic
    print len(dataMat)
    dataMat.shape = -1, testMusic
    dataMat = np.absolute(dataMat)
    dataMat = np.mat(dataMat)
    print dataMat.shape
    meanVals = mean(dataMat, axis = 1) #axis = 0 mean(list) axis = 1 mean(line)
    meanVals = np.array(meanVals)
    fw = open(PROJECT_DOC_DIR + u"meanVals.txt", "w")
    dataToWrite = []
    for i in xrange(len(meanVals)):
        dataToWrite.append(str(meanVals[i][0]) + "\t")
    fw.writelines(dataToWrite)
    fw.close()
    meanVals = np.mat(meanVals)
#    meanRemoved = dataMat - meanVals
    meanRemoved = dataMat
    covMat = np.cov(meanRemoved, rowvar = 0)
    print covMat.shape
    eigVals, eigVects = np.linalg.eig(np.mat(covMat))
    print eigVals.shape
    print eigVects.shape
    eigValIndex = argsort(eigVals)
    eigValIndex = eigValIndex[:-(topNFeat + 1):-1]
    topEigVects = eigVects[:,eigValIndex]
    sumEig = 0
    for i in xrange(len(eigVals)):
        sumEig += eigVals[i]
    partition = 0
    printForOnce = False
    for i in xrange(len(eigVals)):
        partition += eigVals[i]
        ratio = float(partition) / float(sumEig)
        if ratio < 0.95:
            print "%d, %f" % (i, ratio)
        else:
            if printForOnce == False:
                print "%d, %f" % (i, ratio)
                printForOnce = True
    F = dataMat * eigVects
    del dataMat
    del covMat
    del eigVals
    del eigVects
    del meanRemoved
    print F.shape
    F = np.array(F)
    dataToWrite = []
    fw = open(PROJECT_DOC_DIR + PROJECT_PCASPACE_DIR, "w")
    for line in xrange(len(F)):
        for i in xrange(len(F[line])):
            dataToWrite.append(str(F[line][i]) + "\t")
    fw.writelines(dataToWrite)
    fw.close()

def sigmoid(inX):
    return 1.0 / (1 + exp(-inX))

def logRegres(class0, class1, tag = "", draw3D = False):
    global TOTAL_MUSIC
    global MFCC_TOP_FEAT
    global WLT_TOP_FEAT
    global BEAT_TOP_FEAT
    global MFCCPara_dimension
    global wltPara_dimension
    global beats_dimension
    global frequencyPara_dimension
    trainPartition = 0.7
    PATH = PROJECT_RAWFEAT_DIR
    
    trainC0Num = 0
    trainC1Num = 0
    totalC0Num = 0
    totalC1Num = 0

    # generate trainClassData and trainClassLabels
    totalClassLabels = []
    totalClassData = []
    fr = open(PROJECT_DOC_DIR + PATH, "r")
    for line in fr:
        line = line.strip().split("\t")
#        print len(line)
#        return
        tmp = []
        if "@" in line[0]:
            genre = line[0].strip().split("@")
            if class0 == genre[1]:
                totalClassLabels.append(0)
                for index in xrange(wltPara_dimension):
                    line[MFCCPara_dimension + 1 + index] = float( \
                        line[MFCCPara_dimension + 1 + index]) / 3.0
                line[frequencyPara_dimension + 1] = float(line[frequencyPara_dimension + 1]) / 1
                line[frequencyPara_dimension + 2] = float(line[frequencyPara_dimension + 2]) / 1000
                line[frequencyPara_dimension + 3] = float(line[frequencyPara_dimension + 3]) / 1000
                line[frequencyPara_dimension + 4] = float(line[frequencyPara_dimension + 4]) / 1000
                line[frequencyPara_dimension + 5] = float(line[frequencyPara_dimension + 5]) / 1000
                line[frequencyPara_dimension + 6] = float(line[frequencyPara_dimension + 6]) / 1000
                line[frequencyPara_dimension + 7] = float(line[frequencyPara_dimension + 7]) * 5
                line[frequencyPara_dimension + 8] = float(line[frequencyPara_dimension + 8]) * 5
                tmp.extend(line[1 : MFCCPara_dimension + 81])
                tmp.extend(line[MFCCPara_dimension + 442 : MFCCPara_dimension + 522])
                tmp.extend(line[frequencyPara_dimension + 1 : ])
#                tmp.extend(line[1 : ])
#                print line[1:21]
#                print line[MFCCPara_dimension + 1:MFCCPara_dimension+21]
#                print line[frequencyPara_dimension + 1:]
#                return

                totalClassData.append(tmp)
                totalC0Num += 1
            if class1 == genre[1]:
                totalClassLabels.append(1)
                for index in xrange(wltPara_dimension):
                    line[MFCCPara_dimension + 1 + index] = float( \
                        line[MFCCPara_dimension + 1 + index]) / 3.0
                line[frequencyPara_dimension + 1] = float(line[frequencyPara_dimension + 1]) / 1
                line[frequencyPara_dimension + 2] = float(line[frequencyPara_dimension + 2]) / 1000
                line[frequencyPara_dimension + 3] = float(line[frequencyPara_dimension + 3]) / 1000
                line[frequencyPara_dimension + 4] = float(line[frequencyPara_dimension + 4]) / 1000
                line[frequencyPara_dimension + 5] = float(line[frequencyPara_dimension + 5]) / 1000
                line[frequencyPara_dimension + 6] = float(line[frequencyPara_dimension + 6]) / 1000
                line[frequencyPara_dimension + 7] = float(line[frequencyPara_dimension + 7]) * 5
                line[frequencyPara_dimension + 8] = float(line[frequencyPara_dimension + 8]) * 5
                tmp.extend(line[1 : MFCCPara_dimension + 81])
                tmp.extend(line[MFCCPara_dimension + 442 : MFCCPara_dimension + 522])
                tmp.extend(line[frequencyPara_dimension + 1 : ])
#                tmp.extend(line[1 : ])
                
                totalClassData.append(tmp)
                totalC1Num += 1
    for i in xrange(len(totalClassData)):
        for j in xrange(len(totalClassData[i])):
            totalClassData[i][j] = float(totalClassData[i][j])
    fr.close()
#    return

    trainC0Num = int(trainPartition * totalC0Num)
    trainC1Num = int(trainPartition * totalC1Num)

    accuracy = 0
    testIter = 20
    for testIt in xrange(testIter):
        trainC0Index = rnd.sample(range(totalC0Num), trainC0Num)
        trainC1Index = rnd.sample(range(totalC1Num), trainC1Num)
        testC0Index = [val for val in range(totalC0Num) if val not in trainC0Index]
        testC1Index = [val for val in range(totalC1Num) if val not in trainC1Index]
        for i in xrange(len(trainC1Index)):
            trainC1Index[i] += totalC0Num
        for i in xrange(len(testC1Index)):
            testC1Index[i] += totalC0Num

        print "class0 number: ", totalC0Num
        print "class1 number: ", totalC1Num
        print "train class0 number: ", trainC0Num
        print "train class1 number: ", trainC1Num
        print "test class0 number: ", totalC0Num - trainC0Num
        print "test class1 number: ", totalC1Num - trainC1Num
    #    print "trainC0Index: ", trainC0Index
    #    print "trainC1Index: ", trainC1Index
    #    print "testC0Index: ", testC0Index
    #    print "testC1Index: ", testC1Index

        tmpData = totalClassData
        tmpLabels = totalClassLabels
        trainClassData = []
        trainClassLabels = []
        testClassData = []
        testClassLabels = []
        for index in xrange(len(tmpData)):
            if index in trainC0Index:
                trainClassData.append(tmpData[index])
                trainClassLabels.append(tmpLabels[index])
            if index in trainC1Index:
                trainClassData.append(tmpData[index])
                trainClassLabels.append(tmpLabels[index])
            if (index not in trainC0Index) and (index not in trainC1Index):
                testClassData.append(tmpData[index])
                testClassLabels.append(tmpLabels[index])

        # logRegres
        w0 = []
        trainClassData = mat(trainClassData)
        trainClassLabels = mat(trainClassLabels).T
        print "trainClassData.shape: ", trainClassData.shape
        print "trainClassLabels.shape: ", trainClassLabels.shape
        m, n = trainClassData.shape
        alpha = 1
        maxIter = 5000
        weights = np.ones((n, 1))
        for it in xrange(maxIter):
            w0.append(weights[100])
            prediction = sigmoid(trainClassData * weights)
            error = trainClassLabels - prediction
            weights = weights + alpha * trainClassData.T * error
            
#        pl.figure(figsize=(8,4))
#        pl.subplot(111)
#        pl.plot(np.arange(len(w0)), np.array(w0), alpha = 0.7)
#        pl.show()

        # test
        error = 0
        for index in xrange(len(testClassData)):
            prediction = sigmoid(np.array(testClassData[index], np.float) * weights)
            print prediction
    #        print prediction
            if prediction >= 0.5:
                prediction = 1
            else:
                prediction = 0
            if prediction != testClassLabels[index]:
                print testClassData[index][:5]
                error += 1
        accuracy += (1.0 - float(error) / float(len(testClassData))) * 100
        print "error: ", error
        print "accuracy: %f%%" % (accuracy / (testIt + 1))

        if (draw3D):    
            #draw 3D
            fr = open(PROJECT_DOC_DIR + PATH, "r")

            if tag == "instruments":
                acoustic = []
                piano = []
                relaxing = []
                hiphop = []
                rock = []
                for line in fr:
                    line = line.strip().split("\t")
                    if "acoustic" in line[0]:
                        acoustic.append(line[1 : ])
                    if "piano" in line[0]:
                        piano.append(line[1 : ])
                    if "relaxing" in line[0]:
                        relaxing.append(line[1 : ])
                    if "hiphop" in line[0]:
                        hiphop.append(line[1 : ])
                    if "rock" in line[0]:
                        rock.append(line[1 : ])
                acoustic = np.array(acoustic, np.float)
                acoustic = acoustic.T
                piano = np.array(piano, np.float)
                piano = piano.T
                relaxing = np.array(relaxing, np.float)
                relaxing = relaxing.T
                hiphop = np.array(hiphop, np.float)
                hiphop = hiphop.T
                rock = np.array(rock, np.float)
                rock = rock.T

        #        print weights
                mark = ['s','o','^','v','>','<','d','p','h','8','+','*']
                fig=plt.figure()
                ax = Axes3D(fig)
                ax.scatter(acoustic[0], acoustic[1], acoustic[2],  marker = mark[0], color = np.array([[0,0,1,0.5]]))
                ax.scatter(piano[0], piano[1], piano[2],  marker = mark[1], color = np.array([[0,1,0, 0.5]]))
                ax.scatter(hiphop[0], hiphop[1], hiphop[2],  marker = mark[6], color = np.array([[1,0,0, 0.5]]))
                ax.scatter(relaxing[0], relaxing[1], relaxing[2],  marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.5]]))
                ax.scatter(rock[0], rock[1], rock[2],  marker = mark[11], color = np.array([[0,0,0, 0.5]]))

                ax.set_xlabel("MFCC")
                ax.set_ylabel("wlt")
                ax.set_zlabel("beat")
                LRx = np.arange(1,60000,750)
                LRy = np.arange(1,160000,2000)
                LRx, LRy = np.meshgrid(LRx, LRy)
                weights = np.array(weights)
                LRz = (-1 * weights[0] * LRx - weights[1] * LRy) / weights[2]

                ax.plot_surface(LRx, LRy, LRz, alpha = 0.1)
                plt.show()

            if tag == "all":
                alternative = []
                blues = []
                electronic = []
                folkcountry = []
                funksoulrnb = []
                jazz = []
                pop = []
                raphiphop = []
                rock = []
                for line in fr:
                    line = line.strip().split("\t")
                    if "alternative" in line[0]:
                        alternative.append(line[1 : ])
                    if "blues" in line[0]:
                        blues.append(line[1 : ])
                    if "electronic" in line[0]:
                        electronic.append(line[1 : ])
                    if "folkcountry" in line[0]:
                        folkcountry.append(line[1 : ])
                    if "funksoulrnb" in line[0]:
                        funksoulrnb.append(line[1 : ])
                    if "jazz" in line[0]:
                        jazz.append(line[1 : ])
                    if "pop" in line[0]:
                        pop.append(line[1 : ])
                    if "raphiphop" in line[0]:
                        raphiphop.append(line[1 : ])
                    if "rock" in line[0]:
                        rock.append(line[1 : ])
                
                alternative = np.array(alternative, np.float)
                alternative = alternative.T
                blues = np.array(blues, np.float)
                blues = blues.T
                electronic = np.array(electronic, np.float)
                electronic = electronic.T
                folkcountry = np.array(folkcountry, np.float)
                folkcountry = folkcountry.T
                funksoulrnb = np.array(funksoulrnb, np.float)
                funksoulrnb = funksoulrnb.T
                pop = np.array(pop, np.float)
                pop = pop.T
                raphiphop = np.array(raphiphop, np.float)
                raphiphop = raphiphop.T
                rock = np.array(rock, np.float)
                rock = rock.T
                jazz = np.array(alternative, np.float)
                jazz = alternative.T

                print weights
                mark = ['s','o','^','v','>','<','d','p','h','8','+','*']
                fig=plt.figure()
                ax = Axes3D(fig)
                ax = fig.add_subplot(111, projection = "3d")
                ax.scatter(alternative[0], alternative[1], alternative[2], marker = mark[0], color = np.array([[0,0,1,0.3]]), label = str("alternative"))
                ax.scatter(blues[0], blues[1], blues[2], marker = mark[1], color = np.array([[0,1,0, 0.3]]), label = str("blues"))
                ax.scatter(electronic[0], electronic[1], electronic[2], marker = mark[6], color = np.array([[1,0,0, 0.3]]), label = str("electronic"))
                ax.scatter(folkcountry[0], folkcountry[1], folkcountry[2], marker = mark[3], color = np.array([[0.7,0.3,0.2, 0.3]]), label = str("folkcountry"))
                ax.scatter(funksoulrnb[0], funksoulrnb[1], funksoulrnb[2], marker = mark[11], color = np.array([[0.3,0.1,0.9, 0.3]]), label = str("funksoulrnb"))
                ax.scatter(jazz[0], jazz[1], jazz[2], marker = mark[11], color = np.array([[0.65,0.145,0.8, 0.3]]), label = str("jazz"))
                ax.scatter(pop[0], pop[1], pop[2], marker = mark[11], color = np.array([[0.2,0.5789,0.67, 0.3]]), label = str("pop"))
                ax.scatter(raphiphop[0], raphiphop[1], raphiphop[2], marker = mark[11], color = np.array([[0.57,0.2,0.3, 0.3]]), label = str("raphiphop"))
                ax.scatter(rock[0], rock[1], rock[2], marker = mark[11], color = np.array([[0,0,0, 0.3]]), label = str("rock"))

                ax.set_xlabel("MFCC")
                ax.set_ylabel("wlt")
                ax.set_zlabel("beat")
                LRx = np.arange(1,60000,750)
                LRy = np.arange(1,160000,2000)
                LRx, LRy = np.meshgrid(LRx, LRy)
                weights = np.array(weights)
                LRz = (-1 * weights[0] * LRx - weights[1] * LRy) / weights[2]

                ax.plot_surface(LRx, LRy, LRz, alpha = 0.1)
                plt.show()

def sklearnSVM(class0, class1):
    global TOTAL_MUSIC
    global MFCC_TOP_FEAT
    global WLT_TOP_FEAT
    global BEAT_TOP_FEAT
    global MFCCPara_dimension
    global wltPara_dimension
    global beats_dimension
    global frequencyPara_dimension
    trainPartition = 0.9
    PATH = PROJECT_RAWFEAT_DIR
#    PATH = PROJECT_PCADATA_DIR
    
    trainC0Num = 0
    trainC1Num = 0
    totalC0Num = 0
    totalC1Num = 0

    # generate trainClassData and trainClassLabels
    totalClassLabels = []
    totalClassData = []
    fr = open(PROJECT_DOC_DIR + PATH, "r")
    for line in fr:
        line = line.strip().split("\t")
#        print len(line)
#        return
        tmp = []
        if "@" in line[0]:
#            print line[0]
            genre = line[0].strip().split("@")
            if class0 == genre[1]:
                totalClassLabels.append(0)
                
                for index in xrange(wltPara_dimension):
                    line[MFCCPara_dimension + 1 + index] = float( \
                        line[MFCCPara_dimension + 1 + index]) / 2.0
                line[frequencyPara_dimension + 1] = float(line[frequencyPara_dimension + 1]) / 1
                line[frequencyPara_dimension + 2] = float(line[frequencyPara_dimension + 2]) / 1000
                line[frequencyPara_dimension + 3] = float(line[frequencyPara_dimension + 3]) / 1000
                line[frequencyPara_dimension + 4] = float(line[frequencyPara_dimension + 4]) / 1000
                line[frequencyPara_dimension + 5] = float(line[frequencyPara_dimension + 5]) / 1000
                line[frequencyPara_dimension + 6] = float(line[frequencyPara_dimension + 6]) / 1000
                line[frequencyPara_dimension + 7] = float(line[frequencyPara_dimension + 7]) * 5
                line[frequencyPara_dimension + 8] = float(line[frequencyPara_dimension + 8]) * 5
#                tmp.extend(line[1 : MFCCPara_dimension + 51])
#                tmp.extend(line[MFCCPara_dimension + 1025 : MFCCPara_dimension + 1025 + 50])
#                tmp.extend(line[frequencyPara_dimension + 1 : ])
                
                tmp.extend(line[1 : ])

                '''
                sumS = 0
                for i in xrange(441):
                    sumS = sumS + float(line[MFCCPara_dimension + i])
                print sumS
                S = 0
                for i in xrange(441):
                    S += float(line[MFCCPara_dimension + i])
                    if S / sumS < 0.95:
                        print i, S / sumS
                    else: return
                '''
                

#                print line[1:21]
#                print line[MFCCPara_dimension + 1:MFCCPara_dimension+21]
#                print line[frequencyPara_dimension + 1:]
#                return

                totalClassData.append(tmp)
                totalC0Num += 1
            if class1 == genre[1]:
                totalClassLabels.append(1)
                
                for index in xrange(wltPara_dimension):
                    line[MFCCPara_dimension + 1 + index] = float( \
                        line[MFCCPara_dimension + 1 + index]) / 2.0
                line[frequencyPara_dimension + 1] = float(line[frequencyPara_dimension + 1]) / 1
                line[frequencyPara_dimension + 2] = float(line[frequencyPara_dimension + 2]) / 1000
                line[frequencyPara_dimension + 3] = float(line[frequencyPara_dimension + 3]) / 1000
                line[frequencyPara_dimension + 4] = float(line[frequencyPara_dimension + 4]) / 1000
                line[frequencyPara_dimension + 5] = float(line[frequencyPara_dimension + 5]) / 1000
                line[frequencyPara_dimension + 6] = float(line[frequencyPara_dimension + 6]) / 1000
                line[frequencyPara_dimension + 7] = float(line[frequencyPara_dimension + 7]) * 5
                line[frequencyPara_dimension + 8] = float(line[frequencyPara_dimension + 8]) * 5
#                tmp.extend(line[1 : MFCCPara_dimension + 51])
#                tmp.extend(line[MFCCPara_dimension + 1025 : MFCCPara_dimension + 1025 + 50])
#                tmp.extend(line[frequencyPara_dimension + 1 : ])
                
                tmp.extend(line[1 : ])
                
                totalClassData.append(tmp)
                totalC1Num += 1
    for i in xrange(len(totalClassData)):
        for j in xrange(len(totalClassData[i])):
            totalClassData[i][j] = float(totalClassData[i][j])
    fr.close()

    trainC0Num = int(trainPartition * totalC0Num)
    trainC1Num = int(trainPartition * totalC1Num)

    accuracy = 0
    testIter = 20
    for testIt in xrange(testIter):
        trainC0Index = rnd.sample(range(totalC0Num), trainC0Num)
        trainC1Index = rnd.sample(range(totalC1Num), trainC1Num)
        testC0Index = [val for val in range(totalC0Num) if val not in trainC0Index]
        testC1Index = [val for val in range(totalC1Num) if val not in trainC1Index]
        for i in xrange(len(trainC1Index)):
            trainC1Index[i] += totalC0Num
        for i in xrange(len(testC1Index)):
            testC1Index[i] += totalC0Num

        print "class0 number: ", totalC0Num
        print "class1 number: ", totalC1Num
        print "train class0 number: ", trainC0Num
        print "train class1 number: ", trainC1Num
        print "test class0 number: ", totalC0Num - trainC0Num
        print "test class1 number: ", totalC1Num - trainC1Num
    #    print "trainC0Index: ", trainC0Index
    #    print "trainC1Index: ", trainC1Index
    #    print "testC0Index: ", testC0Index
    #    print "testC1Index: ", testC1Index

        tmpData = totalClassData
        tmpLabels = totalClassLabels
        trainClassData = []
        trainClassLabels = []
        testClassData = []
        testClassLabels = []
        for index in xrange(len(tmpData)):
            if index in trainC0Index:
                trainClassData.append(tmpData[index])
                trainClassLabels.append(tmpLabels[index])
            if index in trainC1Index:
                trainClassData.append(tmpData[index])
                trainClassLabels.append(tmpLabels[index])
            if (index not in trainC0Index) and (index not in trainC1Index):
                testClassData.append(tmpData[index])
                testClassLabels.append(tmpLabels[index])

        # SVM
        from sklearn import svm
        trainClassData = np.mat(trainClassData)
        trainClassLabels = np.array(trainClassLabels)
#        trainClassLabels = np.mat(trainClassLabels)
        print "trainClassData.shape: ", trainClassData.shape
        print "trainClassLabels.shape: ", trainClassLabels.shape
        clf = svm.SVC(C = 1, tol = 1e-10, kernel = "linear")
        clf.fit(trainClassData, trainClassLabels)

#        prediction = clf.predict(testClassData)
#        print prediction
        
        # test
        error = 0
        for index in xrange(len(testClassData)):
            tmp = testClassData[index]
            tmp = np.array(tmp)
#            print tmp
#            print tmp.shape
            prediction = clf.predict(tmp)
#            print prediction
    #        print prediction
            if prediction >= 0.5:
                prediction = 1
            else:
                prediction = 0
            if prediction != testClassLabels[index]:
#                print testClassData[index][:5]
                error += 1
        accuracy += (1.0 - float(error) / float(len(testClassData))) * 100
        print "error: ", error
        print "accuracy: %f%%" % (accuracy / (testIt + 1))
        

        
    
    
                

    




















def find(word): # test
    fr = open(PROJECT_DOC_DIR + PROJECT_RAWFEAT_DIR, "r")
    count = 0
    for line in fr:
        count += 1
        if word in line.strip().split("\t")[0]:
            print line.strip().split("\t")[0]
    
def createWave():
    data = ones(1000000)
    for i in xrange(len(data) - 1):
	data[i + 1] = data[i] * 1.00001 * math.pow(-1, (i % 2))
    return data

def showData(fileName, start, end):
    global PROJECT_DIR
    fileName = fileName.decode("UTF-8")
    path = PROJECT_DIR + fileName
    f = wave.open(path, "rb")
    # 读取格式信息
    # (nchannels, sampwidth, framerate, nframes, comptype, compname)
    params = f.getparams()
    nchannels, sampwidth, framerate, nframes, comtype, comname = params[:6]
    print (nchannels)
    print (sampwidth)
    print (framerate)
    print (nframes)
    print (comtype)
    print (comname)
    # 读取波形数据
    str_data = f.readframes(nframes)
    str_data = np.fromstring(str_data, dtype = np.short)
    f.close()
    print(str_data[start : end])
    return nchannels, sampwidth, framerate, nframes, comtype, comname

def loadMusic(fileName, param = "null"):
    global PROJECT_DIR
    fileName = fileName.decode("UTF-8")
    path = PROJECT_DIR + fileName
    f = wave.open(path, "rb")
    # 读取格式信息
    # (nchannels, sampwidth, framerate, nframes, comptype, compname)
    params = f.getparams()
    nchannels, sampwidth, framerate, nframes = params[:4]
    # 读取波形数据
    str_data = f.readframes(nframes)
    f.close()
    #将波形数据转换为数组
    wave_data = np.fromstring(str_data, dtype = np.short)
    if (param == "w"):
        fw = open(pathTxt, "w")
        i = 0
        while (i < len(wave_data)):
            temp = (str(wave_data[i]))
            temp += " "
            fw.write(temp)
            i += 1
        print "Write finish"
        fw.close()
    return nchannels, sampwidth, framerate, nframes, len(wave_data), wave_data

def writeWAVFile(data, channels, fileName):
    global PROJECT_DIR
    global params
#    global data

#    getData(fileName)
    path = PROJECT_DIR + fileName + u"_custom.wav"
    if os.path.exists(path):
        os.remove(path)
        print "same file detected, delete"
    print "data length in WRITE WAVE: ", (len(data))
#    dataToWrite = np.arange(0, len(data))
#    for i in xrange(len(data)):
#        dataToWrite[i] = data[i]
    fw = wave.open(path, "w")
#    nframes = len(data) / 2
    params[0] = channels
    params[3] = int(len(data) / channels)
    fw.setparams(params)
    fw.writeframes(data.tostring())
    fw.close()

def readTxt():
#    path = u"H:\ccc.txt"
    path = PROJECT_DIR + u"MP3Frames.txt"
    fr = open(path, "rb")
    for line in fr:
        print line
    fr.close()

def Hilbert(fileName):
    global data
    global PROJECT_DIR
    
    unicodeFlag = isinstance(fileName, unicode)
    #人为输入文件名
    if (unicodeFlag == False):
        path, fileName = getPath(fileName)
    else:
        #读取程序返回的unicode编码的路径
        path = fileName
        fileName = fileName.strip().split("\\")[-1]
        
    data = getData(fileName)
    ldata = np.zeros(len(data) / 2, dtype = np.short)
    for i in xrange(len(ldata)):
        ldata[i] = data[i * 2]
    t1 = time.time()
    Hdata = fftpack.hilbert(ldata)
    t2 = time.time()
    print "hilbert takes %fs." % (t2 - t1)
    for i in xrange(len(ldata)):
        Hdata[i] = sqrt(int(ldata[i]) * int(ldata[i]) + int(Hdata[i]) * int(Hdata[i]))
    writeWAVFile(Hdata)

def draw(level):
    if level == 1:
        fr = open("H:\\var_L.txt", "r")
        num_L = []
        for line in fr:
            line = line.replace("\n", "")
            line = string.atof(line)
            num_L.append(line)
        fr = open("H:\\var_R.txt", "r")
        num_R = []
        for line in fr:
            line = line.replace("\n", "")
            line = string.atof(line)
            num_R.append(line)
        fr.close()
    
        pl.subplot(111)
        pl.plot(num_L, alpha = 0.8)
        pl.plot(num_R, color = "red", alpha = 0.8)
        pl.xlabel("")
        pl.ylabel("var")
        pl.show()

    if level == 2:
        fr = open("H:\\varVar_L.txt", "r")
        num_L = fr.readline()
        num_L = num_L.strip().split("\t")
#        fr = open("H:\\varVar_R.txt", "r")
#        num_R = []
#        for line in fr:
#            line = line.replace("\n", "")
#            line = string.atof(line)
#            num_R.append(line)
        fr.close()

        x = np.linspace(0, 5, len(num_L))

        pl.figure(figsize=(8,4))
        pl.subplot(111)
        ylim(0, 400)
        pl.plot(x, num_L, alpha = 1, label = "Stability Vector")
        legend(loc = "upper right")
#        pl.plot(num_R, color = "red", alpha = 0.8)
        pl.xlabel("Time(s)")
        pl.ylabel("Stability")
        pl.show()

    if level == 2.1:
        data = getData("test_stable.wav")
        data.shape = -1, 2
        data = data.T
        data = data[0]

        x = np.linspace(0, 5, len(data))

        pl.figure(figsize=(8,4))
        pl.subplot(111)
        ylim(-15000, 15000)
        pl.plot(x, data, alpha = 1, label = "Music Signal")
        legend(loc = "upper right")
        pl.xlabel("Time(s)")
        pl.ylabel("Amplitude")
        pl.show()

    if level == 2.2:
#        data = getData(PROJECT_TEST_DIR + u"100101@806_Excuse Me_highPass(0.9).wav")
#        data = getData(PROJECT_TEST_DIR + u"100101@806_Excuse Me_highPass(0.9)_beatFilter(False).wav")
        data = getData(PROJECT_TEST_DIR + u"100101@806_Excuse Me_highPass(0.9)_beatFilter(False)_countBeats(3528-False).wav")
#        data = getData(PROJECT_TEST_DIR + u"100191@Interlude_highPass(0.9).wav")
#        data = getData(PROJECT_TEST_DIR + u"100191@Interlude_highPass(0.9)_beatFilter(False).wav")
        data.shape = -1, 2
        data = data.T
        data = data[0]
        data = data[7938000:9040500]
#        data = data[529200:1631700]

        x = np.linspace(0, 25, len(data))

        pl.figure(figsize=(8,4))
        pl.subplot(111)
#        ylim(-32768, 32768)
#        ylim(-10000, 10000)
        ylim(0, 32768)
#        pl.plot(x, data, alpha = 1, label = "Music Signal")
#        pl.plot(x, data, alpha = 1, label = "Signals in Beat-area")
        pl.plot(x, data, alpha = 1, label = "beats")
        legend(loc = "upper right")
        pl.xlabel("Time(s)")
        pl.ylabel("Amplitude")
        pl.show()

    if level == 3:
        fr = open("H:\\HvarVar_L.txt", "r")
        num_L = []
        for line in fr:
            line = line.replace("\n", "")
            line = string.atof(line)
            num_L.append(line)
        fr = open("H:\\HvarVar_R.txt", "r")
        num_R = []
        for line in fr:
            line = line.replace("\n", "")
            line = string.atof(line)
            num_R.append(line)
        fr.close()
    
        pl.subplot(111)
        pl.plot(num_L, alpha = 0.8)
        pl.plot(num_R, color = "red", alpha = 0.8)
        pl.xlabel("")
        pl.ylabel("HvarVar")
        pl.show()


def outputDwtMusic(fileName, channel, level):
    data = getData(fileName)
    dwtData = np.zeros(0)
    tmpData = np.zeros(0)
    if channel == 2:
        data.shape = -1, 2
        data = data.T
        data = data[0]
    tempName = fileName.strip().split(".")[0]

    tmpData = data
    tmpData.shape = -1, len(tmpData)
    Clen = 0
    for dwtLevel in xrange(level):
        for i in xrange(len(tmpData)):
            CA, CD = dwt(tmpData[i], "db4")
            CA = np.short(CA)
            CD = np.short(CD)
            Clen = len(CA)
            dwtData = np.append(dwtData, CA)
            dwtData = np.append(dwtData, CD)
        dwtData.shape = -1, Clen
        tmpData = dwtData
        dwtData = np.zeros(0)
    dwtData = tmpData
    print dwtData.shape
    dwtExtData = np.zeros(0)
    tmpExtData = np.zeros(0)
    
    for i in xrange(len(dwtData)):
        for lev in xrange(pow(2, level)):
            tmpExtData = np.append(tmpExtData, dwtData[i])
        tmpExtData.shape = -1, Clen
        tmpExtData = tmpExtData.T
        tmpExtData = np.concatenate(tmpExtData)
        dwtExtData = np.append(dwtExtData, tmpExtData)
        tmpExtData = np.zeros(0)
    dwtExtData.shape = -1, Clen * pow(2, level)
    
        
    freqRange = 22100 / pow(2, level)
    for i in xrange(pow(2, level)):
#        if abs(freqRange * i - 80) < 2 or \
#           abs(freqRange * i - 690) < 2 or \
#           abs(freqRange * i - 1000) < 2 or \
#           abs(freqRange * i - 4500) < 2 or \
#           abs(freqRange * i - 16000) < 2 or \
#           abs(freqRange * i - 18000) < 2 :
        for lev in xrange(pow(2, level)):
            tmpExtData = np.append(tmpExtData, dwtData[i])
        tmpExtData.shape = -1, Clen
        tmpExtData = tmpExtData.T
        tmpExtData = np.concatenate(tmpExtData)
        num = str(freqRange * i) + "-" + str(freqRange * (i + 1))
        writeWAVFile(short(tmpExtData), 1, tempName.decode("gb18030") + num)
        tmpExtData = np.zeros(0)
    



    
