import matplotlib.pylab as plt
import numpy as np
from scipy.linalg import svd
from scipy.fft import fft
import pywt
import read
import sys

plt.figure(dpi=50)

def noiseFilterSVD(data):
  resSize = int(np.sqrt(len(data)))
  arr = data[:resSize * resSize]
  res = np.array(arr).reshape(resSize, resSize)
  U, Sig, Vh = svd(res)
  SigList = list(Sig)
  RESERVE = 2
  for i in range(RESERVE, len(SigList)):
    SigList[i] = 0.0
  SigNew = np.mat(np.diag(SigList))
  res = np.array(U * SigNew * Vh)
  return np.squeeze(res.reshape(resSize * resSize, -1), axis = (1)).tolist()

def noiseFilterAvg(data):
  Wid = 300
  res = []
  avg = np.average(data[:Wid])
  for i in range(Wid, len(data)):
    avg += data[i] / Wid
    res.append(avg)
    avg -= data[i - Wid] / Wid
  return res

def noiseFilterWavelet(data):
  typ = 'db8'
  wavelet = pywt.Wavelet(typ)
  maxlev = pywt.dwt_max_level(len(data), wavelet.dec_len)
  threshold = 0.2
  coeffs = pywt.wavedec(data, typ, level = maxlev)
  for i in range(len(coeffs)):
    coeffs[i] = pywt.threshold(coeffs[i], threshold * max(coeffs[i]))
  res = pywt.waverec(coeffs, typ)
  return res.tolist()


def dealbyDistance(data, distance, func): # len(data) must equal to distance
  point = 0
  res = []
  for i in range(1,len(distance)):
    if distance[i] != distance[i-1]:
      res = res + (func(data[point: i - 1]))
      point = i
  res = res + (func(data[point:]))
  return res

def displaySpectrum(x, t): # 显示语谱图
  data = np.array(x)
  data = data * 1.0 / max(data)

  plt.subplot(2,2,1+2*t)
  plt.plot(data)
  plt.title("wave")

  sr = 4
  ft = fft(data)
  framesize = int(200 * sr)
  magnitude = np.absolute(ft)  # 对fft的结果直接取模（取绝对值），得到幅度magnitude
  frequency = np.linspace(0, sr, len(magnitude))  # (0, 16000, 121632)
  plt.subplot(2,2,2+2*t)
  plt.specgram(data, NFFT = framesize, Fs = sr, window = np.hanning(M = framesize))
  plt.title("yupu")

def getPercentage(x, y):
  return np.abs((x - y) / y)

def adjustAvg(data, avg):
  sum = 0
  cnt = 0
  for i in range(len(data)):
    if getPercentage(data[i], avg) < 0.01:
      sum += data[i]
      cnt += 1
  return sum / cnt
  

def getDefectList(func, data, dis):
  data = func(data)
  avg = np.average(data)
  for i in range(0, 3):
    avg = adjustAvg(data, avg)
  print(avg)
  res = {}
  for i in range(len(data)):
    if getPercentage(data[i], avg) >= 0.01:
      res[dis[i]] = max(res.get(dis[i], 0), getPercentage(data[i], avg))
  return res

def removeStartEnd(val):
  todo = []
  for i in range(len(val)):
    last = val["Distance"][i]
  for i in range(len(val)):
    if (1 <= val["Distance"][i] and val["Distance"][i] <= 4) or (last - 5 <= val["Distance"][i] and val["Distance"][i] <= last):
      todo.append(i)
  return val.drop(todo, axis = 0)




def getDefectFromOneFile(obj, trans):
  columnName = ['A', 'B', 'C', 'D', 'E', 'F']
  # displaySpectrum(obj['D'].values.tolist(), 0)
  obj = removeStartEnd(obj)
  obj.reset_index(drop = True, inplace = True)
  # clearData = dealbyDistance(obj['A'], obj['Distance'], noiseFilterSVD)
  res = {}
  for i in columnName:
    res[i] = getDefectList(noiseFilterAvg, obj[i] if trans else obj[i].T.T, obj['Distance'])
  return res
  # print(ansDict)
  # print(len(ansDict))
  # displaySpectrum(clearData, 1)
  # plt.show()

def getAllDefect():
  fileList = ["20210517221359_0.kwd", "20210517221542_1.kwd", "20210517221730_0.kwd", "20210517221935_1.kwd", "20210517222122_0.kwd", "20210517222309_1.kwd", "20210517222453_0.kwd", "20210517222636_1.kwd", "20210517222818_0.kwd", "20210517222936_1.kwd"]
  res = []
  for i in range(len(fileList)):
   obj = read.getFromFile(fileList[i])
   res.append(getDefectFromOneFile(obj, i % 2 == 1))
  return res

ori_stdout = sys.stdout
f = open("./defect.txt", "w")
sys.stdout = f

print(getAllDefect())

sys.stdout = ori_stdout
f.close()