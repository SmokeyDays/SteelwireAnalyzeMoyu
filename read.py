from ast import Str
import pandas as pd
pd.__version__

def getFromFile(fileName):
  res = pd.read_csv("./data/" + fileName)
  res.columns = ['A', 'B', 'C', 'D', 'E', 'F', 'Direction', 'Distance', 'Time']
  return res