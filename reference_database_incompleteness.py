#___REFERENCE DATABASE INCOMPLETENESS EXPERIMENT____

from sklearn.model_selection import train_test_split
from google.colab import drive
drive.mount("/content/drive")
import xml.etree.ElementTree as ET
!pip install xmltodict
!pip install Bio
import json
import xmltodict
from Bio.Blast import NCBIXML
import pandas as pd
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import pyplot as plt
import _pickle as cPickle
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt

gbif_diversity = pd.read_excel("gbif.diversity.xlsx")
gbif_diversity.to_csv ("gbif_diversity.csv",
                  index = None,
                  header=True)
gbif_diversity = pd.DataFrame(pd.read_csv("gbif_diversity.csv"))

otu_to_sample = pd.read_csv('otu.sample.pivot.txt',  delim_whitespace=True, header=None, on_bad_lines='warn', index_col = 0)
otu_to_sample = otu_to_sample.transpose()
otu_to_sample.sort_values("OTU", axis = 0, ascending=True, inplace = True, na_position = "first")

final_df = pd.merge(gbif_diversity, otu_to_sample, left_on='sample',right_on='OTU')
final_df = final_df.drop("OTU", axis=1)
final_df = final_df.drop("sample", axis=1)
original_header = set()
final_df = final_df.dropna()

df = pd.read_csv('blast_output_yay.csv')
pct = list(df["Bitscore"])
pct = sorted(pct)

ranges = [y for y in range(len(pct))]
plt.plot(ranges, pct)
plt.plot(ranges, pct)

gbif_diversity = pd.read_excel("gbif.diversity.xlsx")
gbif_diversity.to_csv ("gbif_diversity.csv",
                  index = None,
                  header=True)
gbif_diversity = pd.DataFrame(pd.read_csv("gbif_diversity.csv"))

otu_to_sample = pd.read_csv('otu.sample.pivot.txt',  delim_whitespace=True, header=None, on_bad_lines='warn', index_col = 0)
otu_to_sample = otu_to_sample.transpose()
otu_to_sample.sort_values("OTU", axis = 0, ascending=True, inplace = True, na_position = "first")

print(gbif_diversity)
print(otu_to_sample)

final_df = pd.merge(gbif_diversity, otu_to_sample, left_on='sample',right_on='OTU')
final_df = final_df.drop("OTU", axis=1)
final_df = final_df.drop("sample", axis=1)
final_df = final_df.dropna()
final_df

all_otus = set(final_df.columns[3:])

#find quartiles
import numpy

values = df["Bitscore"]
quartiles = numpy.percentile(values, [15, 85])
print(quartiles)

#get OTU splits

results = []
names = set()
i=0
OTU_split = [[], [], [], [], []]
for index, row in df.iterrows():
  i+=1
  if row["Query"] in all_otus:
    if row["Bitscore"] <= quartiles[0]:
      OTU_split[0].append(row["Query"])
    if row["Bitscore"] > quartiles[0]:
      OTU_split[1].append(row["Query"])
    if row["Bitscore"] <= quartiles[1]:
      OTU_split[2].append(row["Query"])
    if row["Bitscore"] > quartiles[1]:
      OTU_split[3].append(row["Query"])
    OTU_split[4].append(row["Query"])
print(OTU_split)
print([len(x) for x in OTU_split])

#train models, store data

def trainModel(final_df, use_otus, name, num, x_train, x_test, y_train, y_test):
  drop_otus = list(set(final_df.columns[3:]) - set(use_otus))
  x_test = x_test.drop(columns=drop_otus)
  x_train = x_train.drop(columns=drop_otus)
  model = RandomForestRegressor(n_estimators=1000, min_samples_leaf = 10, max_features = int(len(x_train.columns)/3), random_state=17)
  model.fit(x_train, y_train)
  res = evaluate(name, x_test, y_test, model)
  print(num+1, res)
  with open("/content/drive/MyDrive/aar_results_really_final/" + name, 'wb') as f:
    cPickle.dump(model, f)
    print(name)
  return [model, x_test, y_test, res]

def evaluate(name, x_test, y_true, model_info):
  rf = model_info[0]
  y_pred = [x for x in rf.predict(x_test)]
  res = dict()
  res["mean absolute error"] = metrics.mean_absolute_error(y_true, y_pred)
  res["mean squared error"] = metrics.mean_squared_error(y_true, y_pred)
  res["root mean squared error"] = metrics.mean_squared_error(y_true, y_pred, squared=False)
  res["mean absolute percent error"] = metrics.mean_absolute_percentage_error(y_true, y_pred)
  res["explained variance score"] = metrics.explained_variance_score(y_true, y_pred)
  res["max error"] = metrics.max_error(y_true, y_pred)
  res["median error"] = metrics.median_absolute_error(y_true, y_pred)
  res["r2"] = metrics.r2_score(y_true, y_pred)
  return res

data = dict()
for j in range(50):
  data[j] = dict()
  x_train, x_test, y_train, y_test = train_test_split(final_df.drop(['GBIF richness (family)', 'GBIF richness (species)', 'Kreft & Jetz index'],axis='columns'), final_df["Kreft & Jetz index"], test_size = 0.2)
  for i in range(len(OTU_split)):
    split = OTU_split[i]
    data[j][i] = dict()
    data[j][i] = trainModel(final_df, split, "trial"+str(j)+"model"+str(i), i+1, x_train, x_test, y_train, y_test)
  cPickle.dump(data, open("/content/drive/MyDrive/aar_results_really_final/saved_dictionary.pkl", "wb"))

#data analysis - original ("good results")
import _pickle as cPickle
from google.colab import drive
drive.mount("/content/drive")

final_values = []

with open('/content/drive/MyDrive/aar_results_really_final/saved_dictionary.pkl', 'rb') as handle:
    smaller = cPickle.load(handle)