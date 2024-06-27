import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import VotingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from xgboost import XGBRegressor
from sklearn.neural_network import MLPRegressor

gbif_diversity = pd.read_excel("gbif.diversity.xlsx")
gbif_diversity.to_csv ("gbif_diversity.csv",
                  index = None,
                  header=True)
gbif_diversity = pd.DataFrame(pd.read_csv("gbif_diversity.csv"))

otu_to_sample = pd.read_csv('otu.sample.pivot.txt',  delim_whitespace=True, header=None, on_bad_lines='warn', index_col = 0)
otu_to_sample = otu_to_sample.transpose()
otu_to_sample
otu_to_sample.sort_values("OTU", axis = 0, ascending=True, inplace = True, na_position = "first")

print(gbif_diversity)
print(otu_to_sample)

final_df = pd.merge(gbif_diversity, otu_to_sample, left_on='sample',right_on='OTU')
final_df = final_df.drop("OTU", axis=1)
final_df = final_df.drop("sample", axis=1)
final_df = final_df.dropna()
final_df

import numpy as np
from sklearn.decomposition import PCA
import pandas as py

pca = PCA(n_components=100)
principal_components = pca.fit_transform(final_df.iloc[:, 3:])
print(pca.explained_variance_ratio_.cumsum())
principal_df = pd.DataFrame(data = principal_components)
principal_df.insert(0, "num", [i for i in range(len(principal_df))], True)
print(principal_df)

outputs = final_df.iloc[:, 0:3]
outputs.insert(0, "num", [i for i in range(len(outputs))], True)
print(outputs)

final_df_really = pd.merge(outputs, principal_df, left_on="num",right_on="num")
final_df_really

from sklearn.ensemble import StackingRegressor
from sklearn.linear_model import LogisticRegression
from google.colab import drive
import _pickle as cPickle

drive.mount("/content/drive",force_remount =True)

data = dict()

for i in range(100):
  x_train, x_test, y_train, y_test = train_test_split(final_df_really.drop(['num', 'GBIF richness (family)', 'GBIF richness (species)', 'Kreft & Jetz index'],axis='columns'), final_df["Kreft & Jetz index"], test_size = 0.2)

  model = RandomForestRegressor(n_estimators=10000, min_samples_leaf = 5, max_features = int(len(x_train.columns)/3), random_state=17)
  model.fit(x_train, y_train)

  y_pred = [x for x in model.predict(x_test)]
  y_true = y_test

  from scipy.stats import spearmanr
  coef, p = spearmanr(y_true, y_pred)

  print('Mean Absolute Error (MAE):', metrics.mean_absolute_error(y_true, y_pred))
  print('Mean Squared Error (MSE):', metrics.mean_squared_error(y_true, y_pred))
  print('Root Mean Squared Error (RMSE):', metrics.mean_squared_error(y_true, y_pred, squared=False))
  print('Mean Absolute Percentage Error (MAPE):', metrics.mean_absolute_percentage_error(y_true, y_pred))
  print('Explained Variance Score:', metrics.explained_variance_score(y_true, y_pred))
  print('Max Error:', metrics.max_error(y_true, y_pred))
  print('Median Absolute Error:', metrics.median_absolute_error(y_true, y_pred))
  print('R^2:', metrics.r2_score(y_true, y_pred))
  print("Spearman:", coef)
  print("p-value:", p)
  y_true_df = pd.DataFrame(y_true)
  y_true_df.insert(1, "Predictions", y_pred, True)
  y_true_df.to_csv("results_test.csv")

  with open("/content/drive/MyDrive/aar_results_ideal_model_2/" + str(i), 'wb') as f:
    cPickle.dump(model, f)
  with open("/content/drive/MyDrive/aar_results_ideal_model_2/results_test" + str(i) + ".csv", 'wb') as f:
    y_true_df.to_csv(f)
  print(str(i))

  data[i] = [x_train, x_test, y_train, y_test, model, y_true, y_pred, coef, p]
cPickle.dump(data, open("/content/drive/MyDrive/aar_results_ideal_model_2/saved_dictionary.pkl", "wb"))

from sklearn.ensemble import StackingRegressor
from sklearn.linear_model import LogisticRegression
from google.colab import drive
import _pickle as cPickle

drive.mount("/content/drive")
with open('/content/drive/MyDrive/aar_results_ideal_model_2/saved_dictionary.pkl', 'rb') as handle:
    data = cPickle.load(handle)

ordered = [i+1 for i in range(67)]

import numpy as np
from scipy.stats import spearmanr

def generate_permutations(num_permutations, num_elements):
    permutations = []
    for _ in range(num_permutations):
        permutation = np.random.permutation(num_elements) + 1  # Adding 1 to shift range from 0-66 to 1-67
        permutations.append(permutation)
    return permutations

num_permutations = 100
num_elements = 67
permutations = generate_permutations(num_permutations, num_elements)

arr2 = []

for permutation in permutations:
  coef, p = spearmanr(ordered, permutation)
  arr2.append(coef)

sorted_data2 = np.sort(arr2)

# Calculate the CDF values
cdf = np.arange(1, len(sorted_data2) + 1) / len(sorted_data2)

# Plot the CDF
plt.plot(sorted_data2, cdf, marker='.', color='red', linestyle='none')
plt.xlabel('Value')
plt.ylabel('Cumulative Probability')
plt.title('Cumulative Distribution Function')
plt.grid(True)
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from statistics import mean

array = []
for i in range(100):
  array.append(data[i][7])
print(mean(array))
print(median(array))
print(max(array))

import numpy as np
import matplotlib.pyplot as plt

# Generate an array of 100 elements ranging from -1 to 1
# Sort the data
sorted_data = np.sort(array)

# Calculate the CDF values
cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
cdf = np.arange(1, len(sorted_data2) + 1) / len(sorted_data2)

# Plot the CDF
plt.plot(sorted_data, cdf, marker='.', color='orange', linestyle='none')
plt.plot(sorted_data2, cdf, marker='.', color='green', linestyle='none')

plt.xlabel('Value')
plt.ylabel('Cumulative Probability')
plt.title('Cumulative Distribution Function')
plt.grid(True)
plt.show()

