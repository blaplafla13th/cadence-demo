# coding=utf-8
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''Utility functions for GIP.

(1) normalization: MinMax Normalizer
(2) renormalization: Recover the data from normalzied data
(3) rounding: Handlecategorical variables after imputation
(4) rmse_loss: Evaluate imputed data in terms of RMSE
(5) binary_sampler: sample binary random variables
(6) uniform_sampler: sample uniform random variables
(7) sample_batch_index: sample random batch index
'''

# Necessary packages
import numpy as np
import pandas as pd

def normalization (data, parameters=None):
  '''Normalize data in [0, nm_alleles - 1] range.

  Args:
    - data: original data

  Returns:
    - norm_data: normalized data
    - norm_parameters: min_val, max_val for each feature for renormalization
  '''

  # Parameters
  _, dim = data.shape
  norm_data = data.copy()

  if parameters is None:

    # MixMax normalization
    min_val = np.zeros(dim)
    max_val = np.zeros(dim)

    # For each dimension
    for i in range(dim):
      min_val[i] = np.nanmin(norm_data[:,i])
      norm_data[:,i] = norm_data[:,i] - np.nanmin(norm_data[:,i])
      max_val[i] = np.nanmax(norm_data[:,i])
      norm_data[:,i] = norm_data[:,i] / (np.nanmax(norm_data[:,i]) + 1e-6)

      # Return norm_parameters for renormalization
    norm_parameters = {'min_val': min_val,
                       'max_val': max_val}

  else:
    min_val = parameters['min_val']
    max_val = parameters['max_val']

    # For each dimension
    for i in range(dim):
      norm_data[:,i] = norm_data[:,i] - min_val
      norm_data[:,i] = norm_data[:,i] / (max_val + 1e-6)

    norm_parameters = parameters

  return norm_data, norm_parameters


def renormalization (norm_data, norm_parameters):
  '''Renormalize data from [0, 1] range to the original range.

  Args:
    - norm_data: normalized data
    - norm_parameters: min_val, max_val for each feature for renormalization

  \Returns:
    - renorm_data: renormalized original data
  '''

  min_val = norm_parameters['min_val']
  max_val = norm_parameters['max_val']

  _, dim = norm_data.shape
  renorm_data = norm_data.copy()
  try:
    for i in range(dim):
      renorm_data[:,i] = renorm_data[:,i] * (max_val[i] + 1e-6)
      renorm_data[:,i] = renorm_data[:,i] + min_val[i]
  except:
    for i in range(dim):
      renorm_data[:,i] = renorm_data[:,i] * (max_val + 1e-6)
      renorm_data[:,i] = renorm_data[:,i] + min_val

  return renorm_data


def rounding (imputed_data, data_x):
  '''Round imputed data for categorical variables.

  Args:
    - imputed_data: imputed data
    - data_x: original data with missing values

  Returns:
    - rounded_data: rounded imputed data
  '''

  _, dim = data_x.shape
  rounded_data = imputed_data.copy()

  for i in range(dim):
    temp = data_x[~np.isnan(data_x[:, i]), i]
    # Only for the categorical variable
    if len(np.unique(temp)) < 20:
      rounded_data[:, i] = np.round(rounded_data[:, i])

  return rounded_data


def rmse_loss (ori_data, imputed_data, data_m):
  '''Compute RMSE loss between ori_data and imputed_data

  Args:
    - ori_data: original data without missing values
    - imputed_data: imputed data
    - data_m: indicator matrix for missingness

  Returns:
    - rmse: Root Mean Squared Error
  '''

  ori_data, norm_parameters = normalization(ori_data)
  imputed_data, _ = normalization(imputed_data, norm_parameters)

  # Only for missing values
  nominator = np.sum(((1-data_m) * ori_data - (1-data_m) * imputed_data)**2)
  denominator = np.sum(1-data_m)

  rmse = np.sqrt(nominator/float(denominator))

  return rmse

def binary_sampler(p, rows, cols):
  '''Sample binary random variables.

  Args:
    - p: probability of 1
    - rows: the number of rows
    - cols: the number of columns

  Returns:
    - binary_random_matrix: generated binary random matrix.
  '''
  np.random.seed(7)
  unif_random_matrix = np.random.uniform(0., 1., size = [rows, cols])
  binary_random_matrix = 1*(unif_random_matrix < p)
  return binary_random_matrix


def uniform_sampler(low, high, rows, cols):
  '''Sample uniform random variables.

  Args:
    - low: low limit
    - high: high limit
    - rows: the number of rows
    - cols: the number of columns

  Returns:
    - uniform_random_matrix: generated uniform random matrix.
  '''
  np.random.seed(7)
  return np.random.uniform(low, high, size = [rows, cols])


def sample_batch_index(total, batch_size):
  '''Sample index of the mini-batch.

  Args:
    - total: total number of samples
    - batch_size: batch size

  Returns:
    - batch_idx: batch index
  '''
  np.random.seed(7)
  total_idx = np.random.permutation(total)
  batch_idx = total_idx[:batch_size]
  return batch_idx

def get_data(X, miss_rate):
  # Parameters
  no, dim = X.shape

  # Introduce missing data
  data_m = binary_sampler(1-miss_rate, no, dim)
  miss_data_x = X.copy()
  miss_data_x[data_m == 0] = ".|."
  return X, miss_data_x

def save_hap(data, output_path):
  cols = list(range(data.shape[1]))
  df = pd.DataFrame(data = data, columns = cols)
  df.to_csv(output_path, header = False, index = False, sep ="\t", mode='a')