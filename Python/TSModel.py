# Example: python3 TSModel.py '../DataSets/TemperatureData.csv' 4206 1100 1
import sys

from pandas import read_csv
from pandas import DataFrame
from pandas import concat
import matplotlib.pyplot as plt
from statsmodels.tsa.ar_model import AutoReg
from sklearn.metrics import mean_squared_error
from math import sqrt
from numpy import count_nonzero
import time

series = read_csv(sys.argv[1], header=None, index_col=None)

values = DataFrame(series.values)
dataframe = concat([values.shift(1), values], axis=1)
dataframe.columns = ['t-1', 't+1']

N=int(sys.argv[2])
Lag=int(sys.argv[3])
ss=int(sys.argv[4])

X = series.values
X = X[1:len(X):ss]
train, test = X[1:N], X[N:]
t0=time.time()
model = AutoReg(train, lags=Lag-1)
model_fit = model.fit()
print('Running time for model computation:')
print(time.time()-t0)

print('Coefficients: %s' % model_fit.params)
print('Number of nonzero coefficients:')
print(count_nonzero(model_fit.params))
predictions = model_fit.predict(start=len(train), end=len(train)+len(test)-1, dynamic=False)
rmse = sqrt(mean_squared_error(test, predictions))
print('Test RMSE: %.3f' % rmse)
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(test)
ax1.plot(predictions, color='red')
ax2.stem(model_fit.params)
plt.show()
