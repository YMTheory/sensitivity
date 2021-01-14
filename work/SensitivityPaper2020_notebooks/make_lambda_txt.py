import pandas as pd

lambda_h5 = "/Users/czyz1/lc-nexouser/workdir/lambda/20_12_02_DNN1_024/critical_lambda_eres_9_resolution_0.008.h5"

df = pd.read_hdf(lambda_h5)
df.head()