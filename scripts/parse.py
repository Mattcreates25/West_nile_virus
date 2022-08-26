# import the panda library
import pandas as pd

#read in the metadata csv
df = pd.read_csv('results/all_metadata.csv')

#check the format
df.head(15)

#create a variable with split data
new_metadata = df['strain'].str.split('.', n=1, expand=True)

#rename the columns
new_meta = new_metadata.rename(columns={0:'strain',1:'name'})

#export as a csv
new_meta.to_csv('newmeta.csv',index=False)


