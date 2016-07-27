```
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
%matplotlib inline
```
#read in data
```
motif_scores_df = pd.read_csv('/home/shs038/veh_kla/motif_scores_C57BL6J.tsv', sep='\t')
motif_scores_df.index = motif_scores_df['ID'].values
del motif_scores_df['ID']
```
#vehicle scores
```
veh_df=motif_scores_df[motif_scores_df['Factors'].str.contains('c57bl6_atac_veh')]
del veh_df['Factors']
del veh_df['chr']
```
#KLA scores 
```
kla_df=motif_scores_df[motif_scores_df['Factors'].str.contains('c57bl6_atac_kla')]
del kla_df['Factors']
del kla_df['chr']
```
#normalized by max 
```
max_veh=veh_df.max(axis=0)
max_kla=kla_df.max(axis=0)
veh_normalized_df=veh_df/max_veh
kla_normalized_df=kla_df/max_kla
```
#remove negative value and NaN
```
veh_normalized_df[veh_normalized_df<0]=0
kla_normalized_df[kla_normalized_df<0]=0
veh_normalized_plot=np.nan_to_num(veh_normalized_df)
kla_normalized_plot=np.nan_to_num(kla_normalized_df)
```
#plot normalized motif scores
```
sns.distplot(veh_normalized_plot[veh_normalized_plot!=0])
plt.ylabel('Frequency')
plt.xlabel('motif scores')
plt.title('veh')
plt.show()
sns.distplot(kla_normalized_plot[kla_normalized_plot!=0])
plt.ylabel('Frequency')
plt.xlabel('motif scores')
plt.title('kla')
```
#function to calculate correlation coefficient 
```
def find_correlation(df):
    '''
    input: a dataframe contains all motifs scores
    output: a dataframe contains correlation coefficient of each motif pai
    '''
    motifs = df.columns.values
    correlation_df=np.zeros((df.shape[1],df.shape[1]),dtype=np.float)
    correlation_df=pd.DataFrame(correlation_df, columns=motifs)
    correlation_df.index=motifs
    for i in range(df.shape[1]-1):
        for j in range(i+1,df.shape[1]):
            motif_paris=df.iloc[:,[i,j]]#select two moitfs
            #remove data that two motfis do not co-occur
            motif_paris=motif_paris[motif_paris.ix[:,0]!=0]
            motif_paris=motif_paris[motif_paris.ix[:,1]!=0]
            #calculate correlation
            coef=np.corrcoef(motif_paris.ix[:,0],motif_paris.ix[:,1])
            correlation_df.ix[i,j]=coef[0,1]
    #reshape dataframe
    Pairs=[]
    Correlation=[]
    #loop in part of count data that contain meaningful correlation
    for i in range (correlation_df.shape[1]-1):
        for j in range (i+1,correlation_df.shape[1]):
            #put motif pair and correlation into the empty list
            motif_pairs=(motifs[i],motifs[j])
            Pairs.append(motif_pairs)
            Correlation.append(correlation_df.ix[i,j])
    #reshape the dataframe
    reshaped_df=pd.DataFrame({'Correlation': Correlation}, index=Pairs)
    return reshaped_df
```
#calculate correlation coefficient 
```
veh_correlation=find_correlation(veh_normalized_df)
kla_correlation=find_correlation(kla_normalized_df)
```
#compare correlation coefficient of motif paris under veh and kla treatment
```
correlation_df=pd.concat([veh_correlation, kla_correlation], axis=1)
correlation_df.columns = ['veh', 'kla']
correlation_df['vel-kla']=correlation_df['veh']-correlation_df['kla']
```
#plot difference
```
sns.distplot(correlation_df['vel-kla'])
plt.ylabel('Frequency')
plt.xlabel('vel-kla')
plt.title('Motif Correlation Difference under vel and kal') 
```
#show significant difference
```
correlation_df[(correlation_df['vel-kla'].abs() > 0.2) ].sort('vel-kla', ascending=False)
```
