import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#reading the csv file into a dataframe
df = pd.read_csv('SalmonTPM.wideSamples.csv')
df = df.transpose()
#deleting the first row of indexes
new_header = df.iloc[0] #grab the first row for the header
df = df[1:] #take the data less the header row
df.columns = new_header #set the header row as the df header
#adding a column of the indexes of the row that represent the samples
df['Sample'] = df.index
#changing the view of the dataframe to be melt
df_new = pd.melt(df,id_vars =['Sample'],value_vars = ['ADARB1','ADARB2','ADAR','ADARB2-AS1'])
#reading a csv with a table that represent our two grups
df_sra = pd.read_csv('SraRunTableTestSyndrome.csv')
#selecting only the 2 relevant columns
df_sra_new = df_sra.loc[:,["Run","analysis_group"]]
#changing the name of the column
df_sra_new = df_sra_new.rename(columns={'Run': 'Sample'})
#merging the two tables to 1 by the Sanple column
merge_df = pd.merge(df_new,df_sra_new,on='Sample')
#saving the merge dataframe in a csv file
merge_df.to_csv("mergeTable", index = False)
#creating a boxplot
ax = sns.boxplot(x="analysis_group", y="value", data=merge_df)
g = sns.FacetGrid(merge_df, col="GeneSymbol")
g.map(sns.boxplot, "analysis_group", "value", palette="bright")
sns.color_palette("pastel")
plt.show()
