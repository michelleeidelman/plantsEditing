import pandas as pd
from scipy.stats import mannwhitneyu

df = pd.read_csv("mergeTable")
grouped = df.groupby("analysis_group")
control = grouped.get_group("Control")
case = grouped.get_group("Case")

#mannwhitneyu test:
mannwhitneyuResult = mannwhitneyu(control['value'], case['value'])
print(str(mannwhitneyuResult))

a=3