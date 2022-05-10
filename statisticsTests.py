import pandas as pd
from scipy.stats import mannwhitneyu

# seperating the data by the two groups - control and case
df = pd.read_csv("mergeTable")
grouped = df.groupby("analysis_group")
control = grouped.get_group("Control")
case = grouped.get_group("Case")

# seperating the data by genes:
# for the control group:
grouped1 = control.groupby("GeneSymbol")
ADARB1_1 = grouped1.get_group("ADARB1")
ADARB2_1 = grouped1.get_group("ADARB2")
ADAR_1 = grouped1.get_group("ADAR")
ADARB2_AS1_1 = grouped1.get_group("ADARB2-AS1")

# for the case group:
grouped2 = case.groupby("GeneSymbol")
ADARB1_2 = grouped2.get_group("ADARB1")
ADARB2_2 = grouped2.get_group("ADARB2")
ADAR_2 = grouped2.get_group("ADAR")
ADARB2_AS1_2 = grouped2.get_group("ADARB2-AS1")


listOfControlGenes = [ADARB1_1, ADARB2_1, ADAR_1, ADARB2_AS1_1]
listOfCaseGenes = [ADARB1_2, ADARB2_2, ADAR_2, ADARB2_AS1_2]
listOfGenes = ["ADARB1", "ADARB2","ADAR","ADARB2_AS1"]
listOfTests = ["mannwhitneyu"] * 4
listOfGroups = ["control-case"] * 4



statistics = []
pval =[]

# mannwhitneyu test:
for i in range(len(listOfControlGenes)):
    mannwhitneyuResult = mannwhitneyu(listOfControlGenes[i]['value'], listOfCaseGenes[i]['value'])
    statistics.append(mannwhitneyuResult[0])
    pval.append(mannwhitneyuResult[1])

data = {'Gene':listOfGenes,'Test':listOfTests, 'groups':listOfGroups,'statistic':statistics, 'Pvalue':pval}

# creating a new dataframe with the tests results
new_df = pd.DataFrame(data, columns = ['Gene', 'Test', 'groups', 'statistic', 'Pvalue'])
new_df.to_csv('statisticAnalysis',index=False)

a=3