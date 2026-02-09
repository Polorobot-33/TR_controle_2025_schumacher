import pandas as pd

raw_data = pd.read_csv("comp_figures/results.csv", sep=",")
raw_fail_mask = (raw_data["nb_iter"] >= 1000) | (raw_data["Time"] >= 20)
raw_success_data = raw_data.drop(raw_data[raw_fail_mask].index)

def arange_data(table) :
    return pd.pivot_table(table, values=["Time","nb_iter","calc_time"], index=["case_name","model_name","collision","collocation","nh"])

data = arange_data(raw_success_data)

#print(data.head())

# averages
print("Mean calculation time by collision model :")
print(data.groupby("collision")["calc_time"].agg("mean"))


print("\nMean number of iterations by collision model :")
print(data.groupby("collision")["nb_iter"].agg("mean"))


print("\nMean final time by collision model :")
print(data.groupby("collision")["Time"].agg("mean"))

print("Estimated number of failed optimizations by collision model :")
print(arange_data(raw_data[raw_fail_mask]).groupby("collision")["Time"].count())