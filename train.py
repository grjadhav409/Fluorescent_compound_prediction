
Models={}

def train(X,Y):

    return df # accuracy , prec, recall, MCC 

for model in models["key"]:
    for i in desc_dfs:
        X = i.drop[["target"]]
        Y = i[[target]]
        score_df = model(X,Y)
        score_df.save_csv(PATH)