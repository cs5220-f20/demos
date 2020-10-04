import sys
import pandas as pd
import statsmodels.api as sm
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')

def load_pingpong(fname):
    data = pd.read_csv(fname, sep=" ", header=None)
    data.columns = ["size", "time"]
    return data

def fit_hockney(data):
    X = data['size']
    X = sm.add_constant(X)
    y = data['time']
    model = sm.OLS(y, X)
    results = model.fit()
    alpha = results.params['const']
    beta = results.params['size']
    return alpha, beta

def main(fnames):
    if len(fnames) > 1:
        in_fname = fnames[0]
        out_fname = fnames[1]
    else:
        base_fname = fnames[0]
        in_fname = "{0}.txt".format(base_fname)
        out_fname = "{0}.svg".format(base_fname)
    data = load_pingpong(in_fname)
    alpha, beta = fit_hockney(data)
    print("alpha = {0}, beta = {1}".format(alpha, beta))
    data.plot.scatter(x="size", y="time")
    plt.plot(data['size'], alpha+beta*data['size'])
    plt.savefig(out_fname, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])
