from scripts.mutation import *
import plotly.express as px


def mut():
    df = pd.read_csv("../mut_test/sweep.csv")
    df = df.drop(["jkey_dec", "jlock_dec"], axis=1)
    tm = pivot_matrix(df, "tau")
    dtm = np.sum(np.abs(np.diff(tm, axis=0)[::2]), axis=0)
    return dtm


im = mut()
px.imshow(im, zmin=0).show("browser")
