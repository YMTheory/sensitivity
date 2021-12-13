import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LSQUnivariateSpline

def FindIntersectionByQuadraticInterpolation( xvals, yvals, SplineFunction ):

    if len(xvals)!=len(yvals):
        print('ERROR: need same length arrays for ' +\
            'quadratic interpolation')
        raise ValueError
    print('Xvals:')
    print(xvals)
    print('Yvals:')
    print(yvals)

    # First, select only lambda values on the upward slope, so we find
    # the upper limit
    mask = np.zeros(len(yvals),dtype=bool)
    mask[1:] = ( yvals[1:] - yvals[:-1] ) > 0.
    # Next, select only values near the critical lambda threshold (~2.7)
    mask = mask&(yvals>0.5)&(yvals<6.)

    xfit = np.linspace(0.,30.,600)

    if True:
        try:
            p = np.polyfit( xvals, yvals, 2)
            print('Quadratic fit: {}'.format(p))
            yfit = p[0]*xfit**2 + p[1]*xfit + p[2]
        except np.RankWarning:
            p = np.polyfit( xvals, yvals, 1)
            print('Linear fit: {}'.format(p))
            yfit = p[0]*xfit + p[1]
        yspline = SplineFunction(xfit)
        crossing_idx = np.where( (yfit - yspline)>0. )[0][0]
        crossing = xfit[crossing_idx]
    else:
        yfit = np.zeros(len(xfit))
        crossing_idx = -1
        crossing = -1.

    return xfit, yfit, crossing, crossing_idx

if __name__ == '__main__':

    X = []
    for i in range(122):
        X.append(i * .25)

    colors = ['b','r','k','c','g','m','y']
    labels = ['lt=0.25 yr', 'lt=0.5 yr', 'lt=1 yr', 'lt=2 yr', 'lt=5 yr', 'lt=10 yr', "Wilks' Approximation"]
    # labels = ['lt=0.25 yr', 'lt=0.5 yr', "Wilks' Approximation"]
    # labels = ['lt=1 yr', 'lt=2 yr', "Wilks' Approximation"]
    # labels = ['lt=5 yr', 'lt=10 yr', "Wilks' Approximation"]

    for index, D in enumerate(['024']):
        Y=None
        # folder = "/Users/czyz1/lc-nexouser/workdir/lambda/21_01_21_DNN1_{}/".format(D)
        folder = "/Users/czyz1/lc-nexouser/workdir/lambda/21_03_22_DNN1_{}/".format(D)
        # for res in ['0.008', '0.011', '0.015', '0.018']:
        # for res in ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016']:
        for lt in ['0.25', '0.5', '1.0', '2.0', '5.0', '10.0']:
        # for lt in ['0.25', '0.5']:
        # for lt in ['1.0', '2.0']:
        # for lt in ['5.0', '10.0']:
            file = folder + "lambdas_lt={}.txt".format(lt)
            Y = np.genfromtxt(file)
            # file = folder + "lambdas_res={}.txt".format(res)
            # Ynew = np.genfromtxt(file)
            # if Y is not None:
            #     Y += Ynew
            # else:
            #     Y = Ynew
            # print(Y)
            # spline_xn = np.array([.75, 2.25, 4.5, 7, 10, 15, 21, 30])  # defines the locations of the knots
            spline_xn = np.array([.75, 2.5, 4.5, 10, 15, 21, 30])  # defines the locations of the knots

            # spline_xn = np.array([1, 3, 6, 10, 15, 21, 30])  # defines the locations of the knots

            # spline_xn = np.array([1., 2, 3, 4, 5., 7., 10., 12.5, 15, 17.5, 20., 22.5, 25, 27.5, 30])  # defines the
            SplineFunc = LSQUnivariateSpline(X, Y, t=spline_xn, k=3)
            xfit, yfit, crossing, crossing_idx = FindIntersectionByQuadraticInterpolation(X,Y, SplineFunc)
            # if index == 0:
            plt.plot(xfit, SplineFunc(xfit),label='Fit, LT (years) = '+lt)
            # elif index == 1:
            #     plt.plot(xfit, SplineFunc(xfit), '-.', label=labels[index] + '_' + res)
            # plt.plot(X,Y,label='Raw, LT (years) = '+lt)
        # plt.plot(xfit,yfit)
        plt.ylim([.5, 4])
        plt.xlim([0, 10])
        # np.savetxt(folder + "lambdas_res=0.008.txt", Y)
    plt.plot(X, [2.706 for i in range(122)],label=labels[-1])
    plt.legend( loc='lower right')
    plt.savefig('/Users/czyz1/lc-nexouser/output/plots/21_03_01/lambda_resolution=0.008.png', dpi=200, bbox_inches='tight')
    plt.show()