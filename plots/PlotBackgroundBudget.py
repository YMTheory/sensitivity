#!/usr/bin/env python

"""PlotBackgroundBudget.py: create figures of background budget

Input is from DB Excel file. Plots are generated grouping by Isotope, Material, or Component."""

import os.path
from math import floor, log10

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def SqrtSumSq(x):
    x = np.array(x.tolist())
    return np.sqrt(np.sum(np.multiply(x, x)))


def make_plot(groupby, filename, xlimit=5.0e-6, ):
    # grouping with custom aggregration
    df2 = df.groupby(groupby + ["CV?"]).agg(
        {'C.V.': np.sum, 'Error': SqrtSumSq, 'Limit 90% C.L.': np.sum, 'CV?': all})

    # add a variable for sorting that considers both CV and 90CL appropriately
    df2['value'] = df2['C.V.'] * df2['CV?'] + df2['Limit 90% C.L.'] * (~df2['CV?'])
    df2.sort_values('value', ascending=True, inplace=True)
    print(df2)

    # loop over all the rows to fill the lists used for plotting
    label = []
    value = []
    err = []
    xuplims = []
    xlimit_arrow = 10 ** floor(log10(xlimit)) * 10
    for index, row in df2.iterrows():
        if row['Limit 90% C.L.'] < 1e-6:
            continue

        label.append(' '.join(index[:len(groupby)]))

        # the counts are in 10 years and 3 tonne so divide by 30 to get the cts/tonne/year
        if row['CV?']:
            xuplims.append(0)
            value.append(row['C.V.'] / 30.)
            err.append(row['Error'] / 30.)
        else:
            xuplims.append(1)
            value.append(row['Limit 90% C.L.'] / 30.)
            err.append(value[-1] - xlimit_arrow)

    # now do the plotting
    nn = len(label)
    y = np.arange(0, 2 * nn, 2)
    yerr = np.ones(nn) * 0.8

    fig, ax0 = plt.subplots()
    ax0.errorbar(value, y, xerr=err, yerr=yerr, xuplims=xuplims, fmt=',')
    ax0.set_yticks(y)
    ax0.set_yticklabels(label)
    ax0.set_xscale('log')
    ax0.set_xlim(xlimit)
    ax0.set_xlabel("cts/ROI/tonne/year")
    plt.subplots_adjust(left=0.32, )
    # plt.show()

    plt.savefig(filename, dpi=200)
    return


df = pd.read_excel('../tables/Summary_v73_2016-09-09_0nu.xlsx', sheetname='SS_ExpectedCounts', header=0,
                   skiprows=4, skip_footer=7, parse_cols="A:C,AI:AK", )

# this is for printing without breaking into multiple lines
pd.set_option('display.expand_frame_repr', False)

# remove bb0n and bb2n rows
df = df[~df['Isotope'].isin(['bb2n', 'bb0n'])]

# Add a column that tracks whether the CV>0 or not
df['CV?'] = pd.Series(df['C.V.'] > 0, index=df.index)
print(df)

# Plot by Material, Isotope, and Component separately
make_plot(['Material'], os.path.expanduser("~/Scratch/BackgroundBudgetByMaterial.png"))
make_plot(['Isotope'], os.path.expanduser("~/Scratch/BackgroundBudgetByIsotope.png"), 5.0e-4)
make_plot(['Component'], os.path.expanduser("~/Scratch/BackgroundBudgetByComponent.png"), 5.0e-6)

# This is by Material & Isotope
make_plot(["Material", "Isotope"], os.path.expanduser("~/Scratch/BackgroundBudget.png"), 5.0e-6)
