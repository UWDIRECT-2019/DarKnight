import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns


def vector_magnitude(data):
    """
    A function used to compute the average and std of magnitude of path vectors.
    """
    a = []
    for i in range(len(data)):
        s = 0
        for j in range(data.shape[1]):
            s += (data.iloc[i][j])**2
        s = np.sqrt(s)
        a.append(s)
    aveg = np.average(a)
    std = np.std(a)
    print ('The average magnitude is:',aveg)
    print ('The std magnitude is:',std)
    #return aveg,std

def vector_angle(rct,prd):
    """
    A function used to compute the average and std of angle of path vectors
    """
    #u = []
    #d = []
    angle = []
    for i in range(len(rct)):
        up = 0
        rm = 0
        pm = 0
        for j in range(rct.shape[1]):
            up += rct.iloc[i][j] * prd.iloc[i][j]  #numerator
            rm += (rct.iloc[i][j])**2  # the magnitude of reactant vector
            pm += (prd.iloc[i][j])**2  # the magnitude of product vector
        #u.append(up)
        rm = np.sqrt(rm)
        pm = np.sqrt(pm)
        cos = up/(rm*pm)
        a = math.degrees(math.acos(cos))
        #d.append(rm*pm)
        angle.append(a)
    aveg = np.average(angle)
    std = np.std(angle)
    print('The average angle is:',aveg)
    print('The std angle is:',std)

def difference(lact,lprd):
    """
    Function utilized to calculate the difference between
    the actual and predicted products latent vectors
    """
    d = []
    for i in range(len(lact)):
        s = 0
        for j in range(lact.shape[1]):
            s += (lact.iloc[i][j] - lprd.iloc[i][j])**2
        s = np.sqrt(s)
        d.append(s)
    return d

def distribution(data,rct,prd):
    """
    Function utilized to show the magnitude and angle of each reation
    """
    a = []
    for i in range(len(data)):
        s = 0
        for j in range(data.shape[1]):
            s += (data.iloc[i][j])**2
        s = np.sqrt(s)
        a.append(s)

    angle = []
    for i in range(len(rct)):
        up = 0
        rm = 0
        pm = 0
        for j in range(rct.shape[1]):
            up += rct.iloc[i][j] * prd.iloc[i][j]  #numerator
            rm += (rct.iloc[i][j])**2  # the magnitude of reactant vector
            pm += (prd.iloc[i][j])**2  # the magnitude of product vector
        #u.append(up)
        rm = np.sqrt(rm)
        pm = np.sqrt(pm)
        cos = up/(rm*pm)
        b = math.degrees(math.acos(cos))
        #d.append(rm*pm)
        angle.append(b)

    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot2grid((1, 2), (0, 0))
    ax.hist(a, color = 'blue', edgecolor = 'black')
    plt.title('magnitude distribution')
    plt.xlabel('magnitude')
    plt.ylabel('numbers of reatction')

    ax1 = plt.subplot2grid((1, 2), (0, 1))
    ax1.hist(angle, color = 'blue', edgecolor = 'black')
    plt.title('angle distribution')
    plt.xlabel('angle')
    plt.ylabel('numbers of reatction')
