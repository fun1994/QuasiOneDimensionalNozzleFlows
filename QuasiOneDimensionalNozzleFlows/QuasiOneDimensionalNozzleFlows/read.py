# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 13:48:13 2024

@author: HFKJ059
"""

import numpy as np


def read_1d(path, filename):
    with open("./" + path + "/" + filename + ".txt", "r") as file:
        data = file.read()
    data = data.split()
    for i in range(len(data)):
        data[i] = float(data[i])
    data = np.array(data)
    return data

def read_2d(path, filename):
    data = []
    with open("./" + path + "/" + filename + ".txt", "r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            data_temp = line.split()
            data.append(data_temp)
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = float(data[i][j])
    data = np.array(data)
    return data

def read(path):
    x = read_1d(path, "x")
    t = read_1d(path, "time")
    rho = read_2d(path, "rho")
    V = read_2d(path, "V")
    T = read_2d(path, "T")
    p = read_2d(path, "p")
    Ma = read_2d(path, "Ma")
    m = read_2d(path, "m")
    U1 = read_2d(path, "U1")
    U2 = read_2d(path, "U2")
    U3 = read_2d(path, "U3")
    return x, t, rho, V, T, p, Ma, m, U1, U2, U3


x, t, rho, V, T, p, Ma, m, U1, U2, U3 = read("purelySubsonicIsentropicNozzleFlow")
