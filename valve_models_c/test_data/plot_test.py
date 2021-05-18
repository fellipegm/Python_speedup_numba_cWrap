#!/usr/bin/env python3
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import os


dir = os.getcwd()
os.chdir(dir+"/test_data")

orig_data = np.genfromtxt("dados.csv",delimiter=",")
sim_data = np.genfromtxt("simulation.csv",delimiter=",")

plt.figure()
plt.plot(orig_data[:,0], orig_data[:,4], 'k')
plt.plot(sim_data[:,0], sim_data[:,2], 'r')
plt.legend(['matlab', 'c'])
plt.title('x')

plt.figure()
plt.plot(orig_data[:,0], orig_data[:,2], 'k')
plt.plot(sim_data[:,0], sim_data[:,1], 'r')
plt.legend(['matlab', 'c'])
plt.title('P')
plt.show()