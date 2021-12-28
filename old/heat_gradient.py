import numpy as np
import pandas as pd
from statistics import mean
import matplotlib.pyplot as plt

class GradientFromSPH:

    def __init__(self, path, tag, radius_lim):
        self.df = pd.read_csv(path, skiprows=2)
        self.df = self.df[self.df['tag'] == tag][self.df['radius'] <= radius_lim]
        self.df = self.df.sort_values(by=['radius'])
        self.min_rad = min(self.df['radius'])
        self.max_rad = max(self.df['radius'])

    def gravity(self, radius, surface_g=9.8, radius_earth=6378 * 1000):
        return surface_g * (1 - (radius / radius_earth))

    def temperature_gradient(self, num_samples=1000):
        inc = (self.max_rad - self.min_rad) / num_samples
        points = list(np.arange(self.min_rad, self.max_rad + inc, inc))
        means = []
        for index, p in enumerate(points):
            if index != 0:
                df = self.df[self.df['radius'] >= points[index - 1]][self.df['radius'] <= p]
                if len(df['temperature']) > 0:
                    m = mean(df['temperature'])
                else:
                    m = means[-1]
                means.append(m)
        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_subplot(111)
        # ax.scatter(
        #     self.df['radius'], self.df['temperature'], s=1, color='blue'
        # )
        # ax.plot(
        #     points[1:], means, linewidth=2.0, color='red'
        # )
        # ax.scatter(
        #     points[1:], means, s=10, color='black'
        # )
        # plt.show()
        return points[1:], means

    def mean_dT_dz_from_avg(self, num_samples=1000):
        points, temps = self.temperature_gradient(num_samples=num_samples)
        s = list(zip(points, temps))
        dT_dz = []
        for index, i in enumerate(s):
            p, t = i
            if index > 0:
                p_prev, t_prev = s[index - 1]
                dT_dz.append((t - t_prev) / (p - p_prev))
        # print(mean(dT_dz))
        # fig = plt.figure(figsize=(16, 9))
        # ax = fig.add_subplot(111)
        # ax.plot(
        #     points[1:],
        #     dT_dz
        # )
        # plt.show()
        # return points[1:], dT_dz
        return mean(dT_dz)

