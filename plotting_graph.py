import matplotlib.pyplot as plt
import numpy as np

choice = int(input("Book's code (1) or problem's code? (2): "))

filename = "wave.data" if choice == 1 else "ch4prob13.txt"

with open(filename, "r") as f:
    nums = [[float(num) for num in line.split()] for line in f.readlines()]

x_axis = np.array([point[0] for point in nums])
y_axis = np.array([point[1] for point in nums])

plt.plot(x_axis, y_axis)

plt.xlim([-1.0, 6.0]) if choice == 1 else plt.xlim([-1.0, 6.0])  # [-1.0, 6.0]

plt.show()
