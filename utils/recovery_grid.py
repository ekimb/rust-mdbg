import matplotlib.pyplot as plt
import numpy
  
# ---DENSITY PLOT---  
x = list(numpy.arange(0.001, 0.0055, 0.0005))
y = [1.49, 22.68, 73.11, 93.20, 96.60, 97.00, 97.18, 97.33, 97.42]
plt.plot(x, y)
plt.xlabel('Density')
plt.ylabel('K-min-mer recovery (%)')
plt.show()
# ---K-MINMER LENGTH PLOT---  
x = list(range(10, 51))
y = [97.92,
    97.85,
    97.79,
    97.74,
    97.69,
    97.63,
    97.58,
    97.53,
    97.47,
    97.42,
    97.37,
    97.31,
    97.26,
    97.20,
    97.14,
    97.07,
    97.00,
    96.93,
    96.84,
    96.73,
    96.60,
    96.43,
    96.22,
    95.92,
    95.55,
    95.14,
    94.64,
    93.96,
    93.09,
    91.93,
    90.49,
    88.65,
    86.45,
    83.91,
    80.94,
    77.60,
    73.82,
    69.88,
    65.35,
    60.61,
    55.80]
plt.plot(x, y)
plt.xlabel('K-min-mer length')
plt.ylabel('K-min-mer recovery (%)')
plt.show()