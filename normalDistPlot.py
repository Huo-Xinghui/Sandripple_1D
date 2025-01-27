import numpy as np

import matplotlib.pyplot as plt

# Parameters for the normal distribution
mean = 300
std_dev = 100

# Generate x values
x = np.linspace(50, 600, 40)

# Calculate the normal distribution values
y = (1 / (std_dev * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std_dev) ** 2)

# Plot the normal distribution curve
plt.plot(x, y, label='Normal Distribution (mean=300, std=100)')
plt.xlabel('x')
plt.ylabel('Probability Density')
plt.title('Normal Distribution Curve')
plt.legend()
plt.grid(True)
plt.show()