import numpy as np
import matplotlib.pyplot as plt

# Generate some sample data
data = np.random.rand(5, 5)

# Create heatmap
fig, ax = plt.subplots()
im = ax.imshow(data)

# Add colorbar
cbar = ax.figure.colorbar(im, ax=ax)

# Set plot title and labels
ax.set_title("Heatmap of 2D array values")
ax.set_xlabel("X-axis label")
ax.set_ylabel("Y-axis label")

plt.show()

