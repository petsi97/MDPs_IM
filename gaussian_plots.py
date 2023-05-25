import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math

def plot_gaussian():
    mu = 0
    variance = 1
    sigma = math.sqrt(variance)

    for i in range(20):
        mu = 100
        sigma = 2*mu*(i+1)/20
        x = np.linspace(0, mu + 3*sigma, 100)
        plt.plot(x, stats.norm.pdf(x, mu, sigma))
        
    plt.show()

    
if __name__ == "__main__":        
    plot_gaussian()