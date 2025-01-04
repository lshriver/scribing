import pennylane as qml
from pennylane import numpy as np
from matplotlib import pyplot as plt

def plot_opt_surface(costs, params, cost_func):

    n = 200
    theta1 = np.linspace(-2*np.pi, 2*np.pi, n)
    theta2 = np.linspace(-2*np.pi, 2*np.pi, n)
    # Grid of elements
    X, Y = np.meshgrid(theta1, theta2)
    # Cost function on grid
    z = []
    for y in Y:
        for x in X[0]:
            z.append(cost_func([x, y[0]]))
    Z = np.reshape(z, (n, n))

    #Plotting the cost function of grid
    fig = plt.figure(figsize = (10, 10))
    ax = fig.add_subplot(projection='3d')
    ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10, edgecolor='grey', linewidth=0.5)

    # Plotting the optimization result
    param1, param2 = np.split(params, 2, axis=1)
    ax.scatter(param1, param2, costs, c='red', s = 20, marker = 'o')

    # Set labels
    ax.set_xlabel(r'$\theta_1$')
    ax.set_ylabel(r'$\theta_2$')
    ax.set_title(r'Optimization of Quantum Circuit Using Adam Optimizer')
    ax.legend(('Grid of cost function', 'Optimization results'))  