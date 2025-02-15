import numpy as np

def random_modular_graph(n, c, p, r):
    """
    Build a random modular graph, given number of modules, and link density.
    
    Parameters:
    n (int): Number of nodes
    c (int): Number of clusters/modules
    p (float): Overall probability of attachment
    r (float): Proportion of links within modules
    
    Returns:
    tuple: (adjacency matrix, modules to which the nodes are assigned)
    """
    
    z = round(p * (n - 1))  # z = p(n - 1) - Erdos-Renyi average degree

    # Assign nodes to modules
    modules = [list(range(round((k - 1) * n / c), round(k * n / c))) for k in range(1, c + 1)]

    adj = np.zeros((n, n), dtype=int)  # Initialize adjacency matrix

    for i in range(n):
        for j in range(i + 1, n):
            module_i = (i * c) // n  # The module to which i belongs
            module_j = (j * c) // n  # The module to which j belongs
            
            if module_i == module_j:
                # Probability of attachment within module
                if np.random.rand() <= r * z / (n / c - 1):
                    adj[i, j] = 1
                    adj[j, i] = 1
            else:
                # Probability of attachment across modules
                if np.random.rand() <= z * (1 - r) / (n - n / c):
                    adj[i, j] = 1
                    adj[j, i] = 1

    return adj, modules

# Example usage
n = 100  # Number of nodes
c = 5    # Number of modules
p = 0.1  # Overall probability of attachment
r = 0.5  # Proportion of links within modules

adjacency_matrix, modules = random_modular_graph(n, c, p, r)
print("Adjacency Matrix:\n", adjacency_matrix)
print("Modules:\n", modules)