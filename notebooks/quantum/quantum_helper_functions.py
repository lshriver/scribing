from qiskit import transpile
from qiskit_aer import Aer, AerSimulator, AerJob

def execute(circuits, backend=None, shots=1024, **kwargs):
    if backend is None:
        backend = AerSimulator()
        
    # Transpile the circuits for the backend
    transpiled_circuits = transpile(circuits, backend)
    
    # Run the circuits on the backend
    job = backend.run(transpiled_circuits, shots=shots, **kwargs)
    
    return job

from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt

def run_circuit(qc,simulator='statevector_simulator',shots=1,hist=True):
    # Tell Qiskit how to simulate our circuit
    backend = Aer.get_backend(simulator)
    
    # Execute the qc
    results = execute(qc,backend, shots=shots).result().get_counts()
    
    # Plot the results
    return plot_histogram(results, figsize=(18,4)) if hist else results             

# Calculating the angle that represents a certain probability
from math import asin, sqrt

def prob_to_angle(prob):
    """"
    Converts a given P(psi) value into an equivalent theta value
    """
    return 2*asin(sqrt(prob))