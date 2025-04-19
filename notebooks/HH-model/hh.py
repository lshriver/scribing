import numpy as np
import matplotlib.pyplot as plt

def alpha_n(V):
    """Rate constant alpha_n(V) for gating variable n."""
    return (0.01 * (10.0 - (V + 65.0))) / (np.exp((10.0 - (V + 65.0)) / 10.0) - 1.0)

def beta_n(V):
    """Rate constant beta_n(V) for gating variable n."""
    return 0.125 * np.exp(-(V + 65.0) / 80.0)

def alpha_m(V):
    """Rate constant alpha_m(V) for gating variable m."""
    return (0.1 * (25.0 - (V + 65.0))) / (np.exp((25.0 - (V + 65.0)) / 10.0) - 1.0)

def beta_m(V):
    """Rate constant beta_m(V) for gating variable m."""
    return 4.0 * np.exp(-(V + 65.0) / 18.0)

def alpha_h(V):
    """Rate constant alpha_h(V) for gating variable h."""
    return 0.07 * np.exp(-(V + 65.0) / 20.0)

def beta_h(V):
    """Rate constant beta_h(V) for gating variable h."""
    return 1.0 / (np.exp((30.0 - (V + 65.0)) / 10.0) + 1.0)

def hodgkin_huxley_sim(t_start=-30.0, t_end=200.0, dt=0.01):
    """Simulate the classic Hodgkin-Huxley model with multiple external current steps and record ionic currents and capacitance.
    
    External current steps:
      - 5 to 25 ms: depolarizing pulse (10 µA/cm² with noise)
      - 50 to 70 ms: hyperpolarizing pulse (-5 µA/cm²)
      - 100 to 120 ms: depolarizing pulse (7 µA/cm² with noise)
      - 150 to 170 ms: hyperpolarizing pulse (-3 µA/cm²)
      - 170 to 200 ms: long step held at 5 µA/cm²
      
    Args:
        t_start (float): start time in ms.
        t_end (float): end time in ms.
        dt (float): timestep in ms.
        
    Returns:
        t (np.ndarray): time vector.
        V_trace (np.ndarray): membrane potential over time (mV).
        n_trace, m_trace, h_trace (np.ndarray): gating variables over time.
        I_ext_trace (np.ndarray): external current over time.
        IK_trace, INa_trace, IL_trace (np.ndarray): individual ionic currents (µA/cm²).
        I_ion_trace (np.ndarray): net ionic current (µA/cm²).
        Cm_trace (np.ndarray): membrane capacitance over time (µF/cm²) (constant here).
    """
    # Conductances (mS/cm²)
    gNa = 120.0
    gK = 36.0
    gL = 0.3

    # Reversal potentials (mV)
    ENa = 50.0    
    EK = -77.0    
    EL = -54.4    
    
    # Membrane capacitance (µF/cm²); assumed constant here
    Cm = 1.0

    # Time vector
    t = np.arange(t_start, t_end + dt, dt)
    nSteps = len(t)
    
    # Allocate buffers for storing results
    V_trace = np.zeros(nSteps)
    n_trace = np.zeros(nSteps)
    m_trace = np.zeros(nSteps)
    h_trace = np.zeros(nSteps)
    I_ext_trace = np.zeros(nSteps)
    
    # Ionic current traces
    IK_trace = np.zeros(nSteps)
    INa_trace = np.zeros(nSteps)
    IL_trace = np.zeros(nSteps)
    I_ion_trace = np.zeros(nSteps)
    
    # Membrane capacitance trace (constant)
    Cm_trace = np.ones(nSteps) * Cm
    
    # Initial conditions
    V = -65.0   # membrane potential (mV)
    n_gate = alpha_n(V) / (alpha_n(V) + beta_n(V))  # steady-state n
    m_gate = alpha_m(V) / (alpha_m(V) + beta_m(V))  # steady-state m
    h_gate = alpha_h(V) / (alpha_h(V) + beta_h(V))  # steady-state h
    
    # Run the simulation loop
    for i, tt in enumerate(t):
        # Store current state
        V_trace[i] = V
        n_trace[i] = n_gate
        m_trace[i] = m_gate
        h_trace[i] = h_gate
        
        # External current steps:
        if (tt >= 5 and tt < 25):
            I_ext = 10 + 2 * np.random.randn()   # depolarizing pulse with noise
        elif (tt >= 50 and tt < 70):
            I_ext = -5                           # hyperpolarizing pulse
        elif (tt >= 100 and tt < 120):
            I_ext = 7 + 1 * np.random.randn()      # depolarizing pulse with noise
        elif (tt >= 150 and tt < 170):
            I_ext = -3                           # hyperpolarizing pulse
        elif (tt >= 170 and tt < 200):
            I_ext = 5                            # long step at 5 µA/cm²
        else:
            I_ext = 0.0
        I_ext_trace[i] = I_ext
        
        # Compute rate constants
        an = alpha_n(V)
        bn = beta_n(V)
        am = alpha_m(V)
        bm = beta_m(V)
        ah = alpha_h(V)
        bh = beta_h(V)
        
        # Update gating variables (Euler method)
        dn_dt = an * (1.0 - n_gate) - bn * n_gate
        dm_dt = am * (1.0 - m_gate) - bm * m_gate
        dh_dt = ah * (1.0 - h_gate) - bh * h_gate
        
        n_gate += dn_dt * dt
        m_gate += dm_dt * dt
        h_gate += dh_dt * dt
        
        # Compute ionic currents (µA/cm²)
        IK = gK * (n_gate**4) * (V - EK)
        INa = gNa * (m_gate**3) * h_gate * (V - ENa)
        IL = gL * (V - EL)
        I_ion = IK + INa + IL
        
        # Store ionic currents
        IK_trace[i] = IK
        INa_trace[i] = INa
        IL_trace[i] = IL
        I_ion_trace[i] = I_ion
        
        # Update membrane potential using Euler's method: dV/dt = (I_ext - I_ion)/Cm
        V += dt * (I_ext - I_ion) / Cm
    
    return t, V_trace, n_trace, m_trace, h_trace, I_ext_trace, IK_trace, INa_trace, IL_trace, I_ion_trace, Cm_trace

if __name__ == "__main__":
    # Run the simulation from -30 ms to 200 ms
    (t, V, n, m, h, I_ext,
     IK, INa, IL, I_ion, Cm) = hodgkin_huxley_sim(t_start=-30.0, t_end=200.0, dt=0.01)
    
    # Subselect data for t >= 0 so plots start at t = 0
    idx = t >= 0
    t_plot = t[idx]
    V_plot = V[idx]
    I_ext_plot = I_ext[idx]
    n_plot = n[idx]
    m_plot = m[idx]
    h_plot = h[idx]
    IK_plot = IK[idx]
    INa_plot = INa[idx]
    IL_plot = IL[idx]
    I_ion_plot = I_ion[idx]
    Cm_plot = Cm[idx]
    
    # First Figure: Original plots (Voltage, External Current, and gating variables)
    plt.figure(figsize=(10, 18))
    
    # Membrane Potential (dark purple)
    plt.subplot(5, 1, 1)
    plt.plot(t_plot, V_plot, label="V (mV)", color="#4B0082")
    plt.title("Hodgkin-Huxley Model Simulation")
    plt.ylabel("Membrane Potential (mV)")
    plt.grid(True)
    plt.legend()
    
    # External Current (bright light blue)
    plt.subplot(5, 1, 2)
    plt.plot(t_plot, I_ext_plot, label="I_ext (µA/cm²)", color="#87CEFA")
    plt.ylabel("External Current (µA/cm²)")
    plt.grid(True)
    plt.legend()
    
    # Potassium channel activation (n) in orange
    plt.subplot(5, 1, 3)
    plt.plot(t_plot, n_plot, label="n (K⁺ activation)", color="#FFA500")
    plt.ylabel("n")
    plt.grid(True)
    plt.legend()
    
    # Sodium channel activation (m) in blue
    plt.subplot(5, 1, 4)
    plt.plot(t_plot, m_plot, label="m (Na⁺ activation)", color="b")
    plt.ylabel("m")
    plt.grid(True)
    plt.legend()
    
    # Sodium channel inactivation (h) in green
    plt.subplot(5, 1, 5)
    plt.plot(t_plot, h_plot, label="h (Na⁺ inactivation)", color="#008000")
    plt.xlabel("Time (ms)")
    plt.ylabel("h")
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.show()
    
    # Second Figure: Ionic currents and membrane capacitance
    plt.figure(figsize=(10, 18))
    
    # Membrane Capacitance (should be constant)
    plt.subplot(4, 1, 1)
    plt.plot(t_plot, Cm_plot, label="Cm (µF/cm²)", color="k")
    plt.title("Membrane Capacitance and Ionic Currents")
    plt.ylabel("Cm (µF/cm²)")
    plt.grid(True)
    plt.legend()
    
    # Net Ionic Current (I_ion)
    plt.subplot(4, 1, 2)
    plt.plot(t_plot, I_ion_plot, label="I_ion (net ionic current)", color="m")
    plt.ylabel("I_ion (µA/cm²)")
    plt.grid(True)
    plt.legend()
    
    # Individual Ionic Currents: IK, INa, and IL
    plt.subplot(4, 1, 3)
    plt.plot(t_plot, IK_plot, label="IK (K⁺ current)", color="#FF4500")  # Orange-red
    plt.plot(t_plot, INa_plot, label="INa (Na⁺ current)", color="#1E90FF")  # Dodger blue
    plt.plot(t_plot, IL_plot, label="IL (Leak current)", color="#32CD32")  # Lime green
    plt.ylabel("Ionic Currents (µA/cm²)")
    plt.grid(True)
    plt.legend()
    
    # In case you want to compare net vs. individual, overlay net current
    plt.subplot(4, 1, 4)
    plt.plot(t_plot, I_ion_plot, label="I_ion (net ionic)", color="m", linestyle="--")
    plt.plot(t_plot, IK_plot, label="IK", color="#FF4500", alpha=0.7)
    plt.plot(t_plot, INa_plot, label="INa", color="#1E90FF", alpha=0.7)
    plt.plot(t_plot, IL_plot, label="IL", color="#32CD32", alpha=0.7)
    plt.xlabel("Time (ms)")
    plt.ylabel("Currents (µA/cm²)")
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.show()
