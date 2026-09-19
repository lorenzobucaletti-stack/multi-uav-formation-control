# Distributed Formation Control of Nonlinear Multi-Agent Systems

[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Distributed cooperative control architecture for a heterogeneous multi-agent fleet of three planar autonomous robots required to track dynamic reference trajectories (linear, circular, sinusoidal) while maintaining a rigid equilateral-triangular formation.

---

## Architecture Overview

Direct linear distributed consensus is unfeasible on heterogeneous nonlinear systems. The proposed architecture solves this through a two-stage hierarchical framework:

1. **Local Input-Output Feedback Linearization (IOFL):**
   - Each agent is governed by distinct nonlinear affine drift dynamics and input matrices:
     - **Agent 1:** Uncoupled dynamics with quadratic aerodynamic drag terms.
     - **Agent 2:** State-dependent diagonal decoupling matrix governed by trigonometric bounds.
     - **Agent 3:** Complex cross-axis asymmetric couplings.
   - The IOFL control law analytically cancels the physical nonlinearities, mapping each physical robot onto a uniform, linear double-integrator virtual plant ($\ddot{p}_i = \nu_i$).
   - Since the vector relative degree is strictly $[2, 2]^T$ for all agents (summing to the state dimension 4), internal zero dynamics are completely absent.

2. **Distributed Cooperative Formation Layer:**
   Two distinct multi-agent cooperative methodologies were designed, implemented, and benchmarked:
   - Dynamic Output-Feedback via Loop-Shaping (Frequency Domain)
   - Distributed State-Variable Feedback with Luenberger State Estimation (Time Domain / LQR)

---

## Comparative Methodologies

### 1. Loop-Shaping Approach (Output-Feedback)
- **Concept:** Uses only relative position measurements ($\epsilon_i^p$) across an undirected communication graph with leader pinning.
- **Compensator Synthesis:** Homogenized plants $G(s) = 1/s^2$ have an intrinsic $-180^\circ$ phase. To avoid high-frequency noise amplification caused by a pure PD controller, a Lead Compensator was designed:

$$K_{\text{Lead}}(s) = K_p + \frac{K_d s}{1 + \tau_{\text{lead}}s}$$

- **Stability & Margins:** Synthesized via modal decomposition of the pinned Laplacian matrix $L_p = L + \gamma \Pi_1$. Guaranteed a Phase Margin $\ge 55^\circ$ across all pinned Laplacian eigenvalues.
- **Feedforward Acceleration Injection:** To cancel tracking delays during accelerating trajectories (circular/sinusoidal), a reference acceleration feedforward term $\ddot{r}(t)$ was integrated.
- **Robustness:** The sensitivity functions ($S(s)$ and W(s)) ensure optimal attenuation of low-frequency disturbances ($5\text{ rad/s}$) and high-frequency sensor noise ($100\text{ rad/s}$), preventing actuator saturation.

### 2. State-Variable Feedback (SVFB / LQR + Observer)
- **Concept:** Optimal full-state consensus protocol formulated in the time domain.
- **LQR ARE Solution:** Solves the continuous-time Algebraic Riccati Equation to establish optimal feedback gains balancing tracking error penalties ($Q$) and actuation efforts ($R$).
- **Luenberger Observer:** Because velocity states are not measured directly, dual Riccati-tuned observers reconstruct the full state vector in real time from noisy position measurements.
- **Performance:** Achieves fast convergence and near-zero tracking errors in nominal conditions, but shows heightened sensitivity and saturation chatter under persistent external sinusoidal disturbances.

---

## Experimental Benchmark Summary

| Metric / Scenario | State-Feedback (LQR + Observer) | Loop-Shaping (Lead Compensator + FF) |
| :--- | :--- | :--- |
| **Feedback Information** | Full state $[p_x, v_x, p_y, v_y]^T$ (via Observer) | Relative output position only $[p_x, p_y]^T$ |
| **Convergence Speed** | Faster under nominal conditions ($t_{\text{conv}} \approx 2.49\text{ s}$) | Smooth transient ($t_{\text{conv}} \approx 4.3\text{ s}$) |
| **Persistent Disturbance Rejection** | Degraded (persistent oscillations and tracking offset) | Robust (active sensitivity shaping filters disturbances) |
| **Actuator Saturation ($\pm 8\text{ N}$)** | Frequently reaches saturation bounds on sharp turns | Smooth control action, strictly respects input bounds |
| **Sensor Noise Sensitivity** | Requires accurate noise covariance scaling ($R_o$) | Inherently filtered via high-frequency pole ($\omega_p = 50\text{ rad/s}$) |

---

## Repository Structure

```text
├── docs/
│   └── part_II_CPS_Group8.pdf     # Full technical conference-style paper
├── src/
│   ├── loop_shaping/
│   │   ├── lp_new_proj.m          # Pure PD vs Lead Compensator benchmark
│   │   └── lp_new_proj_ff.m       # Lead Compensator with Feedforward dynamics
│   └── state_feedback/
│       └── sf_proj.m              # LQR distributed tracking + Luenberger observer
└── README.md
