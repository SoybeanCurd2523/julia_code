using JuMP
import Ipopt
import Plots

h_0 = 1                      # Initial height
v_0 = 0                      # Initial velocity
m_0 = 1.0                    # Initial mass
m_T = 0.6                    # Final mass
g_0 = 1                      # Gravity at the surface
h_c = 500                    # Used for drag
c = 0.5 * sqrt(g_0 * h_0)    # Thrust-to-fuel mass
D_c = 0.5 * 620 * m_0 / g_0  # Drag scaling
u_t_max = 3.5 * g_0 * m_0    # Maximum thrust
T_max = 0.2                  # Number of seconds
T = 1_000 # 1000             # Number of time steps

# \Delta
Δt = 0.2 / T;                # Time per discretized step

model = Model(Ipopt.Optimizer)

# It is good practice for nonlinear programs to always provide a starting solution for each variable.
@variable(model, x_v[1:T] >= 0 , start = v_0)
@variable(model, x_h[1:T] >= 0 , start = h_0)
@variable(model, x_m[1:T] >= m_T, start = m_0)
@variable(model, 0 <= u_t[1:T] <= u_t_max, start = 0)
