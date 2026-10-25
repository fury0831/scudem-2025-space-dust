# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Fengrui Zhang, Yiye Huang, Zening Wang
"""Historical notebook algorithms, preserved for comparison.

These functions retain known numerical and endpoint-label limitations.
Importing this module runs no experiment. The modern CLI uses guarded RK4.
"""
import time
import numpy as np

RHO_P = 3000.0
C_D = 1.0
L_V = 8.0e6
LAMBDA = 0.1
G = 9.81
RHO_0 = 1.225
H = 8500.0

def particle_model(t, y):
    h, v, m = y
    if m <= 0:
        return np.array([0, 0, 0])
    try:
        rho_atm = RHO_0 * np.exp(-h / H)
        r = (3.0 * m / (4.0 * np.pi * RHO_P))**(1/3.0)
        A = np.pi * r**2
    except (ValueError, OverflowError):
        return np.array([0, 0, 0])
    dh_dt = -v
    F_drag = 0.5 * rho_atm * C_D * A * v**2
    dv_dt = G - (F_drag / m)
    dm_dt = - (LAMBDA * rho_atm * A * v**3) / (2.0 * L_V)
    return np.array([dh_dt, dv_dt, dm_dt])

def solve_rk4_trajectory(f, y0, t_span, dt):
    t = np.arange(t_span[0], t_span[1] + dt, dt)
    n_steps = len(t)
    y_history = np.zeros((n_steps, len(y0)))
    y_history[0] = y0
    y = y0.copy()
    
    for i in range(n_steps - 1):
        t_i = t[i]
        k1 = f(t_i, y)
        k2 = f(t_i + 0.5 * dt, y + 0.5 * dt * k1)
        k3 = f(t_i + 0.5 * dt, y + 0.5 * dt * k2)
        k4 = f(t_i + dt, y + dt * k3)
        y = y + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        
        # 检查终止条件
        if y[0] <= 0 or y[2] <= 0: # h <= 0 or m <= 0
            y_history[i+1:] = y
            t = t[:i+2]
            break
        
        y_history[i+1] = y
        
    return t, y_history.T

def solve_rk4_manual(f, y0, t_span, dt=0.1):
    t = np.arange(t_span[0], t_span[1] + dt, dt)
    n_steps = len(t)
    y = y0.copy()
    for i in range(n_steps - 1):
        t_i = t[i]
        k1 = f(t_i, y)
        k2 = f(t_i + 0.5 * dt, y + 0.5 * dt * k1)
        k3 = f(t_i + 0.5 * dt, y + 0.5 * dt * k2)
        k4 = f(t_i + dt, y + dt * k3)
        y = y + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        if y[0] <= 0 or y[2] <= 0:
            break
    return y

def run_survivability_grid_high_res_v():
    # --- 定义网格 ---
    GRID_SIZE_R = 5    # X轴 (半径) - 5 个点就够了
    GRID_SIZE_V = 50   # Y轴 (速度) - 50 个点，高分辨率
    
    r0_values_mm = np.linspace(0.05, 0.5, GRID_SIZE_R)
    v0_values_kms = np.linspace(11.2, 30.0, GRID_SIZE_V)
    
    survival_grid = np.zeros((GRID_SIZE_V, GRID_SIZE_R))
    
    h0 = 100.0 * 1e3
    t_span = [0, 100]
    
    total_sims = GRID_SIZE_R * GRID_SIZE_V
    print(f"--- 开始 {GRID_SIZE_R}x{GRID_SIZE_V} = {total_sims} 次高分辨率模拟 ---")
    start_time = time.time()
    
    sim_count = 0
    # --- 循环遍历网格 ---
    for i, r0_mm in enumerate(r0_values_mm):
        for j, v0_kms in enumerate(v0_values_kms):
            
            m0 = RHO_P * (4.0/3.0) * np.pi * (r0_mm * 1e-3)**3
            y0 = np.array([h0, v0_kms * 1e3, m0])
            
            y_final = solve_rk4_manual(particle_model, y0, t_span)
            m_final = max(0, y_final[2])
            survival_perc = (m_final / m0) * 100
            survival_grid[j, i] = survival_perc
            
            sim_count += 1
            if sim_count % 50 == 0:
                print(f"  ...已完成 {sim_count}/{total_sims} 次模拟")
                
    end_time = time.time()
    print(f"--- 模拟完成，总用时: {end_time - start_time:.2f} 秒 ---")
    
    return r0_values_mm, v0_values_kms, survival_grid

def plot_historical_trajectory(t, y_results, r0_mm, v0_km_s, dt):
    """Original cell-4 plot layout; returns a figure instead of calling show."""
    import matplotlib.pyplot as plt
    h = y_results[0] / 1e3
    v = y_results[1] / 1e3
    m = y_results[2]
    m_perc = (m / m[0]) * 100
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15), sharex=True)

    ax1.plot(t, h)
    ax1.set_ylabel("Altitude (km)")
    ax1.set_title(f"修正后的模拟: r0={r0_mm}mm, v0={v0_km_s}km/s, dt={dt}s")
    ax1.grid(True)

    ax2.plot(t, v)
    ax2.set_ylabel("Velocity (km/s)")
    ax2.grid(True)

    ax3.plot(t, m_perc)
    ax3.set_ylabel("Mass Remaining (%)")
    ax3.set_xlabel("Time (s)")
    ax3.set_ylim(-10, 110) # 允许我们看到它变为负值
    ax3.grid(True)

    plt.tight_layout()
    return fig

def plot_historical_grid(r_vals, v_vals, s_grid):
    """Original cell-6 plot; its survival label is historically misleading."""
    import matplotlib.pyplot as plt
    plt.figure(figsize=(11, 8))
    X, Y = np.meshgrid(r_vals, v_vals)

    # 使用 'contourf' (填充等高线) 会使过渡更平滑
    plt.contourf(X, Y, s_grid, levels=50, cmap='viridis')

    plt.colorbar(label='幸存质量百分比 (%)')
    plt.xlabel('初始半径 (mm)')
    plt.ylabel('初始速度 (km/s)')
    plt.title('宇宙尘埃幸存图 (高分辨率速度轴)')
    return plt.gcf()
