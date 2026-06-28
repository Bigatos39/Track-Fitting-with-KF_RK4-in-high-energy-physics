# Track Fitting with Kalman Filter + Runge-Kutta 4 in High Energy Physics

A nonlinear charged particle track fitting framework combining **Kalman Filter (KF)** and **4th-order Runge-Kutta (RK4)** propagation for a 10-detector spectrometer system with magnetic deflection.

---

## Physics Overview

An electron beam (e⁻, 1 GeV/c) traverses a 10-detector system:

```
DetIn1→DetIn2→DetIn3→DetIn4→DetIn5 ══[B=6.54T]══ DetOut1→DetOut2→DetOut3→DetOut4→DetOut5
|←————————— drift (300mm) ————————→|←— dipole —→|←————————— drift (415mm) ————————→|
                                        θ ≈ 45°
```

- **Detector material**: Xenon (Xe) — introduces multiple scattering
- **Magnetic field**: B = 6.54 T, fieldLength = 400 mm, rotation Y30°
- **Bending angle**: θ ≈ 45°, ρ ≈ 510 mm

---

## Method

### State Vector (4D)

```
x = (y, ty, x, tx)
```

- **y, ty**: transverse coordinate and slope — unaffected by B field (B = By)
- **x, tx**: bending coordinate and slope — deflected by Lorentz force

### RK4 Equations of Motion

| Region  | dy/ds | dty/ds | dx/ds | dtx/ds |
| ------- | ----- | ------ | ----- | ------ |
| Drift   | ty    | 0      | tx    | 0      |
| B field | ty    | 0      | tx    | +1/rho |

### Kalman Filter Steps

**Predict:**

```
x_pred   = RK4(x, h, rho)
Cov_pred = F * Cov * F^T + Q
```

where F is the Jacobian matrix:

```
F = | 1  h_y  0   0  |
    | 0   1   0   0  |
    | 0   0   1  h_x |
    | 0   0   0   1  |
```

**Update:**

```
residual = z - H * x_pred
S        = H * Cov_pred * H^T + R
K        = Cov_pred * H^T * S^(-1)
x        = x_pred + K * residual
Cov      = (I - K*H) * Cov_pred * (I - K*H)^T + K*R*K^T   [Joseph form]
```

where measurement matrix H selects y and x:

```
H = | 1  0  0  0 |
    | 0  0  1  0 |
```

### Multiple Scattering (Q Matrix)

The process noise is measured by fitting the scattering angle distribution in Xe with a 4-Gaussian mixture. The core sigma `sigma_ms = 0.01631 rad` is used in Q:

```
Q = | h^2 * s2   h * s2     0          0      |
    | h * s2     s2          0          0      |
    | 0          0          h^2 * s2   h * s2  |
    | 0          0          h * s2     s2      |

where s2 = sigma_ms^2
```

### Momentum Reconstruction

After fitting, momentum is reconstructed from the bending angle between DetIn5 and DetOut1:

```
delta_tx = tx_out - tx_in
p [GeV/c] = 0.3 * B [T] * L_mag [m] / |delta_tx|
```

---

## Results

### Track Fitting Validation (e⁻, 1 GeV/c, Xe)

| Diagnostic | Result           | Expected  |
| ---------- | ---------------- | --------- |
| Pull y     | N(0, 0.9997)     | N(0, 1) ✓ |
| Pull x     | N(0, 1.001)      | N(0, 1) ✓ |
| Chi2 y     | mean = 9.998     | ~10 ✓     |
| Chi2 x     | mean = 10.03     | ~10 ✓     |
| p-value    | flat, mean=0.499 | uniform ✓ |

### Momentum Measurement

| Quantity               | Value       |
| ---------------------- | ----------- |
| p mean                 | 0.977 GeV/c |
| True momentum          | 1.0 GeV/c   |
| Bias                   | -2.3%       |
| Resolution sigma(dp/p) | 7.0%        |

---

## Project Structure

```
.
├── simulation/
│   ├── 10Detector_nonlinear_vacuum.g4bl   # G4BL simulation — vacuum
│   ├── 10Detector_nonlinear_xe.g4bl       # G4BL simulation — Xenon
│   └── convert_to_root.C                  # Convert ASCII → ROOT TTree
├── src/
│   ├── kf_rk4.C                           # Main KF + RK4 macro (ROOT)
│   └── multiple_scattering.C              # Fit scattering angle distribution
└── results/
    ├── pull_chi2_pval.pdf                  # Tracking validation plots
    ├── momentum.pdf                        # Momentum reconstruction
    └── multipleScattering_fitting.pdf      # 4-Gaussian fit of theta distribution
```

---

## Requirements

- [ROOT](https://root.cern) >= 6.x (C++17)
- [G4Beamline](http://g4beamline.muonsinc.com) (for simulation)

---

## Usage

```bash
# 1. Run G4BL simulation
g4bl simulation/10Detector_nonlinear_xe.g4bl

# 2. Convert ASCII output to ROOT
root -l simulation/convert_to_root.C

# 3. Fit multiple scattering angle distribution
root -l src/multiple_scattering.C

# 4. Run KF + RK4 track fitting
root -l src/kf_rk4.C
```

---

## Key Parameters

```cpp
// Physics
p_MeV      = 1000.0    // beam momentum [MeV/c]
B_T        = 6.54      // dipole field [T]
rho        = 510.2     // bending radius [mm]

// Geometry
h_drift_in  = 300.0   // DetIn spacing [mm]
h_B         = 1006.0  // arc length DetIn5→DetOut1 [mm]
h_drift_out = 415.0   // DetOut spacing [mm]
h_x_B       = 483.0   // Delta-x DetIn5→DetOut1 [mm]
h_x_out     = 303.0   // Delta-x between DetOut [mm]

// Measurement
sigma_det  = 0.1      // detector resolution [mm]
sigma_ms   = 0.01631  // multiple scattering core sigma [rad] — fitted from Xe data
```

---

## Author

**Vi Nguyen** — Graduate thesis  
Ho Chi Minh City University of Science
