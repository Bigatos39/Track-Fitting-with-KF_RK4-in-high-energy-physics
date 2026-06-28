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

$$\mathbf{x} = (y,\ t_y,\ x,\ t_x)^T$$

- **y, ty**: transverse coordinate and slope — unaffected by B field (B = Bŷ)
- **x, tx**: bending coordinate and slope — deflected by Lorentz force

### RK4 Equations of Motion

| Region  | dy/ds | dty/ds | dx/ds | dtx/ds |
| ------- | ----- | ------ | ----- | ------ |
| Drift   | ty    | 0      | tx    | 0      |
| B field | ty    | 0      | tx    | +1/ρ   |

### Kalman Filter Steps

**Predict:**
$$\mathbf{x}_{\text{pred}} = \text{RK4}(\mathbf{x}, h, \rho)$$
$$P_{\text{pred}} = F P F^T + Q$$

**Update:**
$$K = P_{\text{pred}} H^T (H P_{\text{pred}} H^T + R)^{-1}$$
$$\mathbf{x} = \mathbf{x}_{\text{pred}} + K(z - H\mathbf{x}_{\text{pred}})$$

### Multiple Scattering (Q matrix)

The process noise is measured by fitting the scattering angle distribution in Xe with a 4-Gaussian mixture. The core sigma σ_ms = 0.01631 rad is used in Q:

$$Q = \begin{pmatrix} h^2\sigma_{ms}^2 & h\sigma_{ms}^2 & 0 & 0 \\ h\sigma_{ms}^2 & \sigma_{ms}^2 & 0 & 0 \\ 0 & 0 & h^2\sigma_{ms}^2 & h\sigma_{ms}^2 \\ 0 & 0 & h\sigma_{ms}^2 & \sigma_{ms}^2 \end{pmatrix}$$

---

## Results

### Track Fitting Validation

| Diagnostic | Result                    | Expected     |
| ---------- | ------------------------- | ------------ |
| Pull y     | N(0, 0.9997)              | N(0, 1) ✓    |
| Pull x     | N(0, 1.001)               | N(0, 1) ✓    |
| χ² y       | mean = 9.998              | ~ndof = 10 ✓ |
| χ² x       | mean = 10.03              | ~ndof = 10 ✓ |
| p-value    | flat [0, 1], mean = 0.499 | uniform ✓    |

### Momentum Measurement

Momentum is reconstructed from the bending angle:

$$p = \frac{0.3 \cdot B \cdot L_{\text{mag}}}{|\Delta t_x|}$$

| Quantity           | Value                         |
| ------------------ | ----------------------------- |
| p mean             | 0.977 GeV/c (true: 1.0 GeV/c) |
| Resolution σ(Δp/p) | 7.0%                          |

---

## Files

| File                        | Description                                         |
| --------------------------- | --------------------------------------------------- |
| `10Detector_nonlinear.g4bl` | G4Beamline simulation — 10 detectors, dipole magnet |
| `convert_to_root.C`         | Convert G4BL ASCII output to ROOT TTree             |
| `kf_rk4.C`                  | Main KF + RK4 track fitting (ROOT macro)            |
| `Total_KF_4D_xe.pdf`        | Result plots — pull, chi2, p-value                  |

---

## Requirements

- [ROOT](https://root.cern) ≥ 6.x
- [G4Beamline](http://g4beamline.muonsinc.com) (for simulation)

---

## Usage

```bash
# 1. Run simulation
g4bl 10Detector_nonlinear.g4bl

# 2. Convert to ROOT
root -l convert_to_root.C

# 3. Run KF track fitting
root -l kf_rk4.C
```

---

## Key Parameters

```cpp
// Physics
p    = 1000 MeV/c   // beam momentum
B    = 6.54 T        // dipole field
rho  = 510 mm        // bending radius

// Geometry
h_drift_in  = 300 mm  // spacing between DetIn
h_B         = 1006 mm // arc length DetIn5→DetOut1
h_drift_out = 415 mm  // arc length between DetOut

// Measurement
sigma_det = 0.1 mm    // detector resolution
sigma_ms  = 0.01631 rad // multiple scattering (Xe, fitted)
```

---

## Author

**Vi Nguyen** — Graduate thesis project  
Ho Chi Minh City University of Science
