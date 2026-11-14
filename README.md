# Space Mission Analysis: Earth-5797 Bivoj

**Preliminary asteroid rendezvous design through impulsive orbit transfers, heliocentric trajectory construction, and patched-conic analysis.**

MATLAB implementation of a mission-analysis study to asteroid **5797 Bivoj**, developed for the *Space Mission Analysis Laboratory* at **Politecnico di Milano**. The analysis connects three design problems: transfer to an eccentric Earth parking orbit, construction of a heliocentric rendezvous trajectory, and characterization of the departure and arrival hyperbolae.

The selected heliocentric solution has a reported **68.24-day flight time** and a **7.193 km/s endpoint velocity-matching cost**. Resolving Earth escape and asteroid arrival gives a **7.603 km/s departure-and-capture budget**, excluding the initial transfer to the parking orbit.

![Heliocentric Earth-Bivoj transfer, showing the departure and rendezvous geometry.](figura_s2.3.png)

**Author:** Aleksandar Franchini Petrovic  
**Degree program:** Aerospace Engineering  
**Academic year:** 2025-2026  
**Academic advisors:** Prof. Mauro Massari and Prof. Fabio Ferrari

[Technical report](Space_Mission_Analysis_Laboratory.pdf) · [Mission model](#mission-model) · [Numerical methods](#numerical-methods) · [Reference results](#reference-results) · [Execution](#execution) · [Source files](#source-files)

## Mission model

The trajectory is represented by Keplerian conic arcs under Newtonian two-body dynamics. Maneuvers are instantaneous velocity discontinuities; the interplanetary and planetocentric descriptions are connected through their hyperbolic excess velocities.

| Analysis segment | Central body | Reference frame | Design objective |
| :--- | :--- | :--- | :--- |
| Initial orbit to parking orbit | Earth | Geocentric equatorial inertial | Compare a standard four-impulse sequence with a two-impulse near-intersection transfer |
| Interplanetary cruise | Sun | Heliocentric ecliptic | Construct an elliptic transfer satisfying the adopted rendezvous phasing |
| Earth departure | Earth | Geocentric equatorial inertial | Match the outgoing excess-velocity vector while reducing injection cost |
| Asteroid arrival | Bivoj | Asteroid-centered orbital plane | Estimate the impulse required to match the asteroid and enter a prescribed circular orbit |

The Earth parking orbit is eccentric: its semi-major axis is approximately **46,750 km**, with eccentricity **0.7507** and inclination **39.65°**. It is obtained from the prescribed Cartesian state rather than assumed circular.

### Units and notation

| Quantity | Code notation | Unit / convention |
| :--- | :--- | :--- |
| Semi-major axis | `a` | km; negative for hyperbolic trajectories |
| Eccentricity | `e` | Dimensionless |
| Inclination | `i` | rad |
| Right ascension of the ascending node | `OM` | rad |
| Argument of periapsis | `om` | rad |
| True anomaly | `th` | rad |
| Position and velocity | `r`, `v` | km and km/s; three-component column vectors |
| Gravitational parameter | `mu` | km³/s² |
| Maneuver cost | `DV` | km/s |
| Time of flight | `DT`, `tT` | s; converted to hours or days for reporting |

Angles are in radians unless explicitly stated otherwise. Calls to `rotx` and `rotz` convert angles to degrees; the ecliptic-to-equatorial rotation uses a prescribed obliquity of **23.45°**. Suffixes `i`, `f`, `t`, `p`, and `h` identify initial, final, transfer, parking, and hyperbolic orbits. In the heliocentric analysis, `e` and `a` identify Earth and the asteroid, while `1` and `2` identify departure and arrival.

## Numerical methods

### 1. Geocentric transfer to the parking orbit

Two maneuver architectures are evaluated for the same initial and target orbital geometries.

**Standard sequence.** The baseline separates the geometry changes into an orbital-plane rotation, an apsidal rotation, and an apoapsis-to-periapsis bitangent transfer. The last operation requires two impulses, giving four impulses in total. Velocity increments are obtained analytically, and coast times are evaluated through the eccentric anomaly and Kepler's equation.

**Near-single-impulse sequence.** The alternative exploits the proximity of the two orbits near their plane-intersection line. The main-impulse point on the final orbit is determined by solving

$$
\widehat{\mathbf h}_i\cdot\mathbf r_f(\theta_f)=0,
$$

where $\widehat{\mathbf h}_i$ is the initial orbit's unit angular-momentum vector. The selected root defines a transfer orbit constrained to retain the initial orientation and apoapsis radius:

$$
i_t=i_i,\qquad \Omega_t=\Omega_i,\qquad \omega_t=\omega_i,\qquad r_{a,t}=r_{a,i}.
$$

A small tangential correction at the common apoapsis establishes the intersection; the subsequent vector impulse inserts the spacecraft into the parking orbit. Despite its name, this architecture contains **two finite impulses**. For the reference case, the correction is approximately **4 m/s**, compared with **1.817 km/s** for the main maneuver.

![Near-single-impulse transfer between the initial geocentric orbit and the parking orbit.](figura_s1.6.png)

### 2. Heliocentric transfer construction and rendezvous phasing

For a prescribed Earth departure anomaly $\theta_{e,1}$, the transfer is parameterized by the asteroid arrival anomaly $\theta_{a,2}$ and transfer argument of periapsis $\omega_t$. The two endpoint position vectors define the transfer plane. After transforming those vectors to the corresponding perifocal frame, `direct.m` reconstructs the conic from

$$
e_t=\frac{r_1-r_2}{r_2\cos\theta_{t,2}-r_1\cos\theta_{t,1}},
\qquad
a_t=\frac{r_1\left(1+e_t\cos\theta_{t,1}\right)}{1-e_t^2},
$$

with $r_j=\lVert\mathbf r_j\rVert$. The design branch is elliptic, with $a_t>0$ and $0<e_t<1$. This is a geometric conic construction with a separate timing calculation, rather than a prescribed-time Lambert formulation.

The heliocentric endpoint cost is

$$
J(\theta_{a,2},\omega_t)=
\left\lVert\mathbf v_{t,1}-\mathbf v_{e,1}\right\rVert+
\left\lVert\mathbf v_{a,2}-\mathbf v_{t,2}\right\rVert.
$$

The rendezvous condition requires the spacecraft and asteroid flight times to agree:

$$
\Delta T(\theta_{a,2},\omega_t)=
\Delta t_t(\theta_{a,2},\omega_t)-\Delta t_a(\theta_{a,2})=0.
$$

The implemented search proceeds as follows:

1. Minimize the geometric Earth-asteroid orbital separation using `fminsearch` and adopt the resulting configuration as the phasing reference.
2. Propagate the asteroid backward to the prescribed Earth departure anomaly by solving Kepler's equation with `fsolve`.
3. Inspect the eccentricity and cost curves to choose an elliptic transfer branch and an initial guess for $\omega_t$.
4. For each trial arrival anomaly, minimize $J$ over $\omega_t$ with `fminsearch`; `transferTime.m` returns the resulting flight time.
5. Minimize the absolute timing residual over $\theta_{a,2}$ with an outer `fminsearch`, then recover the transfer elements and endpoint states.
6. Compare the recorded solutions for discrete Earth departure anomalies. The selected sample is $\theta_{e,1}=340^\circ$.

This is a **nested local search**. The code screens elliptic admissibility graphically and uses unconstrained minimization; it does not establish a global optimum or solve the general equality-constrained minimum-cost problem. The departure-anomaly comparison plotted by the script contains recorded numerical results rather than an automated parameter sweep.

### 3. Patched-conic departure and arrival

**Earth departure.** The heliocentric departure velocity difference supplies the geocentric excess-velocity vector after rotation into the equatorial frame. Its magnitude determines the departure hyperbola's semi-major axis:

$$
\mathbf v_{\infty,E}=R_x(23.45^\circ)
\left(\mathbf v_{t,1}-\mathbf v_{e,1}\right),
\qquad
a_h=-\frac{\mu_E}{v_{\infty,E}^2}.
$$

For a prescribed parking-orbit injection anomaly, a three-equation nonlinear system enforces the conic radius, the hyperbolic asymptote condition, and the angular separation between the injection radius and outgoing asymptote. `fsolve` returns the hyperbolic eccentricity and the injection and asymptotic anomalies. The impulse is evaluated from the full velocity-vector difference at the common injection position. The recorded discrete search selects an injection anomaly of **177.5°** on the parking orbit.

**Asteroid arrival.** With the prescribed mass $M_a=6.7021\times10^{10}$ kg, the model uses $\mu_a=GM_a\approx4.4732\times10^{-9}$ km³/s². For a circular capture orbit of radius $r_p$, the arrival impulse is

$$
\Delta v_{\mathrm{cap}}=
\sqrt{v_{\infty,a}^2+\frac{2\mu_a}{r_p}}
-\sqrt{\frac{\mu_a}{r_p}}.
$$

As $\mu_a/r_p\to0$, the turning angle approaches zero and $\Delta v_{\mathrm{cap}}\to v_{\infty,a}$. Within this model, arrival is dominated by active relative-velocity cancellation. The reference visualization uses **$r_p=2.5$ km**, measured from the asteroid's center.

## Reference results

The following values summarize the accompanying report and the numerical reference inputs embedded in the scripts. Values are rounded for presentation.

### Parking-orbit transfer

| Architecture | Impulses | Total Δv [km/s] | Transfer-arc time [h] |
| :--- | ---: | ---: | ---: |
| Standard maneuver sequence | 4 | 4.483 | 13.75 |
| Near-single-impulse sequence | 2 | 1.821 | 2.79 |

The alternative reduces the reported transfer cost by approximately **59.4%**. The times refer to the transfer arcs accumulated by the scripts; they exclude the initial coast to the first maneuver, subsequent parking-orbit phasing, and launch operations.

### Selected heliocentric transfer

| Parameter | Reference value |
| :--- | ---: |
| Earth departure true anomaly | 340.00° |
| Asteroid arrival true anomaly | 43.51° |
| Transfer semi-major axis | 204.81 × 10⁶ km |
| Transfer eccentricity | 0.2813 |
| Transfer inclination | 2.37° |
| Transfer RAAN | 261.05° |
| Transfer argument of periapsis | 180.07° |
| Flight time | 68.24 days |
| Heliocentric departure velocity difference | 3.967 km/s |
| Heliocentric arrival velocity difference | 3.226 km/s |
| Total heliocentric endpoint cost | 7.193 km/s |

### Departure-and-capture budget

| Contribution | Δv [km/s] |
| :--- | ---: |
| Parking-orbit injection into the Earth departure hyperbola | 4.377 |
| Asteroid arrival and circular-orbit insertion | 3.226 |
| **Patched-conic departure and capture** | **7.603** |

The heliocentric endpoint cost and patched-conic budget describe alternative accounting levels and must not be added together. Including the selected **1.821 km/s** initial parking-orbit transfer gives an arithmetic subtotal of **9.424 km/s**, before launch, corrections, margins, and any additional phasing maneuvers. This subtotal does not constitute an integrated, epoch-synchronized mission budget.

## Execution

### Requirements

- **MATLAB R2021a or later** for the `Name=Value` plotting syntax used in the scripts.
- **Optimization Toolbox** for [`fsolve`](https://www.mathworks.com/help/optim/ug/fsolve.html). The local minimizations use MATLAB's `fminsearch`.
- Degree-based [`rotx`](https://www.mathworks.com/help/phased/ref/rotx.html) and [`rotz`](https://www.mathworks.com/help/phased/ref/rotz.html) functions on the MATLAB path, for example from Phased Array System Toolbox. Their convention must match the documented active rotation matrices.

All numerical scenario inputs are embedded in the source. No runtime ephemeris download or external data service is required.

### Run an analysis

Keep the supplied `.m` files together in the repository root, alongside this README, the report, and the original PNG figures. Open that directory in MATLAB and check function availability:

```matlab
addpath(pwd);
which fsolve -all
which rotx -all
which rotz -all
```

Run one entry-point script in a clean workspace. For example, to compute the heliocentric transfer:

```matlab
clearvars;
close all;
clc;

section2;

fprintf('Heliocentric endpoint cost: %.3f km/s\n', DV);
fprintf('Transfer time: %.2f days\n', tT / 86400);
```

The scripts create figures and leave their calculated variables in the base workspace. They reuse variable names and figure numbers, so inspect or save each result before running another entry point. The heliocentric script includes nested minimizations inside its diagnostic sampling loops and can take longer than the other analyses.

| Entry point | Calculation | Principal workspace outputs |
| :--- | :--- | :--- |
| [`section1_project1.m`](section1_project1.m) | Standard four-impulse parking-orbit transfer | `DV1`-`DV4`, `DV`, `DT1`-`DT3`, `DT` |
| [`section1_project2.m`](section1_project2.m) | Near-single-impulse parking-orbit transfer | `DV1`, `DV2`, `DV`, `DT`, `a_t`, `e_t`, `thI_f` |
| [`section2.m`](section2.m) | Heliocentric transfer and rendezvous phasing | `DV`, `tT`, `a_t`, `e_t`, `i_t`, `OM_t`, `om_t`, `th2_a` |
| [`section3_project1.m`](section3_project1.m) | Earth departure hyperbola | `DV`, `a_h`, `e_h`, `i_h`, `OM_h`, `om_h`, `rP`, `del`, `DEL` |
| [`section3_project2.m`](section3_project2.m) | Asteroid arrival and circular capture | `DV`, `a_h`, `e_h`, `rP`, `vP_h`, `vP_c`, `del`, `DEL` |

The scripts are independently executable with their embedded reference inputs. They do not form an automated data pipeline: the departure excess-velocity vector in `section3_project1.m` and arrival excess speed in `section3_project2.m` are stored numerical inputs. After changing the heliocentric solution, update both before interpreting the patched-conic results. Likewise, orbital data are duplicated in `section2.m` and `transferTime.m` and must remain consistent.

## Source files

| Function | Responsibility |
| :--- | :--- |
| [`par2car.m`](par2car.m) | Keplerian elements to inertial Cartesian position and velocity; supports elliptic and hyperbolic conics |
| [`car2par.m`](car2par.m) | Cartesian state to classical orbital elements for nondegenerate orbital geometries |
| [`por2for.m`](por2for.m) | Elliptic-orbit reconstruction from a position, prescribed orientation, and apoapsis radius |
| [`changeOrbitalPlane.m`](changeOrbitalPlane.m) | Plane-change geometry and impulse with preserved orbit shape |
| [`changePericenterArg.m`](changePericenterArg.m) | Apsidal-rotation impulse and the two orbit-intersection anomaly pairs |
| [`bitangentTransfer.m`](bitangentTransfer.m) | Signed tangential velocity increments and half-period flight time for the selected apsidal transfer |
| [`direct.m`](direct.m) | Transfer-conic elements from two position vectors and a prescribed argument of periapsis |
| [`transferTime.m`](transferTime.m) | Scenario-specific inner cost minimization and flight-time evaluation |
| [`timeOfFlight.m`](timeOfFlight.m) | Forward elliptic time of flight, including a single anomaly wrap |
| [`orbitProp.m`](orbitProp.m) | Forward or backward elliptic anomaly propagation through Kepler's equation |
| [`plotOrbitPlane.m`](plotOrbitPlane.m) | Planar conic sampling for visualization |
| [`plotOrbitSpace.m`](plotOrbitSpace.m) | Inertial three-dimensional conic sampling for visualization |

The standard state-conversion interface is `[r, v] = par2car(...)`. Passing the optional final argument `1` reverses the output order to `[v, r]`; the transfer scripts use this mode to construct velocity-valued anonymous functions. Angles returned by the helpers are in radians. `bitangentTransfer.m` returns signed tangential increments; use their absolute magnitudes when assembling budgets for other transfer configurations.

## Modeling scope and interpretation

- **Assumed reference phasing.** The study associates the minimum-distance orbital geometry with **9 January 2092**, then propagates fixed Keplerian elements relative to that configuration. The code does not ingest epoch-tagged ephemerides or verify simultaneous occupancy of those geometric minimum-distance points. Calendar dates inferred from this assumption are illustrative rather than validated launch windows.
- **Local numerical solutions.** Solver initial guesses are chosen through diagnostic plots, and the displayed departure-anomaly comparisons use stored samples. Solver defaults, branch selection, and the timing residual affect reproducibility. For a modified scenario, check convergence, elliptic admissibility, endpoint-position agreement, and the final timing mismatch explicitly.
- **Idealized dynamics and propulsion.** The analysis excludes third-body perturbations within each conic arc, nonspherical gravity, solar radiation pressure, finite-duration burns, thrust and mass constraints, and operational margins. The patched-conic calculation uses asymptotic velocity matching without explicit finite sphere-of-influence boundary integration.
- **Small-body proximity operations.** The circular capture orbit is a point-mass reference construction. Its stability under solar perturbations, the asteroid's irregular gravity field, rotation, and shape is not assessed. The arrival analysis estimates velocity matching, not the feasibility of sustained orbital operations.
- **Function domains.** `timeOfFlight.m` and `orbitProp.m` are elliptic routines. `car2par.m` uses classical elements and is singular for circular or equatorial cases; `direct.m` assumes non-collinear endpoints and a non-equatorial transfer plane. These helpers should be used within their stated geometric domains.

The accompanying report contains the derivations, intermediate orbital configurations, diagnostic plots, and tabulated reference solutions. The repository documents a preliminary academic trajectory-design study; spacecraft sizing, launch-system selection, and flight-dynamics validation require additional analysis.
