# Space Mission Analysis Laboratory

See the [Final Report](Space_Mission_Analysis_Laboratory.pdf) for the complete derivations, numerical procedure, and results.

Preliminary mission design for rendezvous with the near-Earth asteroid **5797 Bivoj**, developed in MATLAB at Politecnico di Milano. The analysis connects transfer to the parking orbit, direct Earth–Bivoj transfer, and Earth-departure and asteroid-arrival hyperbolic trajectories under the **two-body approximation**, **impulsive maneuvers**, and **patched-conics approximation**.

## Transfer to the Parking Orbit

A standard sequence of orbital-plane change, change in the argument of periapsis, and orbital-shape change is compared with a **near-single-impulse maneuver**. The alternative exploits the near-intersection of the initial and final orbits: a small apoapsis correction produces a transfer orbit passing through $\mathbf{r}_I$, while preserving the initial orbital plane, orientation, and apoapsis radius. The main impulse at $\mathbf{r}_I$ completes the transfer, reducing $\Delta v$ from **4.483 to 1.821 km/s** and shortening the transfer time. A direct-transfer search with freely varying endpoint anomalies and argument of periapsis provides a numerical cross-check.

## Direct Transfer

The transfer geometry is parameterized by $\theta_{e,1}$, $\theta_{a,2}$, and $\omega_t$. The endpoint position vectors and $\omega_t$ determine the transfer orbit; Keplerian propagation from an assumed minimum-distance configuration couples Earth's and Bivoj's departure anomalies. For prescribed departure conditions, the numerical procedure combines local minimization of $\Delta v(\theta_{a,2},\omega_t)$ with respect to $\omega_t$ and root-finding on the time-delay function along the resulting branch, enforcing $\Delta T=0$. Repeating the procedure over discrete departure configurations identifies the lowest-cost solution within the explored family.

## Hyperbolic Trajectories

The **Earth-departure hyperbola** is reconstructed from $\mathbf{v}_\infty$ and $\theta_{p,H}$ through a nonlinear system for $e_h$, $\theta_{h,H}$, and $\theta_\infty$. Discrete variation of the maneuver point determines the most favorable injection configuration.

For the **asteroid-arrival hyperbola**, an asymptotic analysis in the limit $\mu/r_p\to0$ establishes negligible gravitational deflection and $\Delta v\to v_\infty$. Capture cost becomes essentially independent of the periapsis radius over the operational interval considered, allowing its selection through other operational requirements. Arrival is therefore interpreted as a **deep-space rendezvous maneuver**, requiring active reduction of the spacecraft's relative velocity to a near-zero value.

**Author:** Aleksandar Franchini Petrovic · Aerospace Engineering · Academic Year 2025–26
