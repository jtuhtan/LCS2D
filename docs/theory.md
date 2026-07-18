# Theory & numerics

A compact but complete account of what LCS2D computes and how. Notation
follows Shadden et al. (2005).

## The flow map

A velocity field **v**(x) advects a fluid particle according to

    dx/dt = v(x),    x(t0) = x0 .

The **flow map** Φᵀ(x0) sends an initial position x0 to its position after
integration time T. All LCS diagnostics are properties of this map: they
ask *how strongly neighboring initial positions are pulled apart*.

Because the fields here are steady, only the integration duration matters
(not the start time), and trajectories coincide with streamlines — the
interesting stretching happens near hyperbolic (saddle-type) stagnation
points and the separatrices connecting them.

## FTLE — finite-time Lyapunov exponent

Linearizing the flow map about x0 gives the deformation-gradient
(Jacobian) matrix J = ∇Φᵀ(x0). The **right Cauchy–Green strain tensor**

    C = Jᵀ J

is symmetric positive semi-definite; its largest eigenvalue λmax(C) is the
squared factor by which an infinitesimal disc around x0 is stretched along
its most-stretching direction. The FTLE is the average exponential rate:

    σ(x0) = ln λmax(C) / (2 |T|) .

Ridges of σ computed **forward** in time (T > 0) mark *repelling* material
lines; ridges of σ computed **backward** (T < 0) mark *attracting* ones.
Together they form the transport skeleton of the flow.

### The auxiliary-grid discretization

J is approximated by central differences over four auxiliary particles
placed a small distance δ from each seed (up/down/left/right — the
"secondary seeds" of the legacy MATLAB code). After advecting all five
points:

    J = [ (x⁺ − x⁻)/2δ   (xᵘ − xᵈ)/2δ ]
        [ (y⁺ − y⁻)/2δ   (yᵘ − yᵈ)/2δ ]

where superscripts denote the particles started right/left/up/down of the
seed. λmax of the 2×2 tensor C is evaluated in closed form. Because δ is
independent of the seed spacing, ridges sharper than the seed grid are
still detected — the classic advantage of the auxiliary grid over
differencing neighboring seeds.

Choose δ well below the seed spacing (default: `res/100`). The computed
field should be *insensitive* to δ; if it is not, δ is too large (the
difference is no longer local) or absurdly small (floating-point
cancellation — position differences approach machine epsilon).

## FSLE — finite-size Lyapunov exponent

The FSLE asks the converse question: not "how far do particles separate in
a fixed time" but **"how long until they separate by a fixed factor?"**
With initial separation d0 = δ and threshold df = α·d0, the first time τ
at which *any* of the four auxiliary particles is at least df from the
primary gives

    λ(x0) = ln(α) / τ .

Seeds that never reach the threshold within `t_max` are assigned λ = 0
(and τ = NaN in the output). Unlike the FTLE offset, d0 is a *physical*
scale — sweeping it probes stretching at different scales (Aurell et
al. 1997; d'Ovidio et al. 2004).

## Numerical integration

Trajectories are integrated with the classical 4th-order **Runge–Kutta**
scheme (local error O(dt⁵)); `--integrator euler` selects forward Euler,
which is what the legacy MATLAB code used, for comparison only. The
requested time is hit exactly by rounding to a whole number of steps and
adjusting dt accordingly.

The default step obeys a CFL-style bound of dt = min(dx, dy)/max|v| for
FTLE — at most one grid cell per step. With RK4 the step is limited by
accuracy, not stability, and this default was verified on the
double-gyre benchmark to match a 4×-finer step to ~10⁻⁴ relative error
(halve `--dt` once on your own data as the same self-check). FSLE and
plain advection use the more conservative 0.5·min(dx, dy)/max|v|: the
FSLE separation time τ is quantized in units of dt, so its step size is
a measurement-resolution choice, not just an accuracy one.

**Velocity sampling** is bilinear on the regular grid: exact for linear
fields, C⁰ across cell edges, and cheap enough to evaluate four times per
step for RK4. Points outside the domain (or in NaN-masked cells) return
zero velocity, so escaping particles freeze at the boundary. Interpret
ridges that hug the domain edge with this rule in mind, and use `--margin`
to keep seeds away from unreliable borders.

**Backward time** is handled by integrating with a negative dt — no field
manipulation needed, and forward/backward runs use identical machinery.

## Verification

The test suite (and `lcs2d selftest`) checks the implementation against
flows with exact answers:

* **Uniform strain** (u = ax, v = −ay): FTLE = a everywhere; FSLE crossing
  time τ = ln(α)/a. The engine matches to ~10⁻¹² (FTLE) and to the
  time-discretization limit (FSLE).
* **Rigid rotation** (u = −y, v = x): closed orbits return to their start
  after 2π; RK4 error ~10⁻⁷ at dt = 0.01.
* **Double gyre**: the forward-FTLE ridge lies on the separatrix x = 1.
* The NumPy and Numba backends must agree to machine precision.

## References

* Shadden, S.C., Lekien, F., Marsden, J.E. (2005). Definition and
  properties of Lagrangian coherent structures from finite-time Lyapunov
  exponents in two-dimensional aperiodic flows. *Physica D* 212, 271–304.
* Haller, G. (2015). Lagrangian coherent structures. *Annual Review of
  Fluid Mechanics* 47, 137–162.
* Aurell, E., Boffetta, G., Crisanti, A., Paladin, G., Vulpiani, A.
  (1997). Predictability in the large: an extension of the concept of
  Lyapunov exponent. *J. Phys. A* 30, 1–26.
* d'Ovidio, F., Fernández, V., Hernández-García, E., López, C. (2004).
  Mixing structures in the Mediterranean Sea from finite-size Lyapunov
  exponents. *Geophys. Res. Lett.* 31, L17203.
* Tuhtan, J.A. (2012). [Investigating Lagrangian Coherent Structures in
  Environmental Flows](https://www.academia.edu/5643321/Investigating_Lagrangian_Coherent_Structures_in_Environmental_Flows)
  — background for the original MATLAB implementation.
