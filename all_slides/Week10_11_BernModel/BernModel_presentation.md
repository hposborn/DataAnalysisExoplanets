---
title: "Project 6 - Exoplanet Imaging"
---
[contents]

---

[P "Demographics of exoplanets"]

- On a per-system (or per-planet) basis, exosolar planetary systems appear extremely diverse[+]
- But, across populations, multiple trends have been identified in the underlying exoplanet distribution[+]

-v-

### Compact multi-planet systems are common

[i=kepler_orrery_kruse.gif]

- As many as 67\% of stars host multiple planets within 200d.[+]

-v-

### A "valley" in radius around 1.8Re

[i=Radius-valley.jpg]

-v-

### Giant planet occurrence peak at the snow line

[i=Fernandes_giant_planet_occurrence.png]

-v-

### Giant planet occurrence vs metallicity

[ism=Buchhave_metallicity_occurrence.png]

- Occurrence scales steeply as $f_{\text{giant}} \propto 10^{2[\text{Fe/H}]}$, supporting core accretion [+]
- Low-mass planets show weak or negligible host-star metallicity dependence [+]
<!--
-v-

### Giant planet–metallicity correlation

[ism=lopez_fortney_mass_Z_giants.jpg]-->

-v-

## Linking to formation

- All of these observations must be linked to underlying physical processes happening to planets, including:
  - Formation / core accretion[+]
  - Migration[+]
  - Evolution (e.g. evaporation, etc)[+]

*The question is therefore*:[+]
1) Can we model planet populations from **first principles**.[+]
2) Can this exercise match the distribution of observed planets[+]
3) Can the theoretical populations predict demographics[+]

This is the goal of **Planetary Population Synthesis**[+]

---

[P "Modelling planetary processes"]

### Global model of Planet Formation

A planetary population synthesis model must track:
1. The structure and evolution of the protoplanetary **gas** disk[+]
2. The structure and evolution of the disk of **solids** (dust, pebbles, planetesimals)[+]
3. The **accretion of solids** leading to the growth of the planetary solid core[+]
4. The **accretion of H/He** leading to the growth of the planetary gaseous envelope[+]
5. **Orbital migration** resulting from the exchange of angular momentum[+]
6. **N-body interactions** amongst (proto)planets[+]
7. The **interior structure** of the planets.[+]

-v-

### Global model of Planet Formation

[i=BernModelConsiderations.png]

-v-

### Global model of Planet Formation

* Modelling the full system in 3D would be impossibly slow, so each physical interaction is parameterised in 1D based on detailed physics & simulations, or simplification assumptions. E.g.:
  * Protoplanetary disc physics can be parameterised as equations of a 1D axisymmetric disk
  * Planets can be assumed spherically symmetric with hydrostatic interiors
* Key mass/energy/angular momentum parameters are tracked & conserved

<!--
### Architecture of the Generation III Bern Model
* **Global end-to-end coupling ([Alibert et al. 2005](https://arxiv.org/abs/astro-ph/0412444), [Emsenhuber et al. 2021a](https://arxiv.org/abs/2007.05561)):**
  * Simultaneously integrates **10–100 lunar-mass embryos** ($0.01 M_\oplus$) per disc from $t=0$ to $10 \text{Gyr}$ [+]
* **Three coupled temporal & physical regimes:**
  * **Disc phase ($0–10~\text{Myr}$):** 1D axisymmetric gas/planetesimal disc + Type I/II orbital migration + N-body interactions [+]
  * **Transition phase ($\sim 10–20~\text{Myr}$):** Disc dispersal, giant impact phase, and dynamical instabilities [+]
  * **Evolution phase ($20~\text{Myr}–10~\text{Gyr}$):** Long-term thermodynamic cooling, envelope contraction, and atmospheric photoevaporation [+]
* Each part of physics is approximated to 1D (i.e. distance from star) tracking key mass/energy/angular momentum conservations.-->

-v-

### Disc Physics
* We need to track how the surface density of gas & solids evolves in a disc (and, relatedly, how temperature, pressure, scale height, etc change).
* The Bern model uses **Viscous gas disc evolution** ([Lynden-Bell & Pringle 1974](https://ui.adsabs.harvard.edu/abs/1974MNRAS.168..603L)) modified by viscosity parameterised as $\alpha$.
  - Initial conditions include a magnetospheric cavity close to the star ($R_{\rm in}$), then a power-law increase (with power $p_g \approx -1$) followed by a drop-off at $R_{\rm out}$, with the dust-to-gas ratio $f{\rm dg}$ determining its composition.
  - For statistical power, disc parameters ($R_{\rm in}$, $R_{\rm out}$, $p_g$, $f{\rm dg}$) can be sampled from realistic observed distributions.

<!--* The Bern model uses:[+]
  * **Viscous gas disc evolution** ([Lynden-Bell & Pringle 1974](https://ui.adsabs.harvard.edu/abs/1974MNRAS.168..603L)) - 1D vertical-averaged $\alpha$-disc solving $\frac{\partial \Sigma_g}{\partial t} = \frac{3}{r}\frac{\partial}{\partial r}\left[r^{1/2}\frac{\partial}{\partial r}\left(\nu \Sigma_g r^{1/2}\right)\right] - \dot{\Sigma}_{\text{photo}}$ with $\nu = \alpha c_s H$ [+]
  * **Photoevaporative dispersal** ([Clarke et al. 2001](https://arxiv.org/abs/astro-ph/0101292), [Mordasini et al. 2012](https://arxiv.org/abs/1201.1036)) - Internal EUV wind + external FUV wind from birth cluster stars carve gaps and truncate outer discs within $1–10~\text{Myr}$ [+]
  * **Solid disc physics** ([Thiabaud et al. 2014](https://arxiv.org/abs/1403.5240)) - Planetesimal surface density $\Sigma_s$ evolves dynamically via aerodynamic gas drag and planetary stirring [+]
  * **Condensation fronts** - Stoichiometric equilibrium condensation tracks refractory/volatile icelines ($\text{H}_2\text{O}, \text{CO}, \text{CO}_2, \text{CH}_4$) to set bulk planet composition [+]-->

-v-

### Formation Physics (Solid Accretion)
* For **planetesimal accretion**, core groth occurs via the accretion of small planetestimals, $M_c = \Omega\Sigma_p R^2_{\rm capture} F_G$. 
* Hence, accretion n this **oligarchic growth regime** can be increased by: [+]
  - Close-in planets (faster Keplerian frequency $\Omega$)[+]
  - Higher surface density of planetesimals ($\Sigma_p$)[+]
  - Larger capture radius $R_{\rm capture}$ (from larger mass)[+]
  - More gravitational focussing ($F_G$, informed by planetesimal velocities).[+]
* Following all dust particles would be computationally impossible, so a few dozen **planetary embryos** are included as seeds[+]
* Pebble accretion progresses differently, and requires protoplanets to have grown large enough (100s km at 1AU) such that pebbles "settle" onto the bodies rather than approach ballistically.[+]
* For **gas accretion**, the planetary interior structure is required...[+]

-v-

### Planetary Interiors (and gas accretion)

* Numerical integration (in 1D) of spherically symmetric structure equations derived from mass conservation, hydrostatic equilibrium, energy conservation, & energy transport (see last week's lecture notes!).[+]
* In the initial **attached phase** ($M_p<10M_\oplus$), the outer boundary of a planet transitions directly to the disc.[+]
  * But radiative cooling (and contraction) of the planet results in an empty shell around the planet which must be filled by nebular gas... [+]
  * Hence the rate of growth is regulated by the rate of the envelope's (Kelvin-Helmholtz) cooling.[+]
* Runaway cooling at $M_p\sim10M_\oplus$ results in a **detached phase**, which is not limited by cooling.[+]
  * The Bondi accretion rate is now limited by the feeding zone, and up to 90% of nearby gas mass is accreted[+]
* **Gap opening** at high masses leads to a reduction in accretion.

-v-

### Orbital migration

* Masses interact gravitationally with the gaseous disc, exchanging angular momentum.
* In most cases, the planet loses angular momentum and moving inwards.
* **Type I** and **Type II** migration then occurs depending on whether the planet has opened a gap or not.
  - In some special locations (e.g. condensation fronts with rapid surface density changes), the type I migration may act outwards, concentrating material at ice lines, hence potentially trapping migrating planets.

* Planets and protoplanets can also undergo **mutual interactions** with each other
* This is typically modelled via **N-body simulations**.

<!--* **Core accretion & envelope structure ([Pollack et al. 1996](https://ui.adsabs.harvard.edu/abs/1996Icar..124...62P), [Fortier et al. 2013](https://arxiv.org/abs/1210.4009)):**
  * Planetesimal accretion in the oligarchic regime accounting for envelope-enhanced capture radii [+]
  * 1D hydrostatic interior equations solved at every timestep; runaway gas accretion triggers once $M_{\text{core}} \gtrsim M_{\text{crit}}$ ($\sim 10~M_\oplus$) [+]
* **Orbital migration & N-body dynamics ([Paardekooper et al. 2011](https://arxiv.org/abs/1007.3116), [Dittkrist et al. 2014](https://arxiv.org/abs/1402.5969)):**
  * **Type I (low-mass):** Lindblad + non-isothermal corotation torques create outward migration zones ("convergence traps") [+]
  * **Type II (gap-opening):** Viscous-coupled migration driven by disc angular momentum transport [+]
  * **Mercury symplectic integrator** tracks mutual gravitational scattering, resonant trapping, and giant impacts [+]-->

<!--
-v-

### Subsection 3: Evolution Physics
* **Long-term interior cooling ([Mordasini et al. 2012b](https://arxiv.org/abs/1201.1037)):**
  * Post-disc planets contract and cool over Gyr timescales via Kelvin-Helmholtz luminosity, governed by grain/molecular opacities [+]
  * Radiogenic decay ($^{40}\text{K}, ^{232}\text{Th}, ^{235}\text{U}, ^{238}\text{U}$) provides internal heating for rocky/icy cores [+]
* **Atmospheric escape ([Jin et al. 2014](https://arxiv.org/abs/1407.2375)):**
  * XUV-driven hydrodynamic photoevaporation strips $\text{H/He}$ envelopes from close-in sub-Neptunes, carving the **radius valley** and **Neptune desert** [+]
* **Tidal evolution & bloating:**
  * Star-planet tidal dissipation circularizes eccentric orbits and drives tidal inspiral inside the Roche limit [+]-->

-v-

[i=Bernmodel_confusogram.png]

-v-

### Initial conditions

- Initial conditions can be drawn from prior distributions to ensure every formation simulations has slightly varying starting points.[+]
- This typically follows variation seen in real-world parameters, including:[+]
  - Metallicity ($[Fe/H]=-0.02\pm0.22$) and corresponding dust-to-gas ratio ($f_{\rm dg}\sim0.02$})[+]
  - Initial disc mass ($\log_10{M_D/M_\odot} = -2 \pm 0.5$)[+]
  - Disc lifetime (following observed curves with $2.5\pm1.5$ Myr).[+]
- In addition, **planetary embryos** are placed in the disc (10-100 moon-mass protoplanets).[+]
- Usually, planetary population syntheses are performed for a fixed stellar mass of $1.0M_\odot$ (although M-dwarf simulations have also been performed)[+]

-v-

[i=BernModelConsiderations.png]

---

[P "Results of the Bern Model"]

[i=mordasini_bern_model.mov]

-v-

### Synthetic Mass–Period ($M–P$) Diagram
* **Distinct demographic loci ([Emsenhuber et al. 2021b](https://arxiv.org/abs/2007.05562)):**
  * **Inner Super-Earth / Sub-Neptune pile-up ($P \sim 3–100~\text{d}, M \sim 1–20~M_\oplus$):** Bounded at the inner disc edge (magnetospheric cavity, $\sim 0.03–0.1~\text{AU}$) [+]
  * **Cold Jupiter peak ($a \sim 1–5~\text{AU}, M \sim 0.3–10~M_J$):** Cores form beyond the water iceline where solid surface density jumps [+]
  * **Hot Jupiters ($P < 10~\text{d}$):** Formed via either early Type II migration or late high-eccentricity tidal migration [+]
  * **Sub-Jovian / Neptune desert ($P < 3~\text{d}, M \sim 20–100~M_\oplus$):** Scarcity caused by runaway accretion bypassing intermediate masses + XUV envelope stripping [+]

-v-

### Observing bias

* In order to match the predicted distribution with the real distribution, we also have to simulate **which planets would be detected**[+]
* For example, RVs require a semi-amplitude above some noise floor ($K\propto M_p\sin{i} M_s^{-2/3} . P^{-1/3} > \frac{\sigma_{\rm RV}}{\sqrt{N_{\rm RV}}}$)[+]
  - Simulation produces planetary inclinations, but we can sample multiple viewing angles to increase samples
* Transits require the combined transit depth above the noise floor [+]
  - Transit depth defined by $\delta = (\frac{R_p}{R_s})^2$[+]
  - Photometric noise $\sigma_p$ (in ppm/hr) converts to total noise as $t_{\rm survey}$: $\frac{\sigma_p}{\sqrt{t_D . N_{\rm tr}}} = \frac{\sigma_p}{\sqrt{t_D . P/t_{\rm survey}}} $[+]
  - Transit probability requires inclination as $b = \frac{a}{R_s}\cos{i} <1$
  - **Mutual inclinations** between planets determine whether we see a multi-system or not).
* Microlensing sensitivity is defined by mass ratio & distance (with detectability peaking at the "Einstein radius")[+]
* Imaging sensitivty is defined by contrast (light reflected/emitted from star as a fraction of the stellar flux) and separation (projected distance from the star: $\rm{sep}/{\rm arcsec} = \frac{a_{\rm pl}/\rm{ AU}}{d_{\rm star}/\rm{pc}}$)[+]

-v-

### Demographic comparisons (HARPS)

[i=BernModel_HARPS_comparison.png]

-v-

### Demographic comparisons (Radius valley)

[i=Burn_2023_RadiusValley.png]

- Suggested radius valley between rocky & water-rich planets ([Burn 2024](https://www.nature.com/articles/s41550-023-02183-7))

-v-

### Demographic comparisons (Metallicity trends)

[i=BernModel_MetallicityTrends.png]

-v-

### Four Emergent System Architectures
* **Classification from NGPPS II ([Emsenhuber et al. 2021b](https://arxiv.org/abs/2007.05562), [Emsenhuber et al. 2023](https://arxiv.org/pdf/2303.00012)):**
  * **Class I (Low-mass / Terrestrial & Super-Earth systems):** Low solid-mass discs ($M_s \lesssim 50~M_\oplus$); in-situ accretion + post-disc giant impacts; high multiplicity [+]
  * **Class II (Migrated Sub-Neptune chains):** Intermediate solid mass ($50–150~M_\oplus$); rapid Type I migration into resonant chains near the inner disc edge [+]
  * **Class III (Mixed Super-Earth + Cold Giant systems):** High solid mass ($>150~M_\oplus$); outer core crosses $M_{\text{crit}}$ to form a cold Jupiter without destabilizing inner planets [+]
  * **Class IV (Dynamically active / Eccentric Giant systems):** Multiple gas giants undergo planet-planet scattering, ejecting companions and leaving eccentric survivors [+]

-v-

### Four Emergent System Architectures

[i=BernModel_four_classes.png]

-v-

### Four Emergent System Architectures

[i=BernModel_four_architectures.png]

-v-

### Applications to M-Dwarfs & Low-Mass Stars
* **Stellar mass scaling study ($0.1–1.0~M_\odot$, [Burn et al. 2021](https://arxiv.org/abs/2105.04596)):**
  * Disc mass scales linearly with stellar mass ($M_{\text{disc}} \propto M_\star$), limiting total solid reservoir [+]
* **Key theoretical predictions:**
  * **Strong suppression of gas giants** below $M_\star \lesssim 0.5~M_\odot$ (slower accretion means cores fail to reach $M_{\text{crit}}$ before disc dispersal) [+]
  * **Abundant compact multi-planet systems** of Earth-to-super-Earth mass planets ($1–5~M_\oplus$) at short periods [+]
  * **High water mass fractions:** inward migration brings ice-rich embryos from beyond the closely packed snowline into the habitable zone [+]
<!--* **Tension with observations ([Schlecker et al. 2022](https://arxiv.org/abs/2205.12971)):**
  * Observed M-dwarf giant planet occurrence (~2–3% in CARMENES/RV) exceeds pure planetesimal-accretion predictions, motivating hybrid pebble accretion [+]-->

-v-

### Additional Physics & Ongoing Extensions
* **Hybrid Pebble + Planetesimal Accretion ([Voelkel et al. 2021](https://arxiv.org/abs/2004.03492), [Burn & Mordasini 2024](https://arxiv.org/abs/2410.00093)):**[+]
  * Drifting mm/cm pebbles accelerate outer core growth around low-mass stars and massive discs, resolving cold-giant formation timescales [+]
* **MHD Wind-Driven Disc Evolution ([Weder et al. 2023](https://arxiv.org/abs/2304.12380v1)):**[+]
  * Replaces pure turbulent $\alpha$-viscosity with magnetohydrodynamic disc winds, altering surface density profiles and slowing inward Type I migration [+]
* **Envelope Enrichment & Atmospheric Chemistry ([Marboeuf et al. 2014](https://arxiv.org/html/1609.00960v1), [Mordasini et al. 2016](https://iopscience.iop.org/article/10.3847/0004-637X/832/1/41/pdf)):**[+]
  * Ablation of accreted icy planetesimals enriches giant envelopes in C, O, and N, linking **C/O ratios** directly to formation location relative to icelines [+]
* **Binary stars ([Venturini et al 2026](https://arxiv.org/abs/2601.14214))**:[+]
  * Including binary stars truncates discs resulting in less planet mass and predictions for exoplanet demographics in binaries.[+]


---

[P "The project"]

## Observations about the data
- Stored on [DACE](https://dace.unige.ch/populationAnalysis)[+]
- **Three Nominal models** the NGPPS Bern model is split into three according to the number of embryos: NG74, NG75, NG76. [+]
  - You can effectively take all three as independent and comine their outputs[+]
- In addition, there exists models for sub-solar masses[+]
- **Extreme amount of data when including time** - it's best to either take snapshot at a given age, or only look at tracks for a few systems with detailed temporal range.[+]
- **Large number of protoplanets** (especially beyond the ice line) - you may want to cut low-mass planets when visualising[+]
- The "tracks" (when downloaded from DACE) lack a header - see the csv file in the code folder.[+]

-v-

## The example file

- Accesses "snapshots" at 5Gyr for a given model (NG74,NG75,NG76).
- Plots these with positions, mass & compositions.
- Introduces a handful of "Detection functions" determining which planets would be observable around a typical star in:
  - RVs (HARPS), Transits (Kepler)
  - And predictions for: Microlensing (Roman), Imaging (HabEx)
- Shows how to download evolution tracks for a subset of planetary systems (close-in intermediate-mass planets in this case - max of 10 to plot).

-v-

## Project Goals

- Get some understanding of how different types of planets form (e.g. their tracks)
- Develop some understanding of how some
- Do a deep dive into some part of the parameter space (and do some data science and visualisations with the outputs)

-v-

## Project ideas
- Can you find equivalent formation tracks for known multi-planet systems (like 55 Cancri, or GJ876) and what are their evolutionary paths?[+]
- What sort of planets get ejected?[+]
- How fast is inward migration for earthlike (and hot jupiter) planets?[+]
- Which protoplanets ended up furthest from their birthplace (`a_start`), and which stayed in the same place (and why)?[+]
- What changes with stellar mass? (NGM10-NGM14)[+]
- What are the most extreme cases (e.g. mass, radius, eccentricity, multiplicity, etc) in the models and are they represented in reality?[+]

-v-

## My own (ongoing) project:

- The Eccentricity distribution of small planets is difficult to measure, but especially important at long-periods (predicting an orbital period from a single transit requires informed eccentricity prior)
- As eccentricities are known precisely, the Bern model provides a good $e$ prior.

[ism=Eccentricity_from_BernModel.png]
