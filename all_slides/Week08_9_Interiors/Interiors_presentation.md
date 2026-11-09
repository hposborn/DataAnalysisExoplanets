---
title: "Project 4 - Interiors"
---

# Modelling planetary interiors

[contents]

---

[P "The building blocks of planets"]

- Angular momentum conservation during star formation means discs form[+]
- Dust falls to the midplane and can clump together, growing into protoplanets[+]

[i=ALMA-DSHARP.png]

-v-

## Planet formation

- The **elemental abundances** of material in the disc should be **identical to the star**.[+]
  - We see similar elemental abundance ratios (Fe, Si, Mg, etc) across all terrestrial planets.[+]
- Solid dust grains lack H, He, Ne, Ar, etc[+]
- **Dust composition is modified with orbital distance** due to sublimation & condensation of volatiles at ice lines (H$_2$O, but also CO$_2$, CH$_4$, etc).[+]
  - Close-in material is only rocky (Fe, Mg, Si, O)[+]
  - Dust in distant regions is more volatile-dominated (C, N, S, O).[+]

-v-

## Planet formation

- The change in compositions with distance is clearly visible in the difference between rocky planets, icy moons, & comets

[i="SolidMassFractions.png][+]

-v-

## Growth

- Pebble-sized dust clumps together into **planetesimals** - gravitationally bound and not blown around by the gas ($>1$ km, 100kyr).[+]
- Collisions between planetesimals and/or accretion of smaller material from the disc (pebbles) grows these into **protoplanets** ($>1000$km, $\sim1$Myr) and beyond.[+]
  - Protoplanets undergo **differentiation** which sorts planets by density.[+]
  - Typically: Core = free iron (and siderophile elements); mantle = magnesiums silicates; atmosphere = volatiles & H-He[+]
  - Once cooled, the atmosphere may further separate into ice/ocean layers under H-He (e.g. Ice Giants)[+]

[ism="https://upload.wikimedia.org/wikipedia/commons/1/12/MIT_12.007_Geobiology%2C_Spring_2007_-_Differentiation_of_Earth.jpg"][+]

-v-

## Growth in the gas disc

- Planets that grow massive enough ($\sim 10M_\oplus$) inevitably start to accrete nebular gas extremely efficiently - _runaway accretion_ of gas giants.[+]
- Small ($\sim1M_\oplus$) planets can also accrete some gas ($<1$%) from the gas disc producing a primordial atmosphere of H-He.[+]
- After $\sim10$Myr, the gas disc is evaporated by stellar high-energy radiation, stopping any further accretion (either of gas or pebbles)[+]
- Protoplanets can continue to grow via giant impacts (e.g. Earth's moon-forming impact)[+]

-v-

## Migration

- Planets peturb the disc around them, but as there is more dust outside this leads to a deceleration of the planet causing migration inwards.[+]
  - This can be rapid ($\ltapprox 1$Myr), especially for large ($1-20M_\oplus$) planets [+]

- Giant planets ($>0.1M_J$) are large enough to open a gap in the disc, but torques with the gap edge still lead to (slower) inward migration[+]
- Other planets can cause large influence on planetary position:[+]
  - An exterior planet may halt migration (e.g Jupiter-Saturn) [+]
  - Gravitational interactions can also scatter planets in/out after the gas disc is stripped [+]

-v-

## Migration

[i=Migration_Ataiee1_eng.gif]

-v-

## Mass loss

- Planets can also lose material:[+]
  - Close-in planets can lose mass through **evaporation** of gaseous envelope (H/He)[+]
  - **Collisions/giant impacts** can also remove material, especially outer layers (e.g. Mercury's large core?) [+]

-v-

## The fundamentals

- There are four fundamental reservoirs:
  - **Iron**[+]
  - **Silicates** (sometimes considered with Iron for simplicity)[+]
  - **Volatiles** (dominated by water, especially for formation $<5$AU; sometimes assume equal parts ice:rock)[+]
  - **Gas** (H & He; can be assumed to be protosolar ratio, i.e. 27% He by mass)[+]
- **Small planets** ($\lessapprox15M_\oplus$) are dominated by material from **dust**, and:[+]
  - Orbital distance of formation determines **volatile content**[+]
  - Stellar abundances determine iron/rock content (e.g. core-mass fraction, CMF)[+]
- Large planets are dominated by **H-He gas**[+]

---

[P "Modelling interior structures"]

## Observations

- Typically we have (if we are lucky) two direct measurements of the planet:[+]
  - Planetary mass (from e.g. RVs or TTVs)[+]
  - Planetary radius (from e.g. transits)[+]
  - (In addition, planetary orbit/temperature)[+]
- These, combined with other assumptions, can give us useful information on what combinations of the four materials might be present[+]
- Such models are complex but the data is extremely limited, meaning interior structure modelling is always highly **degenerate**/**underconstrained**/**prior-dominated**.[+]

-v-

## Giant planet modelling

- Given the mass is dominated by gas, the simplest assumption is to only model two materials: a gas envelope, and some fraction of rock/iron/volatiles[+]
- Volatiles : rock+iron ratio (ice-to-rock) is assumed 1:1[+]
- Need to assume something about the core-envelope boundary:[+]
  - Fully separated (with some small fraction of metals in the atmosphere $Z_{\rm atm}$)[+]
  - Fuzzy/dilute core, as suggested for Jupiter (Juno observations)[+]
  - Fully mixed[+]

[ism=giant_interior_core_envelope.png][+]

-v-

## Modelling Assumptions

- Modelling is built on the following key differential equations:
  - **Conservation of mass**: $ \frac{\partial m}{\partial r} = 4\pi r^2 \rho $ [+]
  - **Hydrostatic equilibrium**: $\frac{\partial P}{\partial r} = - \rho g $ - Gravity is balanced by pressure gradient.[+]
  - **Energy transport**: $\frac{\partial T}{\partial r} = \frac{\partial P}{\partial r} \frac{T}{P}\nabla_T$ - Hotter layers transport thermal energy to cooler layers[+]
  - **Energy conservation**: $\frac{\partial L}{\partial r} = 4\pi r^2 \rho(\dot{\epsilon} - T\frac{\partial S}{\partial t}) $ - L is the local intrinsic luminosity, & $\frac{\partial S}{\partial t}$ is the rate of change of specific entropy with time. [+]
    - $\dot{\epsilon}$ is nuclear energy (e.g. fusion or decay; $\dot{\epsilon}\sim0$) [+]
-v-

## Thermal effects
- The interiors of planets are hotter than their surfaces[+]
  - Planets also start substantially hotter (i.e. higher local intrinsic luminosity), hence cool & contract with act.[+]
- Thermal transport outwards can occur through radiation, conduction, or, (when these are insufficient), convection.[+]
- The **opacity** of the atmosphere (especially near the blackbody peak for a layer's temperature) determines if radiation is efficient.[+]
  - In the upper atmosphere of giants ($T\lessapprox1200$K) H$_2$O, CH$_4$, and collision-induced absorption by H$_2$.[+]
  - At deeper and hotter layers, metals (e.g. K, Na, TiO) and ions (e.g. $H^{\rm -}$) can cause convective layers[+]
<!--- The temperature gradient, $\nabla_T$, depends on the process that dominates the energy transport - adiabatic (simplest), or convective (which may also mix compositional gradients)[+]-->

-v-

## Equations of State

- In all materials, increasing pressure leads to increased temperature[+]
- In some materials, moving to a different region of P-T space causes fundamental phase change (e.g. solid-liquid-gas)[+]
- Knowing the exact P-T curve and phase changes that a material follows in P-T space is key to understanding planetary interiors.[+]
- For this, we need to know the equations of state; $f\left(p,V,T\right)=0$[+]

[ism=hydrogen_P_T.png][+]
 
-v-

## Equations of State

- Classically, we would measure a material's pressure & volume at different temperatures.[+]
- But giant planet interiors require extreme pressures & temperatures not reachable in the lab.[+]
- Additional complexities include its partially-degenerate fluid state, and interactions of atoms, molecules and ions.[+]
  - Metallic hydrogen should be soluble to rock-forming minerals (Fe, MgO, SiO$_2$, etc) in Jupiter's core.[+]
  - H$_2$O & H/He are also soluble at moderate temperatures/pressures.[+]
- Varying amounts of heavy elements also add complexity - often pure EoSs are "mixed" together[+]
- _Ab initio_ calculations are used to compute EoS at extreme pressures[+]

-v-

## Forward modelling interiors

<!--- Hydrostatic equlibrium-->
**We cannot typically take the observed mass & radius and compute expected composition**[+]
- Instead we need to take some combination of materials (e.g. 2Me iron, 3Me silicates, 4Me water, 0.1Me H-He) and ask _what radius would this planet have?_.[+]
- Modelling the interior requires the following:[+]
  1) Guess a central pressure[+]
  2) Build the core layer outwards following the EoS (and tracking e.g. p, T, V, phase, etc) until the total layer mass is reached[+]
  3) At a new layer, switch to the corresponding EoS[+]
  4) When total mass is reached, calculate final pressure ($P_{\rm surf}$) & radius[+]
  5) Iterate this process until $P_{\rm surf}=0$, producing the expected radius for the given input parameters[+]

-v-

## Constraining mass fractions 

- We really care about the fraction of a planet in each material type [+]
- To get this, we need to sample the forward model for different input parameter prior distributions [+]
- We can then use this to get a **posterior** for each material telling us what the allowed mass fractions are[+]

-v-

## Mass regimes 1) Rocky planets

- Planets with earthlike radii/mass can be assumed to have no H/He.[+]
- This is especially true for highly irradiated planets[+]
- The most simplistic approach would be to assume they are only formed of an iron core and silicate mantle.[+]
  - This approach is linear, allowing a simple _inversion_ to estimate the CMF.[+]
- This has enabled, for example, tests of whether Fe & Mg abundances in exoplanet host stars cause [observable changes in core mass fractions](https://arxiv.org/pdf/2102.12444)[+]

[ism=PlanetInteriorTrends.png][+]

-v-

## Mass regimes 1) Rocky planets

[i=mass_radius_trappist1.jpg]

-v-

## Mass regimes 2) Intermediate planets

- Super-Earths and sub-Neptunes ($1.5<R<5R_\oplus$ & $3<M<30M_\oplus$) are the most complex planets
- They may host all four materials (e.g. water and volatiles)
- With only mass and radius, we have large degeneracies between solutions:
  - A $4R_\oplus$/$10M_\oplus$ planet could be 95% rock+iron and 5% H-He, or 40% rock+iron and 60% water.
- The treatment of volatiles is particularly complex, as the following locations of water produce vastly different EoS:
  - Water completely [absorbed into a high-pressure silicate mantle](https://www.nature.com/articles/s41550-024-02347-z)
  - Water as an [ocean layer below a H-He atmosphere](https://arxiv.org/pdf/2406.12794)
  - Water fully mixed with H-He gas as a miscible supercritical atmosphere (a **Neptini** - like mixed [vodka & water](https://astrobites.org/2026/03/12/mixing-planetary-cocktails/))

-v-

## Mass regimes 2) Intermediate planets

[i=Water_fraction_GJ9827d.png] 

-v-

## Mass regimes 3) Giant planets

- Can be once again treated simplistically as a mix of H-He and "metals"[+]
- Metallicity can be measured in two ways:[+]
- As a **mass ratio** with respect to solar values, denoted by $Z$[+]
  - $\log{Z}=1$ means 10x solar metallicity - atmospere is 95% H-He and ~5% metals[+]
  - $\log{Z}\lessapprox1.85$ given solar metal mass fraction is $Z\approx0.014$[+]
- As a **log number ratio** with respect to hydrogen (normalised to the Sun):[+]
  - $[\text{M}/\text{H}] \equiv \log_{10}\left(\frac{N_{\text{metals}}}{N_{\text{H}}}\right)_{\text{planet}} - \log_{10}\left(\frac{N_{\text{metals}}}{N_{\text{H}}}\right)_{\odot}$[+]
  - An atmosphere with 33% water (by mass) has $M/H=70x$ or $[M/H]=1.84$[+]
  - Mars' CO$_2$ dominated atmosphere has $[M/H]>4$ (>10000x solar).[+]
- Metallicities can also be normalised to stellar (rather than solar) values.[+]

-v-

## Mass regimes 3) Giant planets

- Metallicity measurements can once again allow some testing of hypothetical link between [stellar and planetary abundances](https://iopscience.iop.org/article/10.3847/0004-637X/831/1/64)[+]
- More interestingly, the metal fraction of giant planets clearly decreases with planet mass - largest planets are most star-like.[+]
  - This is expected from runaway accretion):[+]

[ism=lopez_fortney_mass_Z_giants.jpg] [+]

---
[P "Projects"]

## General advice
- Interior stucture models are underconstrained, i.e. we have a lot more complexity in the model than we have constraints from the data. Things like model comparison aren't so useful
- Interior structure models are typically only applicable for their preferred mass range (usually determined by available EoS)[+]
- Giant planets typically derive only "metallicity"[+]
- Sub-Neptunes are much more complex, and may be interesting to see:[+]
  - Water-rich model (e.g. water mass fraction $10<w_{\rm H2O}<40\%$) either with or without H-He gas[+]
  - Gas-only model[+]
  - Potentially both could have the CMF value fixed to either earth-like, or the expected CMF given the stellar Fe & Mg abundances[+]
- Rocky planets are again simple:[+]
  - Could derive the CMF[+]
  - Could see how much water a planet could contain assuming fixed CMF and no H-He gas[+]

-v-

## Other advice
- Only use transiting planets with precise masses[+]
- Hot Jupiters typically do not produce reliable interior structure estimates (high external heating, unknown additional inflation component), so avoid giant planets with $T_{\rm eq}>700$K.[+]

-v-

## Modelling packages

- All of the following will need to be git installed using (e.g.) `git clone https://github.com/cpiaulet/smint.git` (ideally in the `all_code` directory)
- [**smint** (Piaulet)](https://github.com/cpiaulet/smint) - python-based EoS package for rocky and sub-Neptune interiors[+]
- [**ExoMDN** (Baumeister)](https://github.com/philippbaumeister/ExoMDN) - Requires older python versions (due to tensorflow install)[+]
- [**planetsynth** (Muller)](https://github.com/tiny-hippo/planetsynth) - python-based giant planet interior modelling code[+]
- [**GASTLI** (Acuna)](https://gastli.readthedocs.io) - needs a Fortran compiler[+]
- [**PlaNETic** (Egger)](https://github.com/joannegger/plaNETic) - Neural-network-based interior structure modelling (but heavy/slow; not in script)[+]

-v-

## The coding scripts:

- `rocky_interiors.py`:
  - Looks at the seven TRAPPIST-1 planets (core mass fractions)
- `subnept_interiors.py`:
  - Looks at the planets in the TOI-270d system with a water-rich and gas-rich model
- `giant_interiors.py`:
  - Uses **planetsynth** & **GASTLI** to look at HD 148797 and TOI-4504 c giant planets

-v-

## Project ideas:

- Use the exoplanet archive to pick two or three planets with:
  - A radius $2.0\pm0.1 R_\oplus$, but widely varying masses...
  - A mass $18\pm2 M_\oplus$, but varying radii... 
- What are their possible interiors?
- What is the CMF of the super-Mercuries (i.e. rocky planets with densities above 7g/cm$^3$)?[+]
- What mass of "hot Neptune desert" planets (i.e. $2<R<6R_\oplus$, $P<2$d) could be in H$_2$ and/or volatiles?[+]
- Are any planets consistent with 50% rocky interior and 50% water steam? (i.e. no H-He; some steam-worlds are predicted from theory)[+]
- Do rocky, or sub-Neptune planets with masses derived from TTVs (e.g. in multi-planet systems) have different interior structures?[+]
- What are the possible interior structures of sub-Neptunes K2-18b or GJ1214b (and is one more consistent with water)?[+]

---

[P "Extra Resources"]
- [Interior and Evolution of the Giant Planets - Yamila Miguel](https://arxiv.org/html/2302.07656v1)
- [Fundamentals of interior modelling and challenges in the interpretation of observed rocky exoplanets](https://arxiv.org/pdf/2511.10269)
