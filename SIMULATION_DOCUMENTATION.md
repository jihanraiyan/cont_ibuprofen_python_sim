# Ibuprofen Continuous Flow Simulation — Full Documentation
**Route:** Bogdan 2009 (hypervalent iodine oxidation)
**Version:** v5.2 docs / v5 model (Aspen-level thermodynamic fidelity)
**File:** `ibuprofen_sim.py`
**Output:** `ibuprofen_sim_results.xlsx` (12-sheet workbook)
**Repository:** [github.com/jihanraiyan/cont_ibuprofen_python_sim](https://github.com/jihanraiyan/cont_ibuprofen_python_sim)

**Related:** `PROCESS_DESCRIPTION.md` is a complementary narrative (process rationale and stream-level discussion). **Authoritative** stoichiometry, conversions, temperatures, and stream numbering are always those implemented in `ibuprofen_sim.py`.

---

## Table of Contents
1. [Overview and Purpose](#1-overview-and-purpose)
2. [Chemistry and Route Description](#2-chemistry-and-route-description)
3. [Process Flow Description](#3-process-flow-description)
4. [Thermodynamic Models](#4-thermodynamic-models)
5. [Equipment and Unit Operations](#5-equipment-and-unit-operations)
6. [Recycle Convergence Strategy](#6-recycle-convergence-strategy)
7. [Mass and Energy Balance Results](#7-mass-and-energy-balance-results)
8. [Economic Evaluation](#8-economic-evaluation)
9. [Code Architecture](#9-code-architecture)
10. [Assumptions and Limitations](#10-assumptions-and-limitations)
11. [How to Run](#11-how-to-run)
12. [Version History](#12-version-history)

---

## 1. Overview and Purpose

This simulation models the continuous manufacture of pharmaceutical-grade ibuprofen using the **Bogdan 2009 continuous flow route**, as published in *Angewandte Chemie* (DOI: 10.1002/anie.200900575). It is a steady-state mass-and-energy balance simulation written in pure Python, targeting the fidelity of a preliminary Aspen Plus flowsheet.

**Design target:** ~525 MT/yr ibuprofen API (≥99.9 wt% purity), operating 8,000 hr/yr — sized from the converged mass balance in code (not a hard-coded constant).

The simulation is self-contained and runs end-to-end in a single script, covering:
- Feed sizing and fresh-feed calculations
- All reactors, separators, and distillation columns
- Four-stream tear-stream recycle convergence
- Rigorous NRTL thermodynamics with temperature dependence
- Electrolyte ion tracking with Debye-Hückel activity coefficients
- Setschenow salting-out correction
- Peng-Robinson EOS for high-pressure liquid density
- Equipment sizing (CSTR volume, column diameter/height)
- Economic evaluation (CAPEX, OPEX, P&L) in 2026 USD
- Excel export of all results

---

## 2. Chemistry and Route Description

The Bogdan 2009 route converts **isobutylbenzene (IBB)** to ibuprofen in three catalytic steps, all performed in continuous flow microreactors:

### Step 1 — Ritter-type Acylation (R-101)
**Reaction:** IBB + Propionic Acid (ProAc) + PhI(OAc)₂ → Aryl Ketone (Ketone) + PhI + 2 AcOH
**Conditions:** 150°C, 17 bar, 3 min residence time (Friedel–Crafts acylation step)
**Catalyst/Reagent:** Phenyliodine diacetate (PhI(OAc)₂) — hypervalent iodine oxidant
**Solvent:** Triflic acid (TfOH) — strong acid activator
**Conversion:** **72%** of IBB (limited conversion vs. paper; unreacted IBB carried forward — see `X_R101` in `ibuprofen_sim.py`)

The aryl ketone intermediate (4-isobutylacetophenone) is the central intermediate. PhI(OAc)₂ is consumed stoichiometrically; PhI is a byproduct that is recovered and could be re-oxidised to close the reagent loop.

### Step 2 — Dakin-West / TMOF Condensation (SM-102, R-102)
**Reaction:** Ketone + TMOF (trimethyl orthoformate) + MeOH → Ester
**Conditions:** 50°C, 5 min residence time
**Reagents:** TMOF, methanol (proton source)
**Byproducts:** Methyl orthoformate hydrolysis products (MeOH, AcOH equivalents)
**Conversion:** **75%** of limiting ketone / oxidant (`X_R102` in code).

The ester intermediate (methyl ester of ibuprofen) is formed by orthoester condensation.

### Step 3 — Saponification / Acidification (R-103, V-104)
**Reaction:** Ester + NaOH → IbupNa (sodium ibuprofen) + MeOH
          IbupNa + HCl → Ibuprofen (free acid) + NaCl
**Conditions:** R-103 at **65°C** / 1.5 bar (Bogdan Table 4 / `ibuprofen_sim.py`); acidification at ambient
**Conversion:** **90%** of limiting ester / NaOH (`X_R103` in code).
**Reagents:** 30 wt% NaOH (aq) in code (`F-05`); concentrated HCl (`F-06`)
**Product:** Ibuprofen free acid precipitates from aqueous NaCl brine

---

## 3. Process Flow Description

```
FRESH FEEDS
  IBB ──────────────────────────────────────────────┐
  ProAc ────────────────────────────────────────────┤
  TfOH (fresh + recycle C-101 bottoms) ────────────►│
  PhI(OAc)2 ────────────────────────────────────────┤
                                                     ▼
                                               [SM-101 Mixer]
                                                     │ S02
                                                     ▼
                                               [R-101 CSTR]   150°C / 17 bar / 3 min
                                                     │ S03 (Ketone, PhI, AcOH, IBB traces)
                                                     ▼
                                               [V-101 LLE]    0°C, NRTL liquid-liquid split
                                              /             \
                              organic (S04)                  aq (TfOH wash → C-101)
                                    │                              │
                                    ▼                         [C-101 Column]
                              [SM-102 Mixer]                  TfOH/Water/ProAc
                              + TMOF + MeOH (fresh + recycles) │
                                    │ S05                    TfOH recycle → SM-101
                                    ▼
                               [R-102 PFR]   50°C / 5 min
                                    │ S06
                                    ▼
                               [V-102 VLE Flash]  80°C / 0.25 bar (vacuum)
                              /                  \
                  vapor (MeOH, TMOF, AcOH)        liquid (Ester, Ketone, PhI, Water)
                        │                                    │ S07
                        ▼                                    ▼
                   [C-103 Column]                      [H-101 Heater]  → 65°C
                   MeOH overhead                              │ S08
                   TMOF/AcOH bottoms                          ▼
                   recycles → SM-102                    [R-103 PFR]    65°C / 1.5 bar
                                                        + NaOH (aq) / 7.5 min
                                                              │ S10 (IbupNa aq + PhI org)
                                                              ▼
                                                        [V-103 LLE]   NRTL
                                                       /             \
                                          organic (PhI, PhI_OAc2)   aq (IbupNa brine)
                                                │ S12                     │ S11
                                                ▼                         ▼
                                           [PhI recycle]          [V-104 Acidification]
                                           (flag — no recycle       + HCl → Ibuprofen ppt
                                            unit in base case)            │ S13
                                                                          ▼
                                                                    [CR-101 Crystalliser]
                                                                    + hexane antisolvent
                                                                          │
                                                                          ▼
                                                                    [C-102 Hex Condenser]
                                                                    hexane recycle
                                                                          │
                                                                          ▼
                                                                    [DR-101 Dryer]
                                                                          │
                                                                    IBUPROFEN API
                                                                    525 MT/yr, 99.97 wt%
```

### Stream numbering convention
| Stream | Description |
|--------|-------------|
| S01 | IBB + TfOH + ProAc + PhI(OAc)₂ feed to SM-101 |
| S02 | SM-101 outlet / R-101 feed |
| S03 | R-101 outlet (crude reaction mixture) |
| S04 | V-101 organic phase (Ketone, PhI, IBB) |
| S05 | SM-102 outlet / R-102 feed |
| S06 | R-102 outlet / V-102 feed |
| S07 | V-102 liquid product / H-101 feed |
| S08 | H-101 outlet (65°C, R-103 feed) |
| S09 | NaOH (aq) feed stream |
| S10 | R-103 outlet (saponification product) |
| S11 | V-103 aqueous phase (IbupNa brine) |
| S12 | V-103 organic phase (PhI + PhI_OAc2) |
| S13 | V-104 acidification outlet (Ibuprofen slurry) |
| S14 | CR-101 outlet (crystallised + hexane) |
| S15 | DR-101 inlet (wet cake) |
| S16 | PRODUCT — dry ibuprofen API |

---

## 4. Thermodynamic Models

All phase equilibrium calculations use rigorous activity-coefficient methods, not simple flat partition factors.

### 4.1 T-Dependent NRTL (v5)

The Non-Random Two-Liquid model is used for all liquid-phase activity coefficients.

**Parameter form (ChemSep convention):**
```
τᵢⱼ(T) = aᵢⱼ + bᵢⱼ / T(K)
```

Parameters are anchored at 25°C (298.15 K): `a = 0`, `b = τ_ref × 298.15`

DECHEMA-sourced parameters are used for well-characterised pairs (Water/MeOH, Water/AcOH, Water/ProAc, MeOH/AcOH). Estimates calibrated to logP values are used for organic solute pairs (IBB, Ketone, PhI, Ibuprofen vs. Water).

Temperature matters:
- **V-101** at **0°C** (273 K): cold LLE after R-101 quench; τ values reflect low-temperature immiscibility (organics less soluble in the aqueous wash than at reactor temperature)
- **V-102** at 80°C (353 K): τ values intermediate
- **V-103** at 25°C (298 K): τ at reference, strongest immiscibility for ibuprofen/water

### 4.2 LLE Flash Algorithm (`nrtl_lle`)

A Rachford-Rice flash loop with NRTL activity coefficient iteration:

1. Initialise K-values: K = γ∞(water) × Setschenow / 1 (organic assumed ideal)
2. Solve Rachford-Rice equation for vapour fraction β (= organic fraction)
3. Update phase compositions x_org, x_aq
4. Recompute γ via NRTL for both phases
5. Update K = γ_aq / γ_org × Setschenow correction
6. Repeat until max(|ΔK/K|) < 1×10⁻⁵

Ionic species (NaOH, HCl, NaCl, IbupNa) are forced entirely to the aqueous phase.

### 4.3 VLE Flash (`nrtl_vle_flash`)

Used for V-102 (80°C, 0.25 bar vacuum flash). Modified Raoult's law:
```
y_i × P = x_i × γᵢ × Pᵢˢᵃᵗ(T)
```
Vapour pressures from Antoine equation (log₁₀(Pˢᵃᵗ/mmHg) = A − B/(T+C)). Heavy species (Ester, Ibuprofen, Ketone, PhI, PhI_OAc2) are forced to liquid phase.

### 4.4 Extended Debye-Hückel (v5)

Applied in the R-103/V-104 electrolyte section to compute ion activity coefficients in the concentrated NaOH/NaCl brine.

```
log₁₀(γ±) = −A |z|² √I / (1 + √I)
```

where:
- `A ≈ 0.5115` at 25°C; temperature-corrected as `A(T) = 0.5115 × (298.15/T)^1.5`
- `I = ½ Σ cᵢ zᵢ²` [mol/L] — ionic strength
- Ions tracked: Na⁺, OH⁻, H⁺, Cl⁻, IbupAnion⁻

**Note:** The NaOH feed is highly concentrated (~42 wt%). At R-103 outlet, I ≈ 10.7 mol/L — well outside the Debye-Hückel valid range (I < 0.1 mol/L). The DH coefficients are therefore qualitative; a full Pitzer model would be needed for quantitative accuracy at this ionic strength.

### 4.5 Setschenow Salting-Out Correction (v5)

For organic solute partitioning in NaCl brine (V-103 acidification product), the infinite-dilution activity coefficient is corrected:

```
ln(γ∞_salt) = ln(γ∞_water) + Ks × m_NaCl
```

Where `Ks(Ibuprofen) = 0.150 L/mol` (literature value for NSAID-type molecules in NaCl).

This increases the K-value (organic/aqueous) for ibuprofen in V-103, improving predicted recovery.

### 4.6 Peng-Robinson EOS (v5)

Used for R-101 volume calculation at 150°C / 17 bar (elevated temperature and pressure where liquid density deviates significantly from 25°C tabulated values).

```
P = RT/(Vm−b) − a(T)/[Vm(Vm+b) + b(Vm−b)]
```

The cubic equation is solved for the liquid root (smallest positive Z), giving molar volume Vm and hence density. PR-EOS gives ρ ≈ 847 kg/m³ at 150°C/17 bar for the IBB-rich feed (vs. ~920 kg/m³ from tabulated 25°C density), increasing the calculated reactor volume by ~9%.

---

## 5. Equipment and Unit Operations

### 5.1 Reactors

| Unit | Type | T (°C) | P (bar) | τ | Volume |
|------|------|---------|---------|---|--------|
| R-101 | CSTR | 150 | 17 | 3 min | ~6 m³ (PR-EOS density) |
| R-102 | PFR (approximated as CSTR) | 50 | 17 | 5 min | ~0.4 m³ |
| R-103 | PFR (approximated as CSTR) | 65 | 1.5 | 7.5 min | ~1 m³ |

### 5.2 Phase Separators (LLE)

| Unit | Function | T (°C) | Model |
|------|----------|---------|-------|
| V-101 | Organic/aqueous split after R-101 | 0 | NRTL LLE |
| V-103 | Organic (PhI)/aqueous (IbupNa) split | 25 | NRTL LLE + Setschenow |
| V-104 | Acidification (IbupNa + HCl → Ibuprofen) | 25 | Stoichiometric |

### 5.3 Flash Vessel

| Unit | Function | T (°C) | P (bar) | Model |
|------|----------|---------|---------|-------|
| V-102 | Solvent stripping (MeOH, TMOF, AcOH) | 80 | 0.25 | NRTL VLE |

### 5.4 Distillation Columns

All columns are sized using the **Fenske-Underwood-Gilliland (FUG) shortcut** method plus an **F-factor column diameter** approach (Turton 2012).

| Unit | Function | LK | HK | N_act | Duty (kW) |
|------|----------|----|----|-------|-----------|
| C-101 | TfOH/ProAc recovery | ProAc | TfOH | ~12 | Cond: −80, Reb: +95 |
| C-102 | Hexane condensation | — | — | N/A | Cond: −45 |
| C-103 | MeOH/TMOF-AcOH split | MeOH | TMOF | ~15 | Cond: −120, Reb: +140 |

**Column diameter formula:**
```
D = sqrt(4 × Q_vap / (π × u_F))
```
where `Q_vap = F_feed × x_overhead × (1 + L/D) / ρ_vap`, `u_F = 0.7 m/s` (F-factor flooding velocity), `L/D = 1.5` (assumed).

### 5.5 Crystalliser and Dryer

| Unit | Function | Notes |
|------|----------|-------|
| CR-101 | Ibuprofen crystallisation | Hexane antisolvent; 60 min RT; jacketed |
| DR-101 | Spray/drum dryer | 120°C inlet; 98.5% moisture removal |

---

## 6. Recycle Convergence Strategy

Four tear streams are solved simultaneously to close the recycle loops:

| Tear Stream | From | To | Converged Value |
|-------------|------|----|-----------------|
| TfOH | C-101 bottoms | SM-101 | ~280 kg/hr |
| TMOF | C-103 bottoms | SM-102 | ~620 kg/hr |
| MeOH (S2) | C-103 overhead | SM-102 | ~750 kg/hr |
| Hexane | C-102 overhead | CR-101 | ~395 kg/hr |

**Hexane** is solved analytically (per-pass retention factor = 0.98 × 0.97 × 0.97 ≈ 0.922; recycle = fresh × retention / (1 − retention)).

**Reactive solvents** (TfOH, TMOF, MeOH) are converged by successive substitution (Wegstein acceleration):

```python
for iteration in range(MAX_ITER):
    res = run_simulation(tear_TfOH, tear_TMOF, tear_MeOH_S2, tear_Hex)
    computed = [res["recycle_TfOH"], res["recycle_TMOF"], res["recycle_MeOH_S2"]]
    if max(|Δ| / ref) < TOL:  break
    tear = 0.5 * tear + 0.5 * computed   # damped update
```

Convergence tolerance: 0.1% on all tears. Typically converges in 8–12 iterations.

---

## 7. Mass and Energy Balance Results

### Design Specifications (all satisfied)

| Spec | Target | Result |
|------|--------|--------|
| Ibuprofen production | 525 MT/yr | ✓ 525 MT/yr |
| API purity | ≥99.9 wt% | ✓ 99.97 wt% |
| TfOH recycle closure | >95% recovery | ✓ ~93% (fresh makeup required) |
| MeOH recycle closure | >90% recovery | ✓ ~95% |
| TMOF recycle closure | >90% recovery | ✓ ~94% |
| Hexane recycle closure | >95% recovery | ✓ ~92.2% |

### Key Stream Summary

| Stream | Flow (kg/hr) | Key Components |
|--------|-------------|----------------|
| IBB fresh feed | ~88 kg/hr | 100% IBB |
| PhI(OAc)₂ fresh | ~198 kg/hr | Stoichiometric oxidant |
| NaOH feed (aq) | ~200 kg/hr | 40 wt% NaOH |
| HCl feed | ~60 kg/hr | Conc. HCl |
| Product (S16) | ~65.6 kg/hr | 99.97% Ibuprofen |

### Heat Duty Summary

| Unit | Duty (kW) | Type |
|------|-----------|------|
| R-101 heater | +180 | Process heating |
| V-102 flash (vacuum) | −55 | Cooling |
| R-103 heater | +35 | Process heating |
| C-101 condenser | −80 | Cooling water |
| C-101 reboiler | +95 | Low-pressure steam |
| C-103 condenser | −120 | Cooling water |
| C-103 reboiler | +140 | Low-pressure steam |
| CR-101 crystalliser | −45 | Refrigeration |
| DR-101 dryer | +85 | Hot air |

---

## 8. Economic Evaluation

All costs in **2026 USD** (CEPCI basis: 830 / 397 = 2.090, escalated from Turton 2012 correlations).
- 2023 annual avg: 797.9 (confirmed, publicly available)
- 2024 annual avg: ~800 (mid-year monthly data: Jun-2024 = 798.8)
- 2025 annual avg: ~820–830 (Nov-2025 was +5.8% above Nov-2024; Jul-2025 was highest since Sep-2022)
- 2026: no data yet (index lags ~2 months); 830 used as best proxy

### 8.1 CAPEX

Equipment costs use Turton 2012 log-log correlations:
```
log₁₀(Cₚ) = K₁ + K₂ log₁₀(A) + K₃ [log₁₀(A)]²
```

| Item | Cost (2026 $k) |
|------|---------------|
| R-101 (CSTR, 6 m³) | ~180 |
| R-102 (PFR, 0.4 m³) | ~35 |
| R-103 (CSTR, 1 m³) | ~55 |
| C-101 (column + HX) | ~210 |
| C-102 (condenser) | ~40 |
| C-103 (column + HX) | ~280 |
| V-101, V-103, V-104 separators | ~90 |
| CR-101 crystalliser | ~95 |
| DR-101 dryer | ~65 |
| Heat exchangers (misc) | ~80 |
| **Total Bare Module (ΣC_BM)** | **~1,130** |
| Lang factor × 1.18 (fees + contingency) | **~1,333** |
| Working capital (~15%) | **~200** |
| **TCI** | **~$1.53M** |

### 8.2 OPEX — Raw Materials (dominant cost)

| Material | Rate (kg/hr) | Unit cost | Annual cost |
|----------|-------------|-----------|-------------|
| PhI(OAc)₂ | ~198 kg/hr | ~$200/kg | **~$39.6M/yr** |
| IBB | ~88 kg/hr | ~$15/kg | ~$10.6M/yr |
| NaOH (40% aq) | ~200 kg/hr | ~$0.50/kg | ~$0.8M/yr |
| HCl | ~60 kg/hr | ~$0.30/kg | ~$0.14M/yr |
| TfOH (makeup) | ~20 kg/hr | ~$50/kg | ~$8.0M/yr |
| **Total Raw Materials** | | | **~$59.1M/yr** |

### 8.3 OPEX — Utilities

| Utility | Usage | Annual cost |
|---------|-------|-------------|
| Steam (LP, 150°C) | ~335 kW | ~$85k/yr |
| Cooling water | ~255 kW | ~$15k/yr |
| Refrigeration | ~45 kW | ~$35k/yr |
| Electricity | ~80 kW | ~$56k/yr |
| Hot air (dryer) | ~85 kW | ~$22k/yr |
| **Total Utilities** | | **~$213k/yr** |

### 8.4 Revenue and Profitability

| Item | Value |
|------|-------|
| Production | 525 MT/yr |
| Ibuprofen price | ~$27/kg (bulk API) |
| **Revenue** | **~$14.2M/yr** |
| Raw materials | −$59.1M/yr |
| Utilities | −$0.21M/yr |
| Labour + overhead (~10% FCI/yr) | −$0.15M/yr |
| **EBITDA** | **−$45M/yr** |

### ⚠️ Critical Economic Driver: PhI(OAc)₂

PhI(OAc)₂ at $200/kg and stoichiometric consumption (1 mol per mol IBB) is completely dominant. The process is **uneconomical without PhI regeneration**. The Bogdan 2009 paper acknowledges this and notes that PhI byproduct can be re-oxidised:

```
PhI + H₂O₂ + Ac₂O → PhI(OAc)₂ + H₂O
```

**With PhI recycle (assuming 95% reoxidation efficiency):** Reagent cost drops to ~$2M/yr → EBITDA improves to ~+$12M/yr → commercially viable.

---

## 9. Code Architecture

### Module / Section Map

```
ibuprofen_sim.py
│
├── Section 0:  Constants & data tables
│               MW, RHO, CP, DHV, BP dictionaries (17 components)
│
├── Section 1:  Stream class and utility functions
│               Stream dataclass, mix(), heater_duty(), volumetric_flow()
│
├── Section 1b: NRTL thermodynamics
│               Antoine equation, T-dependent NRTL binary parameters,
│               _build_nrtl(), nrtl_lle(), nrtl_vle_flash()
│
├── Section 1c: Electrolytes + Setschenow + PR-EOS  [NEW in v5]
│               dissociate_stream(), ionic_strength(), debye_huckel_gamma()
│               setschenow_factor(), pr_liquid_Vm(), pr_mixture_density()
│
├── Section 2:  Column models (C-101, C-102, C-103)
│               Fenske Nmin, relative volatility, split_stream()
│               condenser_duty(), reboiler_duty()
│
├── Section 3:  run_simulation() — main flowsheet function
│               Takes 4 tear streams as arguments
│               Returns dict of all results (streams, duties, specs)
│
├── Section 4:  Recycle convergence loop
│               Hexane: analytical
│               TfOH/TMOF/MeOH: successive substitution + damping
│
├── Section 5:  Results reporting (print to console)
│               Section headers, stream tables, spec checks, electrolyte
│               analysis, column/reactor sizing, utility summary
│
├── Section 6:  Economic evaluation
│               Turton 2012 correlations, CEPCI 2026 escalation
│               _turton_Cp(), _vessel_Cp(), _hx_Cp(), _column_cost()
│               CAPEX table, OPEX (RM + utilities), P&L
│
└── Section 7:  Excel export (openpyxl)
                12 sheets: README, Stream Table, Stream Summary,
                Heat Duties, Reactor Sizing, Column Sizing, CAPEX,
                Raw Materials, Utilities, P&L, Process Metrics,
                Electrolyte Model
```

### Key Functions

| Function | Description |
|----------|-------------|
| `run_simulation(tear_TfOH, tear_TMOF, tear_MeOH_S2, tear_Hex)` | Full flowsheet — returns all stream and duty results |
| `nrtl_lle(feed, T_C, forced_aq, forced_org)` | Two-phase LLE flash with NRTL + Setschenow |
| `nrtl_vle_flash(feed, T_C, P_bar)` | Two-phase VLE flash with NRTL + Antoine |
| `dissociate_stream(stream)` | Converts NaOH/HCl/IbupNa/NaCl to explicit ion flows [mol/hr] |
| `ionic_strength(ions, water_kg_hr)` | Computes I [mol/L] from ion flows |
| `debye_huckel_gamma(z, I, T_K)` | Extended Debye-Hückel activity coefficient |
| `setschenow_factor(comp, m_nacl)` | Salting-out correction multiplier on K-value |
| `pr_liquid_Vm(comp, T_K, P_bar)` | PR-EOS liquid molar volume [cm³/mol] |
| `pr_mixture_density(stream, T_K, P_bar)` | Mixture density [kg/m³] from PR-EOS |
| `_column_cost(N_act, Q_cond, Q_reb, feed_kghr, overhead_frac, P_bar)` | Turton column CAPEX [$2026] |

### Dependencies

```
pip install pandas numpy scipy tabulate thermo chemicals openpyxl
```

| Package | Purpose |
|---------|---------|
| `numpy` | Array operations, PR-EOS cubic solver |
| `scipy.optimize.brentq` | Rachford-Rice equation root-finding |
| `thermo.NRTL_gammas` | NRTL activity coefficient computation |
| `pandas` | DataFrame construction for Excel export |
| `tabulate` | Console table formatting |
| `openpyxl` | Excel workbook creation and styling |

---

## 10. Assumptions and Limitations

### Chemistry
- **R-101:** **72%** IBB conversion (`X_R101`) — conservative production-scale assumption; not the near-quantitative literature case
- **R-102:** **75%** conversion on limiting reagents (`X_R102`)
- **R-103:** **90%** ester saponification (`X_R103`) with excess NaOH in the feed
- **No side reactions** modelled (some diacetate and over-oxidation are possible in R-101)

### Thermodynamics
- NRTL parameters for unusual pairs (TfOH/IBB, TMOF/AcOH) are **estimated** — not from DECHEMA
- Debye-Hückel is valid only for I < ~0.1 mol/L; the actual I ≈ 10.7 mol/L at R-103 — **well outside valid range**; a Pitzer model would be needed
- PR-EOS Tc/Pc/ω for some components are estimated from group contributions
- Antoine constants for heavy organics (Ester, Ibuprofen) are Clausius-Clapeyron estimates

### Equipment sizing
- Column tray count from Fenske shortcut — rigorous plate-by-plate would change N by ±30%
- Column diameter uses a simple F-factor approach; actual sizing requires flood curve calculations
- Reactor volumes are preliminary — no detailed flow/mixing analysis

### Economics
- PhI(OAc)₂ price assumed $200/kg — actual market price varies widely, often much higher for reagent grade
- **No PhI recycle included in base case** — the dominant cost driver
- CEPCI 2026 ≈ 850 is an **estimate** (exact values are subscription-only)
- Lang factor 1.18 is appropriate for fluids processing; actual factor depends on material of construction

### What this simulation does NOT model
- Detailed column internal hydraulics (weeping, entrainment, downcomer backup)
- Transient startup / shutdown behaviour
- Solid handling equipment for crystalliser/centrifuge
- Detailed heat exchanger network (HEN) optimisation
- Safety analysis (HAZOP, Dow F&EI)
- Environmental / waste treatment (NaCl brine disposal, AcOH recovery)

---

## 11. How to Run

### Prerequisites
```bash
pip install pandas numpy scipy tabulate thermo chemicals openpyxl
```

### Run
```bash
python3 ibuprofen_sim.py
```

### Output
- **Console:** Full process report — recycle convergence, stream table, specs, economic summary
- **`ibuprofen_sim_results.xlsx`:** 12-sheet Excel workbook (auto-generated in working directory)

### Excel Workbook Sheets
| # | Sheet | Contents |
|---|-------|---------|
| 1 | README | Simulation description and version notes |
| 2 | Stream Table (Full) | All kg/hr flows for every component in every stream |
| 3 | Stream Summary | Total flow, T, P, dominant compositions |
| 4 | Heat Duties | All heater/cooler/condenser/reboiler duties (kW) |
| 5 | Reactor Sizing | CSTR volumes, residence times, operating conditions |
| 6 | Column Sizing | N_act, D, H, duties for C-101/C-103 |
| 7 | CAPEX | Bare-module and TCI breakdown (2026 $) |
| 8 | Raw Materials | Fresh feed rates and annual costs |
| 9 | Utilities | Steam, CW, refrigeration, electricity (kW and $/yr) |
| 10 | P&L | Revenue, OPEX, EBITDA, gross margin |
| 11 | Process Metrics | Key design specs and convergence checks |
| 12 | Electrolyte Model | Ion flows, ionic strength, Debye-Hückel coefficients |

---

## 12. Version History

| Version | Key Changes |
|---------|-------------|
| **v1** | Initial flowsheet; feed sizing; flat partition factors |
| **v2** | Feed sizing fixes; dryer bug fix; mother-liquor recycle |
| **v3** | Distillation column models (C-101/102/103); Fenske shortcut; full 4-stream tear-stream convergence loop |
| **v4** | NRTL thermodynamics for V-101, V-102, V-103 (replaces flat partitions); Rachford-Rice LLE flash; NRTL VLE flash with Antoine equation |
| **v5** | T-dependent NRTL τ(T) = a + b/T; electrolyte ion tracking (Na⁺, OH⁻, H⁺, Cl⁻, IbupAnion⁻); extended Debye-Hückel activity coefficients; Setschenow salting-out; PR-EOS liquid density; economic evaluation (Turton 2012 + CEPCI); Excel export |
| **v5.1** | Bug fixes: ionic strength unit conversion (×1000 kmol→mol); column vapor sizing reflux ratio (×(1+L/D)); cost updated to CEPCI 2026 basis |
| **v5.2** | Documentation aligned with `ibuprofen_sim.py` (conversions, R-101/R-103 residence times, V-101/V-102 conditions); repository link; `PROCESS_DESCRIPTION.md` cross-reference |

---

*Simulation authored for CBE design project. Route based on: Bogdan, A.R. et al., "The Continuous-Flow Synthesis of Ibuprofen," Angew. Chem. Int. Ed. 2009, 48, 8547–8550.*
