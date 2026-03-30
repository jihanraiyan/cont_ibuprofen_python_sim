# Continuous Flow Ibuprofen Synthesis — Python Process Simulator

**Route:** Bogdan 2009 (hypervalent iodine oxidation)  
**Model version:** v5 (Aspen-level thermodynamic fidelity)  
**Script:** `ibuprofen_sim.py`  
**Output:** `ibuprofen_sim_results.xlsx` (12-sheet workbook)

> **Related:** `PROCESS_DESCRIPTION.md` covers process rationale and stream-level narrative. Authoritative stoichiometry, conversions, temperatures, and stream numbering are always those implemented in `ibuprofen_sim.py`.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Chemistry](#2-chemistry)
3. [Process Flow](#3-process-flow)
4. [Thermodynamic Models](#4-thermodynamic-models)
5. [Equipment and Unit Operations](#5-equipment-and-unit-operations)
6. [Recycle Convergence](#6-recycle-convergence)
7. [Results](#7-results)
8. [Economics](#8-economics)
9. [Code Architecture](#9-code-architecture)
10. [Assumptions and Limitations](#10-assumptions-and-limitations)
11. [How to Run](#11-how-to-run)
12. [Version History](#12-version-history)

---

## 1. Overview

This simulation models the continuous manufacture of pharmaceutical-grade ibuprofen using the **Bogdan 2009 continuous flow route** (*Angew. Chem. Int. Ed.* 2009, DOI: 10.1002/anie.200900575). It is a steady-state mass-and-energy balance written in pure Python, targeting the fidelity of a preliminary Aspen Plus flowsheet.

**Design target:** ~525 MT/yr ibuprofen API (≥99.9 wt% purity), operating 8,000 hr/yr.

The simulation covers end-to-end:

- Feed sizing and fresh-feed calculations
- All reactors, separators, and distillation columns
- Four-stream tear-stream recycle convergence
- Rigorous T-dependent NRTL thermodynamics (LLE and VLE)
- Electrolyte ion tracking with Debye-Hückel activity coefficients
- Setschenow salting-out correction
- Peng-Robinson EOS for high-pressure liquid density
- Equipment sizing (reactor volumes, column diameter/height)
- Economic evaluation (CAPEX, OPEX, P&L) in 2026 USD
- Excel export of all results (12 sheets)

---

## 2. Chemistry

The Bogdan 2009 route converts **isobutylbenzene (IBB)** to ibuprofen in three catalytic steps, all in continuous flow:

### Step 1 — Friedel-Crafts Acylation (R-101)

```
IBB + ProAc  →  Aryl Ketone + H₂O
           TfOH (cat.), 150°C, 17 bar, 3 min
```

- **Conversion:** 72% of IBB (`X_R101`)
- **Catalyst:** Triflic acid (TfOH) — strong acid activator, fed neat
- By-product water must be removed before Step 2 (hydrolyses PhI(OAc)₂)

### Step 2 — 1,2-Aryl Migration (R-102)

```
Ketone + PhI(OAc)₂ + TMOF + H₂O  →  Ester + PhI + 2 AcOH + MeOH + HCOOCH₃
                                      50°C, 17 bar, 5 min
```

- **Conversion:** 75% of limiting ketone (`X_R102`)
- **TMOF** is a stoichiometric co-reagent (methoxy donor), not just a solvent
- **4-stage PhI(OAc)₂ injection** — bulk concentration (0.44 mol/L) exceeds solubility limit (0.15 mol/L); staged addition keeps per-stage value at ~0.11 mol/L

### Step 3 — Saponification and Acidification (R-103, V-104)

```
Ester + NaOH  →  IbupNa + MeOH         65°C, 1.5 bar, 7.5 min
IbupNa + HCl  →  Ibuprofen + NaCl      ambient
```

- **Conversion:** 90% of limiting ester (`X_R103`)
- NaOH used instead of KOH (3–4× cheaper; identical chemistry)

---

## 3. Process Flow

```
FRESH FEEDS
  IBB ──────────────────────────────────────────────────┐
  ProAc ────────────────────────────────────────────────┤
  TfOH (fresh + recycle from C-101) ───────────────────►│
                                                         ▼
                                                   [SM-101 Mixer]
                                                         │ S-01
                                                         ▼
                                                   [HX-101]  25→150°C
                                                         │ S-02
                                                         ▼
                                                   [R-101 PFR]  150°C / 17 bar / 3 min
                                                         │ S-03
                                                         ▼
                                                   [HX-102]  150→0°C
                                                         │ S-04
                                                         ▼
                                                   [V-101 LLE]  0°C, NRTL
                                              ┌──────────┴──────────┐
                                    organic (S-05)              aqueous (S-06)
                                    Ketone, IBB, TfOH            Water, ProAc
                                          │                           │
                                          │                     [C-101 Column]
                                          │                     TfOH recycle ──► SM-101
                                          │
                                          ▼
                                   [SM-102 Mixer]
                                   + PhI(OAc)₂, TMOF, MeOH
                                   + TMOF/MeOH recycles
                                          │ S-07
                                          ▼
                                   [HX-103]  →50°C
                                          │ S-08
                                          ▼
                                   [R-102 PFR]  50°C / 17 bar / 5 min
                                   (4-stage PhI(OAc)₂ injection)
                                          │ S-09
                                          ▼
                                   [V-102 VLE Flash]  80°C / 0.25 bar
                              ┌──────────┴──────────────┐
                         vapor (S-10)               liquid (S-11)
                         MeOH, TMOF, AcOH            Ester, PhI, Ketone
                              │                           │
                        [C-103 Column]             [SM-103 Mixer]
                        MeOH overhead              + NaOH (aq) F-05
                        TMOF bottoms                     │ S-12
                        recycles → SM-102                ▼
                                                   [R-103 PFR]  65°C / 1.5 bar / 7.5 min
                                                         │ S-13
                                                         ▼
                                                   [HX-104]  65→30°C
                                                         │ S-14
                                                         ▼
                                                   [V-103 LLE]  25°C, NRTL + Setschenow
                                              ┌──────────┴──────────┐
                                    organic (S-15)              aqueous (S-16)
                                    PhI, Ketone, IBB              IbupNa brine
                                    → PhI recovery                    │
                                                                       ▼
                                                               [V-104 Acidification]
                                                               + HCl (F-06)
                                                                       │ S-17
                                                                       ▼
                                                               [LLE Hexane Extraction]
                                                               4 stages, K_D = 8
                                                               + Hexane recycle (C-102)
                                                                       │ S-19
                                                                       ▼
                                                               [Hex-Flash]  70°C / 0.9 bar
                                                               hexane overhead → C-102
                                                                       │
                                                               [CR-101 Crystalliser]
                                                               MeOH solvent switch, 0°C
                                                                       │
                                                               [CF-101 Filter]  87% recovery
                                                                       │
                                                               [DR-101 Dryer]  55°C / 50 mbar
                                                                       │
                                                            IBUPROFEN API (S-27)
                                                            ~525 MT/yr, ≥99.9 wt%
```

### Stream Numbering

| Stream | Description |
|--------|-------------|
| S-01 | SM-101 combined feed |
| S-02 | HX-101 outlet / R-101 feed (150°C) |
| S-03 | R-101 outlet (Ketone, TfOH, H₂O) |
| S-04 | HX-102 outlet (0°C, to V-101) |
| S-05 | V-101 organic phase (Ketone, IBB, TfOH) |
| S-06 | V-101 aqueous phase (Water, ProAc → C-101) |
| S-07 | SM-102 combined feed |
| S-08 | HX-103 outlet / R-102 feed (50°C) |
| S-09 | R-102 outlet |
| S-10 | V-102 vapor overhead (MeOH, TMOF, AcOH → C-103) |
| S-11 | V-102 liquid bottoms (Ester, PhI) |
| S-12 | SM-103 combined feed |
| S-13 | R-103 outlet (IbupNa aq + PhI org) |
| S-14 | HX-104 outlet (30°C, to V-103) |
| S-15 | V-103 organic phase (PhI, Ketone) |
| S-16 | V-103 aqueous phase (IbupNa brine → V-104) |
| S-17 | V-104 acidification outlet (Ibuprofen free acid) |
| S-19 | LLE organic extract (Ibuprofen + Hexane) |
| S-27 | Final product — dry ibuprofen API |

---

## 4. Thermodynamic Models

All phase equilibrium calculations use rigorous activity-coefficient methods, not flat partition factors.

### 4.1 T-Dependent NRTL

The Non-Random Two-Liquid model is used for all liquid-phase activity coefficients, with temperature-dependent interaction parameters:

```
τᵢⱼ(T) = aᵢⱼ + bᵢⱼ / T(K)
```

Parameters anchored at 25°C: `a = 0`, `b = τ_ref × 298.15`. DECHEMA-sourced for standard pairs (Water/MeOH, Water/AcOH, Water/ProAc, MeOH/AcOH); logP-calibrated for organic solutes vs. water.

Temperature evaluation:
- **V-101** at 0°C — cold LLE; τ reflects low-temperature immiscibility
- **V-102** at 80°C — VLE flash
- **V-103** at 25°C — reference temperature, strongest organic/aqueous immiscibility

### 4.2 LLE Flash Algorithm

Rachford-Rice loop with iterative NRTL update:

1. Initialise K from γ∞(water) × Setschenow factor
2. Solve Rachford-Rice for organic fraction β
3. Update phase compositions x_org, x_aq
4. Recompute γ from NRTL for both phases
5. Update K = γ_aq / γ_org × Setschenow
6. Repeat until max(|ΔK/K|) < 1×10⁻⁵

Ionic species (NaOH, HCl, NaCl, IbupNa) are forced to the aqueous phase.

### 4.3 VLE Flash

Modified Raoult's law for V-102:

```
yᵢ × P = xᵢ × γᵢ × Pᵢˢᵃᵗ(T)
```

Vapour pressures from Antoine equation. Heavy species (Ester, Ibuprofen, Ketone, PhI, PhI(OAc)₂) forced to liquid.

### 4.4 Setschenow Salting-Out

Applied in V-103 to correct partition coefficients for NaCl brine:

```
ln(γ∞_salt) = ln(γ∞_water) + Ks × m_NaCl
```

`Ks(Ibuprofen) = 0.150 L/mol` (Sangster 1997).

### 4.5 Extended Debye-Hückel

Ion activity coefficients at R-103 outlet (NaOH/NaCl brine):

```
log₁₀(γ±) = −A|z|²√I / (1 + √I)
```

`A(T) = 0.5115 × (298.15/T)^1.5`. Ions tracked: Na⁺, OH⁻, H⁺, Cl⁻, IbupAnion⁻.

> **Note:** Valid range is I < 0.1 mol/L. Actual I ≈ 10.7 mol/L at R-103 outlet — results are qualitative. A Pitzer model would be needed for quantitative accuracy.

### 4.6 Peng-Robinson EOS

Used for R-101 volume at 150°C / 17 bar. Cubic solved for liquid root Z; gives ρ ≈ 847 kg/m³ vs. ~920 kg/m³ from tabulated 25°C density — increases calculated reactor volume by ~9%.

---

## 5. Equipment and Unit Operations

### Reactors

| Unit | Type | T (°C) | P (bar) | τ (min) | Volume |
|------|------|---------|---------|---------|--------|
| R-101 | PFR | 150 | 17 | 3 | PR-EOS calculated |
| R-102 | PFR (4-stage injection) | 50 | 17 | 5 | Volumetric Q×τ |
| R-103 | PFR | 65 | 1.5 | 7.5 | Volumetric Q×τ |

### Phase Separators

| Unit | Function | T (°C) | Model |
|------|----------|---------|-------|
| V-101 | Organic/aqueous split after R-101 | 0 | NRTL LLE |
| V-103 | PhI organic / IbupNa aqueous split | 25 | NRTL LLE + Setschenow |
| V-104 | Acidification (IbupNa + HCl → Ibuprofen) | 25 | Stoichiometric |

### Flash Vessel

| Unit | Function | T (°C) | P (bar) | Model |
|------|----------|---------|---------|-------|
| V-102 | Solvent stripping (MeOH, TMOF, AcOH) | 80 | 0.25 | NRTL VLE |

### Distillation Columns

Sized using Fenske-Underwood shortcut + F-factor column diameter (Turton 2012).

| Unit | Function | LK | HK | N actual | Q_cond (kW) | Q_reb (kW) |
|------|----------|----|----|----------|-------------|------------|
| C-101 | TfOH/ProAc recovery | ProAc | TfOH | ~12 | −80 | +95 |
| C-102 | Hexane condensation | — | — | condenser only | −45 | — |
| C-103 | MeOH/TMOF-AcOH split | MeOH | TMOF | ~15 | −120 | +140 |

### Crystallisation Train

| Unit | Function | Conditions |
|------|----------|------------|
| Hex-Flash | Hexane evaporation / solvent switch | 70°C, 0.9 bar, 97% hexane removal |
| CR-101 | MeOH crystallisation | 45°C → 0°C, 0.138 kg/kg solubility at 0°C |
| CF-101 | Rotary vacuum filter | 87% crystal recovery |
| DR-101 | Vacuum dryer | 55°C, 50 mbar — hard limit <75°C (mp 77°C) |

---

## 6. Recycle Convergence

Four tear streams closed by damped successive substitution:

| Tear Stream | Source | Destination | SS Value |
|-------------|--------|-------------|----------|
| TfOH | C-101 bottoms | SM-101 | ~280 kg/hr |
| TMOF | C-103 bottoms | SM-102 | ~620 kg/hr |
| MeOH | C-103 overhead | SM-102 | ~750 kg/hr |
| Hexane | C-102 condenser | LLE extractor | ~395 kg/hr |

**Hexane** solved analytically (per-pass retention = 0.98 × 0.97 × 0.97):

```
tear_Hex_SS = F07_Hex × retention / (1 − retention)
```

**Reactive solvents** (TfOH, TMOF, MeOH) by damped successive substitution:

```python
for iteration in range(MAX_ITER):
    res = run_simulation(tear_TfOH, tear_TMOF, tear_MeOH_S2, tear_Hex)
    err = max(|Δ| / ref for each tear)
    if err < TOL:  break
    tear = 0.35 * tear_old + 0.65 * tear_new   # damping factor 0.65
```

Convergence tolerance: 0.01%. Typically converges in 8–12 iterations.

---

## 7. Results

### Design Specifications

| Spec | Target | Status |
|------|--------|--------|
| Annual production | ≥ 500 MT/yr | ✅ ~525 MT/yr |
| API purity | > 99.0 wt% | ✅ 99.97 wt% |
| Residual MeOH | < 3,000 ppm | ✅ |
| V-101 water in organic | < 0.5 wt% | ✅ |
| V-102 AcOH in bottoms | < 1.0 wt% | ✅ |
| NaOH/Ester ratio | ≥ 1.0× | ✅ |
| HCl coverage | ≥ 95% | ✅ |
| C-101 relative volatility | α > 1 | ✅ |
| C-103 relative volatility | α > 1 | ✅ |

### Key Stream Flows (converged)

| Stream | Flow (kg/hr) | Key Components |
|--------|-------------|----------------|
| F-01 IBB | 134.0 | 100% IBB |
| F-04 PhI(OAc)₂ | 327.0 | Stoichiometric oxidant |
| F-05 NaOH (aq) | 24 NaOH + 56 H₂O | 30 wt% NaOH |
| F-06 HCl (aq) | 18 HCl + 162 H₂O | 10 wt% HCl |
| S-27 Product | ~65.6 | 99.97 wt% Ibuprofen |

### Heat Duties

| Unit | Duty (kW) | Utility |
|------|-----------|---------|
| HX-101 feed heater | +180 | MPS steam |
| HX-102 R-101 cooler | −large | Chilled glycol |
| HX-103 SM-102 heater | +small | LPS steam |
| HX-104 R-103 cooler | −small | Cooling water |
| C-101 condenser | −80 | Cooling water |
| C-101 reboiler | +95 | HPS steam |
| C-103 condenser | −120 | Cooling water |
| C-103 reboiler | +140 | LPS steam |
| CR-101 crystalliser | −45 | Chilled glycol |

---

## 8. Economics

All costs in **2026 USD** (CEPCI 830 / 397 = 2.090, from Turton 2012 correlations).

### CAPEX Summary

| Item | Cost ($k) |
|------|----------|
| R-101, R-102, R-103 reactors | ~270 |
| C-101, C-102, C-103 columns | ~530 |
| V-101, V-102, V-103 separators | ~90 |
| Heat exchangers | ~80 |
| CR-101 crystalliser + filter + dryer | ~230 |
| **Total Bare Module (ΣC_BM)** | **~1,200** |
| FCI (× 1.18 contingency + fees) | ~1,416 |
| Working capital (15% FCI) | ~212 |
| **TCI** | **~$1.63M** |

### OPEX

| Category | Annual Cost |
|----------|-------------|
| PhI(OAc)₂ (dominant) | ~$65.4M/yr |
| IBB | ~$2.1M/yr |
| TfOH (makeup only) | ~$0M/yr (recycled) |
| NaOH, HCl | ~$0.2M/yr |
| Utilities (all) | ~$0.2M/yr |
| Labor + overhead | ~$2.1M/yr |
| **Total OPEX** | **~$70M/yr** |

### Revenue and Profitability

| Item | Value |
|------|-------|
| Ibuprofen API (~525 MT/yr × $20/kg) | ~$84M/yr |
| PhI byproduct | ~$1M/yr |
| **Total Revenue** | **~$85M/yr** |
| **EBITDA (base case)** | **~$15M/yr** |

> ⚠️ **Critical driver:** PhI(OAc)₂ dominates raw material cost. The process economics improve dramatically with PhI regeneration: `PhI + H₂O₂ + Ac₂O → PhI(OAc)₂`. With 95% reoxidation efficiency, reagent cost drops to ~$2M/yr and the process becomes clearly profitable.

---

## 9. Code Architecture

```
ibuprofen_sim.py
│
├── Section 0   Constants & data tables
│               MW, RHO, CP, DHV, BP for 17 components
│
├── Section 1   Stream class and utility functions
│               Stream dataclass, mix(), heater_duty(), volumetric_flow()
│
├── Section 1b  NRTL thermodynamics
│               Antoine equation, T-dependent binary parameters,
│               _build_nrtl(), nrtl_lle(), nrtl_vle_flash()
│
├── Section 1c  Electrolytes + Setschenow + PR-EOS
│               dissociate_stream(), ionic_strength()
│               debye_huckel_gamma(), setschenow_factor()
│               pr_liquid_Vm(), pr_mixture_density()
│
├── Section 2   Column models (C-101, C-102, C-103)
│               Fenske Nmin, relative_volatility(), split_stream()
│               condenser_duty(), reboiler_duty()
│
├── Section 3   run_simulation()  ← main flowsheet function
│               Takes 4 tear streams; returns full results dict
│
├── Section 4   Recycle convergence loop
│               Hexane: analytical steady-state
│               TfOH/TMOF/MeOH: damped successive substitution
│
├── Section 5   Console results reporting
│               Stream table, spec checks, column/reactor sizing
│
├── Section 6   Economic evaluation
│               Turton 2012 + CEPCI 2026; CAPEX, OPEX, P&L
│
└── Section 7   Excel export (openpyxl)
                12 sheets — see below
```

### Key Functions

| Function | Description |
|----------|-------------|
| `run_simulation(tear_TfOH, tear_TMOF, tear_MeOH_S2, tear_Hex)` | Full flowsheet pass; returns all stream and duty results |
| `nrtl_lle(feed, T_C, forced_aq, forced_org)` | Two-phase LLE flash with NRTL + Setschenow |
| `nrtl_vle_flash(feed, T_C, P_bar)` | Two-phase VLE flash with NRTL + Antoine |
| `dissociate_stream(stream)` | NaOH/HCl/IbupNa/NaCl → explicit ion flows (mol/hr) |
| `ionic_strength(ions, water_kg_hr)` | Computes I (mol/L) from ion flows |
| `debye_huckel_gamma(z, I, T_K)` | Extended Debye-Hückel γ± |
| `setschenow_factor(comp, m_nacl)` | Salting-out multiplier on K-value |
| `pr_liquid_Vm(comp, T_K, P_bar)` | PR-EOS liquid molar volume (cm³/mol) |
| `pr_mixture_density(stream, T_K, P_bar)` | Mixture density (kg/m³) from PR-EOS |

### Dependencies

```bash
pip install pandas numpy scipy tabulate thermo chemicals openpyxl
```

| Package | Role |
|---------|------|
| `numpy` | Array operations, PR-EOS cubic solver |
| `scipy.optimize.brentq` | Rachford-Rice root-finding |
| `thermo.NRTL_gammas` | NRTL activity coefficient computation |
| `pandas` | DataFrame construction for Excel export |
| `tabulate` | Console table formatting |
| `openpyxl` | Excel workbook creation and styling |

---

## 10. Assumptions and Limitations

### Chemistry
- Conversions (72% / 75% / 90%) are conservative production-scale assumptions vs. near-quantitative lab values
- No side reactions modelled (diacetate formation, over-oxidation possible in R-101)
- TMOF treated as stoichiometric co-reagent — consumes 1 mol per mol ester

### Thermodynamics
- NRTL parameters for uncommon pairs (TfOH/IBB, TMOF/AcOH) are **estimated**, not from DECHEMA
- Debye-Hückel valid for I < 0.1 mol/L; actual I ≈ 10.7 mol/L at R-103 — qualitative only (Pitzer model needed)
- Antoine constants for Ester and Ibuprofen are Clausius-Clapeyron estimates
- PR-EOS Tc/Pc/ω for some components estimated from group contributions

### Equipment
- Column tray count from Fenske shortcut — rigorous plate-by-plate would change N by ±30%
- Reactor volumes are preliminary — no detailed mixing or flow analysis
- No heat exchanger network optimisation

### Economics
- PhI(OAc)₂ price assumed $25/kg reagent grade (actual market price varies significantly)
- Base case does not include PhI regeneration loop — this is the dominant cost driver
- CEPCI 2026 ≈ 830 is an estimate (subscription required for exact values)

### Not modelled
- Column hydraulics (weeping, flooding, downcomer backup)
- Transient startup/shutdown
- Detailed solid handling (centrifuge, conveying)
- HAZOP / safety analysis
- Environmental treatment (NaCl brine, AcOH disposal)

---

## 11. How to Run

### Install dependencies

```bash
pip install pandas numpy scipy tabulate thermo chemicals openpyxl
```

### Run simulation

```bash
python3 ibuprofen_sim.py
```

### Outputs

| Output | Description |
|--------|-------------|
| Console | Recycle convergence table, stream table, spec checks, economic summary |
| `ibuprofen_sim_results.xlsx` | 12-sheet workbook, auto-generated in working directory |

### Excel workbook sheets

| # | Sheet | Contents |
|---|-------|----------|
| 1 | README | Summary and version notes |
| 2 | Stream Table (Full) | All kg/hr flows for every component in every stream |
| 3 | Stream Summary | Total flow, T, P, dominant compositions (top 5) |
| 4 | Heat Duties | All heater/cooler/condenser/reboiler duties (kW, MW, GJ/yr) |
| 5 | Reactor Sizing | Volumes, residence times, PR-EOS density, imperial units |
| 6 | Column Sizing | N_act, D, H, duties, Turton cost for C-101/C-103 |
| 7 | CAPEX | Bare-module costs, FCI breakdown (Turton 2012, CEPCI 2026) |
| 8 | Raw Materials | Fresh feed rates and annual costs |
| 9 | Utilities | Steam, CW, refrigeration, electricity (kW, GJ/yr, $/yr) |
| 10 | P&L | Revenue, OPEX, EBITDA, NPV, payback, unit production cost |
| 11 | Process Metrics | Key specs, step yields, recycle convergence values |
| 12 | Electrolyte Model | Ion flows, ionic strength, Debye-Hückel coefficients |

---

## 12. Version History

| Version | Key Changes |
|---------|-------------|
| v1 | Initial flowsheet; feed sizing; flat partition factors |
| v2 | Feed sizing fixes; dryer bug fix; mother-liquor recycle |
| v3 | Distillation column models (C-101/102/103); Fenske shortcut; full 4-stream tear-stream convergence |
| v4 | NRTL thermodynamics for V-101, V-102, V-103 (replaces flat partitions); Rachford-Rice LLE/VLE flash |
| v5 | T-dependent NRTL τ(T) = a + b/T; electrolyte ion tracking; Debye-Hückel; Setschenow salting-out; PR-EOS density; Turton 2012 economics; Excel export |
| v5.1 | Bug fixes: ionic strength unit conversion; column vapour sizing reflux ratio; CEPCI updated to 2026 basis |
| v5.2 | Documentation aligned with code (conversions, temperatures, stream numbering); GitHub README rewrite |

---

*Route based on: Bogdan, A.R. et al., "The Continuous-Flow Synthesis of Ibuprofen," Angew. Chem. Int. Ed. 2009, 48, 8547–8550.*  
*Scale-up reference: Jolliffe, H.G. & Gerogiorgis, D.I., Chem. Eng. Res. Des. 2015, 97, 175–191.*
