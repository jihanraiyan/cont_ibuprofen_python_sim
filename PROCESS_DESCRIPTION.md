# Continuous Ibuprofen Synthesis — Process Description

**Route:** Bogdan 2009 three-step continuous flow synthesis
**Scale:** ~640 MT/yr ibuprofen API (134 kg/hr IBB feed basis)
**Reference:** Cross-validated against Jolliffe & Gerogiorgis (2015)

---

## 1. Process Overview

Ibuprofen is synthesised from isobutylbenzene (IBB) in three sequential chemical reactions, each carried out in a continuous plug flow reactor (PFR). The route is:

1. **Friedel-Crafts acylation** — IBB + propionic acid → 4'-isobutylacetophenone (Ketone) + H₂O
2. **Aryl migration (Ritter-type rearrangement)** — Ketone + PhI(OAc)₂ + TMOF + H₂O → methyl ester (Ester) + PhI + 2 AcOH + MeOH + HCOOCH₃
3. **Saponification** — Ester + NaOH → ibuprofen sodium salt (IbupNa) + MeOH

The three reactors operate at different temperatures and pressures. Between and after them, a series of separation units remove by-products, recover solvents, and isolate the final API. Four recycle loops are closed at steady state.

This route was chosen because it avoids the Friedel-Crafts aluminium chloride chemistry used in older batch processes, instead using triflic acid (TfOH) as a homogeneous catalyst and hypervalent iodine (PhI(OAc)₂) as a mild oxidant. The result is a continuous, atom-efficient process amenable to production scale.

---

## 2. Raw Material Feeds

| Stream | Component | Rate (kg/hr) | Justification |
|--------|-----------|-------------|---------------|
| F-01 | IBB | 134.0 | Basis stream — sized for ~640 MT/yr API at 72%/75%/90% cascade conversion |
| F-02 | Propionic acid (ProAc) | 81.4 | 1.1 equiv relative to IBB (slight excess to drive Friedel-Crafts to completion) |
| F-03 | Triflic acid (TfOH, neat) | makeup to 749 kg/hr total | 5.0 equiv relative to IBB; fed neat — paper specifies "neat triflic acid", no water |
| F-04 | PhI(OAc)₂ | 327.0 | ~1.015 equiv relative to IBB (slight excess ensures Ketone is limiting, not oxidant) |
| F-04 | TMOF | makeup to 424 kg/hr total | 4.0 equiv relative to IBB; dissolution solvent and stoichiometric co-reagent |
| F-04 | MeOH | makeup to 950 kg/hr total | ~55.9 wt% of step-2 reagent stream; dissolution solvent and stoichiometric co-reagent |
| F-05 | NaOH (30 wt% solution) | 24.0 kg NaOH + 56.0 kg H₂O | 1.1 equiv relative to Ester; slight excess ensures complete saponification |
| F-06 | HCl (10 wt% solution) | 18.0 kg HCl + 162.0 kg H₂O | Sized to cover IbupNa acidification; 1:1 stoichiometric with IbupNa |
| F-07 | Hexane (makeup) | 5.0 | Replaces hexane losses in the LLE recycle loop |

**Why TfOH is fed neat:** Water deactivates PhI(OAc)₂ by hydrolysis. If TfOH were fed as an aqueous solution, the water would carry over into the Reaction 2 feed and degrade the oxidant before it reacts. Neat TfOH also gives a higher effective catalyst concentration in R-101.

**Why 1.1 equivalents for ProAc and NaOH:** A small stoichiometric excess (1.05–1.15×) is standard practice in continuous flow to ensure the desired reagent is never limiting when the feed rate fluctuates slightly. Using 1.96 equiv ProAc (the pre-corrected value) would have required larger solvent recovery columns and increased raw material cost by ~80%.

---

## 3. Section 1 — Friedel-Crafts Acylation

### HX-101: Feed Preheater

S-01 (combined IBB, ProAc, TfOH, TfOH recycle) is heated from 25°C to 150°C before entering R-101.

**Justification:** The Friedel-Crafts acylation requires 150°C to achieve reasonable reaction rate. TfOH-catalysed acylation of electron-rich arenes is fast above 120°C but sluggish at room temperature. Preheating in a dedicated heat exchanger rather than relying on the reactor jacket avoids a large axial temperature gradient in R-101 that would reduce selectivity near the inlet.

### R-101: Friedel-Crafts PFR (150°C, 17 bar, τ = 3 min)

**Reaction:**
```
IBB + ProAc → Ketone + H₂O
      TfOH (cat.)
```

**Why PFR:** A plug-flow reactor is ideal here because the reaction is first-order in IBB and the selectivity to Ketone (vs. di-acylation by-product) is highest at low local concentration of IBB, which a PFR achieves by allowing concentration to fall monotonically along the reactor length. A CSTR operating at the same conversion would produce more di-acylation by-product because the exit concentration (low) exists everywhere in the tank.

**Why 17 bar:** The bubble-point pressure of the IBB/ProAc/TfOH mixture at 150°C is approximately 8–12 bar. Operating at 17 bar ensures the mixture remains fully liquid throughout, avoiding two-phase flow in a continuous PFR which would create slug flow and erratic residence time distribution.

**Why 72% conversion (not 91% as in the paper):** At production scale, pushing to near-complete conversion requires significantly longer residence time, a much larger reactor, and risks thermal runaway from the exothermic acylation. 72% is a deliberate conservative margin that keeps the reactor volume manageable and allows unreacted IBB to be recycled (via the organic phase from V-101) or tolerated in the subsequent step where it is inert.

**Reactor volume:** Calculated using Peng-Robinson EOS mixture density at 150°C/17 bar (not tabulated liquid density) to correctly account for thermal expansion of the liquid mixture at elevated temperature. At 150°C the mixture density drops ~8% versus 25°C, so using room-temperature density would undersize the reactor.

**By-products:** Water (1 mol per mol Ketone) is produced. This water must be removed before Reaction 2 because PhI(OAc)₂ hydrolyses rapidly in the presence of water.

### HX-102: Friedel-Crafts Product Cooler (150°C → 0°C)

The R-101 outlet is quenched to 0°C immediately after the reactor.

**Justification:** Two reasons. First, TfOH continues to catalyse the acylation even after the reagents have been mostly consumed; cooling to 0°C arrests the reaction and prevents over-acylation. Second, the subsequent V-101 liquid-liquid extraction is operated at 0°C where the organic/aqueous partition coefficients are more favourable (lower mutual miscibility, sharper phase split).

### V-101: Liquid-Liquid Extractor (NRTL, 0°C, 17 bar)

The cooled R-101 product is split into an organic phase (S-05) and an aqueous phase (S-06) using NRTL activity coefficient thermodynamics.

**What goes where:**
- **Organic (S-05):** IBB, Ketone, TfOH — all highly hydrophobic or, in the case of TfOH, forced to the organic phase by design
- **Aqueous (S-06):** Water (from Friedel-Crafts reaction) + ProAc (highly water-soluble at low temperature)

**Why remove water here:** If water passes into R-102 with the Ketone, it hydrolyses PhI(OAc)₂ according to PhI(OAc)₂ + H₂O → PhI(OH)(OAc) + AcOH, destroying the oxidant before it can react with the Ketone. The stoichiometry of Reaction 2 requires exactly 1 mol water, and this is provided by the TMOF hydrolysis in situ — any additional water from R-101 would upset this balance.

**Why keep TfOH in the organic phase:** TfOH acts as the Brønsted acid catalyst in both Reaction 1 (Friedel-Crafts) and Reaction 2 (aryl migration). Routing TfOH with the organic phase (S-05 → S-07 → R-102) means it automatically catalyses R-102 without a separate re-addition step. The alternative (routing TfOH to aqueous) would waste catalyst and require a fresh addition.

**Why operate at 17 bar:** The stream enters V-101 already at 17 bar from R-101. Reducing pressure at this point would cause flash vaporisation of the lighter components (TMOF, MeOH if present). Keeping the pressure elevated maintains single-phase liquid feed to V-101 and simplifies separator design.

### C-101: ProAc/Water Recovery Column (ProAc distillate, TfOH bottoms)

The aqueous phase from V-101 (S-06: ProAc + Water) is distilled to recover ProAc overhead and return any trace TfOH from the bottoms back to the R-101 feed.

**Justification:** ProAc (BP 141°C) separates cleanly from Water (BP 100°C) in a simple distillation. The relative volatility is approximately 0.7 (Water more volatile), so the column takes Water overhead and ProAc as bottoms — but in practice, a ProAc/Water azeotrope does not form, so a clean split is achievable. Recovering ProAc reduces fresh ProAc consumption. TfOH (BP 162°C) leaves as bottoms with the ProAc and is recycled directly to SM-101.

---

## 4. Section 2 — Aryl Migration

### SM-102: Mixing of R-102 Feed

The organic phase from V-101 (S-05: Ketone, IBB, TfOH) is combined with F-04 reagents (PhI(OAc)₂, TMOF, MeOH) and any recycled TMOF/MeOH from C-103 to form the R-102 inlet (S-07) at ~10°C.

**Why cool to 10°C before mixing:** PhI(OAc)₂ begins reacting with the Ketone above approximately 25°C. Mixing at 10°C allows a uniform, fully-dissolved feed to be established before the stream is heated to the 50°C reaction temperature. Mixing at reaction temperature would give uneven conversion due to local concentration gradients at the injection points.

### HX-103: R-102 Feed Heater (10°C → 50°C)

**Justification:** Reaction 2 operates at 50°C (paper Table 4: PFR-2 = 50°C). The feed is warmed to this temperature in a dedicated heat exchanger rather than in the reactor, for the same reason as HX-101 — it avoids a large axial thermal gradient in R-102 and ensures the reaction starts at the correct temperature along the full length of the reactor.

### R-102: Aryl Migration PFR, 4-Stage Multi-Injection (50°C, 17 bar, τ = 5 min)

**Reaction:**
```
Ketone + PhI(OAc)₂ + TMOF + H₂O → Ester + PhI + 2 AcOH + MeOH + HCOOCH₃
```

**Why PFR:** The aryl migration is an oxidative rearrangement where selectivity is critical — the Ester intermediate is the desired product, but over-oxidation or side reactions are possible at high local PhI(OAc)₂ concentration. A PFR with declining concentration along its length is more selective than a CSTR.

**Why 17 bar:** Same reason as R-101 — keeps the solvent mixture (TMOF, MeOH) fully liquid at 50°C and prevents flash.

**Why 75% conversion:** Conservative production-scale margin. The paper reports 98% conversion at lab scale (5.84 g/hr), but at 134 kg/hr the heat removal requirements and residence time needed to approach 98% would require a reactor roughly 4× larger. 75% conversion with recycle of unreacted Ketone (via the organic phase from V-103) is the more economical production choice.

**4-stage multi-injection PFR design:** At the production flow rate, adding all 327 kg/hr PhI(OAc)₂ at one inlet point gives a bulk concentration of 0.443 mol/L — approximately 3× above the solubility limit of PhI(OAc)₂ in TMOF/MeOH (0.10–0.15 mol/L). Precipitation of PhI(OAc)₂ in a continuous PFR causes plugging and erratic conversion. The solution is to split the PhI(OAc)₂ feed across 4 equally-spaced injection nozzles along R-102's length:

- Each nozzle delivers 81.8 kg/hr PhI(OAc)₂
- Per-stage inlet concentration: 0.111 mol/L ✓ (below 0.15 mol/L limit)
- By the time the next nozzle injects, the prior dose is 18–25% consumed at 75% overall conversion

This design directly reflects the "DISSOLVE" step shown in the Bogdan/Jolliffe flowsheet, scaled to production requirements. It adds minimal capital cost (4 injection nozzles vs. 1) but eliminates the most critical operability risk in the flowsheet.

**Full stoichiometry (why TMOF and Water are consumed):** TMOF (trimethyl orthoformate, HC(OCH₃)₃) is not merely a solvent in this step — it acts as a methoxy donor in the Ritter-type rearrangement. One mole of TMOF and one mole of water are stoichiometrically consumed per mole of Ester produced, yielding methanol (MeOH) and methyl formate (HCOOCH₃) as additional by-products alongside iodobenzene (PhI) and acetic acid (2× AcOH). Earlier versions of the simulation that omitted TMOF and Water consumption understated the AcOH production and overstated TMOF recovery in C-103.

**Why AcOH must not enter R-103:** R-102 produces 2 mol AcOH per mol Ester (at 75% conversion: ~1,100 mol/hr AcOH, ~66 kg/hr). If this AcOH carries forward to the saponification reactor, it will neutralise NaOH (AcOH + NaOH → NaOAc + H₂O) and consume base that should be saponifying the Ester. This would require roughly 3× more NaOH in F-05 to compensate. The vacuum flash (V-102) is the engineering solution to this problem.

### V-102: Vacuum Flash Drum (80°C, 0.25 bar)

The R-102 outlet is fed to a vacuum flash drum operating at 80°C and 0.25 bar to strip the volatile by-products overhead before saponification.

**What flashes overhead (S-10):** MeOH (BP 65°C), TMOF (BP 103°C), AcOH (BP 118°C), HCOOCH₃ (BP 32°C) — all vaporise readily at 80°C/0.25 bar.

**What stays in liquid bottoms (S-11):** Ester (BP 296°C), PhI (BP 188°C), unreacted Ketone (BP 270°C), IBB (BP 172°C) — all heavy organics with high boiling points that do not flash under these conditions.

**Why 0.25 bar (not atmospheric):** At atmospheric pressure (1 bar), AcOH (BP 118°C) does not flash efficiently at 80°C — its relative volatility over Ester is insufficient. Reducing pressure to 0.25 bar lowers AcOH's effective boiling point to approximately 55°C, ensuring it flashes cleanly and the bottoms AcOH content falls below 1.0 wt% (the saponification feed specification). At 0.30 bar this specification fails (1.06 wt%); at 0.25 bar it passes (0.84 wt%).

**Why this unit operation is not in the paper:** Jolliffe & Gerogiorgis model a lab-scale flowsheet (5.84 g/hr IBB). At lab scale, the absolute AcOH production is tiny and the NaOH excess in the saponification step can absorb it without meaningful cost impact. At production scale (134 kg/hr IBB), sending 66 kg/hr AcOH to R-103 would waste ~49 kg/hr NaOH, adding ~$400k/yr in raw material cost and producing excess sodium acetate waste. V-102 is a net-positive addition at scale.

### C-103: MeOH/TMOF Recovery Column (MeOH distillate, TMOF bottoms)

The V-102 overhead (S-10: MeOH, TMOF, AcOH, HCOOCH₃) is distilled to recover MeOH overhead and TMOF as bottoms, both recycled to SM-102.

**Justification:** TMOF costs approximately $1,200–1,500/MT. Recovering and recycling it from the V-102 overhead reduces fresh TMOF consumption (F-04) to only the makeup quantity needed to replace reaction losses. MeOH (produced by the reaction stoichiometry and from TMOF decomposition) is recovered at the top as it is lighter than TMOF. AcOH exits with the bottoms as a minor impurity — it is non-volatile relative to TMOF under column conditions and is purged from the system with C-103 bottoms.

---

## 5. Section 3 — Saponification

### SM-103: Mixing with NaOH Solution

The V-102 bottoms (S-11: Ester, PhI, unreacted Ketone/IBB) are combined with F-05 (NaOH/Water solution) to form the R-103 feed (S-12) at 40°C.

### R-103: Saponification PFR (65°C, 1.5 bar, τ = 7.5 min)

**Reaction:**
```
Ester + NaOH → IbupNa + MeOH
```

**Why PFR:** Saponification of esters is second-order overall (first-order in each of Ester and NaOH). In a PFR the reactant concentrations are highest at the inlet and fall along the length; this maximises the initial rate and achieves the target conversion efficiently. A CSTR would need to be several times larger to achieve the same 90% conversion.

**Why 65°C (not 70°C as initially coded):** Paper Table 4 explicitly lists PFR-3 at 65°C. The 70°C value was an error introduced when the simulation was first written. At 65°C, the saponification rate constant is still high (typical k ~ 10–50 L/mol/min for methyl esters), and the lower temperature reduces MeOH evaporation and minimises hydrolysis of IBB or PhI side reactions.

**Why 1.5 bar:** The mixture at 65°C contains MeOH (BP 65°C). At atmospheric pressure, MeOH would be at or near its boiling point, potentially causing bubble formation in the PFR and slug flow. Operating at 1.5 bar raises the effective boiling point of MeOH to ~90°C, ensuring fully liquid-phase operation.

**Why 90% conversion:** The saponification is relatively fast and high-conversion is achievable. However, 90% is chosen rather than 99% because the remaining 10% unreacted Ester partitions into the organic phase at V-103 and is effectively purged with the organic waste stream, keeping the purification simple. Pushing to 99%+ conversion would require a much longer reactor and increases the risk of Ester hydrolysis back-reaction under excess NaOH.

**Why NaOH (not KOH as in the paper):** The paper uses KOH following Bogdan 2009's original lab protocol. NaOH is chosen here because: (1) NaOH is approximately 3–4× cheaper than KOH per mole of base; (2) at production scale (~600 mol/hr base requirement) the cost difference is substantial; (3) NaOH saponification of methyl esters proceeds at essentially the same rate as KOH — the hydroxide ion is the active species. The final acidification step produces NaCl (vs. KCl with KOH), which is equally benign as a waste salt.

**NaOH stoichiometry:** F-05 is sized at 1.1 equiv relative to Ester entering R-103, ensuring NaOH is always in slight excess so Ester is the limiting reagent. The simulation checks NaOH/Ester ratio at convergence and flags it if it falls below 1.0×.

**Electrolyte modelling:** The R-103 outlet contains ionic species (Na⁺, IbupA⁻, OH⁻, NaOAc). The simulation computes ionic strength and Debye-Hückel activity coefficients for the key ions. This affects the partition behaviour in the subsequent LLE step (V-103) via the Setschenow salting-out correction.

### HX-104: R-103 Product Cooler (65°C → 30°C)

The saponification product is cooled to 30°C before phase separation.

**Justification:** V-103 LLE is operated at 25–30°C where the hexane/water partition coefficients are most favourable. Cooling also prevents flash vaporisation of MeOH and hexane in the separator.

---

## 6. Section 4 — Phase Separation, Acidification, and LLE

### V-103: Organic/Aqueous Phase Split (NRTL + Setschenow, 25°C, 1 bar)

The cooled R-103 outlet is split into an organic phase (S-15) and an aqueous phase (S-16) by liquid-liquid extraction with hexane.

**What goes where:**
- **Organic (S-15):** PhI, unreacted Ketone, IBB, PhI_OAc2, unreacted Ester — all hydrophobic, low logP compounds that strongly prefer the hexane phase
- **Aqueous (S-16):** IbupNa (ionic salt, fully water-soluble), Water, NaOH, NaOAc, MeOH

**Why hexane:** Hexane is a non-polar solvent with very low water miscibility and a high partition coefficient for the hydrophobic organic by-products (PhI, IBB, Ketone). Hexane is also easily recovered and recycled (low BP 69°C, recovered via C-102 condensation). Alternatives (ethyl acetate, dichloromethane) are either more water-miscible or pose greater environmental/regulatory concerns for pharmaceutical API production.

**Setschenow salting-out correction:** The aqueous phase contains NaCl (from incomplete ionisation) and NaOAc. Dissolved salts reduce the solubility of organic compounds in water (the salting-out effect), increasing their partition into the hexane phase. The simulation corrects the partition coefficients using the Setschenow equation: K_D,corrected = K_D × exp(Kₛ × m_salt), where Kₛ is the Setschenow constant for each compound and m_salt is the salt molality. This is important for Ibuprofen, which has meaningful water solubility at its free acid pKa. Neither IDAES nor DWSIM implements this correction for pharmaceutical intermediates; it is a unique feature of this simulation.

### Acidification (HCl addition → S-17)

The aqueous IbupNa stream (S-16) is mixed with F-06 (dilute HCl) to convert the sodium salt to free ibuprofen acid.

**Reaction:**
```
IbupNa + HCl → Ibuprofen + NaCl
```

**Why acidify before LLE:** IbupNa is ionic (water-soluble, hexane-insoluble). Free ibuprofen acid (pKa ~4.9) is much more hydrophobic (logP ≈ 3.97) and partitions strongly into hexane. Acidifying before hexane extraction therefore achieves a large, efficient single-step separation. Attempting to extract IbupNa directly into hexane would give negligible recovery.

**HCl stoichiometry:** 1.0 equiv HCl relative to IbupNa. Excess HCl would lower the pH further, potentially causing ProAc (present as a trace impurity) to partition into the organic phase and co-extract with ibuprofen.

### LLE: Hexane Extraction of Ibuprofen (4-stage, K_D = 8)

Free ibuprofen from the acidified aqueous stream is extracted into hexane using a 4-stage counter-current liquid-liquid extraction.

**Basis:** Distribution coefficient K_D = 8 (ibuprofen hexane/water at 25°C, pH ~2). With 4 extraction stages and E = K_D × V_org/V_aq ≈ 8–12, fractional recovery reaches 92%.

**Why 4 stages:** Single-stage extraction with K_D = 8 gives ~89% recovery (E/(1+E) with E=8). Four stages of counter-current extraction push this to >99% theoretical recovery; the simulation caps actual recovery at 92% to account for emulsion losses, wetting losses, and non-ideal distribution at production scale.

---

## 7. Section 5 — Crystallisation, Filtration, and Drying

### Hexane Flash (70°C, 0.9 bar)

The ibuprofen-hexane extract (S-19) is heated to 70°C to evaporate 97% of the hexane before crystallisation.

**Justification:** Crystallising from a large volume of hexane would require enormous crystalliser volumes and slow crystal growth. Pre-concentrating by flash evaporation reduces the hexane load to ~3% residual before the solvent switch to MeOH. The hexane vapour (S-21) is condensed and recycled via C-102.

### C-102: Hexane Condenser/Recovery

The hexane overhead from the flash step is condensed and recycled directly to the hexane LLE step.

**Justification:** Hexane loss to atmosphere or effluent is both an environmental and economic concern. The recycle loop (with 5 kg/hr fresh hexane makeup via F-07) maintains the hexane inventory at steady state. This is the fourth and simplest recycle loop in the process.

### CR-101: Crystalliser (Solvent Switch to MeOH, 45°C → 0°C)

MeOH is added to the concentrated ibuprofen (0.6 kg MeOH per kg ibuprofen) to carry out a solvent switch from hexane to MeOH, then the slurry is cooled to 0°C to crystallise ibuprofen.

**Why MeOH (not hexane directly):** Ibuprofen crystallises poorly from hexane — it tends to form oils rather than a filterable solid at the hexane concentrations used here. MeOH gives consistent needle or plate morphology with good filtration characteristics. The ibuprofen solubility in MeOH at 0°C is 0.138 kg/kg, meaning the crystallisation yield is thermodynamically favourable.

**Steady-state mother liquor recycle:** At steady state, the mother liquor (MeOH + dissolved ibuprofen) is recycled back to the crystalliser feed. This means that at steady state, the crystallisation yield is not limited by solubility — essentially all ibuprofen that enters the crystalliser is recovered as solid over multiple recycle passes. The steady-state recovery is therefore set by the filter/centrifuge cake losses (87%), not the crystallisation equilibrium.

### Filter and Dryer (87% yield, 99.7% MeOH removal)

The crystal slurry is filtered and the wet cake dried to remove MeOH to ICH Q3C levels.

**87% filter recovery:** Accounts for cake washing losses (~5%), fine particle entrainment in the mother liquor (~5%), and centrifuge heel (~3%).

**Dryer specification:** The dryer removes 99.7% of residual MeOH from the cake. The ICH Q3C limit for Class 2 solvents (MeOH) in a drug substance is 3,000 ppm. The simulation verifies that the final product (S-27) meets this specification at steady state.

---

## 8. Recycle Loops

Four material recycle streams are closed by damped successive-substitution iteration (convergence within 10 iterations, tolerance 0.01%):

| Loop | Stream | Source | Destination | Purpose |
|------|--------|--------|-------------|---------|
| TfOH recycle | R-TfOH | C-101 bottoms | SM-101 | Recover TfOH catalyst; reduces fresh TfOH makeup |
| TMOF recycle | R-TMOF | C-103 bottoms | SM-102 | Recover TMOF solvent; reduces fresh TMOF makeup |
| MeOH recycle | R-MeOH-S2 | C-103 distillate | SM-102 | Recover MeOH produced in Reaction 2; reduces fresh MeOH makeup |
| Hexane recycle | R-Hex | C-102 condenser | V-103 LLE | Recover hexane extraction solvent; reduces hexane losses to ~5 kg/hr makeup |

**Why recycle TfOH:** TfOH costs approximately $80/kg. At 749 kg/hr total requirement, a zero-recycle process would cost ~$47M/yr in TfOH alone. With recycle, the TfOH fresh feed converges to near-zero once the loop fills — the majority is recycled via C-101.

**Why the recycle loop needs iteration:** Each recycle stream changes the composition of its respective feed mixer, which changes the reactor inlet compositions, which changes the by-product rates, which changes what the separation units produce and recycle. The four loops are coupled (TMOF recycle affects R-102 conversion which affects the AcOH load on V-102 which affects the C-103 overhead composition). Iterative convergence is required rather than direct substitution.

---

## 9. Key Design Decisions vs. Paper

| Parameter | Paper (Jolliffe/Bogdan) | This Design | Reason |
|-----------|------------------------|-------------|--------|
| Scale | 50 kg/yr (lab) | ~640 MT/yr | Production-scale design |
| Reactor type (all) | PFR | PFR | Consistent |
| R-101 temperature | 150°C | 150°C | Consistent |
| R-102 temperature | 50°C | 50°C | Consistent |
| R-103 temperature | 65°C | 65°C | Consistent (corrected from 70°C) |
| Base for saponification | KOH | NaOH | Lower cost; equivalent chemistry |
| Inter-reactor separations | None (telescoped) | V-101 (LLE), V-102 (flash) | Required at production scale: V-101 removes Friedel-Crafts water to protect PhI(OAc)₂; V-102 removes AcOH to protect NaOH stoichiometry |
| PhI(OAc)₂ addition | Single point | 4-stage injection | Bulk concentration (0.44 mol/L) exceeds solubility limit; staged addition keeps per-stage value at 0.11 mol/L |
| R-101 conversion | 91% | 72% | Conservative production-scale margin; recycle of unreacted IBB via organic phase |
| R-102 conversion | 98% | 75% | Conservative production-scale margin |
| R-103 conversion | 99% | 90% | Conservative production-scale margin |
| TfOH feed | Neat | Neat | Consistent; corrected from earlier aqueous formulation |

---

## 10. Process Specifications (Design Basis)

| Specification | Limit | Basis |
|--------------|-------|-------|
| Annual production | ≥ 500 MT/yr | Design target |
| Ibuprofen purity | > 99.0 wt% | ICH Q6A pharmaceutical |
| Residual MeOH | < 3,000 ppm | ICH Q3C Class 2 solvent |
| V-101 water in organic | < 0.5 wt% | PhI(OAc)₂ hydrolysis limit |
| V-102 AcOH in bottoms | < 1.0 wt% | NaOH overconsumption threshold |
| NaOH/Ester ratio | ≥ 1.0 × | Ensure saponification not base-limited |
| HCl/IbupNa coverage | ≥ 95% | Acidification completeness |
| C-101, C-103 relative volatility | α > 1 | Distillation feasibility |

All specifications are met at the converged steady-state operating point.
