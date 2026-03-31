# Continuous-flow ibuprofen synthesis — Python simulation

Steady-state mass-and-energy balance for pharmaceutical-grade **ibuprofen** via the **Bogdan 2009** continuous-flow route (*Angew. Chem.*, hypervalent iodine chemistry). A single script — `ibuprofen_sim.py` — covers feeds, reactors, LLE/VLE flashes, distillation shortcuts, four-stream recycle convergence, electrolyte handling, equipment sizing, economics, and Excel export.

**Repository:** [github.com/jihanraiyan/cont_ibuprofen_python_sim](https://github.com/jihanraiyan/cont_ibuprofen_python_sim)

---

## Quick start

```bash
pip install pandas numpy scipy tabulate thermo chemicals openpyxl
python3 ibuprofen_sim.py
```

The run prints a full report and writes **`ibuprofen_sim_results.xlsx`** (12-sheet workbook) in the working directory.

---

## Documentation

| Resource | What it is |
|----------|------------|
| **[SIMULATION_DOCUMENTATION.md](SIMULATION_DOCUMENTATION.md)** | Full technical manual — NRTL/PR-EOS, unit ops, tears, economics, Excel sheet map |
| **[PROCESS_DESCRIPTION.md](PROCESS_DESCRIPTION.md)** | Narrative process description and design rationale |

**Authoritative** stoichiometry, conversions, temperatures, and stream logic are always those in `ibuprofen_sim.py`.

---

## Model highlights (v5)

- T-dependent NRTL; Rachford–Rice LLE; vacuum VLE with Antoine \(P^\text{sat}\)
- Electrolyte handling (Debye–Hückel, Setschenow) where applicable
- Four-stream recycle (TfOH, TMOF, MeOH, hexane)
- Turton-style CAPEX/OPEX / P&L (basis year in-code)

### Recent model updates

- **V-101:** TfOH forced to the **aqueous** phase → C-101 recovery and recycle to R-101 only (step 2 does not use TfOH in the organic path, per PFD).
- **F-04:** Stoichiometric **water** in the step-2 reagent stream for aryl-migration stoichiometry.
- **F-05 / F-06:** **NaOH** and **HCl** sized **dynamically** from molar flows (base for ester saponification + PhI(OAc)₂ hydrolysis; acid for excess NaOH neutralisation then IbupNa).
- **C-101:** Stricter TfOH overhead loss fraction (`fHK_dist`) with physical basis in comments.

---

## License

[MIT](LICENSE)

---

## Reference

Bogdan, A. R. *et al.* “The Continuous-Flow Synthesis of Ibuprofen.” *Angew. Chem. Int. Ed.* **2009**, *48*, 8547–8550. [DOI: 10.1002/anie.200900575](https://doi.org/10.1002/anie.200900575)

Scale-up context: Jolliffe, H. G. & Gerogiorgis, D. I., *Chem. Eng. Res. Des.* **2015**, *97*, 175–191.
