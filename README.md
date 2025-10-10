# N2BubbleIon Repository Overview

This repository collects the data sets, models, molecular dynamics (MD) inputs, and analysis tools used for the study of nitrogen nanobubbles interacting with hydronium and hydroxide ions in water. The material here enables researchers to reproduce the deep potential (DP) based workflow described in the linked publication and to further explore ion-modulated nanobubble dynamics.

## Dataset and Model Resources

- The dataset for training or testing the DP model is uploaded to [AIS Square](https://www.aissquare.com/datasets/detail?pageType=datasets&name=SCAN_H2O_H3O_OH_N2_Nanobubble&id=261) and [Zenodo](https://zenodo.org/records/14308794).
- The compressed DP model is uploaded to [AIS Square](https://www.aissquare.com/models/detail?pageType=models&name=SCAN_H2O_H3O_OH_N2_Nanobubble&id=256) and [Zenodo](https://zenodo.org/records/14308794).
- In the [INCAR file](https://github.com/Zhang-pchao/N2BubbleIon/tree/main/DP-GEN_Iteration/INCAR), the charge of H₃O⁺ or OH⁻ is controlled by the NELECT keyword, which sets the number of electrons. All bulk subsystems in the DP datasets are neutral because of a homogeneous background charge. As shown in [VASP Wiki](https://www.vasp.at/wiki/index.php/NELECT), if the number of electrons is not compatible with the number derived from the valence and the number of atoms a homogeneous background charge is assumed. If the number of ions specified in the POSCAR file is 0 and NELECT=n, then the energy of a homogeneous electron gas is calculated.

## Software Environment

- vasp_v5.4.4
- deepmd-kit_v2.1.5
- dpgen_v0.11.0
- lammps_v20220623.1

## Repository Structure

The directories in this repository are organized to cover data generation, training, molecular dynamics simulations, and post-processing:

- **Analysis_Scripts/** – Python scripts used to post-process MD trajectories. They cover bubble density calculations (`BubbleDensity_each.py`, `BubbleDensity_final.py`), stress analysis (`Bubble_Stress.py`, `Bubble_eachStress.py`), hydrogen-bond statistics (`BubbleHydrogenBonds_each.py`, `BubbleHydrogenBonds_final.py`, `IonHydrogenBonds.py`), ion–nitrogen coordination analysis (`Many_Ions_Find.py`, `N2_neighbor_ion.py`, `N2_neighbor_ion_each.py`, `Ion_neighbor_N2.py`), cluster turnover tracking (`N2Cluster_Turnover_final.py`), trajectory re-centering (`shift_cluster_center_trj.py`), and virial/kinetic energy conversions (`virial_ke2xyz_trj.py`).
- **DP-GEN_Iteration/** – Input decks for the DP-GEN active learning cycles, including VASP configuration files (`INCAR/INCAR_SCAN_AIMD_E552_C600_e5`, `INCAR/INCAR_SCAN_SP_E552_C600_e5`), example structural templates (`POSCAR/mixIon`, `POSCAR/mix2Ion`), system-specific run settings (`System_index_*` directories with `run.json` files), and the global execution settings (`run_machine.json`).
- **DP_Training/** – Artifacts from training the DP model, such as the trained frozen model (`frozen_model.pb`), learning curve output (`lcurve.out`), the configuration for launching training (`run.json`), and a representative subset of training data under `partial_dataset/interface/` containing raw boxes, coordinates, energies, forces, atom types, and virials (`*.raw`, `set.000`).
- **Molecular_Dynamics/** – LAMMPS inputs for running nanobubble simulations. The directory provides system topologies (`D30L55*.data` files for different ion counts and placements), atomic configurations with ion centers (`*_atomic.data`), and the driving script (`in.bubble`).

Use these descriptions as a starting point to locate the files relevant to your workflow—whether you are re-running the DP-GEN iterations, retraining the model, performing MD simulations, or analyzing trajectory data.

## Related Publication

Hydroxide and Hydronium Ions Modulate the Dynamic Evolution of Nitrogen Nanobubbles in Water. [J. Am. Chem. Soc.](https://pubs.acs.org/doi/10.1021/jacs.4c06641)

```bibtex
@article{Zhang_JAmChemSoc_2024_v146_p19537,
  title        = {{Hydroxide and Hydronium Ions Modulate the Dynamic Evolution of Nitrogen Nanobubbles in Water}},
  author       = {Pengchao Zhang and Changsheng Chen and Muye Feng and Chao Sun and Xuefei Xu},
  year         = 2024,
  journal      = {J. Am. Chem. Soc.},
  volume       = 146,
  number       = 28,
  pages        = {19537--19546},
  doi          = {10.1021/jacs.4c06641},
}
```
