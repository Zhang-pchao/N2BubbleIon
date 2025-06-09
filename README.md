## Dataset and model
  - The dataset for training or testing the DP model is uploaded to [AIS Square](https://www.aissquare.com/datasets/detail?pageType=datasets&name=SCAN_H2O_H3O_OH_N2_Nanobubble&id=261) and [Zenodo](https://zenodo.org/records/14308794).
  - The compressed DP model is uploaded to [AIS Square](https://www.aissquare.com/models/detail?pageType=models&name=SCAN_H2O_H3O_OH_N2_Nanobubble&id=256) and [Zenodo](https://zenodo.org/records/14308794).
  - In the [INCAR file](https://github.com/Zhang-pchao/N2BubbleIon/tree/main/DP-GEN_Iteration/INCAR), the charge of H₃O⁺ or OH⁻ is controlled by the NELECT keyword, which sets the number of electrons. All bulk subsystems in the DP datasets are neutral because of a homogeneous background charge. As shown in [VASP Wiki](https://www.vasp.at/wiki/index.php/NELECT), if the number of electrons is not compatible with the number derived from the valence and the number of atoms a homogeneous background charge is assumed. If the number of ions specified in the POSCAR file is 0 and NELECT=n, then the energy of a homogeneous electron gas is calculated.

## Packages Used
  - vasp_v5.4.4; deepmd-kit_v2.1.5; dpgen_v0.11.0; lammps_v20220623.1

## Paper

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
