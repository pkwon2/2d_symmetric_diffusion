# All-Atom Diffusion with 2D-Templated Motifs, Repeat Proteins, and Symmetry

## Basic Information

This repository contains the modified All-Atom RFdiffusion implementation
used in the following preprint:

**Controlling metal-carbonate phase, form, and function through de novo protein design**

Paul S. Kwon, Xinqi Li, Le Tracy Yu, Harley Pyles, Todd H. Lewis,
Connor Weidle, Andrew J. Borst, Catherine C. Bodinger, Alex Kang,
Hannah Nguyen, Johnny Mendoza, Kenneth D. Carr, Brian Coventry,
Yang Hsia, Zsombor Molnar, Dongsheng Li, Bo Zhang,
Brandi M. Cossairt, Shuai Zhang, Asim K. Bera,
James De Yoreo, and David Baker.

**bioRxiv** (2026).

[Read the preprint](https://doi.org/10.64898/2026.06.10.730916)

### Citation

If you use this code in your research, please cite:

```bibtex
@article{kwon2026controlling,
  title={Controlling metal-carbonate phase, form, and function through de novo protein design},
  author={Kwon, Paul S and Li, Xinqi and Yu, Le Tracy and Pyles, Harley
    and Lewis, Todd H and Weidle, Connor and Borst, Andrew J
    and Bodinger, Catherine C and Kang, Alex and Nguyen, Hannah
    and Mendoza, Johnny and Carr, Kenneth D and Coventry, Brian
    and Hsia, Yang and Molnar, Zsombor and Li, Dongsheng
    and Zhang, Bo and Cossairt, Brandi M and Zhang, Shuai
    and Bera, Asim K and De Yoreo, James and Baker, David},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.06.10.730916},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Overview

This repository extends All-Atom RFdiffusion to support protein design
with 2D-templated motifs, repeat protein architectures, and symmetric
protein assemblies.

The implementation supports:

- **2D-templated motif conditioning:** Conditioning diffusion using

  two-dimensional residue-pair representations to control structural
  relationships between motifs.

- **Repeat protein design:** Generating architectures composed of
  repeated protein units.

- **Cyclic symmetry:** Designing protein assemblies with cyclic
  symmetry, including C2.

- **Dihedral symmetry:** Generating protein assemblies with
  dihedral symmetry, including D3.

- **Combined symmetry:** Supporting combined symmetry specifications,
  such as C3-C2, for designing complex protein architectures.

These methods were developed to explore protein architectures for
biomineralization, including protein arrays and symmetric assemblies
capable of templating inorganic crystal formation.

## Installation and Setup

### 1. Clone the repository

Clone this repository:

```bash
git clone https://github.com/pkwon2/2d_symmetric_diffusion.git
cd 2d_symmetric_diffusion
```

If the repository stores model checkpoints using Git LFS, ensure that
Git LFS is installed and download the checkpoint files:

```bash
git lfs pull
```

### 2. Set environment variables and create directories

The following environment variables simplify the inference commands.

Run these commands from the root of the repository:

```bash
# Repository directory
export RFDIFFUSION_DIR="$(pwd)"
# Model weights
export WEIGHTS_DIR="$RFDIFFUSION_DIR/weights"
# Checkpoint for design generation
export RFDIFFUSION_CKPT_PATH="$WEIGHTS_DIR/BFF_0.20.pt"
# Checkpoint for refinement
export RFDIFFUSION_REFINE_CKPT_PATH="$WEIGHTS_DIR/BFF_3.pt"
# Apptainer environment directory
export ENV_DIR="$RFDIFFUSION_DIR/exec"
# Apptainer image used by all inference examples
export APPTAINER_PATH="$ENV_DIR/SE3nv.sif"
# Output directory
export DESIGN_DIR="$RFDIFFUSION_DIR/outputs"
# Input PDB directory
export INPUT_DIR="$RFDIFFUSION_DIR/inputs"
# Create directories
mkdir -p "$ENV_DIR" "$DESIGN_DIR" "$INPUT_DIR"
```

### 3. Model weights

The design and refinement commands use separate checkpoints:

| Purpose | Environment variable | Checkpoint |
|---|---|---|
| Design generation | `RFDIFFUSION_CKPT_PATH` | `weights/BFF_0.20.pt` |
| Refinement | `RFDIFFUSION_REFINE_CKPT_PATH` | `weights/BFF_3.pt` |

Place the checkpoint files in the repository's `weights/` directory:

```text
2d_symmetric_diffusion/
└── weights/
    ├── BFF_0.20.pt
    └── BFF_3.pt
```

The checkpoint paths are configured by the environment variables in Step 2:

```bash
export RFDIFFUSION_CKPT_PATH="$WEIGHTS_DIR/BFF_0.20.pt"
export RFDIFFUSION_REFINE_CKPT_PATH="$WEIGHTS_DIR/BFF_3.pt"
```

If a checkpoint is not included in the cloned repository or available
through Git LFS, copy it into the `weights/` directory before running
the corresponding command.

Verify that both checkpoints are present:

```bash
ls -lh "$RFDIFFUSION_CKPT_PATH" "$RFDIFFUSION_REFINE_CKPT_PATH"
```

### 4. Download the Apptainer environment

This repository uses an [Apptainer](https://apptainer.org/)
container to provide the dependencies required for All-Atom RFdiffusion.

The container can be downloaded from the RFDpoly repository's
public distribution.

Navigate to the environment directory:

```bash
cd "$ENV_DIR"
```

Download the Apptainer image:

```bash
curl -fL -o "$APPTAINER_PATH" https://files.ipd.uw.edu/pub/2025_RFDpoly/SE3nv.sif
```

Make the downloaded container executable:

```bash
chmod +x "$APPTAINER_PATH"
```

Set the Apptainer path:

```bash
export APPTAINER_PATH="$ENV_DIR/SE3nv.sif"
```

Verify the container:

```bash
ls -lh "$APPTAINER_PATH"
```

Downloading the `.sif` file may take several minutes.

**Note:** The downloaded container provides the software environment.

The modified RFdiffusion implementation and pretrained checkpoints are
provided separately through this repository.

The examples below use NVIDIA GPUs through Apptainer's `--nv` option.

To check the environment before inference:

```bash
test -f "$RFDIFFUSION_CKPT_PATH" || { echo "Missing checkpoint: $RFDIFFUSION_CKPT_PATH"; exit 1; }
test -f "$APPTAINER_PATH" || { echo "Missing container: $APPTAINER_PATH"; exit 1; }
apptainer exec --nv "$APPTAINER_PATH" python --version
```

### 5. Prepare input structures

Some inference examples require input PDB structures defining motif
geometry or symmetry relationships.

The examples below use:

```text
inputs/
├── array_rechain.pdb
└── D3.pdb
```

Place the corresponding PDB files in the `inputs/` directory.

The input structure paths are configured using:

```bash
export INPUT_DIR="$RFDIFFUSION_DIR/inputs"
```

## Example Run Commands

The following examples demonstrate different symmetry specifications
supported by this implementation.

All examples use the All-Atom RFdiffusion configuration. Run the setup commands in your current shell first (they set `APPTAINER_PATH`, `RFDIFFUSION_DIR`, and other required variables):

```bash
--config-name=aa
```

Before running the examples, ensure that the repository paths,
model checkpoint, Apptainer environment, and input PDB structures
have been configured as described above.

Refinement commands are provided after the generation examples for C2,
C3-C2, and D3. They use the `BFF_3.pt` checkpoint configured by
`RFDIFFUSION_REFINE_CKPT_PATH` and require the updated
[refinement mask support](#refinement-mask-support) described in Example 5.
Keep the matching `.trb` file beside each selected input PDB and retain
access to the original motif PDB recorded in it. Each command refines
`design_0.pdb`; change the filename to select another generated design.
The three refined structures are saved beside the input as
`design_0_refined_0.pdb`, `design_0_refined_1.pdb`, and
`design_0_refined_2.pdb`.

### 1. Cyclic Symmetry (C2)

This example generates a protein assembly with C2 cyclic symmetry
using 2D-templated motif conditioning.

The design uses two repeated units with a repeat length of 333 residues.

Create the output directory:

```bash
mkdir -p "$DESIGN_DIR/c2"
```

Run inference:

```bash
apptainer run --nv "$APPTAINER_PATH" \
    "$RFDIFFUSION_DIR/run_inference.py" \
    --config-name=aa \
    inference.ckpt_path="$RFDIFFUSION_CKPT_PATH" \
    inference.model_runner=NRBStyleSelfCond \
    inference.num_designs=1 \
    inference.output_prefix="$DESIGN_DIR/c2/design" \
    contigmap.contigs="['C1-199,20-20,C221-255,26-26 D1-53 B1-199,20-20,B221-255,26-26 F1-53']" \
    inference.two_template=True \
    inference.three_template=True \
    inference.rigid_repeat_motif=False \
    diffuser.T=50 \
    model.symmetrize_repeats=True \
    model.symmsub_k=1 \
    model.main_block=0 \
    model.sym_method=max \
    denoiser.noise_scale_frame=0.05 \
    denoiser.noise_scale_ca=0.05 \
    inference.ij_visible=abcdef \
    inference.motif_only_2d=True \
    preprocess.eye_frames=True \
    inference.supply_motif_seq=True \
    inference.cyclic_protein_indices=False \
    inference.input_pdb="$INPUT_DIR/array_rechain.pdb" \
    inference.n_repeats=2 \
    model.repeat_length=333 \
    model.pseudo_cycle=False \
    inference.align_px0_motif=False
```

**Key parameters:**
| Parameter | Description |
|---|---|
| `inference.n_repeats=2` | Specifies two repeated units. |
| `model.repeat_length=333` | Sets the repeat length to 333 residues. |
| `model.symmetrize_repeats=True` | Enables repeat symmetrization. |
| `inference.motif_only_2d=True` | Enables 2D motif conditioning. |
| `inference.ij_visible=abcdef` | Specifies visible motif information. |
| `model.pseudo_cycle=False` | Disables pseudo-cycle symmetry mode. |
| `inference.input_pdb` | Specifies the input motif structure. |

#### Refine the generated C2 design

After generation completes, refine `c2/design_0.pdb` using
2 repeats of 333 residues.

```bash
export REFINE_INPUT_PDB="$DESIGN_DIR/c2/design_0.pdb"

apptainer run --nv "$APPTAINER_PATH" \
    "$RFDIFFUSION_DIR/run_inference.py" \
    --config-name=aa \
    inference.ckpt_path="$RFDIFFUSION_REFINE_CKPT_PATH" \
    inference.model_runner=NRBStyleSelfCond \
    inference.num_designs=3 \
    inference.two_template=True \
    inference.three_template=True \
    inference.rigid_repeat_motif=False \
    diffuser.T=50 \
    diffuser.so3_type=random \
    diffuser.eucl_type=gaussian \
    inference.refine=True \
    inference.refine_recycles=4 \
    model.symmetrize_repeats=True \
    model.symmsub_k=1 \
    model.main_block=0 \
    model.sym_method=max \
    inference.motif_only_2d=True \
    inference.supply_motif_seq=True \
    inference.cyclic_protein_indices=True \
    inference.input_pdb="$REFINE_INPUT_PDB" \
    inference.ij_visible=abcdef \
    inference.n_repeats=2 \
    model.repeat_length=333 \
    inference.align_px0_motif=False
```

### 2. Combined C3-C2 Symmetry

This example demonstrates combined C3 and C2 symmetry.

This configuration can be useful for designing protein arrays
and larger assemblies with multiple rotational symmetry axes.

The design contains six repeated units, each consisting of
40 residues.

Create the output directory:

```bash
mkdir -p "$DESIGN_DIR/c3-c2"
```

Run inference:

```bash
apptainer run --nv "$APPTAINER_PATH" \
    "$RFDIFFUSION_DIR/run_inference.py" \
    --config-name=aa \
    inference.ckpt_path="$RFDIFFUSION_CKPT_PATH" \
    inference.model_runner=NRBStyleSelfCond \
    inference.num_designs=3 \
    inference.output_prefix="$DESIGN_DIR/c3-c2/design" \
    contigmap.contigs="['40-40 40-40 40-40 40-40 40-40 40-40']" \
    inference.two_template=True \
    inference.three_template=True \
    inference.rigid_repeat_motif=False \
    diffuser.T=50 \
    model.symmetrize_repeats=True \
    model.symmsub_k=1 \
    model.main_block=0 \
    model.sym_method=max \
    denoiser.noise_scale_frame=0.05 \
    denoiser.noise_scale_ca=0.05 \
    inference.motif_only_2d=True \
    preprocess.eye_frames=True \
    inference.supply_motif_seq=True \
    inference.cyclic_protein_indices=True \
    inference.input_pdb="$INPUT_DIR/D3.pdb" \
    inference.n_repeats=6 \
    model.repeat_length=40 \
    model.pseudo_cycle=c3-c2 \
    inference.align_px0_motif=False
```

**Key parameters:**
| Parameter | Description |
|---|---|
| `inference.n_repeats=6` | Specifies six repeated units. |
| `model.repeat_length=40` | Sets each repeat length to 40 residues. |
| `model.pseudo_cycle=c3-c2` | Selects combined C3-C2 symmetry mode. |
| `model.symmetrize_repeats=True` | Enables repeat symmetrization. |
| `inference.num_designs=3` | Generates three designs. |

#### Refine the generated C3-C2 design

After generation completes, refine `c3-c2/design_0.pdb` using
6 repeats of 40 residues.

```bash
export REFINE_INPUT_PDB="$DESIGN_DIR/c3-c2/design_0.pdb"

apptainer run --nv "$APPTAINER_PATH" \
    "$RFDIFFUSION_DIR/run_inference.py" \
    --config-name=aa \
    inference.ckpt_path="$RFDIFFUSION_REFINE_CKPT_PATH" \
    inference.model_runner=NRBStyleSelfCond \
    inference.num_designs=3 \
    inference.two_template=True \
    inference.three_template=True \
    inference.rigid_repeat_motif=False \
    diffuser.T=50 \
    diffuser.so3_type=random \
    diffuser.eucl_type=gaussian \
    inference.refine=True \
    inference.refine_recycles=4 \
    model.symmetrize_repeats=True \
    model.symmsub_k=1 \
    model.main_block=0 \
    model.sym_method=max \
    inference.motif_only_2d=True \
    inference.supply_motif_seq=True \
    inference.cyclic_protein_indices=True \
    inference.input_pdb="$REFINE_INPUT_PDB" \
    model.pseudo_cycle=c3-c2 \
    inference.n_repeats=6 \
    model.repeat_length=40 \
    inference.align_px0_motif=False
```

### 3. Dihedral Symmetry (D3)

This example generates a protein assembly with D3 dihedral symmetry.

D3 symmetry consists of a threefold rotational axis and
three perpendicular twofold rotational axes.

The design contains six repeated units of 40 residues each.

Create the output directory:

```bash
mkdir -p "$DESIGN_DIR/d3"
```

Run inference:

```bash
apptainer run --nv "$APPTAINER_PATH" \
    "$RFDIFFUSION_DIR/run_inference.py" \
    --config-name=aa \
    inference.ckpt_path="$RFDIFFUSION_CKPT_PATH" \
    inference.model_runner=NRBStyleSelfCond \
    inference.num_designs=3 \
    inference.output_prefix="$DESIGN_DIR/d3/design" \
    contigmap.contigs="['40-40 40-40 40-40 40-40 40-40 40-40']" \
    inference.two_template=True \
    inference.three_template=True \
    inference.rigid_repeat_motif=False \
    diffuser.T=50 \
    model.symmetrize_repeats=True \
    model.symmsub_k=1 \
    model.main_block=0 \
    model.sym_method=max \
    denoiser.noise_scale_frame=0.05 \
    denoiser.noise_scale_ca=0.05 \
    inference.motif_only_2d=True \
    preprocess.eye_frames=True \
    inference.supply_motif_seq=True \
    inference.cyclic_protein_indices=True \
    inference.input_pdb="$INPUT_DIR/D3.pdb" \
    inference.n_repeats=6 \
    model.repeat_length=40 \
    model.pseudo_cycle=d3 \
    inference.align_px0_motif=False
```

**Key parameters:**
| Parameter | Description |
|---|---|
| `inference.n_repeats=6` | Specifies six repeated units. |
| `model.repeat_length=40` | Sets each repeat length to 40 residues. |
| `model.pseudo_cycle=d3` | Enables the D3 symmetry mode. |
| `model.symmetrize_repeats=True` | Enables repeat symmetrization in this example. |
| `inference.num_designs=3` | Generates three designs. |

#### Refine the generated D3 design

After generation completes, refine `d3/design_0.pdb` using
6 repeats of 40 residues.

```bash
export REFINE_INPUT_PDB="$DESIGN_DIR/d3/design_0.pdb"

apptainer run --nv "$APPTAINER_PATH" \
    "$RFDIFFUSION_DIR/run_inference.py" \
    --config-name=aa \
    inference.ckpt_path="$RFDIFFUSION_REFINE_CKPT_PATH" \
    inference.model_runner=NRBStyleSelfCond \
    inference.num_designs=3 \
    inference.two_template=True \
    inference.three_template=True \
    inference.rigid_repeat_motif=False \
    diffuser.T=50 \
    diffuser.so3_type=random \
    diffuser.eucl_type=gaussian \
    inference.refine=True \
    inference.refine_recycles=4 \
    model.symmetrize_repeats=True \
    model.symmsub_k=1 \
    model.main_block=0 \
    model.sym_method=max \
    inference.motif_only_2d=True \
    inference.supply_motif_seq=True \
    inference.cyclic_protein_indices=True \
    inference.input_pdb="$REFINE_INPUT_PDB" \
    model.pseudo_cycle=d3 \
    inference.n_repeats=6 \
    model.repeat_length=40 \
    inference.align_px0_motif=False
```

### 4. Motif Scaffolding with D3 Symmetry

This example uses six input PDB motif segments (chains A–F, residues 18–158), with 17 designed residues before and 37 designed residues after each motif. Each repeat contains 195 residues (17 + 141 + 37). It enables 2D motif conditioning and cyclic protein indices under the D3 symmetry configuration.

```bash
mkdir -p "$DESIGN_DIR/d3_motif_scaffolding"

apptainer run --nv "$APPTAINER_PATH" \
    "$RFDIFFUSION_DIR/run_inference.py" \
    --config-name=aa \
    inference.ckpt_path="$RFDIFFUSION_CKPT_PATH" \
    inference.model_runner=NRBStyleSelfCond \
    inference.num_designs=3 \
    inference.output_prefix="$DESIGN_DIR/d3_motif_scaffolding/design" \
    "contigmap.contigs=['17-17,A18-158,37-37 17-17,B18-158,37-37 17-17,C18-158,37-37 17-17,D18-158,37-37 17-17,E18-158,37-37 17-17,F18-158,37-37']" \
    inference.two_template=True \
    inference.three_template=True \
    inference.rigid_repeat_motif=False \
    diffuser.T=50 \
    model.symmetrize_repeats=True \
    model.symmsub_k=1 \
    model.main_block=0 \
    model.sym_method=max \
    denoiser.noise_scale_frame=0.05 \
    denoiser.noise_scale_ca=0.05 \
    inference.motif_only_2d=True \
    preprocess.eye_frames=True \
    inference.supply_motif_seq=True \
    inference.cyclic_protein_indices=True \
    inference.input_pdb="$INPUT_DIR/D3.pdb" \
    inference.ij_visible=abcdef \
    inference.n_repeats=6 \
    model.repeat_length=195 \
    model.pseudo_cycle=d3 \
    'potentials.guiding_potentials=["type:olig_contacts,weight_intra:1,weight_inter:0.06"]' \
    potentials.olig_intra_all=True \
    potentials.olig_inter_all=True \
    potentials.guide_scale=2 \
    potentials.guide_decay=cubic \
    inference.align_px0_motif=False
```

Key distinctions from the preceding D3 example are the explicit fixed motif ranges in `contigmap.contigs`, `model.repeat_length=195` and the `inference.ij_visible=abcdef` settings. The input `inputs/D3.pdb` must include the referenced chain IDs and residue numbers.

### 5. Frame Refinement with D3 Repeat Symmetry

The generation command in Example 4 produces the Cα scaffold. Refine that
scaffold with `BFF_3.pt` to predict the backbone frames, using the separate
`RFDIFFUSION_REFINE_CKPT_PATH` checkpoint variable.

#### Refinement mask support

The refinement command below requires the updated mask condition in
`NRBStyleSelfCond._get_3template_masks()` in `inference/model_runners.py`.
Within the `if self._conf.inference.motif_only_2d:` block, the condition is:

```python
if refine or not self._conf.model.symmetrize_repeats:
```

This routes refinement through the existing full-structure motif-mask
logic, which restores motif indices and pair visibility from the saved
metadata. It allows `model.symmetrize_repeats=True` during refinement.
Generation with repeat symmetrization still uses the repeat-specific
mask branch. Older code with only
`if not self._conf.model.symmetrize_repeats:` raises
`AssertionError: refine not yet implemented for symmetry/repeat` for
this command.

#### Select a generated input

Select one of the PDB files generated by Example 4:

```bash
export REFINE_INPUT_PDB="$DESIGN_DIR/d3_motif_scaffolding/design_0.pdb"
```

To refine a different design, replace `design_0.pdb` with its filename,
such as `design_6.pdb` if it exists from a previous generation run.
Keep the matching `.trb` file beside the selected input PDB, and retain
access to the original motif PDB recorded in that file.
Refinement uses these files to restore the original motif conditioning.
The explicit `inference.ij_visible` setting below is superseded by the
saved refinement metadata. Refinement constructs contigs from the
generated structure, so the original generation contigs are not repeated
in this command.

#### Run refinement

This configuration has been successfully run with the updated model runner:

```bash
apptainer run --nv "$APPTAINER_PATH" \
    "$RFDIFFUSION_DIR/run_inference.py" \
    --config-name=aa \
    inference.ckpt_path="$RFDIFFUSION_REFINE_CKPT_PATH" \
    inference.model_runner=NRBStyleSelfCond \
    inference.num_designs=3 \
    inference.two_template=True \
    inference.three_template=True \
    inference.rigid_repeat_motif=False \
    diffuser.T=50 \
    diffuser.so3_type=random \
    diffuser.eucl_type=gaussian \
    inference.refine=True \
    inference.refine_recycles=4 \
    model.symmetrize_repeats=True \
    model.symmsub_k=1 \
    model.main_block=0 \
    model.sym_method=max \
    inference.motif_only_2d=True \
    inference.supply_motif_seq=True \
    inference.cyclic_protein_indices=True \
    inference.input_pdb="$REFINE_INPUT_PDB" \
    inference.ij_visible=abcdef \
    inference.n_repeats=6 \
    model.repeat_length=195 \
    model.pseudo_cycle=d3 \
    inference.align_px0_motif=False
```

**Key refinement parameters:**

| Parameter | Description |
|---|---|
| `inference.ckpt_path="$RFDIFFUSION_REFINE_CKPT_PATH"` | Selects the `BFF_3.pt` refinement checkpoint. |
| `inference.refine=True` | Enables refinement of the generated input structure. |
| `inference.refine_recycles=4` | Uses four model recycles for refinement. |
| `model.symmetrize_repeats=True` | Ties repeat pair features inside the model. |
| `model.pseudo_cycle=d3` | Selects the D3 repeat-pair symmetry pattern. |
| `inference.n_repeats=6` | Specifies six repeats. |
| `model.repeat_length=195` | Specifies 195 residues per repeat. |
| `inference.cyclic_protein_indices=True` | Uses cyclic protein-index features, matching generation. |

Leave `inference.pseudo_symmetry` out of this command. The model's
`model.pseudo_cycle=d3` setting controls the repeat-pair pattern.
`inference.pseudo_symmetry=d3` is a separate operation that replaces
coordinates with rotations about fixed global axes, which can disrupt
a scaffold generated in a different orientation or repeat order.
Do not explicitly pass `inference.pseudo_symmetry=null`: the checkpoint
override parser attempts to cast it through `NoneType` and raises
`TypeError: NoneType takes no arguments`. Omitting the argument retains
its default `null` value.

The refinement command uses the checkpoint's frame preprocessing defaults;
it does not set `preprocess.eye_frames=True` or
`inference.start_from_input=True`. The refinement path bypasses the
denoiser update, so generation-stage denoiser noise-scale settings are
not included here. Repeat-feature tying does not mathematically guarantee
exact all-atom coordinate symmetry; inspect the refined structures.

Refinement derives its output prefix from the selected input PDB path by
replacing `.pdb` with `_refined`. For `design_6.pdb`, the three outputs are
`design_6_refined_0.pdb`, `design_6_refined_1.pdb`, and
`design_6_refined_2.pdb`, saved beside the input PDB.

## Additional Parameters

The following parameters are used across the inference examples.
| Parameter | Description |
|---|---|
| `--config-name=aa` | Uses the All-Atom RFdiffusion configuration. |
| `inference.ckpt_path` | Path to the design or refinement checkpoint for the current stage. |
| `inference.model_runner=NRBStyleSelfCond` | Uses the self-conditioning inference runner. |
| `diffuser.T=50` | Sets the diffusion process to 50 steps. |
| `inference.two_template=True` | Enables two-template functionality. |
| `inference.three_template=True` | Enables three-template functionality. |
| `inference.rigid_repeat_motif=False` | Disables rigid repeat-motif handling. |
| `model.symmsub_k=1` | Sets the symmetry-subunit parameter. |
| `model.main_block=0` | Selects model block configuration 0. |
| `model.sym_method=max` | Uses max as the symmetry aggregation method. |
| `denoiser.noise_scale_frame=0.05` | Controls frame noise during denoising. |
| `denoiser.noise_scale_ca=0.05` | Controls C-alpha noise during denoising. |
| `inference.motif_only_2d=True` | Enables the 2D motif-conditioning mode. |
| `preprocess.eye_frames=True` | Enables eye-frame preprocessing. |
| `inference.supply_motif_seq=True` | Supplies motif sequence information. |
| `inference.cyclic_protein_indices=False` | Disables cyclic protein-index handling. |
| `inference.align_px0_motif=False` | Disables motif alignment during the specified inference step. |

## Outputs

Generated structures are saved to the directory specified by
`inference.output_prefix`.

For example:

```text
outputs/
├── c2/
├── c3-c2/
├── d3/
└── d3_motif_scaffolding/
```

Each directory contains the inference outputs for its corresponding
symmetry specification. Refinement outputs are saved alongside the
selected input design using a prefix ending in `_refined`.

## Notes

### GPU Requirements

All-Atom RFdiffusion requires an NVIDIA GPU and a compatible CUDA
environment.

The inference commands use Apptainer with the `--nv` flag to expose
NVIDIA GPUs to the container.

### Model Checkpoints

Design generation uses `weights/BFF_0.20.pt` through
`RFDIFFUSION_CKPT_PATH`. Refinement uses `weights/BFF_3.pt` through
`RFDIFFUSION_REFINE_CKPT_PATH`.

Both checkpoint paths are defined relative to the repository directory
by the setup commands.

### Input Structures

The examples require the following PDB files:

```text
inputs/array_rechain.pdb
inputs/D3.pdb
```

Ensure that these files are available before running inference.

### Initial Execution

The first inference run may take longer because certain model
components, including diffusion-related caches, may need
to be initialized.

Subsequent runs may be faster.

### Container Compatibility

The provided Apptainer image is downloaded from the RFDpoly
distribution.

### Associated Preprint

Kwon, P. S., et al. (2026).

**Controlling metal-carbonate phase, form, and function through
de novo protein design.**

**bioRxiv**.

https://doi.org/10.64898/2026.06.10.730916
