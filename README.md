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
# Pretrained checkpoint included in the repository
export RFDIFFUSION_CKPT_PATH="$WEIGHTS_DIR/BFF_0.20.pt"
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

The model checkpoint is located at:

```text
2d_symmetric_diffusion/
└── weights/
    └── BFF_0.20.pt
```

The checkpoint path is automatically configured by the environment
variables defined in Step 2:

```bash
export RFDIFFUSION_CKPT_PATH="$RFDIFFUSION_DIR/weights/BFF_0.20.pt"
```

No additional checkpoint download is required if the full checkpoint
is included in the cloned repository.

Verify that the checkpoint is present:

```bash
ls -lh "$RFDIFFUSION_CKPT_PATH"
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

The modified RFdiffusion implementation and pretrained checkpoint are
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
apptainer run --nv "$APPTAINER_PATH" \\
    "$RFDIFFUSION_DIR/run_inference.py" \\
    --config-name=aa \\
    inference.ckpt_path="$RFDIFFUSION_CKPT_PATH" \\
    inference.model_runner=NRBStyleSelfCond \\
    inference.num_designs=1 \\
    inference.output_prefix="$DESIGN_DIR/c2/design" \\
    contigmap.contigs="['C1-199,20-20,C221-255,26-26 D1-53 B1-199,20-20,B221-255,26-26 F1-53']" \\
    inference.two_template=True \\
    inference.three_template=True \\
    inference.rigid_repeat_motif=False \\
    diffuser.T=50 \\
    model.symmetrize_repeats=True \\
    model.symmsub_k=1 \\
    model.main_block=0 \\
    model.sym_method=max \\
    denoiser.noise_scale_frame=0.05 \\
    denoiser.noise_scale_ca=0.05 \\
    inference.ij_visible=abcdef \\
    inference.motif_only_2d=True \\
    preprocess.eye_frames=True \\
    inference.supply_motif_seq=True \\
    inference.cyclic_protein_indices=False \\
    inference.input_pdb="$INPUT_DIR/array_rechain.pdb" \\
    inference.n_repeats=2 \\
    model.repeat_length=333 \\
    model.pseudo_cycle=False \\
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
apptainer run --nv "$APPTAINER_PATH" \\
    "$RFDIFFUSION_DIR/run_inference.py" \\
    --config-name=aa \\
    inference.ckpt_path="$RFDIFFUSION_CKPT_PATH" \\
    inference.model_runner=NRBStyleSelfCond \\
    inference.num_designs=3 \\
    inference.output_prefix="$DESIGN_DIR/c3-c2/design" \\
    contigmap.contigs="['40-40 40-40 40-40 40-40 40-40 40-40']" \\
    inference.two_template=True \\
    inference.three_template=True \\
    inference.rigid_repeat_motif=False \\
    diffuser.T=50 \\
    model.symmetrize_repeats=True \\
    model.symmsub_k=1 \\
    model.main_block=0 \\
    model.sym_method=max \\
    denoiser.noise_scale_frame=0.05 \\
    denoiser.noise_scale_ca=0.05 \\
    inference.motif_only_2d=True \\
    preprocess.eye_frames=True \\
    inference.supply_motif_seq=True \\
    inference.cyclic_protein_indices=True \\
    inference.input_pdb="$INPUT_DIR/D3.pdb" \\
    inference.n_repeats=6 \\
    model.repeat_length=40 \\
    model.pseudo_cycle=c3-c2 \\
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
apptainer run --nv "$APPTAINER_PATH" \\
    "$RFDIFFUSION_DIR/run_inference.py" \\
    --config-name=aa \\
    inference.ckpt_path="$RFDIFFUSION_CKPT_PATH" \\
    inference.model_runner=NRBStyleSelfCond \\
    inference.num_designs=3 \\
    inference.output_prefix="$DESIGN_DIR/d3/design" \\
    contigmap.contigs="['40-40 40-40 40-40 40-40 40-40 40-40']" \\
    inference.two_template=True \\
    inference.three_template=True \\
    inference.rigid_repeat_motif=False \\
    diffuser.T=50 \\
    model.symmetrize_repeats=True \\
    model.symmsub_k=1 \\
    model.main_block=0 \\
    model.sym_method=max \\
    denoiser.noise_scale_frame=0.05 \\
    denoiser.noise_scale_ca=0.05 \\
    inference.motif_only_2d=True \\
    preprocess.eye_frames=True \\
    inference.supply_motif_seq=True \\
    inference.cyclic_protein_indices=True \\
    inference.input_pdb="$INPUT_DIR/D3.pdb" \\
    inference.n_repeats=6 \\
    model.repeat_length=40 \\
    model.pseudo_cycle=d3 \\
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

### 4. Motif Scaffolding with D3 Symmetry

This example uses six input PDB motif segments (chains A–F, residues 1–158) and appends 37 designed residues to each segment. It enables 2D motif conditioning and cyclic protein indices under the D3 symmetry configuration.

```bash
mkdir -p "$DESIGN_DIR/d3_motif_scaffolding"
apptainer run --nv "$APPTAINER_PATH" \
    "$RFDIFFUSION_DIR/run_inference.py" \
    --config-name=aa \
    inference.ckpt_path="$RFDIFFUSION_CKPT_PATH" \
    inference.model_runner=NRBStyleSelfCond \
    inference.num_designs=3 \
    inference.output_prefix="$DESIGN_DIR/d3_motif_scaffolding/design" \
    contigmap.contigs="['A1-158,37-37 B1-158,37-37 C1-158,37-37 D1-158,37-37 E1-158,37-37 F1-158,37-37']" \
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
    potentials.guiding_potentials=[\"type:olig_contacts,weight_intra:1,weight_inter:0.06\"] \
    potentials.olig_intra_all=True potentials.olig_inter_all=True potentials.guide_scale=2 \
    potentials.guide_decay="cubic" \
    inference.align_px0_motif=False
```

Key distinctions from the preceding D3 example are the explicit fixed motif ranges in `contigmap.contigs`, `model.repeat_length=195` and the `inference.ij_visible=abcdef` settings. The input `inputs/D3.pdb` must include the referenced chain IDs and residue numbers.

## Additional Parameters

The following parameters are shared across the inference examples.
| Parameter | Description |
|---|---|
| `--config-name=aa` | Uses the All-Atom RFdiffusion configuration. |
| `inference.ckpt_path` | Path to the pretrained model checkpoint. |
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
└── d3/
```

Each directory contains the inference outputs for its corresponding
symmetry specification.

## Notes

### GPU Requirements

All-Atom RFdiffusion requires an NVIDIA GPU and a compatible CUDA
environment.

The inference commands use Apptainer with the `--nv` flag to expose
NVIDIA GPUs to the container.

### Model Checkpoint

The provided examples use:

```text
weights/BFF_0.20.pt
```

The checkpoint path is automatically defined relative to the
repository directory.

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
