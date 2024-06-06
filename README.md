# All-Atom Diffusion w/ 2D templated motifs, repeat proteins, and symmetry. 

## Set up 

1. Make a new folder where you would like to place your copy of the repository. Once there, do `git clone https://github.com/w-ahern/rf_diffusion.git`
2. `cd` into the `rf_diffusion` that was just created. Now, do `git checkout aa-repeats-3template`.
3. Now, you need to get the RF2 allatom code in your `rf_diffusion` folder. To do this, do:
```
git submodule init
git submodule update
```
4. `cd` into `RF2-allatom/rf2aa` which should now be full of code. Now do `git checkout 3template_pseudo1`.
5. Stay where you are, do `git pull`.

If this series of steps doesn't work, let DJ know! 


## Running the code
Note that getting backbones from this protocol is a two step process. (1) You will diffuse a backbone w/ the diffusion model (2) You must then refine the backbone with a refinement model. Don't be alarmed if the backbones straight out of diffusion are displayed with dashed lines in the backbone in pymol. This is normal, and will go away once you refine the backbone w/ the refinement model.

### Location of saved model weights (as of Oct 24, 2023). 
The current best checkpoints for T2D templated motif scaffolding all live in here: `/home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions`. 
You will have to make your own copy of the weights for now until I find a better place to put them, as we cannot share files on the digs when running scripts. 

Checkpoint I reccomend at the moment for **protein-only stuff**: `/home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions/train_session2023-09-14_1694727421.6775959/models/BFF_0.20.pt` 

Checkpoints I reccomend at the moment for **protein and ligand stuff**:

- `/mnt/home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions/train_session2023-10-03_1696375068.089569/models/BFF_6.pt`
- `/mnt/home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions/train_session2023-10-03_1696375068.089569/models/BFF_5.pt`

I have tested `BFF_5` slightly more than `BFF_6`, but I would test both and see what you think. 

Checkpoints I reccomend for **refinement after diffusion** (see below for description of refinement): 
- `/home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions/train_session2023-09-22_1695413763.3657782/models/BFF_4.pt`

### Asymmetric motif scaffolding w/ ligands
The following is an example of a script which takes one of Cullen's metalo-hydrolase motifs containing a small molecule, and some discontiguous chunks, and scaffolding them into a monomer. The flags that you the user can change are 
- `inference.num_designs`: num designs you want
- `inference.output_prefix`: prefix you want 
- `contigmap.contigs`: Contigs - each motif chunk in the contigs is represented in the `ij_visible` flag below via a letter. The letter mappings are `first chunk --> a`, `second chunk --> b`, etc. I.e., the position from left to right in the contigs where a chunk is corresponds to its letter in the alphabet for the `ij_visible` flag. See description of this flag below. 
- `denoiser.noise_scale_ca`
- `inference.ij_visible`: Used to specify which chunks of motif are constrained w.r.t every other chunk. The syntax is dash (`-`) separated groups of letters, where each group specifies chunks that are constrained w.r.t each other. E.g., if I have 3 motif chunks and I only care that the first two are rigid w.r.t each other (in their conformation in the input pdb), and the third one could float wherever, I would specify `inference.ij_visible='ab-c'`. If you want them all constrained w.r.t eachother, do `inference.ij_visible='abc'`. And if none of them are constrained w.r.t eachother and all free to float relative to each other, do `inference.ij_visible=a-b-c`. **NOTE**: Importantly, if you have a ligand and want it constrained relative to other protein chunks in your contigs, you must **add an extra letter to ij_visible** corresponding to the ligand. In Cullen's metalo-enzyme example below, this is exemplified by having `inference.ij_visible='abcde'` even though there are only `abcd` worth of chunks in contigs. `e` there is for the ligand being constrained rigidly with respect to the other chunks.
  
- `inference.ligand`: name of ligand in the input pdb. Must be a three-letter name otherwise pdb writing gets messed up! 
  
```
#!/bin/bash 

script="/mnt/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/run_inference.py"

ckpt_3template='/home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions/train_session2023-10-03_1696375068.089569/models/BFF_5.pt'

pdb='/mnt/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/combined_0_mut_with_lig.pdb'

# Partially diffusing 1hk9 with no noise using new sym until find clashing example 
# Then can be used as input to test differences in predictions between new sym and old sym 
prefix="./experiments/test_eye/test_cullen_1"

/software/containers/SE3nv.sif $script --config-name=aa \
inference.model_runner='NRBStyleSelfCond' \
inference.num_designs=15 \
inference.output_prefix=$prefix \
inference.ckpt_path=$ckpt_3template \
contigmap.contigs=[\'20,B63-65,40,B92-96,40,B102-102,40,B149-149,20\'] \
inference.two_template=True \
inference.three_template=True \
inference.cautious=True \
inference.input_pdb=$pdb \
diffuser.T=50 \
denoiser.noise_scale_frame=0.9 \
denoiser.noise_scale_ca=0.9 \
inference.ij_visible='abcde' \
inference.motif_only_2d=True \
preprocess.eye_frames=True \
inference.supply_motif_seq=True \
inference.ligand='EAC' \
```

### Repeat protein motif scaffolding 

You can also make repeat proteins that scaffold a motif, and you can have the multiple motif copies have pre-defined distances/orientations w.r.t each other, or they can be figured out on the fly by the model. 

Below is an example of specifing a repeat-protein diffusion motif scaffolding run, where the multiple motif copies are constrained to be what is in the input pdb file (using the `ij_visible` flag described above). 

The only non-trivial flags that you should care about for now are: 
- `contigmap.contigs`: Can be used to specify entire repeat-protein's worth of contigs, in the case of constraining motif copies w.r.t other copies in the input pdb. However, can also specify an asymmetric unit's worth of contig if running in the mode where the multiple motif copies are not constrained w.r.t each other. The length of the contigs string you provide must be consistent with the `n_repeats` and `repeat_length` flags below. 
- `n_repeats`: Specifies how many repeat units you are simulating/designing.
- `model.repeat_length`: Yes this is redundant and could be calculated automatically, but for now, please put the length of an asymmetric unit in your repeat protein here. This is `L_tot/n_repeats`, where `L_tot` is the length of the entire repeat protein being designed in the simulation.

```
#!/bin/bash 

script="/mnt/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/run_inference.py"

ckpt_3template='/home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions/train_session2023-09-14_1694727421.6775959/models/BFF_0.20.pt'


# The full PDB used as symm template 
pdb='/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/repeat_input.pdb'

# Partially diffusing 1hk9 with no noise using new sym until find clashing example 
# Then can be used as input to test differences in predictions between new sym and old sym 
prefix="./experiments/thuddy/test1"

/software/containers/SE3nv.sif $script --config-name=aa \
inference.model_runner='NRBStyleSelfCond' \
inference.num_designs=25 \
inference.output_prefix=$prefix \
inference.ckpt_path=$ckpt_3template \
inference.n_repeats=3 \
contigmap.contigs=[\'60,A53-66,60,A98-111,60,A143-156\'] \
inference.two_template=True \
inference.three_template=True \
inference.rigid_repeat_motif=False \
inference.cautious=True \
inference.input_pdb=$pdb \
diffuser.T=50 \
model.symmetrize_repeats=True \
model.repeat_length=74 \
model.symmsub_k=1 \
model.main_block=0 \
model.sym_method='max' \
denoiser.noise_scale_frame=0.5 \
denoiser.noise_scale_ca=0.5 \
inference.ij_visible='abc' \
inference.motif_only_2d=True \
preprocess.eye_frames=True \
inference.supply_motif_seq=True \
```
### new features
#### stop at early steps
Bcov found it useful to stop diffusion traj if protein complex interfaces are not formed well after certain steps, because the continous success will low rate anyway. We add the feature to stop the diffusion traj after `inference.check_rmsd_step` if `inference.contig_rmsd` is still larger than cutoff. LA find it useful when doing partial diffusion for binder/enzymes. LA does not recommend it for motif grafting, because you high chance need some terminal motif change to allow your key rotamers to be respected.

#### examples
`experiment1_enzyme_motif_grafting.sh`: command example for motif grafting
`experiment2_free_diffuse_pseudocycle_sm.sh`: command example for pseudocycle diffusion
`experiment3_enzyme_partial_diffusion.sh`: command example for enzyme partial diffusion

note: currently, the last string cannot be contig, for example, you can do [\'A1-100,5,A199-200,1-1\'], but not [\'A1-100,5,A199-200\']. This is due to how the contig read function writes.
note2: one key feature found to interfere with protein quality is `noise_scale_ca`, while scaffold diversity increases with large noise, protein quality drop as well. Recommended range [0-0.5], can test [0,0.01,0.1,0.2...0.5] for design

