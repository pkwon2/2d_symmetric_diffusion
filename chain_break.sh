
#!/bin/bash 

script='/mnt/home/lhtran/software/rf_diffusion_pseudocycle/rf_diffusion/run_inference.py'
ckpt='/mnt/home/lhtran/software/rf_diffusion_repeat_dev/train_session2023-10-03_1696375068.089569/models/BFF_7.pt'
pdb='/home/lhtran/projects/sm_binders/ef2_binder/EF2_0696_MET.pdb'
pref='/net/scratch/lhtran/repeat_dev/chain_break_debug/EF2_160_C4_dimmer'

/software/containers/SE3nv.sif $script --config-name=aa \
    inference.model_runner='NRBStyleSelfCond' \
    contigmap.contigs=[\'160-160\'] \
    inference.output_prefix=$pref \
    model.symmetrize_repeats=True \
    model.repeat_length=40 \
    model.symmsub_k=3 \
    inference.num_designs=20 \
    inference.input_pdb=$pdb \
    inference.ckpt_path=$ckpt \
    model.sym_method='max' \
    model.copy_main_block_template=True \
    model.main_block=1 \
    model.pseudo_cycle=True \
    inference.pseudocycle_break=80 \
    diffuser.T=50 \
    inference.ligand='EF2' \
    inference.ij_visible='a' \
    inference.two_template=True \
    inference.three_template=True \
    inference.motif_only_2d=True \
    preprocess.eye_frames=True \
    inference.supply_motif_seq=True \
    inference.pseudo_symmetry='c4' \
    inference.align_px0_motif=False \
    inference.n_repeats=4 2>&1 | tee ./EF2_160_C4_dimmer.log
