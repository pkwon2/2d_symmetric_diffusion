#!/bin/bash

script='/mnt/home/lhtran/software/rf_diffusion_pseudocycle/rf_diffusion/run_inference.py'
ckpt='/mnt/home/lhtran/software/rf_diffusion_repeat_dev/train_session2023-10-03_1696375068.089569/models/BFF_7.pt'
#pdb='/home/lhtran/projects/sm_binders/ef2_binder/EF2_0696_MET.pdb'
pdb='/home/lhtran/projects/sm_binders/est_binder/round_1/EST_0330_0001_relax_fix_genpot_Met.pdb'
pref='/net/scratch/lhtran/repeat_dev/chain_break/EST_120_C4_dimmer_T_30_break_sym'
#pref='/net/scratch/lhtran/repeat_dev/chain_break/EF2_160_C2_dimmer_user_T_25_break_sym'


/software/containers/SE3nv.sif $script --config-name=aa \
    inference.model_runner='NRBStyleSelfCond' \
    contigmap.contigs=[\'120-120\'] \
    inference.output_prefix=$pref \
    model.symmetrize_repeats=True \
    model.repeat_length=30 \
    model.symmsub_k=3 \
    inference.num_designs=30 \
    inference.input_pdb=$pdb \
    inference.ckpt_path=$ckpt \
    model.sym_method='max' \
    model.copy_main_block_template=True \
    model.main_block=1 \
    model.pseudo_cycle=True \
    inference.pseudocycle_break=60 \
    diffuser.T=50 \
    inference.ligand='EST' \
    inference.ij_visible='abcd' \
    inference.two_template=True \
    inference.three_template=True \
    inference.motif_only_2d=True \
    preprocess.eye_frames=True \
    inference.supply_motif_seq=True \
    inference.pseudo_symmetry='c4' \
    inference.align_px0_motif=False \
    inference.T_break_sym=30\
    inference.n_repeats=4 2>&1 | tee ./EST_120_C4_dimmer_break_sym.log
