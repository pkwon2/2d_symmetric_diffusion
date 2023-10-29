#!/bin/bash 

script='/mnt/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/run_inference.py'
ckpt='/mnt/home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions/train_session2023-10-03_1696375068.089569/models/BFF_7.pt'
pdb='/mnt/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/EF2_0696_MET.pdb'
pref='./experiments/test_pseudo/c4_asu50_EF2'

/software/containers/SE3nv.sif $script --config-name=aa \
    inference.model_runner='NRBStyleSelfCond' \
    contigmap.contigs=[\'200-200\'] \
    inference.output_prefix=$pref \
    model.symmetrize_repeats=True \
    model.repeat_length=50 \
    model.symmsub_k=3 \
    inference.num_designs=20 \
    inference.input_pdb=$pdb \
    inference.ckpt_path=$ckpt \
    model.sym_method='max' \
    model.copy_main_block_template=True \
    model.main_block=1 \
    model.pseudo_cycle=True \
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
    inference.n_repeats=4
