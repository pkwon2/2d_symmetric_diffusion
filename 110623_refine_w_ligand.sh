#!/bin/bash 
#SBATCH -c 1
#SBATCH --mem=12g
#SBATCH -p gpu
#SBATCH --gres=gpu:a4000:1 
#SBATCH -t 60
#SBATCH -J dife_rep_epoch5

script="/mnt/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/run_inference.py"

ckpt='/home/davidcj/projects/train_2template/rf_diffusion/train_session2023-11-04_1699130567.230787/models/BFF_1.pt'

pdb='/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/experiments/test_FAD/test1_3.pdb'

# Partially diffusing 1hk9 with no noise using new sym until find clashing example 
# Then can be used as input to test differences in predictions between new sym and old sym 

/software/containers/SE3nv.sif $script --config-name=aa \
inference.model_runner='NRBStyleSelfCond' \
inference.num_designs=2 \
inference.ckpt_path=$ckpt \
inference.two_template=True \
inference.three_template=True \
inference.cautious=False \
inference.input_pdb=$pdb \
diffuser.T=50 \
inference.motif_only_2d=True \
inference.supply_motif_seq=True \
diffuser.so3_type='random' \
diffuser.eucl_type='gaussian' \
inference.refine_recycles=4 \
inference.refine=True \
inference.refine_w_ligand=True \
inference.ligand='FAD'
