#!/bin/bash 
#SBATCH -c 1
#SBATCH --mem=12g
#SBATCH -p gpu
#SBATCH --gres=gpu:a4000:1 
#SBATCH -t 60
#SBATCH -J dife_rep_epoch5

script="/mnt/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/run_inference.py"

#ckpt_3template='/home/davidcj/projects/train_2template/rf_diffusion/train_session2023-09-06_1694049198.4710052/models/BFF_12.pt'
ckpt_3template='/home/davidcj/projects/train_2template/rf_diffusion/train_session2023-09-14_1694727421.6775959/models/BFF_0.10.pt'

pdb='/home/davidcj/projects/train_2template/rf_diffusion/train_session2023-09-13_1694634209.8062925/training_pdbs/epoch_1_80_diff_pdb_aa_t_200_true.pdb'

# Partially diffusing 1hk9 with no noise using new sym until find clashing example 
# Then can be used as input to test differences in predictions between new sym and old sym 
prefix="./experiments/test_eye/test_epoch1_80_w_seq_and_center_finetune" 

/software/containers/SE3nv.sif $script --config-name=aa \
inference.model_runner='NRBStyleSelfCond' \
inference.num_designs=15 \
inference.output_prefix=$prefix \
inference.ckpt_path=$ckpt_3template \
contigmap.contigs=[\'142,A165-168,8,A177-180,7,A188-193,46\'] \
inference.two_template=True \
inference.three_template=True \
inference.cautious=True \
inference.input_pdb=$pdb \
diffuser.T=50 \
denoiser.noise_scale_frame=0.5 \
denoiser.noise_scale_ca=0.5 \
inference.ij_visible='abc' \
inference.motif_only_2d=True \
preprocess.eye_frames=True \
inference.supply_motif_seq=True \
