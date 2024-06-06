#!/bin/bash 
#SBATCH -c 1
#SBATCH --mem=12g
#SBATCH -p gpu
#SBATCH --gres=gpu:a4000:1 
#SBATCH -t 60
#SBATCH -J dife_rep_epoch5

scripts="./run_inference.py"
ckpt_3template='/projects/ml/aa_template/checkpoints/train_session2023-10-03_1696375068.089569/models/BFF_5.pt'

pdb='./experiments/o2_102601_connect_dj.pdb'

# Partially diffusing 1hk9 with no noise using new sym until find clashing example 
# Then can be used as input to test differences in predictions between new sym and old sym 
prefix="./experiments/test_eye/test_lauko_w_ligand_BFF5" 

/software/containers/SE3nv.sif $scripts --config-name=aa \
inference.model_runner='NRBStyleSelfCond' \
inference.num_designs=2 \
inference.output_prefix=$prefix \
inference.ckpt_path=$ckpt_3template \
contigmap.contigs=[\'50-50,A1-4,50-50,B5-9,50-50\'] \
inference.two_template=True \
inference.three_template=True \
inference.cautious=True \
inference.input_pdb=$pdb \
diffuser.T=50 \
denoiser.noise_scale_frame=0.9 \
denoiser.noise_scale_ca=0.9 \
inference.ij_visible='abc' \
inference.motif_only_2d=True \
preprocess.eye_frames=True \
inference.supply_motif_seq=True \
inference.ligand='mu2' \
inference.align_motif=True \
inference.contig_rmsd=1000 #contig_rmsd as filter is not recommended for motig grafting jobs, just set a high number so it always pass. 
