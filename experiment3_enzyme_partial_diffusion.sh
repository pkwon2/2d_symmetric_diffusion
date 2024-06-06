#!/bin/bash 

script='./run_inference.py'
ckpt='/projects/ml/aa_template/checkpoints/train_session2023-11-12_1699834651.0304408/models/BFF_1.pt'
pdb='./experiments/1CGZ_stateA_8B3C_B_ERD_switch_clean_suplig_fix_1CGZ_stateA.pdb'
pref='./experiments/test_enzyme_diffuse/CHS'

/software/containers/SE3nv.sif $script --config-name=aa \
    inference.model_runner='NRBStyleSelfCond' \
	contigmap.contigs=[\'A1-189,5-8,A196-334,2,A336-387,1-1\'] \
	model.symmetrize_repeats=False \
	inference.num_designs=2 \
	inference.input_pdb=$pdb \
	inference.ckpt_path=$ckpt \
	diffuser.T=50 \
	inference.ligand='PYC' \
	inference.ij_visible='abcd' \
	inference.two_template=True \
	inference.three_template=True \
	inference.motif_only_2d=True \
	preprocess.eye_frames=True \
	inference.supply_motif_seq=True \
	denoiser.noise_scale_ca=0 \
	inference.output_prefix=$pref \
	inference.contig_rmsd=0.25 \
	inference.check_rmsd_step=10
