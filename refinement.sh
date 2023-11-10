script='/mnt/home/lhtran/software/rf_diffusion_pseudocycle/rf_diffusion/run_inference.py'
#/projects/ml/aa_template/checkpoints/refine/train_session2023-11-04_1699130567.230787/BFF_3.pt
ckpt='/global/cfs/cdirs/m4129/users/lhtran/rf_diffusion_pseudocycle_chkpt/refinement/BFF_1.pt'
ckpt='/projects/ml/aa_template/checkpoints/refine/train_session2023-11-04_1699130567.230787/BFF_3.pt'
#pdb='/pscratch/sd/l/lhtran/projects/dig_cid/diffusion/0172/BFF_12/sym2_90aa_69.pdb'
pdb='//home/lhtran/projects/rf_diffusion_pseudocycle/sym4_50aa_0.pdb'
pdb='/mnt/home/lhtran/software/rf_diffusion_pseudocycle/DIG_0122_ref2015_PCsym3_70aa_T50_bbf13_22.pdb'
pref='/home/lhtran/projects/rf_diffusion_pseudocycle/'

# Partially diffusing 1hk9 with no noise using new sym until find clashing example 
# Then can be used as input to test differences in predictions between new sym and old sym 

/software/containers/SE3nv.sif $script --config-name=aa \
inference.model_runner='NRBStyleSelfCond' \
inference.output_prefix=$pref \
inference.num_designs=2 \
inference.ckpt_path=$ckpt \
inference.two_template=True \
inference.three_template=True \
inference.cautious=False \
inference.input_pdb=$pdb \
diffuser.T=50 \
inference.pseudo_symmetry='c1' \
inference.pseudocycle_break=210 \
inference.ij_visible='a' \
inference.motif_only_2d=True \
inference.supply_motif_seq=True \
diffuser.so3_type='random' \
diffuser.eucl_type='gaussian' \
inference.refine_recycles=4 \
inference.refine=True \
inference.refine_w_ligand=True \
inference.ligand='DIG'

