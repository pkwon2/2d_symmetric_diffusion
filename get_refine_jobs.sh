#!/bin/bash

/software/containers/SE3nv.sif util/make_refine_jobs.py --outdir /home/linnaan/software/rf_diffusion_pseudocycle_ligand/rf_diffusion/experiments/test_pseudo/ --run_script '/home/linnaan/software/rf_diffusion_pseudocycle_ligand/rf_diffusion/run_inference.py' --ckpt_refine '/projects/ml/aa_template/checkpoints/refine/train_session2023-09-22_1695413763.3657782/models/BFF_4.pt' --n_refine 2 --n_per_job 100


