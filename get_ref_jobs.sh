#!/bin/bash

outdir='/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/experiments/test_FAD/'
script='/mnt/home/davidcj/projects/rf_diffusion_2template/rf_diffusion/run_inference.py'
ckpt='/home/davidcj/projects/train_2template/rf_diffusion/train_session2023-11-04_1699130567.230787/models/BFF_1.pt'

/software/containers/SE3nv.sif ./make_refine_jobs.py --outdir $outdir --run_script $script --ckpt_refine $ckpt --n_refine 2 --n_per_job 100 --ligand 'FAD' --refine_w_ligand


