import os,sys,glob,json,copy

from detect_left import is_left_handed 
import argparse

def make_refine_jobs(outdir,
                     run_script,
                     n_refine=2,
                     n_per_job=100,
                     ckpt_refine='/home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions/train_session2023-09-22_1695413763.3657782/models/BFF_4.pt',
                     ligand=None,
                     refine_w_ligand=False):
    """
    Creates refinement jobs for a folder of outputs
    """

    pdbs = glob.glob(outdir+'/*.pdb')

    BASE_INF_DICT = {'model_runner':'NRBStyleSelfCond',
                                'num_designs':n_refine,
                                'ckpt_path':ckpt_refine,
                                'two_template':True,
                                'three_template':True,
                                'cautious':True,
                                'motif_only_2d':True,
                                'supply_motif_seq':True,
                                'refine_recycles':4,
                                'refine':True,
                                'refine_w_ligand':refine_w_ligand}

    if ligand:
        assert refine_w_ligand, 'detected refine_w_ligand False but also detected supplied ligand is not None'
        BASE_INF_DICT['ligand'] = ligand
    if refine_w_ligand:
        assert ligand is not None, 'detected refine_w_ligand True but did not recieve specific ligand.'

    BASE_DIFF_DICT = {'T':50,
                    'so3_type':'random',
                    'eucl_type':'gaussian'}


    BASE_REFINE_DICT = {'inference':BASE_INF_DICT,
                        'diffuser':BASE_DIFF_DICT
                        }


    jobs = []

    for i,d in enumerate(pdbs):

        is_left,_ = is_left_handed(d)

        if True: ### NOTE: USED TO BE IF IS_LEFT BUT FXN IS TOO STINGENT AT THE MOMENT

            cur_dict = BASE_REFINE_DICT.copy()

            cur_dict['inference']['input_pdb'] = d

            jobs.append(copy.deepcopy(cur_dict))


    # write jsons
    json_outdir = os.path.join(outdir,'refine/')
    if not os.path.exists(json_outdir):
        os.makedirs(json_outdir)

    JOBS_PER_JSON = n_per_job
    for i in range(0,len(jobs),JOBS_PER_JSON):
        json_out = os.path.join(json_outdir, f'job_{i}.json')
        with open(json_out, 'w') as f:
            json.dump(jobs[i:i+JOBS_PER_JSON], f)



    # write job to execute refinement w/ sourcing jsons
    jsons = glob.glob(os.path.join(json_outdir, '*.json'))
    with open(os.path.join(json_outdir, 'refine_tasks.txt'), 'w') as fp:
        for j in jsons:
            line = f'/software/containers/SE3nv.sif {run_script} '
            line += f'--config-name aa inference.json_args="{j}" '
            fp.write(line+'\n')

    print('All done. Wrote jobs to '+os.path.join(json_outdir, 'refine_tasks.txt'))


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--outdir', '-o', type=str, required=True,
                        help='Output directory for refinement jobs')

    parser.add_argument('--run_script', '-r', type=str, required=True,
                        help='Path to your copy of run_inference.py')

    parser.add_argument('--ckpt_refine', '-c', type=str, required=True,
                        help='Path to your own copy of /home/davidcj/projects/train_2template/rf_diffusion/good_training_sessions/train_session2023-09-22_1695413763.3657782/models/BFF_4.pt')

    parser.add_argument('--n_refine', '-n', type=int, default=2,
                        help='Number of refinement runs per input') 

    parser.add_argument('--n_per_job', '-p', type=int, default=100,
                        help='Number of refinement runs per gpu job') 

    parser.add_argument('--ligand', '-l', type=str, default=None,
                        help='ligand in the pdbs being refinend')

    parser.add_argument('--refine_w_ligand', '-ref_w_lig', action='store_true', default=False,
                        help='If True, refine with ligand present. Make sure you use a model that was trained to refine in the presence of ligands')

    args = parser.parse_args()

    args.outdir = os.path.abspath(args.outdir)

    make_refine_jobs(args.outdir, 
                     args.run_script, 
                     args.n_refine, 
                     args.n_per_job, 
                     args.ckpt_refine,
                     args.ligand,
                     args.refine_w_ligand)


if __name__ == '__main__':
    main()

