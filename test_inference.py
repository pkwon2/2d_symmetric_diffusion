import sys
import os
import glob
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'RF2-allatom'))

import unittest
from unittest import mock
import subprocess
from pathlib import Path
from inspect import signature
from io import StringIO

import hydra
from hydra import compose, initialize
from hydra.core.hydra_config import HydraConfig
from icecream import ic
import torch
import numpy as np

import test_utils
import run_inference
from functools import partial
from rf2aa import tensor_util
import rf2aa.chemical
from rf2aa.RoseTTAFoldModel import RoseTTAFoldModule
import rf2aa.loss
import inference.utils
ic.configureOutput(includeContext=True)

import shutil

from rf2aa import util as rf2aa_util

import ipdb
#import traceback

REWRITE = False
def infer(overrides):
    conf = construct_conf(overrides)
    run_inference.main(conf)
    p = Path(conf.inference.output_prefix + '_0.pdb')
    return p, conf

def construct_conf(overrides):
    initialize(version_base=None, config_path="config/inference", job_name="test_app")
    conf = compose(config_name='aa_small.yaml', overrides=overrides, return_hydra_config=True)
    # This is necessary so that when the model_runner is picking up the overrides, it finds them set on HydraConfig.
    HydraConfig.instance().set_config(conf)
    conf = compose(config_name='aa_small.yaml', overrides=overrides)
    return conf

def get_trb(conf):
    path = conf.inference.output_prefix + '_0.trb'
    return np.load(path,allow_pickle=True)


def clean_up():
    #os.replace(old,dest+os.path.basename(old))
    shutil.rmtree("tmp/")
if os.path.exists("tmp/"):
    clean_up()

# class TestRegression(unittest.TestCase): #TODO

#     def setUp(self) -> None:
#         # Some other test is leaving a global hydra initialized, so we clear it here.
#         if hydra.core.global_hydra.GlobalHydra().is_initialized():
#             hydra.core.global_hydra.GlobalHydra().clear()
#         return super().setUp()

#     def tearDown(self):
#         hydra.core.global_hydra.GlobalHydra.instance().clear()
    
#     # Example regression test.
#     def test_t2(self):
#         run_inference.make_deterministic()
#         pdb, _ = infer([
#             'diffuser.T=3',
#             'inference.num_designs=1',
#             'inference.output_prefix=tmp/test_2',
#             "contigmap.contigs=['1,A10-11,1']",
#             "contigmap.length=4-4",
#             'inference.design_startnum=0',
#             'inference.two_template=True',
#             'inference.three_template=True',
#             'inference.motif_only_2d=True',
#             'model.main_block=1',
#             'inference.align_px0_motif=False',
#             "inference.ij_visible='a'",
#             'inference.supply_motif_seq=True'
#         ])
        
#         pdb_contents = inference.utils.parse_pdb(pdb)
#         cmp = partial(tensor_util.cmp, atol=5e-2, rtol=0)
#         test_utils.assert_matches_golden(self, 'T2', pdb_contents, rewrite=REWRITE, custom_comparator=cmp)


class TestInference(unittest.TestCase):

    def tearDown(self):
        hydra.core.global_hydra.GlobalHydra.instance().clear()
    
    # Test that the motif remains fixed throughout inference.
    # def test_motif_remains_fixed(self):
    #     T = 2
    #     conf = construct_conf([
    #         f'diffuser.T={T}',
    #         'inference.num_designs=1',
    #         'inference.output_prefix=tmp/test_7',
    #         'inference.two_template=True',
    #         'inference.three_template=True',
    #         'inference.motif_only_2d=True',
    #         'model.main_block=1',
    #         'inference.align_px0_motif=False',
    #         "preprocess.eye_frames=True",
    #         "inference.ij_visible=a",
    #         'inference.supply_motif_seq=False'
    #     ])

    #     func_sig = signature(RoseTTAFoldModule.forward)
    #     fake_forward = mock.patch.object(RoseTTAFoldModule, "forward", autospec=True)

    #     def side_effect(self, *args, **kwargs):
    #         ic("mock forward", type(self), side_effect.call_count)
    #         side_effect.call_count += 1
    #         return fake_forward.temp_original(self, *args, **kwargs)
    #     side_effect.call_count = 0

    #     with fake_forward as mock_forward:
    #         mock_forward.side_effect = side_effect
    #         run_inference.main(conf)

    #         mapped_calls = []
    #         for args, kwargs in mock_forward.call_args_list:
    #             args = (None,) + args[1:]
    #             argument_binding = func_sig.bind(*args, **kwargs)
    #             argument_map = argument_binding.arguments
    #             argument_map = tensor_util.cpu(argument_map)
    #             mapped_calls.append(argument_map)
        
    #     is_motif = 1
    #     def constant(mapped_call):
    #         c = {}
    #         c['xyz'] = mapped_call['xyz'][0,is_motif]
    #         is_sidechain_torsion = torch.ones(3*rf2aa.chemical.NTOTALDOFS).bool()
    #         is_sidechain_torsion[0:2] = False
    #         is_sidechain_torsion[3:5] = False
    #         c['alpha'] = mapped_call['alpha_t'][0,0,is_motif]
    #         # Remove backbone torsions
    #         c['alpha'][~is_sidechain_torsion] = torch.nan
    #         c['sctors'] = mapped_call['sctors'][0, is_motif]
    #         # Remove backbone torsions
    #         c['sctors'][0:2] = torch.nan
    #         return c
        
    #     constants = []
    #     for mapped_call in mapped_calls:
    #         constants.append(constant(mapped_call))
    #     #     ipdb.set_trace()
    #     ipdb.set_trace()
    #     self.assertEqual(len(constants), T)
    #     cmp = partial(tensor_util.cmp, atol=1e-9, rtol=1e-4)
    #     for i in range(1, T):
    #         test_utils.assertEqual(self, cmp, constants[0], constants[i])
    
    # def test_motif_fixed_in_output(self):
    #     output_pdb, conf = infer([
    #         f'diffuser.T={T}',
    #         'inference.num_designs=1',
    #         'inference.output_prefix=tmp/test_7',
    #         'inference.two_template=True',
    #         'inference.three_template=True',
    #         'inference.motif_only_2d=True',
    #         'model.main_block=1',
    #         'inference.align_px0_motif=False',
    #         "preprocess.eye_frames=True",
    #         "inference.ij_visible=a",
    #         'inference.supply_motif_seq=False'
    #     ])

    #     input_feats = inference.utils.parse_pdb(conf.inference.input_pdb)
    #     output_feats = inference.utils.parse_pdb(output_pdb)

    #     trb = get_trb(conf)
    #     is_motif = torch.tensor(trb['con_hal_idx0'])
    #     is_motif_ref = torch.tensor(trb['con_ref_idx0'])
    #     n_motif = len(is_motif)
        

    #     input_motif_xyz = input_feats['xyz'][is_motif_ref]
    #     output_motif_xyz = output_feats['xyz'][is_motif]
    #     atom_mask = input_feats['mask'][is_motif]
    #     self.assertEqual(n_motif, 2)

    #     # Backbone only
    #     backbone_atom_mask = torch.zeros((n_motif, 14)).bool()
    #     backbone_atom_mask[:,:3] = True
    #     backbone_rmsd = rf2aa.loss.calc_crd_rmsd(
    #             torch.tensor(input_motif_xyz)[None],
    #             torch.tensor(output_motif_xyz)[None],
    #             backbone_atom_mask[None])
    #     # The motif gets rotated and translated, so the accuracy is somewhat limited
    #     # due to the precision of coordinates in a PDB file.
    #     self.assertLess(backbone_rmsd, 0.04) #TODO verify 0.04 is a reasonable threshold

    #     # All atoms
    #     rmsd = rf2aa.loss.calc_crd_rmsd(
    #             torch.tensor(input_motif_xyz)[None],
    #             torch.tensor(output_motif_xyz)[None],
    #             torch.tensor(atom_mask)[None])
    #     self.assertLess(rmsd, 0.02)

#############
#
#############

    def test_diffuse_bb_wo_ligand_wo_motif(self):
        T = 2
        L = 4
        output_pdb, conf = infer([
            f'diffuser.T={T}',
            'inference.num_designs=1',
            'inference.output_prefix=tmp/test_diffuse_bb_wo_ligand_wo_motif',
            'inference.two_template=True',
            'inference.three_template=True',
            'inference.motif_only_2d=True',
            'model.main_block=1',
            'inference.align_px0_motif=False',
            "preprocess.eye_frames=True",
            "inference.ij_visible=a",
            'inference.supply_motif_seq=False',
            f"contigmap.contigs=['{L}-{L}']"
        ])
        print(output_pdb)
        input_feats = inference.utils.parse_pdb(conf.inference.input_pdb,parse_hetatom=True)
        output_feats = inference.utils.parse_pdb(output_pdb, parse_hetatom=True)

        trb = get_trb(conf)
        is_motif = torch.tensor(trb['con_hal_idx0'])
        is_motif_ref = torch.tensor(trb['con_ref_idx0'])
        n_motif = len(is_motif)
        cmp = partial(tensor_util.cmp, atol=1e-9, rtol=1e-4)
        #test length of output
        diff_len  = len(output_feats['seq'])# - len(input_feats['xyz'])
        test_utils.assertEqual(self,cmp, diff_len, L)

    def test_diffuse_bb_wo_ligand_wo_motif_symmetry(self):
        T = 3
        L = 40 #10xc4
        c_sym = 'c4'
        output_pdb, conf = infer([
            f'diffuser.T={T}',
            'inference.num_designs=1',
            'inference.output_prefix=tmp/test_diffuse_bb_wo_ligand_wo_motif_symmetry',
            'inference.two_template=True',
            'inference.three_template=True',
            'inference.motif_only_2d=True',
            'model.main_block=1',
            'inference.align_px0_motif=False',
            "preprocess.eye_frames=True",
            "inference.ij_visible=a",
            'inference.supply_motif_seq=False',
            f"contigmap.contigs=['{L}-{L}']",
            f"inference.pseudo_symmetry={c_sym}"
        ])
            #
        input_feats = inference.utils.parse_pdb(conf.inference.input_pdb)#,parse_hetatom=True)
        output_feats = inference.utils.parse_pdb(output_pdb)#, parse_hetatom=True)

        trb = get_trb(conf)
        is_motif = torch.tensor(trb['con_hal_idx0'])
        is_motif_ref = torch.tensor(trb['con_ref_idx0'])
        n_motif = len(is_motif)
        cmp = partial(tensor_util.cmp, atol=1e-9, rtol=1e-4)

        diff_len  = len(output_feats['seq'])
        test_utils.assertEqual(self,cmp, diff_len, L)
        
        n_sym = int(c_sym[1:])
        unit_len = diff_len//n_sym
        ref_unit_xyz = torch.tensor(output_feats['xyz'][:unit_len,:1]) #only Ca
        ref_unit_xyz = ref_unit_xyz.reshape(unit_len,3)
        #ipdb.set_trace()

        for i in range(1,n_sym):
            unit_xyz = torch.tensor(output_feats['xyz'][i*unit_len:(i+1)*unit_len,:1])
            unit_xyz = unit_xyz.reshape(unit_len,3)
            unit_rmsd = float(rf2aa_util.kabsch(ref_unit_xyz, unit_xyz)[0])
            #ipdb.set_trace()
            self.assertLess(unit_rmsd, 0.1) #0.1 is high but ok for now, kabsh is not perfect anw

    def test_diffuse_bb_wo_ligand_wo_motif_repeat(self):
        T = 3
        L = 40 #10xc4
        n_repeat = 4#c_sym = 'c3'
        output_pdb, conf = infer([
            f'diffuser.T={T}',
            'inference.num_designs=1',
            'inference.output_prefix=tmp/test_diffuse_bb_wo_ligand_wo_motif_repeat',
            'inference.two_template=True',
            'inference.three_template=True',
            'inference.motif_only_2d=True',
            'model.main_block=1',
            'inference.align_px0_motif=False',
            "preprocess.eye_frames=True",
            "inference.ij_visible=a",
            'inference.supply_motif_seq=False',
            f"contigmap.contigs=['9,A501-501,9,A501-501,9,A501-501,9,A501-501']",
            f"inference.n_repeats={n_repeat}"
        ])
            #
        input_feats = inference.utils.parse_pdb(conf.inference.input_pdb)#,parse_hetatom=True)
        output_feats = inference.utils.parse_pdb(output_pdb)#, parse_hetatom=True)

        trb = get_trb(conf)
        is_motif = torch.tensor(trb['con_hal_idx0'])
        is_motif_ref = torch.tensor(trb['con_ref_idx0'])
        n_motif = len(is_motif)
        cmp = partial(tensor_util.cmp, atol=1e-9, rtol=1e-4)
        ipdb.set_trace()
        diff_len  = len(output_feats['seq'])
        test_utils.assertEqual(self,cmp, diff_len, L) #not sure what else to test here, seq is all A's

    def test_diffuse_bb_wo_ligand_with_motif_repeat(self):
        T = 3
        L = 40 #10xc4
        n_repeat = 4#c_sym = 'c3'
        output_pdb, conf = infer([
            f'diffuser.T={T}',
            'inference.num_designs=1',
            'inference.output_prefix=tmp/test_diffuse_bb_wo_ligand_with_motif_repeat',
            'inference.two_template=True',
            'inference.three_template=True',
            'inference.motif_only_2d=True',
            'model.main_block=1',
            'inference.align_px0_motif=False',
            "preprocess.eye_frames=True",
            "inference.ij_visible=a",
            'inference.supply_motif_seq=True',
            f"contigmap.contigs=['9,A501-501,9,A501-501,9,A501-501,9,A501-501']",
            f"inference.n_repeats={n_repeat}"
        ])
            #
        input_feats = inference.utils.parse_pdb(conf.inference.input_pdb)#,parse_hetatom=True)
        output_feats = inference.utils.parse_pdb(output_pdb)#, parse_hetatom=True)

        trb = get_trb(conf)
        is_motif = torch.tensor(trb['con_hal_idx0'])
        is_motif_ref = torch.tensor(trb['con_ref_idx0'])
        n_motif = len(is_motif)
        cmp = partial(tensor_util.cmp, atol=1e-9, rtol=1e-4)
        ipdb.set_trace()
        diff_len  = len(output_feats['seq'])
        test_utils.assertEqual(self,cmp, diff_len, L)
        test_utils.assertEqual(self,cmp, output_feats['seq'][is_motif], input_feats['seq'][is_motif_ref])

    def test_diffuse_bb_with_ligand(self):
        T = 3
        L = 40 #10xc4
        output_pdb, conf = infer([
            f'diffuser.T={T}',
            'inference.num_designs=1',
            'inference.output_prefix=tmp/test_diffuse_bb_with_ligand',
            'inference.two_template=True',
            'inference.three_template=True',
            'inference.motif_only_2d=True',
            'model.main_block=1',
            'inference.align_px0_motif=False',
            "preprocess.eye_frames=True",
            "inference.ij_visible=a",
            'inference.supply_motif_seq=True',
            "inference.ligand=LG1"
        ])

        input_feats = inference.utils.parse_pdb(conf.inference.input_pdb,parse_hetatom=True)
        output_feats = inference.utils.parse_pdb(output_pdb,parse_hetatom=True)

        trb = get_trb(conf)
        is_motif = torch.tensor(trb['con_hal_idx0'])
        is_motif_ref = torch.tensor(trb['con_ref_idx0'])
        n_motif = len(is_motif)
        #ipdb.set_trace()
        cmp = partial(tensor_util.cmp, atol=1e-9, rtol=1e-4)

        #test keep specified input motif seq
        #test_utils.assertEqual(self,cmp,""
        assert(n_motif == 1) # per aa_small.yaml
        test_utils.assertEqual(self,cmp, output_feats['seq'][is_motif], input_feats['seq'][is_motif_ref])
        test_utils.assertEqual(self,cmp, len(output_feats['xyz_het']),  len(input_feats['xyz_het']))
        lig_xyz_input = torch.tensor(input_feats['xyz_het'])
        lig_xyz_output = torch.tensor(output_feats['xyz_het'])
        lig_rmsd = float(rf2aa_util.kabsch(lig_xyz_input, lig_xyz_output)[0])
        self.assertLess(lig_rmsd, 0.1) #0.1 is high but ok for now, kabsh is not perfect anw

    def test_diffuse_bb_with_ligand_wo_motif_symmetry(self):
        T = 3
        L = 40 #10xc4
        c_sym = 'c4'
        output_pdb, conf = infer([
            f'diffuser.T={T}',
            'inference.num_designs=1',
            'inference.output_prefix=tmp/test_diffuse_bb_with_ligand_wo_motif_symmetry',
            'inference.two_template=True',
            'inference.three_template=True',
            'inference.motif_only_2d=True',
            'model.main_block=1',
            'inference.align_px0_motif=False',
            "preprocess.eye_frames=True",
            "inference.ij_visible=a",
            'inference.supply_motif_seq=True',
            "inference.ligand=LG1",
            f"contigmap.contigs=['{L}-{L}']",
            f"inference.pseudo_symmetry={c_sym}"
        ])
        input_feats = inference.utils.parse_pdb(conf.inference.input_pdb,parse_hetatom=True)
        output_feats = inference.utils.parse_pdb(output_pdb,parse_hetatom=True)

        trb = get_trb(conf)
        is_motif = torch.tensor(trb['con_hal_idx0'])
        is_motif_ref = torch.tensor(trb['con_ref_idx0'])
        n_motif = len(is_motif)
        #ipdb.set_trace()
        cmp = partial(tensor_util.cmp, atol=1e-9, rtol=1e-4)

        test_utils.assertEqual(self,cmp, output_feats['seq'][is_motif], input_feats['seq'][is_motif_ref])
        test_utils.assertEqual(self,cmp, len(output_feats['xyz_het']),  len(input_feats['xyz_het']))
        lig_xyz_input = torch.tensor(input_feats['xyz_het'])
        lig_xyz_output = torch.tensor(output_feats['xyz_het'])
        lig_rmsd = float(rf2aa_util.kabsch(lig_xyz_input, lig_xyz_output)[0])
        self.assertLess(lig_rmsd, 0.1) #0.1 is high but ok for now, kabsh is not perfect anw

        diff_len  = len(output_feats['xyz'])
        test_utils.assertEqual(self,cmp, diff_len, L) #not sure what else to test here, seq is all A's

        n_sym = int(c_sym[1:])
        unit_len = diff_len//n_sym
        ref_unit_xyz = torch.tensor(output_feats['xyz'][:unit_len,:1]) #only Ca
        ref_unit_xyz = ref_unit_xyz.reshape(unit_len,3)
        #ipdb.set_trace()

        for i in range(1,n_sym):
            unit_xyz = torch.tensor(output_feats['xyz'][i*unit_len:(i+1)*unit_len,:1])
            unit_xyz = unit_xyz.reshape(unit_len,3)
            unit_rmsd = float(rf2aa_util.kabsch(ref_unit_xyz, unit_xyz)[0])
            self.assertLess(unit_rmsd, 0.1) #0.1 is high but ok for now, kabsh is not perfect anw



    def test_diffuse_bb_with_ligand_with_motif_symmetry(self):
        T = 3
        L = 40 #10xc4
        c_sym = 'c4'
        #L = 4 # ['1,A518-518,1'] aa_small.yaml
        n_repeat=4
        output_pdb, conf = infer([
            f'diffuser.T={T}',
            'inference.num_designs=1',
            'inference.output_prefix=tmp/test_diffuse_bb_with_ligand_with_motif_symmetry',
            'inference.two_template=True',
            'inference.three_template=True',
            'inference.motif_only_2d=True',
            'model.main_block=1',
            'inference.align_px0_motif=False',
            "preprocess.eye_frames=True",
            "inference.ij_visible=a",
            'inference.supply_motif_seq=True',
            "inference.ligand=LG1",
            f"contigmap.contigs=['9,A518-518,9,A518-518,9,A518-518,9,A518-518']", #A580=D
            f"inference.n_repeats={n_repeat}",
            f"inference.pseudo_symmetry={c_sym}"
        ])

        input_feats = inference.utils.parse_pdb(conf.inference.input_pdb,parse_hetatom=True)
        output_feats = inference.utils.parse_pdb(output_pdb,parse_hetatom=True)

        trb = get_trb(conf)
        is_motif = torch.tensor(trb['con_hal_idx0'])
        is_motif_ref = torch.tensor(trb['con_ref_idx0'])
        n_motif = len(is_motif)
        #ipdb.set_trace()
        cmp = partial(tensor_util.cmp, atol=1e-9, rtol=1e-4)

        test_utils.assertEqual(self,cmp, output_feats['seq'][is_motif], input_feats['seq'][is_motif_ref])
        test_utils.assertEqual(self,cmp, len(output_feats['xyz_het']),  len(input_feats['xyz_het']))
        lig_xyz_input = torch.tensor(input_feats['xyz_het'])
        lig_xyz_output = torch.tensor(output_feats['xyz_het'])
        lig_rmsd = float(rf2aa_util.kabsch(lig_xyz_input, lig_xyz_output)[0])
        self.assertLess(lig_rmsd, 0.1)
        diff_len  = len(output_feats['xyz'])
        test_utils.assertEqual(self,cmp, diff_len, L) #not sure what else to test here, seq is all A's

        n_sym = int(c_sym[1:])
        unit_len = diff_len//n_sym
        ref_unit_xyz = torch.tensor(output_feats['xyz'][:unit_len,:1]) #only Ca
        ref_unit_xyz = ref_unit_xyz.reshape(unit_len,3)
        #ipdb.set_trace()

        for i in range(1,n_sym):
            unit_xyz = torch.tensor(output_feats['xyz'][i*unit_len:(i+1)*unit_len,:1])
            unit_xyz = unit_xyz.reshape(unit_len,3)
            unit_rmsd = float(rf2aa_util.kabsch(ref_unit_xyz, unit_xyz)[0])
            self.assertLess(unit_rmsd, 0.1) #0.1 is high but ok for now, kabsh is not perfect anw

    def test_deteriministic(self): #complicated/obscure diffusion contig test
        run_inference.make_deterministic()
        T = 3
        L = 40 #10xc4
        c_sym = 'c4'
        #L = 4 # ['1,A518-518,1'] aa_small.yaml
        n_repeat=4
        
        output_pdb, conf = infer([
            f'diffuser.T={T}',
            'inference.num_designs=1',
            'inference.output_prefix=tmp/test_diffuse_bb_with_ligand_with_motif_symmetry',
            'inference.two_template=True',
            'inference.three_template=True',
            'inference.motif_only_2d=True',
            'model.main_block=1',
            'inference.align_px0_motif=False',
            "preprocess.eye_frames=True",
            "inference.ij_visible=a",
            'inference.supply_motif_seq=True',
            "inference.ligand=LG1",
            f"contigmap.contigs=['9,A518-518,9,A518-518,9,A518-518,9,A518-518']", #A580=D
            f"inference.n_repeats={n_repeat}",
            f"inference.pseudo_symmetry={c_sym}"
        ])

        input_feats = inference.utils.parse_pdb(conf.inference.input_pdb,parse_hetatom=True)
        output_feats = inference.utils.parse_pdb(output_pdb,parse_hetatom=True)

        trb = get_trb(conf)
        is_motif = torch.tensor(trb['con_hal_idx0'])
        is_motif_ref = torch.tensor(trb['con_ref_idx0'])
        n_motif = len(is_motif)
        #ipdb.set_trace()
        cmp = partial(tensor_util.cmp, atol=5e-2, rtol=0)
        ipdb.set_trace()
        test_utils.assert_matches_golden(self, 'deterministic_diffusion', output_feats, rewrite=True, custom_comparator=cmp)


## refinement movement of backbone + ligand
## deterministic test result in consistent pdb output


if __name__ == '__main__':
        #ipdb.set_trace()
        clean_up()
        unittest.main()
