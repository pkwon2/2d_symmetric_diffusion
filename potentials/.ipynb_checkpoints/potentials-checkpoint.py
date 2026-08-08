import torch
from icecream import ic 
import numpy as np 
from util import generate_Cbeta
from icecream import ic
import rf2aa

class Potential:
    '''
        Interface class that defines the functions a potential must implement
    '''

    def compute(self, seq, xyz):
        '''
            Given the current sequence and structure of the model prediction, return the current
            potential as a PyTorch tensor with a single entry

            Args:
                seq (torch.tensor, size: [L,?]:    The current sequence of the sample.
                                                     TODO: determine whether this is one hot or an 
                                                     integer representation
                xyz (torch.tensor, size: [L,27,3]: The current coordinates of the sample
            
            Returns:
                potential (torch.tensor, size: [1]): A potential whose value will be MAXIMIZED
                                                     by taking a step along it's gradient
        '''
        raise NotImplementedError('Potential compute function was not overwritten')
    

class Zstretch(Potential):
    """
    Potential for stretching points out along the z axis: 
    """

    def __init__(self, weight=1, scale=1.2):
        self.weight = weight
        self.scale = scale

    def compute(self, seq, xyz):
        """
        Computes gradient of how far each point is from being scaled by self.scale. 
        """
        ca = xyz[:,1] # [L,3]
        z = ca[:,-1] 

        z_scaled = z*self.scale
        delta = z_scaled - z

        return torch.sum(delta)*self.weight
    

class funnel(Potential):
    """
    Applies a potential encouraging funnel shape along some axis

    Generic formulation will be P = f(sqrt(x^2+y^2))*g(z)*v, so we are somehow scaling with (a) the distance from the axis and (b) the z coordinate
    
    """

    def __init__(self, weight=1, rad_cut=15, cut_zscale=30, zcut=10):
        self.weight = weight
        self.rad_cut = rad_cut
        self.cut_zscale = cut_zscale
        self.zcut = zcut

    def compute(self, seq, xyz):

        CA = xyz[:,1] # [L,3]

        x,y,z = CA[:,0], CA[:,1], CA[:,2]

        # loss = -torch.abs(z)*torch.sqrt(x**2 + y**2 + 1e-6)
        loss = self.cone_loss(x,y,z)

        return loss.sum()*self.weight
    
    def cone_loss(self, x,y,z):
        """
        """
        cutoff = self.rad_cut
        radius = torch.sqrt(x**2 + y**2 + z**2)

        loss = torch.zeros_like(x)

        mask = z <= self.zcut 
        ic(mask.shape)
        ic(mask.sum())
        loss[mask]  = -z[mask]*torch.sqrt(x[mask]**2 + y[mask]**2 + 1e-6)
        loss[~mask] = -z[~mask]*self.cut_zscale 

        ic(torch.norm(loss[mask]))
        ic(torch.norm(loss[~mask]))

        # negative sign is important here
        # due to mask definition and loss definition
        return -loss 


class monomer_ROG(Potential):
    '''
        Radius of Gyration potential for encouraging monomer compactness

        Written by DJ and refactored into a class by NRB
    '''

    def __init__(self, weight=1, min_dist=15):

        self.weight   = weight
        self.min_dist = min_dist

    def compute(self, seq, xyz):
        Ca = xyz[:,1] # [L,3]

        centroid = torch.mean(Ca, dim=0, keepdim=True) # [1,3]

        dgram = torch.cdist(Ca[None,...].contiguous(), centroid[None,...].contiguous(), p=2) # [1,L,1,3]

        dgram = torch.maximum(self.min_dist * torch.ones_like(dgram.squeeze(0)), dgram.squeeze(0)) # [L,1,3]

        rad_of_gyration = torch.sqrt( torch.sum(torch.square(dgram)) / Ca.shape[0] ) # [1]

        return -1 * self.weight * rad_of_gyration

class binder_ROG(Potential):
    '''
        Radius of Gyration potential for encouraging binder compactness

        Author: NRB
    '''

    def __init__(self, binderlen, weight=1, min_dist=15):

        self.binderlen = binderlen
        self.min_dist  = min_dist
        self.weight    = weight

    def compute(self, seq, xyz):
        
        # Only look at binder residues
        Ca = xyz[:self.binderlen,1] # [Lb,3]

        centroid = torch.mean(Ca, dim=0, keepdim=True) # [1,3]

        # cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca[None,...].contiguous(), centroid[None,...].contiguous(), p=2) # [1,Lb,1,3]

        dgram = torch.maximum(self.min_dist * torch.ones_like(dgram.squeeze(0)), dgram.squeeze(0)) # [Lb,1,3]

        rad_of_gyration = torch.sqrt( torch.sum(torch.square(dgram)) / Ca.shape[0] ) # [1]

        return -1 * self.weight * rad_of_gyration


class dimer_ROG(Potential):
    '''
        Radius of Gyration potential for encouraging compactness of both monomers when designing dimers

        Author: PV
    '''

    def __init__(self, binderlen, weight=1, min_dist=15):

        self.binderlen = binderlen
        self.min_dist  = min_dist
        self.weight    = weight

    def compute(self, seq, xyz):

        # Only look at monomer 1 residues
        Ca_m1 = xyz[:self.binderlen,1] # [Lb,3]
        
        # Only look at monomer 2 residues
        Ca_m2 = xyz[self.binderlen:,1] # [Lb,3]

        centroid_m1 = torch.mean(Ca_m1, dim=0, keepdim=True) # [1,3]
        centroid_m2 = torch.mean(Ca_m1, dim=0, keepdim=True) # [1,3]

        # cdist needs a batch dimension - NRB
        #This calculates RoG for Monomer 1
        dgram_m1 = torch.cdist(Ca_m1[None,...].contiguous(), centroid_m1[None,...].contiguous(), p=2) # [1,Lb,1,3]
        dgram_m1 = torch.maximum(self.min_dist * torch.ones_like(dgram_m1.squeeze(0)), dgram_m1.squeeze(0)) # [Lb,1,3]
        rad_of_gyration_m1 = torch.sqrt( torch.sum(torch.square(dgram_m1)) / Ca_m1.shape[0] ) # [1]

        # cdist needs a batch dimension - NRB
        #This calculates RoG for Monomer 2
        dgram_m2 = torch.cdist(Ca_m2[None,...].contiguous(), centroid_m2[None,...].contiguous(), p=2) # [1,Lb,1,3]
        dgram_m2 = torch.maximum(self.min_dist * torch.ones_like(dgram_m2.squeeze(0)), dgram_m2.squeeze(0)) # [Lb,1,3]
        rad_of_gyration_m2 = torch.sqrt( torch.sum(torch.square(dgram_m2)) / Ca_m2.shape[0] ) # [1]

        #Potential value is the average of both radii of gyration (is avg. the best way to do this?)
        return -1 * self.weight * (rad_of_gyration_m1 + rad_of_gyration_m2)/2

class binder_ncontacts(Potential):
    '''
        Differentiable way to maximise number of contacts within a protein
        
        Motivation is given here: https://www.plumed.org/doc-v2.7/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html

        Author: PV
    '''

    def __init__(self, binderlen, weight=1, r_0=8, d_0=4):

        self.binderlen = binderlen
        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0

    def compute(self, seq, xyz):

        # Only look at binder Ca residues
        Ca = xyz[:self.binderlen,1] # [Lb,3]
        
        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca[None,...].contiguous(), Ca[None,...].contiguous(), p=2) # [1,Lb,Lb]
        divide_by_r_0 = (dgram - self.d_0) / self.r_0
        numerator = torch.pow(divide_by_r_0,6)
        denominator = torch.pow(divide_by_r_0,12)
        binder_ncontacts = (1 - numerator) / (1 - denominator)
        
        print("BINDER CONTACTS:", binder_ncontacts.sum())
        #Potential value is the average of both radii of gyration (is avg. the best way to do this?)
        return self.weight * binder_ncontacts.sum()

    
class dimer_ncontacts(Potential):

    '''
        Differentiable way to maximise number of contacts for two individual monomers in a dimer
        
        Motivation is given here: https://www.plumed.org/doc-v2.7/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html

        Author: PV
    '''


    def __init__(self, binderlen, weight=1, r_0=8, d_0=4):

        self.binderlen = binderlen
        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0

    def compute(self, seq, xyz):

        # Only look at binder Ca residues
        Ca = xyz[:self.binderlen,1] # [Lb,3]
        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca[None,...].contiguous(), Ca[None,...].contiguous(), p=2) # [1,Lb,Lb]
        divide_by_r_0 = (dgram - self.d_0) / self.r_0
        numerator = torch.pow(divide_by_r_0,6)
        denominator = torch.pow(divide_by_r_0,12)
        binder_ncontacts = (1 - numerator) / (1 - denominator)
        #Potential is the sum of values in the tensor
        binder_ncontacts = binder_ncontacts.sum()

        # Only look at target Ca residues
        Ca = xyz[self.binderlen:,1] # [Lb,3]
        dgram = torch.cdist(Ca[None,...].contiguous(), Ca[None,...].contiguous(), p=2) # [1,Lb,Lb]
        divide_by_r_0 = (dgram - self.d_0) / self.r_0
        numerator = torch.pow(divide_by_r_0,6)
        denominator = torch.pow(divide_by_r_0,12)
        target_ncontacts = (1 - numerator) / (1 - denominator)
        #Potential is the sum of values in the tensor
        target_ncontacts = target_ncontacts.sum()
        
        print("DIMER NCONTACTS:", (binder_ncontacts+target_ncontacts)/2)
        #Returns average of n contacts withiin monomer 1 and monomer 2
        return self.weight * (binder_ncontacts+target_ncontacts)/2

class interface_ncontacts(Potential):

    '''
        Differentiable way to maximise number of contacts between binder and target
        
        Motivation is given here: https://www.plumed.org/doc-v2.7/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html

        Author: PV
    '''


    def __init__(self, binderlen, weight=1, r_0=8, d_0=6):

        self.binderlen = binderlen
        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0

    def compute(self, seq, xyz):

        # Extract binder Ca residues
        Ca_b = xyz[:self.binderlen,1] # [Lb,3]

        # Extract target Ca residues
        Ca_t = xyz[self.binderlen:,1] # [Lt,3]

        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca_b[None,...].contiguous(), Ca_t[None,...].contiguous(), p=2) # [1,Lb,Lt]
        divide_by_r_0 = (dgram - self.d_0) / self.r_0
        numerator = torch.pow(divide_by_r_0,6)
        denominator = torch.pow(divide_by_r_0,12)
        interface_ncontacts = (1 - numerator) / (1 - denominator)
        #Potential is the sum of values in the tensor
        interface_ncontacts = interface_ncontacts.sum()

        print("INTERFACE CONTACTS:", interface_ncontacts.sum())

        return self.weight * interface_ncontacts


class avoid_X(Potential):
    """
    Avoids the X axis 
    """

    def __init__(self, weight=1, alpha=1, max_penalty=10):
        self.weight = weight
        self.alpha = alpha
        self.max_penalty = max_penalty

    def compute(self, seq, xyz):
        ca = xyz[:,1] # [L,3]

        # squared distance from x axis
        sq_xdist   = torch.sqrt(ca[:,1]**2 + ca[:,2]**2)
        # penalty decays as 1/d**2, clamp at 10
        penalty = torch.clamp(self.alpha / (sq_xdist), max=self.max_penalty)

        return -self.weight * penalty.sum()



class monomer_contacts(Potential):
    '''
        Differentiable way to maximise number of contacts within a protein

        Motivation is given here: https://www.plumed.org/doc-v2.7/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html
        Author: PV

        NOTE: This function sometimes produces NaN's -- added check in reverse diffusion for nan grads
    '''

    def __init__(self, weight=1, r_0=8, d_0=2, eps=1e-6):

        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0
        self.eps       = eps

    def compute(self, seq, xyz):

        Ca = xyz[:,1] # [L,3]

        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca[None,...].contiguous(), Ca[None,...].contiguous(), p=2) # [1,Lb,Lb]
        divide_by_r_0 = (dgram - self.d_0) / self.r_0
        numerator = torch.pow(divide_by_r_0,6)
        denominator = torch.pow(divide_by_r_0,12)

        ncontacts = (1 - numerator) / ((1 - denominator))


        #Potential value is the average of both radii of gyration (is avg. the best way to do this?)
        return self.weight * ncontacts.sum()


class ligand_ncontacts(Potential):

    '''
        Differentiable way to maximise number of contacts between binder and target

        Motivation is given here: https://www.plumed.org/doc-v2.7/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html

        Author: PV
    '''


    def __init__(self, weight=1, r_0=8, d_0=4):

        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0

    def compute(self, seq, xyz):

        is_atom = rf2aa.util.is_atom(torch.argmax(seq,dim=1)).cpu()

        # Extract ligand Ca residues
        Ca_l = xyz[is_atom,1] # [Ll,3]

        # Extract binder Ca residues
        Ca_b = xyz[~is_atom,1] # [Lb,3]


        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca_b[None,...].contiguous(), Ca_l[None,...].contiguous(), p=2) # [1,Ll,Lb]
        ligand_ncontacts = -1 * contact_energy(dgram, self.r_0, self.d_0)
        #Potential is the sum of values in the tensor
        ligand_ncontacts = ligand_ncontacts.sum()
        print("LIGAND CONTACTS:", ligand_ncontacts.sum())

        return self.weight * ligand_ncontacts

def make_contact_matrix(nchain, contact_string=None):
    """
    Calculate a matrix of inter/intra chain contact indicators
    
    Parameters:
        nchain (int, required): How many chains are in this design 
        
        contact_str (str, optional): String denoting how to define contacts, comma delimited between pairs of chains
            '!' denotes repulsive, '&' denotes attractive
    """
    alphabet   = [a for a in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
    letter2num = {a:i for i,a in enumerate(alphabet)}
    
    contacts   = np.zeros((nchain,nchain))
    written    = np.zeros((nchain,nchain))
    
    contact_list = contact_string.split(',') 
    for c in contact_list:
        if not len(c) == 3:
            raise SyntaxError('Invalid contact(s) specification')

        i,j = letter2num[c[0]],letter2num[c[2]]
        symbol = c[1]
        
        # denote contacting/repulsive
        assert symbol in ['!','&']
        if symbol == '!':
            contacts[i,j] = -1
            contacts[j,i] = -1
        else:
            contacts[i,j] = 1
            contacts[j,i] = 1
            
    return contacts 


class olig_contacts(Potential):
    """
    Applies PV's num contacts potential within/between chains in symmetric oligomers 

    Author: DJ 
    """

    def __init__(self, 
                 contact_matrix, 
                 weight_intra=1, 
                 weight_inter=1,
                 r_0=8, d_0=2):
        """
        Parameters:
            chain_lengths (list, required): List of chain lengths, length is (Nchains)

            contact_matrix (torch.tensor/np.array, required): 
                square matrix of shape (Nchains,Nchains) whose (i,j) enry represents 
                attractive (1), repulsive (-1), or non-existent (0) contact potentials 
                between chains in the complex

            weight (int/float, optional): Scaling/weighting factor
        """
        print('This is chain contact matrix you are using')
        ic(contact_matrix)
        self.contact_matrix = contact_matrix
        self.weight_intra = weight_intra 
        self.weight_inter = weight_inter 
        self.r_0 = r_0
        self.d_0 = d_0

        # check contact matrix only contains valid entries 
        assert all([i in [-1,0,1] for i in contact_matrix.flatten()]), 'Contact matrix must contain only 0, 1, or -1 in entries'
        # assert the matrix is square and symmetric 
        shape = contact_matrix.shape 
        assert len(shape) == 2 
        assert shape[0] == shape[1]
        for i in range(shape[0]):
            for j in range(shape[1]):
                assert contact_matrix[i,j] == contact_matrix[j,i]
        self.nchain=shape[0]

         
    #   self._compute_chain_indices()

    # def _compute_chain_indices(self):
    #     # make list of shape [i,N] for indices of each chain in total length
    #     indices = []
    #     start   = 0
    #     for l in self.chain_lengths:
    #         indices.append(torch.arange(start,start+l))
    #         start += l
    #     self.indices = indices 

    def _get_idx(self,i,L):
        """
        Returns the zero-indexed indices of the residues in chain i
        """
        assert L%self.nchain == 0
        Lchain = L//self.nchain
        return i*Lchain + torch.arange(Lchain)


    def compute(self, seq, xyz):
        """
        Iterate through the contact matrix, compute contact potentials between chains that need it,
        and negate contacts for any 
        """
        L = len(seq.squeeze())

        all_contacts = 0
        start = 0
        for i in range(self.nchain):
            for j in range(self.nchain):
                # only compute for upper triangle, disregard zeros in contact matrix 
                if (i <= j) and (self.contact_matrix[i,j] != 0):

                    # get the indices for these two chains 
                    idx_i = self._get_idx(i,L)
                    idx_j = self._get_idx(j,L)

                    Ca_i = xyz[idx_i,1]  # slice out crds for this chain 
                    Ca_j = xyz[idx_j,1]  # slice out crds for that chain 
                    dgram           = torch.cdist(Ca_i[None,...].contiguous(), Ca_j[None,...].contiguous(), p=2) # [1,Lb,Lb]

                    divide_by_r_0   = (dgram - self.d_0) / self.r_0
                    numerator       = torch.pow(divide_by_r_0,6)
                    denominator     = torch.pow(divide_by_r_0,12)
                    ncontacts       = (1 - numerator) / (1 - denominator)

                    # weight, don't double count intra 
                    scalar = (i==j)*self.weight_intra/2 + (i!=j)*self.weight_inter

                    #                 contacts              attr/repuls          relative weights 
                    all_contacts += ncontacts.sum() * self.contact_matrix[i,j] * scalar 

        return all_contacts 
                    

class olig_intra_contacts(Potential):
    """
    Applies PV's num contacts potential for each chain individually in an oligomer design 

    Author: DJ 
    """

    def __init__(self, chain_lengths, weight=1):
        """
        Parameters:

            chain_lengths (list, required): Ordered list of chain lengths 

            weight (int/float, optional): Scaling/weighting factor
        """
        self.chain_lengths = chain_lengths 
        self.weight = weight 


    def compute(self, seq, xyz):
        """
        Computes intra-chain num contacts potential
        """
        assert sum(self.chain_lengths) == len(seq.squeeze), 'given chain lengths do not match total sequence length'

        all_contacts = 0
        start = 0
        for Lc in self.chain_lengths:
            Ca = xyz[start:start+Lc]  # slice out crds for this chain 
            dgram = torch.cdist(Ca[None,...].contiguous(), Ca[None,...].contiguous(), p=2) # [1,Lb,Lb]
            divide_by_r_0 = (dgram - self.d_0) / self.r_0
            numerator = torch.pow(divide_by_r_0,6)
            denominator = torch.pow(divide_by_r_0,12)
            ncontacts = (1 - numerator) / (1 - denominator)

            # add contacts for this chain to all contacts 
            all_contacts += ncontacts.sum()

            # increment the start to be at the next chain 
            start += Lc 


        return self.weight * all_contacts

def get_damped_lj(r_min, r_lin,p1=6,p2=12):
    
    y_at_r_lin = lj(r_lin, r_min, p1, p2)
    ydot_at_r_lin = lj_grad(r_lin, r_min,p1,p2)
    
    def inner(dgram):
        return (dgram < r_lin) * (ydot_at_r_lin * (dgram - r_lin) + y_at_r_lin) + (dgram >= r_lin) * lj(dgram, r_min, p1, p2)
    return inner

def lj(dgram, r_min,p1=6, p2=12):
    return 4 * ((r_min / (2**(1/p1) * dgram))**p2 - (r_min / (2**(1/p1) * dgram))**p1)

def lj_grad(dgram, r_min,p1=6,p2=12):
    return -p2 * r_min**p1*(r_min**p1-dgram**p1) / (dgram**(p2+1))

def mask_expand(mask, n=1):
    mask_out = mask.clone()
    assert mask.ndim == 1
    for i in torch.where(mask)[0]:
        for j in range(i-n, i+n+1):
            if j >= 0 and j < len(mask):
                mask_out[j] = True
    return mask_out

def contact_energy(dgram, d_0, r_0):
    divide_by_r_0 = (dgram - d_0) / r_0
    numerator = torch.pow(divide_by_r_0,6)
    denominator = torch.pow(divide_by_r_0,12)
    
    ncontacts = (1 - numerator) / ((1 - denominator)).float()
    return - ncontacts

def poly_repulse(dgram, r, slope, p=1):
    a = slope / (p * r**(p-1))

    #ic(a)
    #ic(torch.abs(r - dgram)**p * slope)
    return (dgram < r) * a * torch.abs(r - dgram)**p * slope

#def only_top_n(dgram


class substrate_contacts(Potential):
    '''
    Implicitly models a ligand with an attractive-repulsive potential.
    '''

    def __init__(self, weight=1, r_0=8, d_0=2, s=1, eps=1e-6, rep_r_0=5, rep_s=2, rep_r_min=1):

        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0
        self.eps       = eps
        ic(rep_r_0, rep_s)
        
        # motif frame coordinates
        # NOTE: these probably need to be set after sample_init() call, because the motif sequence position in design must be known
        self.motif_frame = None # [4,3] xyz coordinates from 4 atoms of input motif
        self.motif_mapping = None # list of tuples giving positions of above atoms in design [(resi, atom_idx)]
        self.motif_substrate_atoms = None # xyz coordinates of substrate from input motif
        r_min = 2
        self.energies = []
        self.energies.append(lambda dgram: s * contact_energy(torch.min(dgram, dim=-1)[0], d_0, r_0))
        if rep_r_min:
            self.energies.append(lambda dgram: poly_repulse(torch.min(dgram, dim=-1)[0], rep_r_0, rep_s, p=1.5))
        else:
            self.energies.append(lambda dgram: poly_repulse(dgram, rep_r_0, rep_s, p=1.5))


    def compute(self, seq, xyz):
        
        # First, get random set of atoms
        # This operates on self.xyz_motif, which is assigned to this class in the model runner (for horrible plumbing reasons)
        self._grab_motif_residues(self.xyz_motif)
        
        # for checking affine transformation is corect
        first_distance = torch.sqrt(torch.sqrt(torch.sum(torch.square(self.motif_substrate_atoms[0] - self.motif_frame[0]), dim=-1))) 

        # grab the coordinates of the corresponding atoms in the new frame using mapping
        res = torch.tensor([k[0] for k in self.motif_mapping])
        atoms = torch.tensor([k[1] for k in self.motif_mapping])
        new_frame = xyz[self.diffusion_mask][res,atoms,:]
        # calculate affine transformation matrix and translation vector b/w new frame and motif frame
        A, t = self._recover_affine(self.motif_frame, new_frame)
        # apply affine transformation to substrate atoms
        substrate_atoms = torch.mm(A, self.motif_substrate_atoms.transpose(0,1)).transpose(0,1) + t
        second_distance = torch.sqrt(torch.sqrt(torch.sum(torch.square(new_frame[0] - substrate_atoms[0]), dim=-1)))
        assert abs(first_distance - second_distance) < 0.01, "Alignment seems to be bad" 
        diffusion_mask = mask_expand(self.diffusion_mask, 2)
        Ca = xyz[~diffusion_mask, 1]

        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca[None,...].contiguous(), substrate_atoms.float(), p=2) # [1,Lb,Lb]

        all_energies = []
        for i, energy_fn in enumerate(self.energies):
            energy = energy_fn(dgram)
            ic(i, energy.sum(), energy.min(), energy.max())
            all_energies.append(energy.sum())
        return - self.weight * sum(all_energies)

        #Potential value is the average of both radii of gyration (is avg. the best way to do this?)
        return self.weight * ncontacts.sum()

    def _recover_affine(self,frame1, frame2):
        """
        Uses Simplex Affine Matrix (SAM) formula to recover affine transform between two sets of 4 xyz coordinates
        See: https://www.researchgate.net/publication/332410209_Beginner%27s_guide_to_mapping_simplexes_affinely

        Args: 
        frame1 - 4 coordinates from starting frame [4,3]
        frame2 - 4 coordinates from ending frame [4,3]
        
        Outputs:
        A - affine transformation matrix from frame1->frame2
        t - affine translation vector from frame1->frame2
        """

        l = len(frame1)
        # construct SAM denominator matrix
        B = torch.vstack([frame1.T, torch.ones(l)])
        D = 1.0 / torch.linalg.det(B) # SAM denominator

        M = torch.zeros((3,4), dtype=torch.float64)
        for i, R in enumerate(frame2.T):
            for j in range(l):
                num = torch.vstack([R, B])
                # make SAM numerator matrix
                num = torch.cat((num[:j+1],num[j+2:])) # make numerator matrix
                # calculate SAM entry
                M[i][j] = (-1)**j * D * torch.linalg.det(num)

        A, t = torch.hsplit(M, [l-1])
        t = t.transpose(0,1)
        return A, t

    def _grab_motif_residues(self, xyz) -> None:
        """
        Grabs 4 atoms in the motif.
        Currently random subset of Ca atoms if the motif is >= 4 residues, or else 4 random atoms from a single residue
        """
        idx = torch.arange(self.diffusion_mask.shape[0])
        idx = idx[self.diffusion_mask].float()
        if torch.sum(self.diffusion_mask) >= 4:
            rand_idx = torch.multinomial(idx, 4).long()
            # get Ca atoms
            self.motif_frame = xyz[rand_idx, 1]
            self.motif_mapping = [(i,1) for i in rand_idx]
        else:
            rand_idx = torch.multinomial(idx, 1).long()
            self.motif_frame = xyz[rand_idx[0],:4]
            self.motif_mapping = [(rand_idx, i) for i in range(4)]


class binder_distance_ReLU(Potential):
    '''
        Given the current coordinates of the diffusion trajectory, calculate a potential that is the distance between each residue
        and the closest target residue.

        This potential is meant to encourage the binder to interact with a certain subset of residues on the target that 
        define the binding site.

        Author: NRB
    '''

    def __init__(self, binderlen, hotspot_res, weight=1, min_dist=15, use_Cb=False):

        self.binderlen   = binderlen
        self.hotspot_res = [res + binderlen for res in hotspot_res]
        self.weight      = weight
        self.min_dist    = min_dist
        self.use_Cb      = use_Cb

    def compute(self, seq, xyz):
        binder = xyz[:self.binderlen,:,:] # (Lb,27,3)
        target = xyz[self.hotspot_res,:,:] # (N,27,3)

        if self.use_Cb:
            N  = binder[:,0]
            Ca = binder[:,1]
            C  = binder[:,2]

            Cb = generate_Cbeta(N,Ca,C) # (Lb,3)

            N_t  = target[:,0]
            Ca_t = target[:,1]
            C_t  = target[:,2]

            Cb_t = generate_Cbeta(N_t,Ca_t,C_t) # (N,3)

            dgram = torch.cdist(Cb[None,...], Cb_t[None,...], p=2) # (1,Lb,N)

        else:
            # Use Ca dist for potential

            Ca = binder[:,1] # (Lb,3)

            Ca_t = target[:,1] # (N,3)

            dgram = torch.cdist(Ca[None,...], Ca_t[None,...], p=2) # (1,Lb,N)

        closest_dist = torch.min(dgram.squeeze(0), dim=1)[0] # (Lb)

        # Cap the distance at a minimum value
        min_distance = self.min_dist * torch.ones_like(closest_dist) # (Lb)
        potential    = torch.maximum(min_distance, closest_dist) # (Lb)

        # torch.Tensor.backward() requires the potential to be a single value
        potential    = torch.sum(potential, dim=-1)
        
        return -1 * self.weight * potential

class binder_any_ReLU(Potential):
    '''
        Given the current coordinates of the diffusion trajectory, calculate a potential that is the minimum distance between
        ANY residue and the closest target residue.

        In contrast to binder_distance_ReLU this potential will only penalize a pose if all of the binder residues are outside
        of a certain distance from the target residues.

        Author: NRB
    '''

    def __init__(self, binderlen, hotspot_res, weight=1, min_dist=15, use_Cb=False):

        self.binderlen   = binderlen
        self.hotspot_res = [res + binderlen for res in hotspot_res]
        self.weight      = weight
        self.min_dist    = min_dist
        self.use_Cb      = use_Cb

    def compute(self, seq, xyz):
        binder = xyz[:self.binderlen,:,:] # (Lb,27,3)
        target = xyz[self.hotspot_res,:,:] # (N,27,3)

        if use_Cb:
            N  = binder[:,0]
            Ca = binder[:,1]
            C  = binder[:,2]

            Cb = generate_Cbeta(N,Ca,C) # (Lb,3)

            N_t  = target[:,0]
            Ca_t = target[:,1]
            C_t  = target[:,2]

            Cb_t = generate_Cbeta(N_t,Ca_t,C_t) # (N,3)

            dgram = torch.cdist(Cb[None,...], Cb_t[None,...], p=2) # (1,Lb,N)

        else:
            # Use Ca dist for potential

            Ca = binder[:,1] # (Lb,3)

            Ca_t = target[:,1] # (N,3)

            dgram = torch.cdist(Ca[None,...], Ca_t[None,...], p=2) # (1,Lb,N)


        closest_dist = torch.min(dgram.squeeze(0)) # (1)

        potential    = torch.maximum(min_dist, closest_dist) # (1)

        return -1 * self.weight * potential

# Dictionary of types of potentials indexed by name of potential. Used by PotentialManager.

def _kabsch_rigid_transform(X, Y, eps=1e-8):
    """
    Fit the proper rigid-body transformation

        Y ~= X @ R.T + t

    using the differentiable Kabsch/SVD algorithm.

    Parameters
    ----------
    X, Y : torch.Tensor, shape [N, 3]
        Corresponding coordinates.
    eps : float
        Numerical stability value.

    Returns
    -------
    R : torch.Tensor, shape [3, 3]
        Proper rotation matrix.
    t : torch.Tensor, shape [3]
        Translation vector.
    rmsd : torch.Tensor, scalar
        RMSD after superposition.
    """
    if X.ndim != 2 or Y.ndim != 2 or X.shape != Y.shape:
        raise ValueError(
            f"X and Y must have the same [N, 3] shape; "
            f"received {tuple(X.shape)} and {tuple(Y.shape)}"
        )

    if X.shape[-1] != 3:
        raise ValueError(f"Expected [N, 3] coordinates; got {tuple(X.shape)}")

    X_center = X.mean(dim=0)
    Y_center = Y.mean(dim=0)

    X0 = X - X_center
    Y0 = Y - Y_center

    # Covariance for the column-vector convention:
    #
    #     y = R @ x + t
    #
    covariance = X0.transpose(0, 1) @ Y0

    U, S, Vh = torch.linalg.svd(covariance, full_matrices=False)

    # Without this correction, SVD can return an improper rotation
    # containing a reflection.
    det_value = torch.det(Vh.transpose(0, 1) @ U.transpose(0, 1))

    correction = torch.ones(
        3,
        dtype=X.dtype,
        device=X.device,
    )
    correction[-1] = torch.where(
        det_value < 0,
        -torch.ones_like(det_value),
        torch.ones_like(det_value),
    )

    D = torch.diag(correction)

    R = Vh.transpose(0, 1) @ D @ U.transpose(0, 1)
    t = Y_center - R @ X_center

    X_aligned = X @ R.transpose(0, 1) + t
    rmsd = torch.sqrt(
        torch.mean(torch.sum((X_aligned - Y) ** 2, dim=-1)) + eps
    )

    return R, t, rmsd


def _axis_from_rotation_matrix(R, eps=1e-8):
    """
    Extract the unoriented rotation axis from a 3x3 rotation matrix.

    For a nonzero rotation, the antisymmetric part of R is proportional to

        sin(theta) * axis

    This is appropriate for C3 rotations because theta is approximately
    120 degrees and therefore sin(theta) is far from zero.

    Parameters
    ----------
    R : torch.Tensor, shape [3, 3]

    Returns
    -------
    axis : torch.Tensor, shape [3]
        Unit vector parallel to the rotation axis.
    """
    axis = torch.stack(
        [
            R[2, 1] - R[1, 2],
            R[0, 2] - R[2, 0],
            R[1, 0] - R[0, 1],
        ]
    )

    return axis / torch.clamp(torch.linalg.norm(axis), min=eps)


def _closest_point_on_rotation_axis(R, t, axis, eps=1e-8):
    """
    Find one point on the rigid transformation's rotation axis.

    A point p on the axis obeys

        R @ p + t = p

    and therefore

        (I - R) @ p = t.

    This equation leaves motion along the axis unconstrained. We add

        axis dot p = 0

    to select the point on the axis closest to the coordinate origin.

    Parameters
    ----------
    R : torch.Tensor, shape [3, 3]
    t : torch.Tensor, shape [3]
    axis : torch.Tensor, shape [3]

    Returns
    -------
    point : torch.Tensor, shape [3]
        A point lying on the fitted rotation-axis line.
    """
    identity = torch.eye(3, dtype=R.dtype, device=R.device)

    A = torch.cat(
        [
            identity - R,
            axis.reshape(1, 3),
        ],
        dim=0,
    )

    b = torch.cat(
        [
            t,
            torch.zeros(1, dtype=t.dtype, device=t.device),
        ],
        dim=0,
    )

    # pinv uses an SVD-based least-squares solution and remains
    # differentiable with respect to the coordinates.
    point = torch.linalg.pinv(A, rtol=eps) @ b

    return point


def _fit_c3_axis(chain_1, chain_2, chain_3, eps=1e-8):
    """
    Fit the C3 symmetry axis of three corresponding protein chains.

    The function independently fits the transformations

        chain_1 -> chain_2
        chain_2 -> chain_3
        chain_3 -> chain_1

    and combines their inferred rotation axes and axis points.

    Parameters
    ----------
    chain_1, chain_2, chain_3 : torch.Tensor, shape [Lchain, 3]
        Corresponding C-alpha coordinates. Residue i in each chain must
        represent the corresponding symmetry-related residue.
    eps : float

    Returns
    -------
    axis_point : torch.Tensor, shape [3]
        Point on the fitted C3 axis.
    axis_direction : torch.Tensor, shape [3]
        Unit vector parallel to the fitted C3 axis.
    fit_rmsd : torch.Tensor, scalar
        Mean cyclic rigid-body fitting RMSD.
    c3_angle_error : torch.Tensor, scalar
        Error in the fitted rotation angles relative to a C3 rotation.
    """
    transformations = []

    for source, target in (
        (chain_1, chain_2),
        (chain_2, chain_3),
        (chain_3, chain_1),
    ):
        R, t, rmsd = _kabsch_rigid_transform(source, target, eps=eps)
        axis = _axis_from_rotation_matrix(R, eps=eps)
        transformations.append((R, t, axis, rmsd))

    # A symmetry axis is an unoriented line, so +axis and -axis are
    # equivalent. Align all signs to the first estimate before averaging.
    reference_axis = transformations[0][2]

    aligned_axes = []
    axis_points = []
    rmsds = []
    angle_errors = []

    target_cosine = chain_1.new_tensor(-0.5)  # cos(120 degrees)

    for R, t, axis, rmsd in transformations:
        sign = torch.where(
            torch.dot(axis, reference_axis) < 0,
            -torch.ones((), dtype=axis.dtype, device=axis.device),
            torch.ones((), dtype=axis.dtype, device=axis.device),
        )

        aligned_axis = sign * axis
        aligned_axes.append(aligned_axis)

        axis_point = _closest_point_on_rotation_axis(
            R,
            t,
            aligned_axis,
            eps=eps,
        )
        axis_points.append(axis_point)
        rmsds.append(rmsd)

        # For a proper rotation:
        #
        #     trace(R) = 1 + 2*cos(theta)
        #
        fitted_cosine = (torch.trace(R) - 1.0) / 2.0
        fitted_cosine = torch.clamp(fitted_cosine, -1.0, 1.0)

        angle_errors.append((fitted_cosine - target_cosine) ** 2)

    axis_direction = torch.stack(aligned_axes).mean(dim=0)
    axis_direction = axis_direction / torch.clamp(
        torch.linalg.norm(axis_direction),
        min=eps,
    )

    # Each fitted point is defined as the closest point on its estimated
    # axis to the global coordinate origin. Averaging these estimates gives
    # a stable point on the consensus axis.
    axis_point = torch.stack(axis_points).mean(dim=0)

    # Project the averaged point into the plane perpendicular to the final
    # consensus direction, maintaining the closest-to-origin convention.
    axis_point = axis_point - torch.dot(
        axis_point,
        axis_direction,
    ) * axis_direction

    fit_rmsd = torch.stack(rmsds).mean()
    c3_angle_error = torch.stack(angle_errors).mean()

    return axis_point, axis_direction, fit_rmsd, c3_angle_error


def _squared_distance_between_lines(
    point_1,
    direction_1,
    point_2,
    direction_2,
    eps=1e-8,
):
    """
    Differentiable squared shortest distance between two infinite 3D lines.

    Line 1:
        point_1 + s * direction_1

    Line 2:
        point_2 + t * direction_2

    This formulation also remains well behaved when the lines are nearly
    parallel.
    """
    direction_1 = direction_1 / torch.clamp(
        torch.linalg.norm(direction_1),
        min=eps,
    )
    direction_2 = direction_2 / torch.clamp(
        torch.linalg.norm(direction_2),
        min=eps,
    )

    delta = point_2 - point_1

    a = torch.dot(direction_1, direction_1)
    b = torch.dot(direction_1, direction_2)
    c = torch.dot(direction_2, direction_2)
    d = torch.dot(direction_1, delta)
    e = torch.dot(direction_2, delta)

    denominator = a * c - b * b

    safe_denominator = torch.clamp(denominator, min=eps)

    s = (c * d - b * e) / safe_denominator
    t = (b * d - a * e) / safe_denominator

    closest_1 = point_1 + s * direction_1
    closest_2 = point_2 + t * direction_2

    nonparallel_distance_sq = torch.sum((closest_1 - closest_2) ** 2)

    # Parallel-line fallback: remove the component of delta along line 1.
    perpendicular_delta = delta - torch.dot(
        delta,
        direction_1,
    ) * direction_1

    parallel_distance_sq = torch.sum(perpendicular_delta ** 2)

    return torch.where(
        denominator > eps,
        nonparallel_distance_sq,
        parallel_distance_sq,
    )


def _line_intersection_center(
    point_1,
    direction_1,
    point_2,
    direction_2,
    eps=1e-8,
):
    """
    Return the midpoint of the closest points on two lines.

    When the axes intersect exactly, this midpoint is their intersection.
    """
    direction_1 = direction_1 / torch.clamp(
        torch.linalg.norm(direction_1),
        min=eps,
    )
    direction_2 = direction_2 / torch.clamp(
        torch.linalg.norm(direction_2),
        min=eps,
    )

    delta = point_2 - point_1

    a = torch.dot(direction_1, direction_1)
    b = torch.dot(direction_1, direction_2)
    c = torch.dot(direction_2, direction_2)
    d = torch.dot(direction_1, delta)
    e = torch.dot(direction_2, delta)

    denominator = torch.clamp(a * c - b * b, min=eps)

    s = (c * d - b * e) / denominator
    t = (b * d - a * e) / denominator

    closest_1 = point_1 + s * direction_1
    closest_2 = point_2 + t * direction_2

    return 0.5 * (closest_1 + closest_2)


class icosahedral_c3_axes_svd(Potential):
    """
    Encourage two fitted C3 symmetry axes to form an icosahedral arrangement.

    Chain organization
    ------------------
    The input is assumed to contain six equal-length chains in this order:

        A, B, C, D, E, F

    Chains A/B/C form the first C3 component.
    Chains D/E/F form the second C3 component.

    The residue ordering must correspond among symmetry-related chains.
    For example, residue i in A, B, and C must represent the same position
    in the repeated subunit.

    Objective
    ---------
    The potential encourages:

    1. A/B/C to be related by approximately 120-degree rotations.
    2. D/E/F to be related by approximately 120-degree rotations.
    3. The two fitted C3 axes to intersect.
    4. The acute angle between the axes to equal the adjacent-face
       C3-C3 angle of an icosahedron:

           cos(theta) = sqrt(5) / 3
           theta ~= 41.8103 degrees

    Because RFdiffusion maximizes potentials, this class returns the
    negative weighted error.
    """

    def __init__(
        self,
        weight=1.0,
        intersection_weight=1.0,
        axis_angle_weight=10.0,
        c3_fit_weight=0.1,
        c3_rotation_weight=1.0,
        center_weight=0.0,
        eps=1e-8,
    ):
        self.weight = weight
        self.intersection_weight = intersection_weight
        self.axis_angle_weight = axis_angle_weight
        self.c3_fit_weight = c3_fit_weight
        self.c3_rotation_weight = c3_rotation_weight
        self.center_weight = center_weight
        self.eps = eps

    def compute(self, seq, xyz):
        """
        Parameters
        ----------
        seq : torch.Tensor
            Unused, but retained for compatibility with Potential.
        xyz : torch.Tensor, shape [L, 27, 3]
            Current all-atom coordinates.

        Returns
        -------
        potential : torch.Tensor, scalar
            Negative total geometric error.
        """
        if xyz.ndim != 3 or xyz.shape[-1] != 3:
            raise ValueError(
                "icosahedral_c3_axes_svd expects xyz with shape "
                f"[L, Natoms, 3]; received {tuple(xyz.shape)}"
            )

        total_length = xyz.shape[0]

        if total_length % 6 != 0:
            raise ValueError(
                "icosahedral_c3_axes_svd requires six equal-length chains. "
                f"Total residue count {total_length} is not divisible by 6."
            )

        chain_length = total_length // 6

        # C-alpha coordinates, shape [6, Lchain, 3].
        ca = xyz[:, 1, :]
        chains = ca.reshape(6, chain_length, 3)

        axis_point_abc, axis_abc, fit_abc, rotation_error_abc = _fit_c3_axis(
            chains[0],
            chains[1],
            chains[2],
            eps=self.eps,
        )

        axis_point_def, axis_def, fit_def, rotation_error_def = _fit_c3_axis(
            chains[3],
            chains[4],
            chains[5],
            eps=self.eps,
        )

        # Axes are unoriented lines. Taking abs(dot) treats +u and -u as
        # the same symmetry axis.
        axis_cosine = torch.abs(torch.dot(axis_abc, axis_def))

        target_axis_cosine = torch.sqrt(
            xyz.new_tensor(5.0)
        ) / xyz.new_tensor(3.0)

        axis_angle_error = (
            axis_cosine - target_axis_cosine
        ) ** 2

        intersection_error = _squared_distance_between_lines(
            axis_point_abc,
            axis_abc,
            axis_point_def,
            axis_def,
            eps=self.eps,
        )

        c3_fit_error = fit_abc ** 2 + fit_def ** 2

        c3_rotation_error = (
            rotation_error_abc + rotation_error_def
        )

        # Optional restraint placing the shared axis intersection at the
        # overall assembly centroid. This does not force the assembly to
        # the global coordinate origin.
        assembly_center = ca.mean(dim=0)

        fitted_intersection = _line_intersection_center(
            axis_point_abc,
            axis_abc,
            axis_point_def,
            axis_def,
            eps=self.eps,
        )

        center_error = torch.sum(
            (fitted_intersection - assembly_center) ** 2
        )

        total_error = (
            self.intersection_weight * intersection_error
            + self.axis_angle_weight * axis_angle_error
            + self.c3_fit_weight * c3_fit_error
            + self.c3_rotation_weight * c3_rotation_error
            + self.center_weight * center_error
        )

        return -self.weight * total_error

# If you implement a new potential you must add it to this dictionary for it to be used by
# the PotentialManager
implemented_potentials = { 'monomer_ROG':          monomer_ROG,
                           'binder_ROG':           binder_ROG,
                           'binder_distance_ReLU': binder_distance_ReLU,
                           'binder_any_ReLU':      binder_any_ReLU,
                           'dimer_ROG':            dimer_ROG,
                           'binder_ncontacts':     binder_ncontacts,
                           'dimer_ncontacts':      dimer_ncontacts,
                           'interface_ncontacts':  interface_ncontacts,
                           'monomer_contacts':     monomer_contacts,
                           'olig_intra_contacts':  olig_intra_contacts,
                           'olig_contacts':        olig_contacts,
                           'icosahedral_c3_axes_svd': icosahedral_c3_axes_svd,
                           
                           'substrate_contacts':   substrate_contacts,
                           'ligand_ncontacts':     ligand_ncontacts,
                           'avoid_X':              avoid_X,
                           'funnel':               funnel,
                           'Zstretch':             Zstretch,}

require_binderlen      = { 'binder_ROG',
                           'binder_distance_ReLU',
                           'binder_any_ReLU',
                           'dimer_ROG',
                           'binder_ncontacts',
                           'dimer_ncontacts',
                           'interface_ncontacts'}

require_hotspot_res    = { 'binder_distance_ReLU',
                           'binder_any_ReLU' }

