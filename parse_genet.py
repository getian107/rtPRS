#!/usr/bin/env python

"""
Parse the reference panel, validation frequency file, posterior effect estimates and the validation and testing datasets.

"""


import os
import numpy as np
from numpy import linalg
import h5py


def parse_ref(ref_file, chrom):
    print('... parse reference file: %s ...' % ref_file)

    ref_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'MAF':[]}
    with open(ref_file) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom:
                ref_dict['CHR'].append(chrom)
                ref_dict['SNP'].append(ll[1])
                ref_dict['BP'].append(int(ll[2]))
                ref_dict['A1'].append(ll[3])
                ref_dict['A2'].append(ll[4])
                ref_dict['MAF'].append(float(ll[5]))

    print('... %d SNPs on chromosome %d read from %s ...' % (len(ref_dict['SNP']), chrom, ref_file))
    return ref_dict


def parse_vld(vld_prefix, vld_phn_cov):
    print('... parse genotypes for the validation dataset: %s ...' % (vld_prefix + '.bed/.bim/.fam'))


    print('... parse phenotypes and covariates for the validation dataset: %s ...' % vld_phn_cov)

    vld_dict = {'ID':[], 'PHN':[], 'COV':[]}
    with open(vld_phn_cov) as ff:
        for line in ff:
            ll = (line.strip()).split()
            vld_dict['ID'].append(ll[1])
            vld_dict['PHN'].append(float(ll[2]))
            vld_dict['COV'].append([float(ii) for ii in ll[3:]])

    phn = np.array(vld_dict['PHN'], ndmin=2).T
    cov = np.column_stack((np.ones(len(phn)), np.array(vld_dict['COV'])))
    lm_res = np.linalg.lstsq(cov, phn, rcond=None)
    cov_beta = lm_res[0]
    #resid = lm_res[1] #this is the sum of squared residuals for the model, not the vector of individual-level residuals as you might expect.

    #phn_mean = np.mean(resid)
    #phn_std = np.std(resid)
    
    resid = phn - cov.dot(cov_beta)
    phn_mean = np.mean(resid)
    phn_std = np.std(resid)

    print('... phenotypes and covariates from %d individuals read from %s ...' % (len(vld_dict['ID']), vld_phn_cov))
    return phn_mean, phn_std, cov_beta


def parse_frq(vld_frq_prefix, chrom):
    print('... parse the frequency of variants for the validation dataset: %s ...' % (vld_frq_prefix + '.frq'))

    ATGC = ['A', 'T', 'G', 'C']
    frq_dict = {'SNP':[], 'A1':[], 'A2':[], 'MAF':[]}
    with open(vld_frq_prefix + '.frq') as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom and ll[2] in ATGC and ll[3] in ATGC:
                frq_dict['SNP'].append(ll[1])
                frq_dict['A1'].append(ll[2])
                frq_dict['A2'].append(ll[3])
                frq_dict['MAF'].append(float(ll[4]))

    print('... %d SNPs on chromosome %d read from %s ...' % (len(frq_dict['SNP']), chrom, vld_frq_prefix + '.frq'))
    return frq_dict


def parse_tst(tst_prefix):
    print('... parse the list of variants for the testing dataset: %s ...' % (tst_prefix + '.bim'))

    tst_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(tst_prefix + '.bim') as ff:
        for line in ff:
            ll = (line.strip()).split()
            tst_dict['SNP'].append(ll[1])
            tst_dict['A1'].append(ll[4])
            tst_dict['A2'].append(ll[5])

    print('... %d SNPs read from %s ...' % (len(tst_dict['SNP']), tst_prefix + '.bim'))
    return tst_dict


def parse_psteff(ref_dict, frq_dict, tst_dict, pst_eff, psi_est, chrom):
    print('... parse posterior effect estimates: %s ...' % pst_eff)
    print('... parse posterior shrinkage estimates: %s ...' % psi_est)

    pst_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(pst_eff) as ff:
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom:
                pst_dict['SNP'].append(ll[1])
                pst_dict['A1'].append(ll[3])
                pst_dict['A2'].append(ll[4])

    print('... %d SNPs on chromosome %d read from %s ...' % (len(pst_dict['SNP']), chrom, pst_eff))


    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    tst_snp = set(zip(tst_dict['SNP'], tst_dict['A1'], tst_dict['A2']))

    ref_snp = set(zip(ref_dict['SNP'], ref_dict['A1'], ref_dict['A2'])) | set(zip(ref_dict['SNP'], ref_dict['A2'], ref_dict['A1'])) | \
              set(zip(ref_dict['SNP'], [mapping[aa] for aa in ref_dict['A1']], [mapping[aa] for aa in ref_dict['A2']])) | \
              set(zip(ref_dict['SNP'], [mapping[aa] for aa in ref_dict['A2']], [mapping[aa] for aa in ref_dict['A1']]))

    pst_snp = set(zip(pst_dict['SNP'], pst_dict['A1'], pst_dict['A2'])) | set(zip(pst_dict['SNP'], pst_dict['A2'], pst_dict['A1'])) | \
              set(zip(pst_dict['SNP'], [mapping[aa] for aa in pst_dict['A1']], [mapping[aa] for aa in pst_dict['A2']])) | \
              set(zip(pst_dict['SNP'], [mapping[aa] for aa in pst_dict['A2']], [mapping[aa] for aa in pst_dict['A1']]))

    if frq_dict != None:
        frq_snp = set(zip(frq_dict['SNP'], frq_dict['A1'], frq_dict['A2'])) | set(zip(frq_dict['SNP'], frq_dict['A2'], frq_dict['A1'])) | \
            set(zip(frq_dict['SNP'], [mapping[aa] for aa in frq_dict['A1']], [mapping[aa] for aa in frq_dict['A2']])) | \
            set(zip(frq_dict['SNP'], [mapping[aa] for aa in frq_dict['A2']], [mapping[aa] for aa in frq_dict['A1']]))

        comm_snp = tst_snp & ref_snp & frq_snp & pst_snp

        print('... %d common SNPs in the reference panel, validation frequency file, posterior effect estimates file, and testing set ...' % len(comm_snp))

    else:
        comm_snp = tst_snp & ref_snp & pst_snp

        print('... %d common SNPs in the reference panel, posterior effect estimates file, and testing set ...' % len(comm_snp))


    if frq_dict != None:
        vld_frq = {}
        for (ii, snp) in enumerate(frq_dict['SNP']):
            a1 = frq_dict['A1'][ii]; a2 = frq_dict['A2'][ii];
            if (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                vld_frq.update({snp: frq_dict['MAF'][ii]})
            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                vld_frq.update({snp: 1-frq_dict['MAF'][ii]})

    pst_beta = {}
    with open(pst_eff) as ff:
        for line in ff:
            ll = (line.strip()).split()
            snp = ll[1]; a1 = ll[3]; a2 = ll[4]
            if (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                pst_beta.update({snp: float(ll[5])})
            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                pst_beta.update({snp: -1*float(ll[5])})

    pst_psi = {}
    with open(psi_est) as ff:
        for line in ff:
            ll = (line.strip()).split()
            snp = ll[0]
            if snp in pst_beta:
                pst_psi.update({snp: float(ll[1])})

    pst_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'MAF':[], 'BETA':[], 'PSI':[], 'FLP':[]}
    for (ii, snp) in enumerate(ref_dict['SNP']):
        if snp in pst_beta:
            pst_dict['SNP'].append(snp)
            pst_dict['CHR'].append(ref_dict['CHR'][ii])
            pst_dict['BP'].append(ref_dict['BP'][ii])
            pst_dict['BETA'].append(pst_beta[snp])
            pst_dict['PSI'].append(pst_psi[snp])

            a1 = ref_dict['A1'][ii]; a2 = ref_dict['A2'][ii]
            if (snp, a1, a2) in comm_snp:
                pst_dict['A1'].append(a1)
                pst_dict['A2'].append(a2)
                pst_dict['FLP'].append(1)
                if frq_dict != None:
                    pst_dict['MAF'].append(vld_frq[snp])
                else:
                    pst_dict['MAF'].append(ref_dict['MAF'][ii])
            elif (snp, a2, a1) in comm_snp:
                pst_dict['A1'].append(a2)
                pst_dict['A2'].append(a1)
                pst_dict['FLP'].append(-1)
                if frq_dict != None:
                    pst_dict['MAF'].append(vld_frq[snp])
                else:
                    pst_dict['MAF'].append(1-ref_dict['MAF'][ii])
            elif (snp, mapping[a1], mapping[a2]) in comm_snp:
                pst_dict['A1'].append(mapping[a1])
                pst_dict['A2'].append(mapping[a2])
                pst_dict['FLP'].append(1)
                if frq_dict != None:
                    pst_dict['MAF'].append(vld_frq[snp])
                else:
                    pst_dict['MAF'].append(ref_dict['MAF'][ii])
            elif (snp, mapping[a2], mapping[a1]) in comm_snp:
                pst_dict['A1'].append(mapping[a2])
                pst_dict['A2'].append(mapping[a1])
                pst_dict['FLP'].append(-1)
                if frq_dict != None:
                    pst_dict['MAF'].append(vld_frq[snp])
                else:
                    pst_dict['MAF'].append(1-ref_dict['MAF'][ii])

    snp_dict = {snp: ii for (ii, snp) in enumerate(tst_dict['SNP'])}
    tst_idx = [snp_dict.get(snp) for snp in pst_dict['SNP']]

    return pst_dict, tst_idx


def parse_ldblk(ldblk_dir, pst_dict, chrom):
    print('... parse reference LD on chromosome %d ...' % chrom)

    snp_dict = {snp: ii for (ii, snp) in enumerate(pst_dict['SNP'])}

    if '1kg' in os.path.basename(ldblk_dir):
        chr_name = ldblk_dir + '/ldblk_1kg_chr' + str(chrom) + '.hdf5'
    elif 'ukbb' in os.path.basename(ldblk_dir):
        chr_name = ldblk_dir + '/ldblk_ukbb_chr' + str(chrom) + '.hdf5'

    hdf_chr = h5py.File(chr_name, 'r')
    n_blk = len(hdf_chr)
    ld_blk = []
    ld_blk = [np.array(hdf_chr['blk_'+str(blk)]['ldblk']) for blk in range(1,n_blk+1)]

    snp_blk = []
    for blk in range(1,n_blk+1):
        snp_blk.append([bb.decode("UTF-8") for bb in list(hdf_chr['blk_'+str(blk)]['snplist'])])

    blk_size = []
    mm = 0
    for blk in range(n_blk):
        idx = [ii for (ii, snp) in enumerate(snp_blk[blk]) if snp in snp_dict]
        blk_size.append(len(idx))
        if idx != []:
            idx_blk = range(mm,mm+len(idx))
            flip = [pst_dict['FLP'][jj] for jj in idx_blk]
            ld_blk[blk] = ld_blk[blk][np.ix_(idx,idx)]*np.outer(flip,flip)

            _, s, v = linalg.svd(ld_blk[blk])
            h = np.dot(v.T, np.dot(np.diag(s), v))
            ld_blk[blk] = (ld_blk[blk]+h)/2

            mm += len(idx)
        else:
            ld_blk[blk] = np.array([])

    return ld_blk, blk_size


