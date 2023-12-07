#!/usr/bin/env python

"""
Online update of posterior SNP effects.

"""


import numpy as np
from numpy import linalg 
from pandas_plink import read_plink1_bin
from pandas_plink import Chunk


def sgd(n_gwas, phn_mean, phn_std, cov_beta, tst_prefix, tst_phn_cov, pst_dict, tst_idx, ld_blk, blk_size, rate, imp, order, chrom):
    print('... parse genotypes for the testing  dataset: %s ...' % (tst_prefix + '.bed/.bim/.fam'))

    tst = read_plink1_bin(bed=(tst_prefix + '.bed'), ref='a0')
    tst = tst[:,tst_idx]
    n_subj, n_snp = tst.shape

    print('... %d individuals and %d SNPs on chromosome %d read from %s ...' % (n_subj, n_snp, chrom, tst_prefix+'.bed'))


    print('... parse phenotypes and covariates for the testing dataset: %s ...' % tst_phn_cov)

    indv_dict = {'ID':[], 'PHN':[], 'COV':[]}
    with open(tst_phn_cov) as ff:
        for line in ff:
            ll = (line.strip()).split()
            indv_dict['ID'].append(ll[1])
            indv_dict['PHN'].append(float(ll[2]))
            indv_dict['COV'].append([float(ii) for ii in ll[3:]])

    phn = np.array(indv_dict['PHN'], ndmin=2).T
    cov = np.column_stack((np.ones(len(phn)), np.array(indv_dict['COV'])))

    iid = list(map(str, indv_dict['ID']))
    tst = tst.sel(sample=iid).reindex(sample=iid)
    n_subj, n_snp = tst.shape
    
    print('... %d individuals matched between the genotype and phenotype files ...' % (n_subj))


    print('... parse allele frequencies and posterior effects ...')

    frq = np.array(pst_dict['MAF'], ndmin=2).T

    beta_new = np.array(pst_dict['BETA'], ndmin=2).T*np.sqrt(2.0*frq*(1.0-frq))
    psi = np.array(pst_dict['PSI'], ndmin=2).T


    if order == '2nd':
        print('... calculate matrix inverse in second-order precedures ...')
        mm = 0
        n_blk = len(blk_size)
        for kk in range(n_blk):
            p = blk_size[kk]
            if p == 0:
                continue
            else:
                idx_blk = range(mm,mm+p)
                ld_blk[kk] = linalg.inv(ld_blk[kk]+np.diag(1.0/psi[idx_blk].T[0]))
                mm += p


    print('... update posterior effects ...')

    prs = np.zeros((n_subj,1))
    for ii in range(n_subj):
        subj_idx = ii+1
        print('... subject %d ...' % subj_idx)

        x = (np.array(tst[ii,:], ndmin=2).T-2*frq)/np.sqrt(2.0*frq*(1.0-frq))
        y_resid = phn[ii] - cov[ii].dot(cov_beta)
        y = (y_resid-phn_mean)/phn_std
        gamma = rate/(n_gwas+ii+1)

        x = np.nan_to_num(x, nan=0)

        prs[ii] = np.matmul(x.T,beta_new)
        beta_old = beta_new

        if order == '1st':
            if imp == 'False':
                beta_new = beta_old-gamma*x*np.matmul(x.T,beta_old)-gamma*(1.0/psi)*beta_old+gamma*y*x
            elif imp == 'True':
                mm = 0
                n_blk = len(blk_size)
                for kk in range(n_blk):
                    p = blk_size[kk]
                    if p == 0:
                        continue
                    else:
                        idx_blk = range(mm,mm+p)
                        x_blk = x[idx_blk]

                        ainv = 1.0/(1.0/gamma+1.0/psi[idx_blk])
                        ainvx = ainv*x_blk

                        slp = 1.0/gamma*(np.diag(ainv.T[0])-np.outer(ainvx,ainvx)/(1.0+np.matmul(x_blk.T,ainvx)))
                        beta_new[idx_blk] = np.matmul(slp, beta_old[idx_blk]+gamma*y*x_blk)
                        mm += p

        if order == '2nd':
            mm = 0
            for kk in range(n_blk):
                p = blk_size[kk]
                if p == 0:
                    continue
                else:
                    idx_blk = range(mm,mm+p)
                    x_blk = x[idx_blk]

                    if imp == 'False':
                        beta_new[idx_blk] = (1.0-gamma)*beta_old[idx_blk]+gamma*np.matmul(ld_blk[kk],y*x_blk)
                    elif imp == 'True':
                        beta_new[idx_blk] = 1.0/(1.0+gamma)*beta_old[idx_blk]+gamma/(1.0+gamma)*np.matmul(ld_blk[kk],y*x_blk)

                    mm += p


    beta_new = beta_new/np.sqrt(2.0*frq*(1.0-frq))

    return indv_dict, prs, beta_new


