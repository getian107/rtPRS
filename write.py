#!/usr/bin/env python

"""
Write updated posterior SNP effects.

"""


def prs(out_dir, indv_dict, prs, rate, imp, order, chrom):
    print('... write PRS ...')

    if imp == 'True':
        filename = out_dir + '_prs_rate%.1f_imp_%s_chr%d.txt' % (rate, order, chrom)
    elif imp == 'False':
        filename = out_dir + '_prs_rate%.1f_exp_%s_chr%d.txt' % (rate, order, chrom)

    with open(filename, 'w') as ff:
        for subj, prs_indv in zip(indv_dict['ID'], prs):
            ff.write('%s\t%.6f\n' % (subj, prs_indv))


def psteff(out_dir, pst_dict, beta_updt, rate, imp, order, chrom):
    print('... write updated posterior SNP effects ...')

    if imp == 'True':
        filename = out_dir + '_psteff_rate%.1f_imp_%s_chr%d.txt' % (rate, order, chrom)
    elif imp == 'False':
        filename = out_dir + '_psteff_rate%.1f_exp_%s_chr%d.txt' % (rate, order, chrom)

    with open(filename, 'w') as ff:
        for snp, bp, a1, a2, beta in zip(pst_dict['SNP'], pst_dict['BP'], pst_dict['A1'], pst_dict['A2'], beta_updt):
            ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))


