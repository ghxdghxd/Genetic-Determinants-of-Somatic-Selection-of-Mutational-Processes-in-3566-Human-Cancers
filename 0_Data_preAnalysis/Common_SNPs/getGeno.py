#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# getGeno.py
# @Author : JT Guo
# @Email  : guojt-4451@163.com
# @Date   : 2018-5-10 09:09:45

import pandas as pd
import glob
import multiprocessing
import sys

def getGeno(info, geno):
    info = pd.read_csv(info, compression="gzip", header=0,
                       keep_default_na=False, dtype="str", sep="\t")
    # geno1 = geno.loc[geno.rsid.isin(info.ID2)]
    # info.ID1[info.ID2.isin(geno1.rsid)]
    geno1 = pd.merge(info, geno, left_on="ID2", right_on="rsid")
    geno1.drop(["ID2", "CEU_MAF", "rsid"], axis=1, inplace=True)
    geno1.rename(columns={"ID1": "ID"}, inplace=True)
    # geno.to_csv(file.replace("new", "rare"), index=False, header=True,
                # compression="gzip", sep="\t", na_rep='NA')
    return(geno1)

cancer=sys.argv[1]

info_files = glob.glob("TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/2CEU_com/*.CEU.COM.AF.gz")

geno = pd.read_csv("TCGA/SNP6_Genotype/6_filterCN_INS_fcGene_geno_thresh_maxProb0.9_info0.9/" + cancer + ".geno.gz",
                   header=0, sep="\t", compression="gzip", dtype="str")

print(cancer)

pool = multiprocessing.Pool(processes=int(12))
result = []
for file in info_files:
    print(file)
    result.append(pool.apply_async(getGeno, (file, geno)))
pool.close()
pool.join()

com_geno = list(map(lambda x: x.get(), q))

com_geno = pd.concat(com_geno, axis = 0)

com_geno.to_csv(cancer + ".geno.gz", index=False, header=True,
            compression="gzip", sep="\t", na_rep='NA')
