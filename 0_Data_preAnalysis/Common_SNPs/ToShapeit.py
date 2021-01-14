#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ToShapeit.py
# @Author : JT Guo (guojt-4451@163.com)
# @Date   : 2018-5-1 14:53:52

import pandas as pd
import glob
import multiprocessing

def to_shapeit(f):
    reader = pd.read_csv(f, compression="gzip", comment="#", sep="\t", iterator=True, header=0)
    loop = True
    chunkSize = 100000
    chunks = []
    while loop:
        try:
            df = reader.get_chunk(chunkSize)
            df.iloc[:, 6:] = df.iloc[:, 6:].replace(
                0, "1,0,0").replace(1, "0,1,0").replace(2, "0,0,1")
            chunks.append(df)
        except StopIteration:
            loop = False
            print("Iteration is stopped")
    out_df = pd.concat(chunks, ignore_index=True)
    out_df.to_csv(f.strip(".gz"), sep="\t", index=None)

if __name__ == '__main__':
    geno_files = glob.glob("TCGA/CEU/*.gz")
    pool = multiprocessing.Pool(processes=int(12))
    for f in geno_files:
        pool.apply_async(to_shapeit, args = (f,))
    pool.close()
    pool.join()
