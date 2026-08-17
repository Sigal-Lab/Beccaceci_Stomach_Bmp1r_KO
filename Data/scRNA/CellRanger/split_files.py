#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 11:56:19 2021

@author: hilmar
"""

import os, shutil

data_folder = "."

ff = os.listdir(data_folder)
ff = [f for f in ff if os.path.splitext(f)[1] == ".gz"]

for f in ff:
    #print(f)
    tmp = f.split("_")
    name = "_".join(tmp[0:(len(tmp)-1)])
    suff_tmp = tmp[len(tmp)-1]

    tmp2 = suff_tmp.split(".")
    suff = ".".join(tmp2[1:len(tmp2)])
    name = name + "_" + tmp2[0]
    print(name, suff)
    
    target_folder = os.path.join(data_folder, name)
    if not os.path.isdir(target_folder):
        os.mkdir(target_folder)

    if suff == "genes.tsv.gz":
        suff = "features.tsv.gz"
        
    shutil.move(os.path.join(data_folder,f), os.path.join(target_folder,suff))
