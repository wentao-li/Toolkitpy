#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:writer.py
@time:1/17/23 3:31 AM
"""
import os

def check_dir(outdir):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    else:
        print(f"{outdir} exists. Exsisting Files will be updated.")
