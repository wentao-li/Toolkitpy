#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:cluster_module.py
@time:1/14/23 10:24 AM
"""

import os


def check_module(module_load, command):
    if module_load is not None:
        #os.system(module_load)
        command = module_load + "\n" + command
    return command 
