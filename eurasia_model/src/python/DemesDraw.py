#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 12:31:38 2023

@author: mkw5910
"""

import demesdraw
import demes

complexModel = demes.load("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/src_PSU_simulations_Projects/complex_demographic_model/input/Demes_Draw_complexModel.yaml")
demesdraw.tubes(complexModel, log_time=True, seed=1234)

demesdraw.utils.separation_heuristic(complexModel)
