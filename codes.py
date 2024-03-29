#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 13:46:47 2019

@author: ryan
"""

def pdgtostring(code):
    if int(code) == 11:
        return 'electron'
    elif int(code) == -11:
        return 'positron'
    elif int(code) == 12:
        return 'nu_e'
    elif int(code) == -12:
        return 'anti_nu_e'
    elif int(code) == 13:
        return 'muon'
    elif int(code) == -13:
        return 'anti_muon'
    elif int(code) == 14:
        return 'nu_mu'
    elif int(code) == -14:
        return 'anti_nu_mu'
    elif int(code) == 2212:
        return 'proton'
    elif int(code) == 2112:
        return 'neutron'
    elif int(code) == 111:
        return 'pi_0'
    elif int(code) == -211:
        return 'pi_minus'
    elif int(code) == 211:
        return 'pi_plus'
    elif int(code) == 22:
        return 'gamma'
    elif int(code) == 311:
        return 'K_naught'
    else:
        return 'unknwnPDG'

category_colors = {
    4: "xkcd:light red",
    5: "xkcd:brick",
    2: "xkcd:cyan",
    21: "xkcd:cerulean",
    3: "xkcd:cobalt",
    31: "xkcd:sky blue",
    1: "xkcd:green",
    10: "xkcd:mint green",
    11: "xkcd:lime green",
    111: "xkcd:goldenrod",
    6: "xkcd:grey",
    0: "xkcd:black"
}

category_labels = {
    1: r"$\nu_e$ CC",
    10: r"$\nu_e$ CC0$\pi$0p",
    11: r"$\nu_e$ CC0$\pi$Np",
    111: r"MiniBooNE LEE",
    2: r"$\nu_{\mu}$ CC",
    21: r"$\nu_{\mu}$ CC $\pi^{0}$",
    3: r"$\nu$ NC",
    31: r"$\nu$ NC $\pi^{0}$",
    4: r"Cosmic",
    5: r"Out. fid. vol.",
    6: r"other",
    0: r"No slice"
}

def categorytostring(cat):
    if cat == 1:
        return "nu_e_other"
    elif cat == 10:
        return "nu_e_cc0pi0p"
    elif cat == 11:
        return "nu_e_cc0pinp"
    elif cat == 2:
        return "nu_mu_other"
    elif cat == 21:
        return "nu_mu_pi0"
    elif cat == 3:
        return "nc"
    elif cat == 31:
        return "nc_pi0"
    elif cat == 4:
        return "cosmic" 
        #https://github.com/ubneutrinos/searchingfornues/blob/master/Selection/AnalysisTools/DefaultAnalysis_tool.cc#L606
    elif cat == 5:
        return "outfv"
    elif cat == 6:
        return "other"
    elif cat == 0:
        return "data"
    else:
        return "unknwnCat"