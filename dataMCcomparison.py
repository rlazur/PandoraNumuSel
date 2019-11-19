#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 18:53:26 2019

@author: ryan
"""
cut_branches = ['nslice', 'crtveto', 'topological_score', 
                'trk_pid_chipr_v', 'trk_pid_chimu_v', 
                'trk_end_[xyz]_v', 'trk_start_[xyz]_v', 
                'trk_score_v', 'trk_distance', 'trk_len_v', 
                'trk_energy_muon_v',
                'reco_nu_vtx_sce_[xyz]', 'reco_nu_vtx_[xyz]',
                'trk_llr_pid_[uvy]_v', 'trk_llr_pid_v', 'trk_llr_pid_score_v',
                'crthitpe', '_closestNuCosmicDist']
MC_branches = ['weightSpline','backtracked_pdg','backtracked_e']
other_branches = ['category','nu_pdg','slpdg',
                  'nu_e', 'trk_energy_muon_v', 
                  'trk_theta_v','true_nu_vtx_[xyz]','ccnc']
all_branches = cut_branches + other_branches

CATS = ['COSMIC','OOFV','NUMUCC','OTHER']


AVx = [-1.55,254.8]
AVy = [-115.53,117.47]
AVz = [0.1, 1036.9]

####################################################################
#####POWER FUNCTIONS IN HELPER SCRIPT
###################################################################
import helperRun3 as helper

scale_NU = helper.scale_NU
scale_EXT = helper.scale_EXT
scale_DIRT = helper.scale_DIRT
POT = helper.POT

import codes

def unique_entries(df,cuts=False):
    return helper.unique_entries(df,cuts=cuts)
  
    

#should only need to do this the first time the script is run
try:
    SETUPDONE
    print("passing setup...")
except:
    print("doing setup...")
    import uproot
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import plotly.graph_objects as go

    #open up most modern samples
    nu_file = uproot.open("ROOTtrees/harddrive/Run3_1107/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run3_reco2_G_reco2.root")
    sel_nu = nu_file["nuselection"]["NeutrinoSelectionFilter"]
    data_file = uproot.open("ROOTtrees/harddrive/Run3_1107/data_bnb_mcc9.1_v08_00_00_25_reco2_G1_beam_good_reco2_1e19.root")
    sel_data = data_file["nuselection"]["NeutrinoSelectionFilter"]
    ext_file = uproot.open("ROOTtrees/harddrive/Run3_1107/data_extbnb_mcc9.1_v08_00_00_25_reco2_G_all_reco2.root")
    sel_ext = ext_file["nuselection"]["NeutrinoSelectionFilter"]
    dirt_file = uproot.open("ROOTtrees/harddrive/Run3_1107/prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run3_reco2_reco2.root")
    sel_dirt = dirt_file["nuselection"]["NeutrinoSelectionFilter"]
SETUPDONE = True

#change LOADDFS to True in terminal to load on non-first running of script
try:
    if LOADDFS:
        JUMPTOEXCEPT
    else:
        print('...not loading dataframes')
except:
    print("loading dataframes...")
    df_data = sel_data.pandas.df(all_branches, flatten=True)
    df_nu = sel_nu.pandas.df(all_branches+MC_branches, flatten=True)
    df_ext = sel_ext.pandas.df(all_branches, flatten=True)
    df_dirt = sel_dirt.pandas.df(all_branches+MC_branches, flatten=True)
    DFs = {'NU':df_nu, 'EXT':df_ext, 'DIRT':df_dirt, 'DATA':df_data}
LOADDFS = False
try:
    #apply calculated columns do dfs
    #chipr/chimu, costheta
    if APPLYCALCOLS:
        JUPPTOEXCEPT
    else:
        print("...not applying calculted columns")
except:
    print("applying calculated columns...")
    helper.apply_calcols(DFs)
APPLYCALCOLS=False
#try:
#    if APPLYCUTS:
#        JUMPTOEXCEPT
#    else:
#        print("...not applying cuts")
#except:
#    print("applying cuts...") 
#    helper.apply_cuts(DFs)
#APPLYCUTS = False
#make categories
#try:
#    if APPLYCATS:
#        JUMPTOEXCEPT
#    else:
#        print("...not applying subcategories")
#except:
#    print("applying subcategories...")
#    helper.apply_subgroups(df_nu)
#APPLYCATS = False


print("...done")

######################################################
### WORKSPACE BELOW FOR NEW FUNCTIONS
#####################################################
# MAKING DATA-MC PLOTS

FVx = [5,251]
FVy = [-110,110]
FVz = [20,986]
QUERY = 'nslice == 1'
QUERY += ' and reco_nu_vtx_sce_x > {} and reco_nu_vtx_sce_x < {}'.format(FVx[0],FVx[1])
QUERY += ' and reco_nu_vtx_sce_y > {} and reco_nu_vtx_sce_y < {}'.format(FVy[0],FVy[1])
QUERY += ' and reco_nu_vtx_sce_z > {} and reco_nu_vtx_sce_z < {}'.format(FVz[0],FVz[1])
QUERY += ' and (crtveto!=1 or crthitpe < 100.) and (_closestNuCosmicDist > 20.)'
QUERY += ' and trk_len_v > 20'
QUERY += ' and topological_score > 0.06'

KIND = 'category' #interaction
VAR = 'trk_len_v'
XLABEL = 'Track Length [cm]'
TITLE = 'MicroBooNE Preliminary {} POT'.format(POT)
BINS = 30
RANGE = (0, 500)
VERSION = '1107_v26'
SAVENAME = "DataMC_{}_{}_{}.jpg".format(VAR,KIND,VERSION)

fig, axes0, axes1 = helper.comp_DATAMC(DFs,
                                        VAR,
                                        cuts = QUERY,
                                        kind = KIND,
                                        take_longest=True,
                                        title = TITLE,
                                        xlabel = XLABEL,
                                        bins = BINS,
                                        range = RANGE
                                        )


fig.savefig("/home/ryan/HEPPA/eLEE/plots/random/{}".format(SAVENAME), pad_inches=0, bbox_inches='tight')
plt.show()

###################################### 
#### Efficiency of the fiducial volume cuts

#fig = plt.figure(figsize=(8,6))
#scale_NU = 0.0107 
#scale_EXT = 0.0718
#scale_DIRT = 0.0342 
#ACCEPTANCE = 'nu_pdg==14 and ccnc==0 and longest and topo06 and len20 and crtveto'
#B = np.linspace(-0.5,20.5,22)
#
#centers,vals,errs = helper.FV_Eff(df_nu,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='k')
#
#plt.xlabel('dist of FV from AV [cm]')
#plt.ylabel(r'$\nu_{\mu}$ CC INC selection efficiency')
#plt.ylim([0,1])
#plt.grid()
#plt.title(r'$\nu_{\mu}$ CC INC Selection Efficiency, BNB overlay',fontsize=25)
#plt.legend(fontsize=12)
#fig.patch.set_facecolor('silver')
#plt.rcParams['savefig.facecolor']='silver'
#plt.savefig("/home/ryan/HEPPA/eLEE_numu/plots/FVstudy/eff_{}.png".format("fveff"))


######################################
#### EFFICIENCY CURVES        
#fig = plt.figure(figsize=(8,6))
#efficiency vs neutrino energy

#only study true NUMUCC events
scale_NU = 0.0107 
scale_EXT = 0.0718
scale_DIRT = 0.0342 
ACCEPTANCE = 'nu_pdg==14 and ccnc==0 and longest'
VAR = 'nu_e'
B = np.linspace(0,3,25)
#B = np.array([0.,0.2,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,2.0,2.5,3.0])

#cuts = 'slice'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='k',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='g',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV1'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='r',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))

#cuts = 'slice and vtxFV and topo06'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='r',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))

#cuts = 'slice and vtxFV and topo10'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='g',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV and topo15'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='b',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV and topo25'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='y',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))

#cuts = 'slice and vtxFV and topo06 and len20'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='b',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV and topo06 and len20 and crtveto'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='g',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV and topo06 and len20'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='r',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV and topo06 and len30'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='y',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV and topo06 and len40'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='b',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))

#cuts = 'slice and vtxFV and topo06 and len20 and tscore80'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='pink',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV and topo06 and len20 and tscore80 and tdist4'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='c',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV and topo06 and len20 and crtveto'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='b',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV1 and topo06 and len20 and crtveto'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='y',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice and vtxFV and topo06 and len20 and contained25'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='g',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
##
#plt.xlabel('{}'.format(VAR))
#plt.xlabel("{}".format(VAR))
#plt.ylabel(r'$\nu_{\mu}$ CC INC selection efficiency')
#plt.ylim([0,1])
#plt.grid()
#plt.title(r'$\nu_{\mu}$ CC INC Selection Efficiency, BNB overlay',fontsize=25)
#plt.legend(fontsize=12)
#fig.patch.set_facecolor('silver')
#plt.rcParams['savefig.facecolor']='silver'
#plt.savefig("/home/ryan/Pictures/HEPPA/MicroBooNE/eLEE/Run3/eff_{}_{}.png".format(VAR,"newFVstudy"))


########################################################
######## PURITY CURVES
#fig = plt.figure(figsize=(8,6))

AVx = [-1.55,254.8]
AVy = [-115.53,117.47]
AVz = [0.1, 1036.9]

 

#efficiency vs neutrino energy
#B = np.linspace(0.5+116-20,0.5+116,21)
ACCEPTANCE = 'nu_pdg==14 and ccnc==0'
VAR = 'reco_nu_vtx_sce_z'
ncms = 75
#B = np.linspace(AVx[0] - 0.5,
#                AVx[0] + ncms + 0.5,
#                ncms+1)
B = np.linspace(AVz[1] - 0.5 - ncms,
                AVz[1] + 0.5,
                ncms+1)
label = 'FVpurity_zhigh'


#cuts = 'slice'
#centers,vals = helper.Pur(DFs,ACCEPTANCE,VAR,cuts,B)
#plt.plot(centers,vals,'o-k',label='{}'.format(cuts))

#cuts = 'slice and vtxFV'
#centers,vals = helper.Pur(DFs,VAR,cuts,B)
#plt.plot(centers,vals,'o-g',label='{}'.format(cuts))
#
#cuts = 'slice and vtxFV1'
#centers,vals = helper.Pur(DFs,VAR,cuts,B)
#plt.plot(centers,vals,'o-r',label='{}'.format(cuts))
#
#cuts = 'slice and vtxFV and topo06 and len20 and crtveto'
#centers,vals = helper.Pur(DFs,VAR,cuts,B)
#plt.plot(centers,vals,'o-b',label='{}'.format(cuts))
#
#cuts = 'slice and vtxFV1 and topo06 and len20 and crtveto'
#centers,vals = helper.Pur(DFs,VAR,cuts,B)
#plt.plot(centers,vals,'o-y',label='{}'.format(cuts))
##
#cuts = 'slice and topo06'
#centers,vals = helper.Pur(DFs,ACCEPTANCE,VAR,cuts,B)
#plt.plot(centers,vals,'o-r',label='{}'.format(cuts))
#
#cuts = 'slice and topo06 and len20'
#centers,vals = helper.Pur(DFs,ACCEPTANCE,VAR,cuts,B)
#plt.plot(centers,vals,'o-y',label='{}'.format(cuts))
#
#cuts = 'slice and topo06 and len20 and crtveto'
#centers,vals = helper.Pur(DFs,ACCEPTANCE,VAR,cuts,B)
#plt.plot(centers,vals,'o-g',label='{}'.format(cuts))

#
#plt.xlabel('{}'.format(VAR))
#plt.xlabel("{} [{}]".format(VAR,'cm'))
#plt.ylabel(r'$\nu_{\mu}$ CC INC selection Purity')
#plt.ylim([0,1])
#plt.grid()
#plt.title(r'$\nu_{\mu}$ CC INC Selection Purity, BNB overlay',fontsize=25)
#plt.legend(fontsize=15)
#fig.patch.set_facecolor('silver')
#plt.rcParams['savefig.facecolor']='silver'
#plt.savefig("/home/ryan/Pictures/HEPPA/MicroBooNE/eLEE/Run3/FVStudy/pur_{}_{}.png".format(VAR,label))

#######################################################
##### CUT FLOW TABLE

#allcuts = ['slice','vtxFV','topo06','len20','crtveto']
#flows = {'fullselection': ['slice','vtxFV','topo06','len20','crtveto']}
#
#df_flow = pd.DataFrame(index=(['NUMUCC']+allcuts), columns=flows)
#df_true = df_nu.query('nu_pdg==14 and OOFV==False and ccnc==0')
##purity calculations
#df_pur = pd.DataFrame(index=(['NUMUCC']+allcuts), columns=flows)
#    
#scale_NU = 0.0107 
#scale_EXT = 0.0718
#scale_DIRT = 0.0342 
#
#def pur_denom(DFs,cuts):
#    tot = 0
#    for mc in DFs:
#        if mc == 'NU':
#            scale = scale_NU
#        elif mc == 'DIRT':
#            scale = scale_DIRT
#        elif mc == 'EXT':
#            scale = scale_EXT
#        else:
#            print(mc)
#            scale = 0
#        tot += unique_entries(DFs[mc].query(cuts))*scale
#    return tot
#
#def pur_num(df_true,cuts):
#    return unique_entries(df_true.query(cuts))*scale_NU
#    
#for flow in flows:
#    cutflow = 'nslice!=-98' #dummy val
#    df_pur[flow]['NUMUCC'] = pur_num(df_true,cuts) / pur_denom(DFs,cutflow)
#    for cut in allcuts:
#        if cut in flows[flow]:
#            cutflow += ' and ' + cut
#            df_pur[flow][cut] = pur_num(df_true,cutflow) / pur_denom(DFs,cutflow)
#        else:
#            df_pur[flow][cut] = 0
#df_pur.to_csv('df_pur.csv')
#
##get a dataframe of pure true, contained numucc events
#
#for flow in flows:
#    cutflow = 'nslice!=-98' #dummy val for smoother code
#    df_flow[flow]['NUMUCC'] = unique_entries(df_true,cuts=cutflow)*scale_NU
#    for cut in allcuts:
#        if cut in flows[flow]:
#            cutflow += " and " + cut
#            df_flow[flow][cut] = unique_entries(df_true,cuts=cutflow)*scale_NU
#        else:
#            df_flow[flow][cut] = 0
#
#df_effs = df_flow.copy()
#for flow in flows:
#    for cut in allcuts:
#        val = df_effs[flow][cut]
#        if val > 0:
#            df_effs[flow][cut] = val*1.0 / df_effs[flow]['NUMUCC']
#    df_effs[flow]['NUMUCC'] = 1
#            
#print(df_flow.replace(0,''))
#print(df_effs.replace(0,''))
#print(df_pur.replace(0,''))
#df_flow.to_csv("df_flow.csv")
#df_effs.to_csv("df_effs.csv")    