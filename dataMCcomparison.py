#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 18:53:26 2019

@author: ryan
"""
cut_branches = ['nslice', 'crtveto', 'topological_score', 'trk_pid_chipr_v', 'trk_pid_chimu_v', 'trk_end_[xyz]_v', 'trk_start_[xyz]_v', 'trk_score_v', 'trk_distance', 'trk_len_v', 'trk_energy_muon_v']
MC_branches = ['weightSpline','backtracked_pdg','backtracked_e']
other_branches = ['category','nu_pdg','slpdg','nu_e', 'trk_energy_muon_v', 'trk_theta_v', 'nu_vtx_[xyz]','true_nu_vtx_[xyz]','ccnc']
all_branches = cut_branches + other_branches

CATS = ['COSMIC','OOFV','NUMUCC','OTHER']


AVx = [-1.55,254.8]
AVy = [-115.53,117.47]
AVz = [0.1, 1036.9]

####################################################################
#####POWER FUNCTIONS IN HELPER SCRIPT
###################################################################
import helperRun3 as helper
import codes
def comp_DATAMC(DFs, VAR, BINEDGES, lw=25, take_longest=True, cuts=False, version=''):
    return helper.comp_DATAMC(DFs,VAR,BINEDGES,lw=lw,take_longest=take_longest,cuts=cuts,version=version)

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
    nu_file = uproot.open("ROOTtrees/harddrive/Run31e19/prodgenie_bnb_nu_uboone_overlay_mcc9.1_run3_G_reco2.root")
    sel_nu = nu_file["nuselection"]["NeutrinoSelectionFilter"]
    data_file = uproot.open("ROOTtrees/harddrive/Run31e19/data_bnb_optfilter_G1_1e19_mcc9.1_reco2.root")
    sel_data = data_file["nuselection"]["NeutrinoSelectionFilter"]
    ext_file = uproot.open("ROOTtrees/harddrive/Run31e19/data_extbnb_mcc9.1_v08_00_00_16_18_run3_G_reco2.root")
    sel_ext = ext_file["nuselection"]["NeutrinoSelectionFilter"]
    dirt_file = uproot.open("ROOTtrees/harddrive/Run31e19/prodgenie_bnb_dirt_overlay_run3_mcc9.1_v08_00_00_18_G_reco2.root")
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
try:
    if APPLYCUTS:
        JUMPTOEXCEPT
    else:
        print("...not applying cuts")
except:
    print("applying cuts...") 
    helper.apply_cuts(DFs)
APPLYCUTS = False
#make categories
try:
    if APPLYCATS:
        JUMPTOEXCEPT
    else:
        print("...not applying subcategories")
except:
    print("applying subcategories...")
    helper.apply_subgroups(df_nu)
APPLYCATS = False


print("...done")
plt.rcParams.update({'font.size': 22})

######################################################
### WORKSPACE BELOW FOR NEW FUNCTIONS
#####################################################
# MAKING DATA-MC PLOTS




#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and topo06', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and crtveto', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and muon', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and topo06 and len20', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and topo06 and len40', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and topo06 and crtveto', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and crtveto', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and contained25', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and topo06 and tscore80 and tdist4 and contained25', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and topo06 and tscore80 and tdist4 and crtveto', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and topo06 and tscore80 and tdist4', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and topo06 and len20 and tscore80 and tdist4 and crtveto', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and vtxFV and topo06 and len20 and crtveto', lw=25, take_longest=True, version='OFFICIAL')
#
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and topo15', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and contained10', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and contained25', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and topo06 and crtveto', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and crtveto', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and contained25', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and topo06 and vtx_y', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and topo06 and contained25 and vtx_y', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and topo06 and contained25 and crtveto', lw=25, take_longest=True, version='OFFICIAL')

#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and len20', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='slice and len40', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'costheta',np.linspace(-1,1,30), cuts='basic', lw=25, take_longest=True, version='OFFICIAL')

#comp_DATAMC(DFs,'trk_len_v',np.linspace(0,225,15), cuts='slice', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'trk_len_v',np.linspace(0,500,15), cuts='slice and vtxFV', lw=25, take_longest=True, version='OFFICIAL500cm')
#comp_DATAMC(DFs,'trk_len_v',np.linspace(0,500,15), cuts='slice and vtxFV and topo06', lw=25, take_longest=True, version='OFFICIAL500cm')
#comp_DATAMC(DFs,'trk_len_v',np.linspace(0,500,15), cuts='slice and vtxFV and topo06 and len20', lw=25, take_longest=True, version='OFFICIAL500cm')
#comp_DATAMC(DFs,'trk_len_v',np.linspace(0,500,15), cuts='slice and vtxFV and topo06 and len20 and crtveto', lw=25, take_longest=True, version='OFFICIAL500cm')
#comp_DATAMC(DFs,'trk_len_v',np.linspace(0,250,15), cuts='slice and vtxFV and topo06 and crtveto and len20', lw=25, take_longest=True, version='OFFICIAL')
#
#comp_DATAMC(DFs,'trk_pid_chipr_v',np.linspace(0,300,30), cuts='slice and vtxFV and topo06', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'trk_pid_chipr_v',np.linspace(0,300,30), cuts='slice and vtxFV', lw=25, take_longest=True, version='OFFICIAL')
#
#
#comp_DATAMC(DFs,'topological_score',np.linspace(0,1,30), cuts='slice and vtxFV', lw=25, take_longest=True, version='OFFICIALpart')
#comp_DATAMC(DFs,'topological_score',np.linspace(0,0.3,30), cuts='slice and vtxFV', lw=25, take_longest=True, version='OFFICIALpart')
#comp_DATAMC(DFs,'topological_score',np.linspace(0,.3,30), cuts='slice and vtxFV and contained25', lw=25, take_longest=True, version='OFFICIALpart')
#comp_DATAMC(DFs,'topological_score',np.linspace(0,.3,30), cuts='slice and crtveto', lw=25, take_longest=True, version='OFFICIALpart')

#comp_DATAMC(DFs,'trk_energy_muon_v',np.linspace(0,2,30), cuts='slice', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'trk_energy_muon_v',np.linspace(0,2,30), cuts='slice and crtveto', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'trk_energy_muon_v',np.linspace(0,2,30), cuts='slice and topo06', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'trk_energy_muon_v',np.linspace(0,2,30), cuts='slice and contained25', lw=25, take_longest=True, version='OFFICIAL')

#comp_DATAMC(DFs,'trk_distance',np.linspace(0,6,30), cuts='slice and vtxFV', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'trk_distance',np.linspace(0,6,30), cuts='slice and vtxFV and topo06', lw=25, take_longest=True, version='OFFICIAL')
#
#comp_DATAMC(DFs,'trk_score_v',np.linspace(0,1,30), cuts='slice and vtxFV and topo06', lw=25, take_longest=True, version='OFFICIAL')
#comp_DATAMC(DFs,'trk_score_v',np.linspace(0,1,30), cuts='slice and vtxFV and topo06', lw=25, take_longest=True, version='OFFICIAL')

######################################
#### EFFICIENCY CURVES        
#fig = plt.figure(figsize=(8,6))
#efficiency vs neutrino energy

#only study true NUMUCC events
scale_NU = 0.0107 
scale_EXT = 0.0718
scale_DIRT = 0.0342 
ACCEPTANCE = 'nu_pdg==14 and ccnc==0 and longest and OOFV==False'
VAR = 'nu_e'
B = np.linspace(0,3,2)
B = np.array([0.,0.2,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,2.0,2.5,3.0])
#B = np.linspace(0,5,30)

cuts = 'slice and vtxFV'
centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='black',label='{}'.format(cuts))
pur = helper.overallPur(DFs,cuts)
print("overall purity of {}: {}".format(cuts, pur))
#
#cuts = 'slice'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='k',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))

#cuts = 'slice and vtxFV'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='y',label='{}'.format(cuts))
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
#cuts = 'slice and vtxFV and topo06 and len20 and contained25'
#centers,vals,errs = helper.Eff(df_nu,VAR,cuts,ACCEPTANCE,B)
#plt.errorbar(centers,vals,yerr=errs,fmt='o-',color='g',label='{}'.format(cuts))
#pur = helper.overallPur(DFs,cuts)
#print("overall purity of {}: {}".format(cuts, pur))
##
#plt.xlabel('{}'.format(VAR))
#plt.xlabel("True Neutrino Energy")
#plt.ylabel(r'$\nu_{\mu}$ CC INC selection efficiency')
#plt.ylim([0,1])
#plt.grid()
#plt.title(r'$\nu_{\mu}$ CC INC Selection Efficiency, BNB overlay',fontsize=25)
#plt.legend(fontsize=12)
#fig.patch.set_facecolor('silver')
#plt.rcParams['savefig.facecolor']='silver'
#plt.savefig("/home/ryan/Pictures/HEPPA/MicroBooNE/eLEE/Run3/eff_{}_{}.png".format(VAR,"finalselection"))


########################################################
######## PURITY CURVES
#fig = plt.figure(figsize=(8,6))
#efficiency vs neutrino energy
B = np.linspace(0,225,20)
VAR = 'trk_len_v'

#cuts = 'slice and vtxFV'
#centers,vals = helper.Pur(DFs,VAR,cuts,B)
#plt.plot(centers,vals,'o-k',label='{}'.format(cuts))
#
#cuts = 'slice and vtxFV and topo06'
#centers,vals = helper.Pur(DFs,VAR,cuts,B)
#plt.plot(centers,vals,'o-r',label='{}'.format(cuts))

#cuts = 'slice and vtxFV and topo06 and len20'
#centers,vals = helper.Pur(DFs,VAR,cuts,B)
#plt.plot(centers,vals,'o-y',label='{}'.format(cuts))
#
#cuts = 'slice and vtxFV and topo06 and len20 and tscore80 and tdist4'
#centers,vals = helper.Pur(DFs,VAR,cuts,B)
#plt.plot(centers,vals,'o-r',label='{}'.format(cuts))
#
#cuts = 'slice and vtxFV and topo06 and len20 and tscore80 and tdist4 and crtveto'
#centers,vals = helper.Pur(DFs,VAR,cuts,B)
#plt.plot(centers,vals,'o-c',label='{}'.format(cuts))
#
#plt.xlabel('{}'.format(VAR))
#plt.xlabel("Track Length [cm]")
#plt.ylabel(r'$\nu_{\mu}$ CC INC selection Purity')
#plt.ylim([0,1])
#plt.grid()
#plt.title(r'$\nu_{\mu}$ CC INC Selection Purity, BNB overlay',fontsize=25)
#plt.legend(fontsize=15)
#fig.patch.set_facecolor('silver')
#plt.rcParams['savefig.facecolor']='silver'
#plt.savefig("/home/ryan/Pictures/HEPPA/MicroBooNE/eLEE/Run3/pur_{}_{}.png".format(VAR,"len"))

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