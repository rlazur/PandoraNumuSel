#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 13:42:05 2019

@author: ryan
"""
import codes

import numpy as np
import matplotlib.pyplot as plt

plt.ioff()
params = {
        'font.size': 20,
        'legend.fontsize': 19,
        'axes.facecolor' : 'white',
        'savefig.facecolor' : 'silver'
        }
plt.rcParams.update(params)
from matplotlib import gridspec
from datetime import datetime

start_time = datetime.now()

#see the POTnormalizeation.pynb script to derive these
RUN = 3 #run3 data
# https://github.com/ubneutrinos/searchingfornues/wiki/NTuple-location
scale_NU = 0.00787 
scale_EXT = 0.0198
scale_DIRT = 0.0173
POT = 5.648e18

AVx = [-1.55,254.8]
AVy = [-115.53, 117.47]
AVz = [0.1, 1036.9]
#other choice
FVx = [AVx[0]+10,AVx[1]-10]
FVy = [AVy[0]+10,AVy[1]-10]
FVz = [AVz[0]+10,AVz[1]-50]
#updated for SCE
FVx1 = [5,251]
FVy1 = [-110,110]
FVz1 = [20,986]

#how far from the edges the containment cut is
CD = 25 #cm

CATS = ['COSMIC','OOFV','NUMUCC','OTHER']
BASIC = 'slice and topo06 and muon'

def searchkeys(key1,sel_file,key2=''):
    keys = sel_file.keys()
    print([key for key in keys if ((key1 in str(key)) and (key2 in str(key)))])
    
def get_histvals(df,var,B):
    '''
    This takes in dataframe, variable (x axis), and binedges
    Returns array of lists of bin vals defined by input params
    '''
    bin_vals = [0]*(len(B)-1)
    for i in range(len(B)-1):
        binmin = B[i]
        binmax = B[i+1]
        #just grab out the stuff in the ith bin
        query = '{} >= @binmin and {} < @binmax'.format(var,var)
        bin_vals[i] = df.query(query)[var].size
    return np.array(bin_vals)

def get_histerrs(df,var,B,scale):
    #the error for scaled histgrams is sqrt(sum(weights**2))
    #get list of arrays of weights per bin
    bin_weights = [0]*(len(B)-1)
    for i in range(len(B)-1):
        binmin = B[i]
        binmax = B[i+1]
        bincut = '{} >= @binmin and {} < @binmax'.format(var,var)
        df_bin = df.query(bincut)
#        if 'weightSpline' in df.keys():
#            bin_weights[i] = np.array(df_bin['weightSpline'])*scale
#        else:
#            bin_weights[i] = np.ones(df_bin.shape[0])*scale
        bin_weights[i] = np.ones(df_bin.shape[0])*scale
    #square weights and sum for each bin
    sw2s = [np.sum(weights**2) for weights in bin_weights]
    errs = np.sqrt(np.array(sw2s))
    return errs
        

def prep_hist(df,scale,VAR,BINEDGES):
    '''
    takes MC dataframes and histogram parameters and returns lists that can be 
    used the pyplot hist function.
    applies cuts to data frame and returns:
        VAR values (list)
        weights for the values (array)
    '''
    B0 = BINEDGES[0]
    B1 = BINEDGES[-1]
    if 'weightSpline' in df.keys():
        df_cut = df.query('{} >= @B0 and {} <= @B1'.format(VAR,VAR))[[VAR,'weightSpline']]
        weights = np.array(df_cut['weightSpline'])*scale
    else:
        df_cut = df.query('{} >= @B0 and {} <= @B1'.format(VAR,VAR))[[VAR]]
        weights = np.ones(df_cut.shape[0])*scale
    xstacked = list(df_cut[VAR])

    return list(xstacked),list(weights)

def comp_DATAMC(DFs, VAR, cuts=False, kind='category', 
                take_longest=True, title='MicroBooNE Preliminary {} POT'.format(POT), xlabel = '',
                **params):
    '''
    This function builds+saves a data-MC comparison histogram
    It needs a dictionary containing all the dataframes made using uproot
    FEATURES:
        option to plot all tracks or just the longest in each event/entry (take_longest)
        option to not implement any cuts at all
        
    NEED TO DO:
        implement spline weights
        
    '''
    try:
        BINEDGES = np.linspace(params['range'][0], params['range'][1], params['bins'] + 1)
    except:
        print("could not find range and/or bins to build user-defined BINEDGES")
        print("will use default BINEDGES of (0,0.1,...,1)")
        BINEDGES = np.linspace(0,1,11)
        
    BINCENTERS = 0.5*(BINEDGES[1:]+BINEDGES[:-1])
    hist_dict = {'MC':{
                    'vals': [], #list of lists of data correlated with WEIGHTS and LABELS (lists)
                    'errs': np.zeros(BINCENTERS.size), #to calculate scaled MC uncertainty on each bin
                    'weights': [], #weights of each component of xstacked (lists)
                    'labels': [], #labels for each component of xstacked (strings)
                    'colors': [] #color codes for each cateogry
                    
                    },
                'EXT':{
                    'vals': [], #list of external data
                    'errs': np.zeros(BINCENTERS.size), #errors on each bin
                    'weights': [], # weights for plt histogram fn
                    'label': [] #label for EXT sample
                    },
                'DIRT':{
                    'vals': [], #list of external data
                    'errs': np.zeros(BINCENTERS.size), #errors on each bin
                    'weights': [], # weights for plt histogram fn
                    'label': [], # label for DIRT sample
                    'color': ["xkcd:dirt brown"] # color for DIRT sample
                        },
                'DATA':{
                    'vals': [], #list of data
                    'errs': np.zeros(BINCENTERS.size), # counting error for data
                    'label': []
                    }
                 }
    for label in DFs:
    #loop through and pull from each dataframe
    #loop will complain if we don't know about one of the labels
    #start by setting the main dataframe (df)
    #taking longest tracks if appropriate
        df = DFs[label].copy()
        if take_longest:
            df = df.query('longest')
        if cuts:
            df = df.query(cuts.replace('basic',BASIC))
        scale = 0
        if label != 'DATA':
        #everything not data is MC ;)
            #there are certain things particular to each datatype
            if label == 'DIRT':
                scale = scale_DIRT
                print('prepping {}...'.format(label))
                x,w = prep_hist(df,scale,VAR,BINEDGES)
                errs = get_histerrs(df,VAR,BINEDGES,scale)
                hist_dict['DIRT']['vals'].append(x)
                hist_dict['DIRT']['weights'].append(w)
                hist_dict['DIRT']['errs'] += errs
                hist_dict['DIRT']['label'].append(label + ': ' + str(round(len(x)*scale)))
            if label == 'EXT':
                scale = scale_EXT
                print('prepping {}...'.format(label))
                x,w = prep_hist(df,scale,VAR,BINEDGES)
                errs = get_histerrs(df,VAR,BINEDGES,scale)
                hist_dict['EXT']['vals'].append(x)
                hist_dict['EXT']['weights'].append(w)
                hist_dict['EXT']['errs'] += errs
                hist_dict['EXT']['label'].append(label + ': ' + str(round(len(x)*scale))) 
            if label == 'NU':
                scale = scale_NU
                if kind == 'category':
                    #get dictionary of each category count
                    unique, counts = np.unique(df['category'], return_counts=True)
                    cat_totals = dict(zip(unique,counts))
                #separate by category
                    for cat in codes.category_labels:
                        if (cat in cat_totals.keys()):
                        #only go through this if the category exists
                            print('prepping {}...'.format(codes.category_labels[cat]))
                            df_cat = df.query('category == {}'.format(cat))
                            #put the cat in the histogram prepper
                            x,w = prep_hist(df_cat,scale,VAR,BINEDGES)
                            errs = get_histerrs(df_cat,VAR,BINEDGES,scale)
                            hist_dict['MC']['vals'].append(x)
                            hist_dict['MC']['weights'].append(w)
                            hist_dict['MC']['errs'] = np.sqrt(errs**2 + hist_dict['MC']['errs']**2)
                            hist_dict['MC']['labels'].append(codes.category_labels[cat] + ': ' + str(round(len(x)*scale)))
                            hist_dict['MC']['colors'].append(codes.category_colors[cat])
                else:
                #don't subserate MC overlay sample
                    print('prepping {}...'.format(label))
                    x,w = prep_hist(df,scale,VAR,BINEDGES)
                    errs = get_histerrs(df,VAR,BINEDGES,scale)
                    hist_dict['vals'].append(x)
                    hist_dict['weights'].append(w)
                    hist_dict['errs'] += errs
                    hist_dict['label'].append(label + ': ' + str(round(len(x)*scale)))
                    hist_dict['colors'].append("xkcd:purple")
    
        elif label == 'DATA':
            #a list of the values of each bin defined by BINEDGES
            hist_dict['DATA']['vals'] = get_histvals(df,VAR,BINEDGES)
            hist_dict['DATA']['errs'] = np.sqrt(np.array(hist_dict['DATA']['vals'])) #counting error
            
    for key in hist_dict:
        print(hist_dict[key]['errs'])
    
    #################################
    ##DATA-MC PLOT
    fig = plt.figure(figsize=(10,8))
    plt.ion()
    fig.patch.set_facecolor('silver')
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    
    axes0 = plt.subplot(gs[0])
    #first build up the MC histogram (MCoverlay + DIRT)
    MC_ystacked_vals,_,_ = axes0.hist(
            hist_dict['MC']['vals'] + hist_dict['DIRT']['vals'],
            label = hist_dict['MC']['labels'] + hist_dict['DIRT']['label'],
            weights = hist_dict['MC']['weights'] + hist_dict['DIRT']['weights'],
            color = hist_dict['MC']['colors'] + hist_dict['DIRT']['color'],
            histtype = 'stepfilled', stacked = True,
            **params
            )
    #return from .hist gives list of lists for each stacked hist
    #we just want the last (cumulative) one
    MC_ystacked_vals = MC_ystacked_vals[-1]
    #then add the hatched EXT hist on top
    EXT_ystacked_vals,_,_ = axes0.hist(
            hist_dict['EXT']['vals'],
            label = 'EXT',
            weights = hist_dict['EXT']['weights'],
            color = 'white',
            hatch = "//",
            bottom = MC_ystacked_vals,
            **params
            )
    
    #make big list of all the plotted values and their weights
    tot_vals = np.concatenate([np.concatenate(hist_dict['MC']['vals']), 
                               np.concatenate(hist_dict['DIRT']['vals']), 
                               np.concatenate(hist_dict['EXT']['vals'])])
    tot_weights = np.concatenate([np.concatenate(hist_dict['MC']['weights']),
                                  np.concatenate(hist_dict['DIRT']['weights']),
                                  np.concatenate(hist_dict['EXT']['weights'])])
    #just a nice outline of the whole thing
    axes0.hist(
            tot_vals,
            weights = tot_weights,
            histtype = 'step',
            edgecolor = 'black',
            **params
            )
    
    #put error bars on the MC hists + EXT
    mc_errs = np.sqrt(np.array(hist_dict['MC']['errs'])**2 + 
                      np.array(hist_dict['DIRT']['errs'])**2 + 
                      np.array(hist_dict['EXT']['errs'])**2)
    bin_sizes = [(BINEDGES[i+1] - BINEDGES[i])/2 for i in range(len(BINCENTERS))]
    axes0.errorbar(
            BINCENTERS,
            list(np.array(MC_ystacked_vals) + np.array(EXT_ystacked_vals)),
            yerr = mc_errs,
            xerr = bin_sizes,
            fmt = 'None', ecolor = 'grey', alpha = 0.7, fillstyle='full'
            )
    #plot data    
    axes0.errorbar(
            BINCENTERS,
            hist_dict['DATA']['vals'],
            yerr = hist_dict['DATA']['errs'],
            xerr = bin_sizes,
            label='data: {}'.format(np.sum(hist_dict['DATA']['vals'])),
            fmt='o', color='k', markersize=5
            )
    
    axes0.set_xlim([BINEDGES[0],BINEDGES[-1]])
    axes0.legend(loc='upper center', fontsize=10, ncol = 3, frameon = False)
    axes0.set_title(title)
    
    ##################################3
    ###Ratio Plot
    axes1 = plt.subplot(gs[1], sharex=axes0)
    #calculate some things
    mc_tot = MC_ystacked_vals + EXT_ystacked_vals
    dat_rat = hist_dict['DATA']['errs'] / hist_dict['DATA']['vals']
    mc_rat = mc_errs / mc_tot    
    ratio_bin_vals = hist_dict['DATA']['vals'] / mc_tot
    ratio_bin_errs = ratio_bin_vals * np.sqrt(dat_rat**2 + mc_rat**2)
    #plot it up
    axes1.errorbar(
            BINCENTERS,
            ratio_bin_vals,
            yerr = ratio_bin_errs,
            xerr = bin_sizes,
            fmt='o',color='k',markersize=5)
    
    ratio_error_mc = np.sqrt(2)*mc_rat
    ratio_error_mc = np.insert(ratio_error_mc,0,ratio_error_mc[0])
    axes1.fill_between(
            BINEDGES,
            1 - ratio_error_mc,
            ratio_error_mc + 1,
            step = "pre",
            color = "grey",
            alpha = 0.5
            )
    
    axes1.axhline(1, linestyle="--", color="k")
    axes1.set_ylim([0.5,1.5])
    axes1.set_ylabel(r'$\frac{BNB}{MC + EXT}$')
    axes1.set_xlabel(xlabel)
    axes1.set_xlim([BINEDGES[0],BINEDGES[-1]])
    
    
    fig.tight_layout()
    
    return fig, axes0, axes1
    
    
def apply_subgroups(df):
    #remember that these dataframes get passed around 
    #like c++ passes around pointers
    #TWO categorization heuristics:
    #1 NUMUCC,cosmic,OOFV,background
#    `outofvd` (for slices that have a neutrino interaction but that interaction has a nu vertex outside of the TPC)
#    `cosmic` (for slices that have a low neutrino purity)
#    `numu` for events that are true numus
#    `nu backgrounds` for events that are not numuCC.
    #cosmic
    #no true numu or nue and there is a reco cosmic
    #reco cosmic is if more than half the hits in the event are overlay
    df['COSMIC'] = df['category'] == 4
    #Fiducial Volume
    FVx = [AVx[0]+10,AVx[1]-10]
    FVy = [AVy[0]+10,AVy[1]-10]
    FVz = [AVz[0]+10,AVz[1]-50]
    xx = np.array(df['true_nu_vtx_x'])
    yy = np.array(df['true_nu_vtx_y'])
    zz = np.array(df['true_nu_vtx_z'])
    oofv = np.zeros(xx.size).astype(bool)
    oofv = oofv|((xx<FVx[0])|(xx>FVx[1]))
    oofv = oofv|((yy<FVy[0])|(yy>FVy[1]))
    oofv = oofv|((zz<FVz[0])|(zz>FVz[1]))
    #no double counting
    #oofv must be truly oofv and not a cosmic
    df['OOFV'] = oofv&(np.invert(df['COSMIC']))
    #numu
    nupdg = np.array(df['nu_pdg']) == 14
    ccnc = np.array(df['ccnc']) == 0
    mu = np.array(df['backtracked_pdg']) == 13
    df['NUMUCC'] = (nupdg&ccnc&mu)&(np.invert(df['COSMIC']|df['OOFV']))
    #it should be NUMUCC if it isn't cosmic or oofv
#    df['NUMUCC'] = df['NUMUCC']&(np.invert(df['COSMIC']&df['OOFV']))
    #group everything else together
    #if it's not NUMUCC, cosmic, or oofv, it's other
    df['OTHER'] = np.invert(df['NUMUCC']|df['COSMIC']|df['OOFV'])
    
    #2 interaction type (coming soon) look at categories branch

def apply_cuts(DFs):        
    for label in DFs:
        df = DFs[label]
        #################################
        ## CUTS
        #################################
        df['slice'] = df['nslice'] == 1
        df['vtxFV'] = (df['reco_nu_vtx_x']>FVx[0])&(df['reco_nu_vtx_x']<FVx[1])&(df['reco_nu_vtx_y']>FVy[0])&(df['reco_nu_vtx_y']<FVy[1])&(df['reco_nu_vtx_z']>FVz[0])&(df['reco_nu_vtx_z']<FVz[1])
        df['vtxFV1'] = (df['reco_nu_vtx_sce_x']>FVx1[0])&(df['reco_nu_vtx_sce_x']<FVx1[1])&(df['reco_nu_vtx_sce_y']>FVy1[0])&(df['reco_nu_vtx_sce_y']<FVy1[1])&(df['reco_nu_vtx_sce_z']>FVz1[0])&(df['reco_nu_vtx_sce_z']<FVz1[1])
        #topo score
        df['notopo'] = df['topological_score'] > 0.00
        df['topo06'] = df['topological_score'] > 0.06
        df['topo10'] = df['topological_score'] > 0.10
        df['topo15'] = df['topological_score'] > 0.15
        df['topo25'] = df['topological_score'] > 0.25
        ###########################################
        #has a muon
        chipr_cut = 60
        chimu_cut = 30
        chirat_cut = 7
        dist_cut = 4
        tscore_cut = 0.8
        chipr = np.array(df['trk_pid_chipr_v'])
        chimu = np.array(df['trk_pid_chimu_v'])
        chirat = np.array(df['chirat'])
        tscore = np.array(df['trk_score_v'])
        dist = np.array(df['trk_distance'])
        #Truth*Truth=Truth, Truth*False=False, False*False=False
        muon = chipr > chipr_cut
        muon &= chimu < chimu_cut
        muon &= chirat > chirat_cut
        muon &= tscore > tscore_cut
        muon &= dist < dist_cut
        df['muon'] = muon
        df['tscore60'] = df['trk_score_v'] > 0.6
        df['tscore70'] = df['trk_score_v'] > 0.7
        df['tscore80'] = df['trk_score_v'] > 0.8
        df['tdist4'] = df['trk_distance'] < 4
        df['tdist3'] = df['trk_distance'] < 3
        df['tdist2'] = df['trk_distance'] < 2
        #sort by trk length (pick longest one later if need be)
        df['longest'] = df['trk_len_v'].groupby("entry").transform(max) == df['trk_len_v']
        
    #    df['muon'] = (chipr>chipr_cut)*(chimu<chimu_cut)*((chipr/chimu)>chirat_cut)*(tscore>tscore_cut)*(dist<dist_cut)
        #########################################
        #muon is contained
        startx = np.array(df['trk_start_x_v'])
        starty = np.array(df['trk_start_y_v'])
        startz = np.array(df['trk_start_z_v'])
        endx = np.array(df['trk_end_x_v'])
        endy = np.array(df['trk_end_y_v'])
        endz = np.array(df['trk_end_z_v'])
        #containment volume
        CD = 10
        CVx = [AVx[0]+CD,AVx[1]-CD]
        CVy = [AVy[0]+CD,AVy[1]-CD]
        CVz = [AVz[0]+CD,AVz[1]-CD]
        contained = np.ones(startx.size).astype(bool)
        contained &= (startx>CVx[0])&(startx<CVx[1])
        contained &= (starty>CVy[0])&(starty<CVy[1])
        contained &= (startz>CVz[0])&(startz<CVz[1])
        contained &= (endx>CVx[0])&(endx<CVx[1])
        contained &= (endy>CVy[0])&(endy<CVy[1])
        contained &= (endz>CVz[0])&(endz<CVz[1])
        df['contained10'] = contained
        #containment volume
        CD = 25
        CVx = [AVx[0]+CD,AVx[1]-CD]
        CVy = [AVy[0]+CD,AVy[1]-CD]
        CVz = [AVz[0]+CD,AVz[1]-CD]
        contained = np.ones(startx.size).astype(bool)
        contained &= (startx>CVx[0])&(startx<CVx[1])
        contained &= (starty>CVy[0])&(starty<CVy[1])
        contained &= (startz>CVz[0])&(startz<CVz[1])
        contained &= (endx>CVx[0])&(endx<CVx[1])
        contained &= (endy>CVy[0])&(endy<CVy[1])
        contained &= (endz>CVz[0])&(endz<CVz[1])
        df['contained25'] = contained
        #############################
        ## experimental y cut
        df['vtx_y'] = (df['reco_nu_vtx_y'] > (AVy[0] + 25))&(df['reco_nu_vtx_y'] < (AVy[1] - 25))
        ############################
        ## fiducial SCE corrected vertex
        
        ###############################
        #length cut
        df['len20'] = df['trk_len_v'] >= 20
        df['len30'] = df['trk_len_v'] >= 30
        df['len40'] = df['trk_len_v'] >= 40
        df['crtveto'] = df['crtveto'] == 0

    print("list of cuts:")
    print(['slice','notopo','topo06','topo15','topo25','muon','contained10','contained25','len','crtveto'])
        
def apply_calcols(DFs):
    for label in DFs:
        df = DFs[label]
        df['costheta'] = np.cos(df['trk_theta_v'])
        chimu = np.array(df['trk_pid_chimu_v'])
        df['trk_pid_chimu_v'] = np.where(chimu==np.inf,99999,chimu)
        df['trk_pid_chimu_v'] = np.where(chimu==-np.inf,-99999,chimu)
        df['chirat'] = df['trk_pid_chipr_v']/df['trk_pid_chimu_v']
        df['longest'] = df['trk_len_v'].groupby("entry").transform(max) == df['trk_len_v']
        
        #there seems ot be a 1 cm offset in the SCE correction in x direction
#        print("correcting reco_nu_vtx_sce_x by 1 cm")
#        df['reco_nu_vtx_sce_x'] = df['reco_nu_vtx_sce_x'] - 1
    
def unique_entries(df,cuts=False,binmin=-99,binmax=99):
    if cuts != False:
        return len(np.unique(df.query(cuts).index.codes[0]))
    else:
        return len(np.unique(df.index.codes[0]))
    
def Eff(df,var,query,acceptance,bin_edges,absval=False,remvecs=False):
#     print(acceptance)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    bins = []
    bin_eff = []
    bin_err = []
    FVy0 = FVy[0]
    
    cum_num = 0
    cum_denom = 0
    
    print(len(np.unique(df.query(query+' and '+acceptance).index.codes[0])))
    for i in range(len(bin_centers)):
        binmin = bin_edges[i]
        binmax = bin_edges[i+1]
        bincut = '{} > @binmin and {} < @binmax'.format(var,var)
        if (absval == True):
            bincut = '({} > {} and {} < {}) or ({} > -{} and {} < -{})'.format(var,binmin,var,binmax,var,binmax,var,binmin)
        if (acceptance != ''): bincut += ' and {}'.format(acceptance)
        df_tmp =  df.query(bincut) # cut on bin range for desired var.
        df_sub = df_tmp.query(query) # apply constraint 
        #print("query/bincut: {}/{}".format(query,bincut))
        if (df_tmp.shape[0] == 0): continue
        #count number of unique entry values
        nentries_sub = len(np.unique(df_sub.index.codes[0]))
        nentries_tmp = len(np.unique(df_tmp.index.codes[0]))
        #print("num/denom: {}/{}={}".format(nentries_sub,nentries_tmp,nentries_sub/nentries_tmp))
        cum_num += nentries_sub
        cum_denom += nentries_tmp
        #print("cumulative num/denom: {}/{}".format(cum_num,cum_denom))
        eff = nentries_sub / float(nentries_tmp)
        #print(eff)
        err = np.sqrt(eff*(1-eff)/nentries_tmp)
        bin_eff.append(eff)
        bin_err.append(err)
        bins.append(bin_centers[i])
        #print 'eff = %.02f @ bin = %.02f'%(eff,bin_centers[i])
    return np.array(bins),np.array(bin_eff),np.array(bin_err)

def FV_Eff(df,acceptance,bin_edges,absval=False,remvecs=False):
#     print(acceptance)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    bins = []
    bin_eff = []
    bin_err = []
    
    for i in range(len(bin_centers)):
        fvcut = bin_centers[i]
        
        df_den = df.query(acceptance)
        fvquery = 'reco_nu_vtx_sce_x > {} and reco_nu_vtx_sce_x < {}'.format(AVx[0]+fvcut,AVx[1]-fvcut)
        fvquery += ' and reco_nu_vtx_sce_y > {} and reco_nu_vtx_sce_y < {}'.format(AVy[0]+fvcut,AVy[1]-fvcut)
        fvquery += ' and reco_nu_vtx_sce_z > {} and reco_nu_vtx_sce_z < {}'.format(AVz[0]+fvcut,AVz[1]-50)
        df_num = df_den.query(fvquery)

        nentries_num = len(np.unique(df_num.index.codes[0]))
        nentries_den = len(np.unique(df_den.index.codes[0]))
        eff = nentries_num / float(nentries_den)
        #print(eff)
        err = np.sqrt(eff*(1-eff)/nentries_den)
        bin_eff.append(eff)
        bin_err.append(err)
        bins.append(bin_centers[i])
        #print 'eff = %.02f @ bin = %.02f'%(eff,bin_centers[i])
    return np.array(bins),np.array(bin_eff),np.array(bin_err)

def Pur(DFs, ACCEPTANCE, var, query, bin_edges):
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    bins = []
    bins_pur = []
    bins_err = []
    df_nu = DFs['NU']
    df_ext = DFs['EXT']
    df_dirt = DFs['DIRT']
    for i in range(len(bin_centers)):
        binmin = bin_edges[i]
        binmax = bin_edges[i+1]
        bincut = '{} > @binmin and {} < @binmax'.format(var,var)
        num_query = ACCEPTANCE+' and '+bincut+' and '+query+' and longest'
        den_query = bincut+' and '+query+' and longest'
        num_scale = np.sum(df_nu.query(num_query)['weightSpline']*scale_NU-1)
        den_scale_nu = np.sum(df_nu.query(den_query)['weightSpline']*scale_NU-1) #only for overlay sample
        den_scale_dirt = np.sum(df_dirt.query(den_query)['weightSpline']*scale_DIRT-1)
        try:
            pur = (unique_entries(df_nu,num_query,binmin,binmax)+num_scale) / ((unique_entries(df_nu,den_query,binmin,binmax)+den_scale_nu) + (unique_entries(df_dirt,den_query,binmin,binmax)+den_scale_dirt) + (unique_entries(df_ext,den_query,binmin,binmax)*scale_EXT))
            #pur = (unique_entries(df_nu,num_query,binmin,binmax)*scale_NU) / ((unique_entries(df_nu,den_query,binmin,binmax)*scale_NU) + (unique_entries(df_dirt,den_query,binmin,binmax)*scale_DIRT) + (unique_entries(df_ext,den_query,binmin,binmax)*scale_EXT))
        except: 
            pur = 0
        bins_pur.append(pur)
        bins.append(bin_centers[i])
    return np.array(bins), np.array(bins_pur)

def overallPur(DFs,cuts):
    df_nu = DFs['NU']
    df_dirt = DFs['DIRT']
    df_ext = DFs['EXT']
    pur = unique_entries(df_nu,cuts + ' and nu_pdg==14 and ccnc==0 and OOFV==False and longest')*scale_NU / (unique_entries(df_nu,cuts)*scale_NU + unique_entries(df_dirt,cuts)*scale_DIRT + unique_entries(df_ext,cuts)*scale_EXT)
    return pur
