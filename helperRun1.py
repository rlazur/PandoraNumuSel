#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 13:42:05 2019

@author: ryan
"""
import numpy as np
import matplotlib.pyplot as plt
#plt.ioff()
plt.rcParams.update({'font.size': 22})
from matplotlib import gridspec
from datetime import datetime

start_time = datetime.now()

#see the POTnormalizeation.pynb script to derive these
scale_NU = 0.038257289303419886 #tor875 / POT
scale_EXT = 0.30488980130260335 #E1DCNT / EXT
scale_DIRT = 0.17585261287706153 #tor875 / POT

AVx = [-1.55,254.8]
AVy = [-115.53, 117.47]
AVz = [0.1, 1036.9]
FVx = [AVx[0]+10,AVx[1]-10]
FVy = [AVy[0]+10,AVy[1]-10]
FVz = [AVz[0]+10,AVz[1]-50]

#how far from the edges the containment cut is
CD = 25 #cm

CATS = ['COSMIC','OOFV','NUCC','OTHER']
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
        query = '{} >= @binmin and {} < @binmax'.format(var,var)
        df_bin = df.query(query)
        if 'weightSpline' in df.keys():
            bin_weights[i] = np.array(df_bin['weightSpline'])*scale
        else:
            bin_weights[i] = np.ones(df_bin.shape[0])*scale
        #square weights
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
        df_cut = df.query('{} >= @B0 and {} < @B1'.format(VAR,VAR))[[VAR,'weightSpline']]
        weights = np.array(df_cut['weightSpline'])*scale
    else:
        df_cut = df.query('{} >= @B0 and {} < @B1'.format(VAR,VAR))[[VAR]]
        weights = np.ones(df_cut.shape[0])*scale
    xstacked = list(df_cut[VAR])
    return list(xstacked),list(weights)

def fig_names(var,cuts,version,time_stamp=False):
    '''
    returns an appropriate save and title name given input
    can optionally put a timestamp on the name
    '''
    if cuts != False:
        cuts = cuts.replace('slice and topo06 and muon','basic')
        if time_stamp:
            ts = str(start_time.month)+str(start_time.day)+str(start_time.year)+'_'+str(start_time.hour)+str(start_time.minute)+str(start_time.second)
        else:
            ts = ''
        if version.upper() == "FINAL" or version.upper() == "OFFICIAL":
            title_name = "{}, {} cuts, Data-MC Comparison".format(var,cuts)
            save_name = "/home/ryan/Pictures/HEPPA/MicroBooNE/eLEE/{}_{}_dataMCcomp{}{}.pdf".format(var,cuts,ts,version)
        else:
            title_name = "{}, {} cuts, Data-MC Comparison {}".format(var,cuts,version)
            save_name = "plots/{}_{}_dataMCcomp{}{}.pdf".format(var,cuts,ts,version)
        
    else:
        if time_stamp:
            ts = str(start_time.month)+str(start_time.day)+str(start_time.year)+'_'+str(start_time.hour)+str(start_time.minute)+str(start_time.second)
        else:
            ts = ''
        if version.upper() == "FINAL" or version.upper() == "OFFICIAL":
            title_name = "{}, {} cuts, Data-MC Comparison".format(var,"No")
            save_name = "/home/ryan/Pictures/HEPPA/MicroBooNE/eLEE/{}_{}_dataMCcomp{}{}.pdf".format(var,"nocuts",ts,version)
        else:
            title_name = "{}, {} cuts, Data-MC Comparison {}".format(var,"No",version)
            save_name = "plots/{}_{}_dataMCcomp{}{}.pdf".format(var,"nocuts",ts,version)

    return save_name,title_name

def comp_DATAMC(DFs, VAR, BINEDGES, lw=25, take_longest=True, cuts=False, version=''):
    '''
    This function builds+saves a data-MC comparison histogram
    It needs a dictionary containing all the dataframes made using uproot
    FEATURES:
        option to plot all tracks or just the longest in each event/entry (take_longest)
        option to not implement any cuts at all
        
    NEED TO DO:
        implement spline weights
        
    '''
    BINCENTERS = 0.5*(BINEDGES[1:]+BINEDGES[:-1])
    hist_dict = {'xstacked': [], #list of lists of data correlated with WEIGHTS and LABELS (lists)
                 'ystacked_vals': np.zeros(BINCENTERS.size), #height of MC data in each bin
                 'LABELS': [], #labels for each component of xstacked (strings)
                 'WEIGHTS': [], #weights of each component of xstacked (lists)
                 'ERRS': np.zeros(BINCENTERS.size) #to calculate scaled MC uncertainty on each bin
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
            if label == 'EXT':
                scale = scale_EXT
            if label == 'NU':
                scale = scale_NU
            if CATS[0] in DFs[label].keys():
            #further categorization of neutrino MC sample possible
            #this categorization is defined by a list defined globally 
            #this case gets treated a bit differently then the above cases
                for cat in CATS:
                    try:
                        print('prepping {}...'.format(cat))
                        df_cat = df.query(cat)
                        #put the cat in the histogram prepper
                        x,w = prep_hist(df_cat,scale,VAR,BINEDGES)
                        errs = get_histerrs(df_cat,VAR,BINEDGES,scale)
                        hist_dict['xstacked'].append(x)
                        hist_dict['WEIGHTS'].append(w)
                        hist_dict['ERRS'] += errs
                        hist_dict['LABELS'].append(cat)
                    except:
                        print('{} not found in {}.\n DANGER: POTENTIAL DOUBLE COUNTING OR UNDERCOUNTING'.format(cat,label))
            elif scale:
            #this will run if it's a dataset without categories
                print('prepping {}...'.format(label))
                x,w = prep_hist(df,scale,VAR,BINEDGES)
                errs = get_histerrs(df,VAR,BINEDGES,scale)
                hist_dict['xstacked'].append(x)
                hist_dict['WEIGHTS'].append(w)
                hist_dict['ERRS'] += errs
                hist_dict['LABELS'].append(label)       
            else:
                print("SCALE UNKNOWN FOR {} SAMPLE".format(label))
                
        elif label == 'DATA':
            #a list of the values of each bin defined by BINEDGES
            data_bin_vals = get_histvals(df,VAR,BINEDGES)
            data_errs = np.sqrt(np.array(data_bin_vals)) #counting error
    #################################
    ##DATA-MC PLOT
    fig = plt.figure(figsize=(13,13))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    save_name,title_name = fig_names(VAR,cuts,version,time_stamp=False)
    axes0 = plt.subplot(gs[0])
    cum_ystacked_vals,_,_ = axes0.hist(hist_dict['xstacked'],bins=BINEDGES,histtype='stepfilled',label=hist_dict['LABELS'],weights=hist_dict['WEIGHTS'],stacked=True,alpha=0.7)
    #bin vals is an array of arrays, each array gives the cumulative bin heights of each hist
    ystacked_vals = cum_ystacked_vals[-1]
    axes0.errorbar(BINCENTERS,data_bin_vals,yerr=data_errs,fmt='o',color='k',markersize=5,label='data')
    axes0.errorbar(BINCENTERS,ystacked_vals,yerr=hist_dict['ERRS'],lw=lw,fmt='_k',ecolor='grey',alpha=0.4,fillstyle='full')
    axes0.legend()
    axes0.set_title(title_name)
    ##################################3
    ###Ratio Plot
    ratio_bin_vals = data_bin_vals / ystacked_vals
    ratio_bin_errs = np.sqrt((hist_dict['ERRS']/ystacked_vals)**2 + (data_errs/data_bin_vals)**2)
    axes1 = plt.subplot(gs[1])
    axes1.errorbar(BINCENTERS,ratio_bin_vals,yerr=ratio_bin_errs,fmt='o',color='k',markersize=5)
    axes1.plot([BINEDGES[0],BINEDGES[-1]],[1,1],'--')
    axes1.set_ylim([0,2])
    axes1.set_ylabel('ratio')
    axes1.set_xlabel(VAR)
    
    fig.savefig(save_name, pad_inches=0, bbox_inches='tight')
    
def apply_subgroups(df):
    #remember that these dataframes get passed around 
    #like c++ passes around pointers
    #TWO categorization heuristics:
    #1 nuCC,cosmic,OOFV,background
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
    nupdg = abs(np.array(df['nu_pdg'])) == 14
    ccnc = np.array(df['ccnc']) == 0
    mu = abs(np.array(df['backtracked_pdg'])) == 13
    df['NUCC'] = (nupdg&ccnc&mu)&(np.invert(df['COSMIC']|df['OOFV']))
    #it should be nucc if it isn't cosmic or oofv
#    df['NUCC'] = df['NUCC']&(np.invert(df['COSMIC']&df['OOFV']))
    #group everything else together
    #if it's not nucc, cosmic, or oofv, it's other
    df['OTHER'] = np.invert(df['NUCC']|df['COSMIC']|df['OOFV'])
    
    #2 interaction type (coming soon) look at categories branch

def apply_cuts(DFs):        
    for label in DFs:
        df = DFs[label]
        #################################
        ## CUTS
        #################################
        df['slice'] = df['nslice'] == 1
        #topo score
        df['notopo'] = df['topological_score'] > 0.00
        df['topo06'] = df['topological_score'] > 0.06
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
        ###############################
        #length cut
        df['len20'] = df['trk_len_v'] >= 20
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
    
def unique_entries(df,cuts=False):
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
    for i in range(len(bin_centers)):
        binmin = bin_edges[i]
        binmax = bin_edges[i+1]
        bincut = '{} > {} and {} < {}'.format(var,binmin,var,binmax)
        if (absval == True):
            bincut = '({} > {} and {} < {}) or ({} > -{} and {} < -{})'.format(var,binmin,var,binmax,var,binmax,var,binmin)
        if (acceptance != ''): bincut += ' and {}'.format(acceptance)
        df_tmp =  df.query(bincut) # cut on bin range for desired var.
        df_sub = df_tmp.query(query) # apply constraint 
        if (df_tmp.shape[0] == 0): continue
        #count number of unique entry values
        nentries_sub = len(np.unique(df_sub.index.codes[0]))
        nentries_tmp = len(np.unique(df_tmp.index.codes[0]))
        eff = nentries_sub / float(nentries_tmp)
        err = np.sqrt(eff*(1-eff)/nentries_tmp)
        bin_eff.append(eff)
        bin_err.append(err)
        bins.append(bin_centers[i])
        #print 'eff = %.02f @ bin = %.02f'%(eff,bin_centers[i])
    return np.array(bins),np.array(bin_eff),np.array(bin_err)

def Pur(DFs, var, query, bin_edges):
    ACCEPTANCE = 'nu_pdg==14 and OOFV==False'
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
        bincut = '{} > {} and {} < {}'.format(var,binmin,var,binmax)
        num_query = ACCEPTANCE+' and '+bincut+' and '+query
        den_query = bincut+' and '+query
        pur = unique_entries(df_nu,num_query)*scale_NU / (unique_entries(df_nu,den_query)*scale_NU + unique_entries(df_dirt,den_query)*scale_DIRT + unique_entries(df_ext,den_query)*scale_EXT)
        bins_pur.append(pur)
        bins.append(bin_centers[i])
    return np.array(bins), np.array(bins_pur)