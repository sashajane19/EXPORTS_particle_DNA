#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sasha Kramer
MBARI
Figure 3 in Kramer et al., 2025
"""
###Plotly treemaps for surface ASVs to different deep pools
import plotly.express as px
import os
import pandas as pd
import numpy as np

os.chdir('/~/18S/treemap/')

#NA surf to deep BULK
nabasv = pd.read_csv(
    'na_surf_asv_bac_bulk.csv', sep=',',
    names=['Values', 'Parents', 'Labels','Deep'])

values = nabasv['Values']
labels = nabasv['Labels']
parents = nabasv['Parents']
deep = nabasv['Deep']

condition = [deep=='None',deep!='None']
choices = ['A','B']
cond = np.select(condition,choices,default='D')

df = pd.DataFrame(
    dict(labels=labels, parents=parents, values=values,deep=deep,cond=cond))

fig = px.treemap(df, path=['parents','labels'], values='values',color='cond',
                  color_discrete_map={'A':'white','B':'chocolate','(?)':'white'})
fig.update_layout(margin = dict(t=1, l=1, r=1, b=1))
fig.update_traces(marker_line_color='grey')
fig.update_layout(uniformtext=dict(minsize=10, mode='hide'))
fig.update_traces(insidetextfont_family='Calibri')
fig.show()

nadasv = pd.read_csv(
    'na_surf_asv_din_bulk.csv', sep=',',
    names=['Values', 'Parents', 'Labels','Deep'])

values = nadasv['Values']
labels = nadasv['Labels']
parents = nadasv['Parents']
deep = nadasv['Deep']

condition = [deep=='None',deep!='None']
choices = ['A','B']
cond = np.select(condition,choices,default='D')

df = pd.DataFrame(
    dict(labels=labels, parents=parents, values=values,deep=deep,cond=cond))

fig = px.treemap(df, path=['parents','labels'], values='values',color='cond',
                  color_discrete_map={'A':'white','B':'salmon','(?)':'white'})
fig.update_layout(margin = dict(t=1, l=1, r=1, b=1))
fig.update_traces(marker_line_color='grey')
fig.update_layout(uniformtext=dict(minsize=10, mode='hide'))
fig.update_traces(insidetextfont_family='Calibri')
fig.show()

naoasv = pd.read_csv(
    'na_surf_asv_och_bulk.csv', sep=',',
    names=['Values', 'Parents', 'Labels','Deep'])

values = naoasv['Values']
labels = naoasv['Labels']
parents = naoasv['Parents']
deep = naoasv['Deep']

condition = [deep=='None',deep!='None']
choices = ['A','B']
cond = np.select(condition,choices,default='D')

df = pd.DataFrame(
    dict(labels=labels, parents=parents, values=values,deep=deep,cond=cond))

fig = px.treemap(df, path=['parents','labels'], values='values',color='cond',
                  color_discrete_map={'A':'white','B':'khaki','(?)':'white'})
fig.update_layout(margin = dict(t=1, l=1, r=1, b=1))
fig.update_traces(marker_line_color='grey')
fig.update_layout(uniformtext=dict(minsize=10, mode='hide'))
fig.update_traces(insidetextfont_family='Calibri')
fig.show()

nacasv = pd.read_csv(
    'na_surf_asv_chl_bulk.csv', sep=',',
    names=['Values', 'Parents', 'Labels','Deep'])

values = nacasv['Values']
labels = nacasv['Labels']
parents = nacasv['Parents']
deep = nacasv['Deep']

condition = [deep=='None',deep!='None']
choices = ['A','B']
cond = np.select(condition,choices,default='D')

df = pd.DataFrame(
    dict(labels=labels, parents=parents, values=values,deep=deep,cond=cond))

fig = px.treemap(df, path=['parents','labels'], values='values',color='cond',
                  color_discrete_map={'A':'white','B':'palegreen','(?)':'white'})
fig.update_layout(margin = dict(t=1, l=1, r=1, b=1))
fig.update_traces(marker_line_color='grey')
fig.update_layout(uniformtext=dict(minsize=10, mode='hide'))
fig.update_traces(insidetextfont_family='Calibri')
fig.show()

nahasv = pd.read_csv(
    'na_surf_asv_hac_bulk.csv', sep=',',
    names=['Values', 'Parents', 'Labels','Deep'])

values = nahasv['Values']
labels = nahasv['Labels']
parents = nahasv['Parents']
deep = nahasv['Deep']

condition = [deep=='None',deep!='None']
choices = ['A','B']
cond = np.select(condition,choices,default='D')

df = pd.DataFrame(
    dict(labels=labels, parents=parents, values=values,deep=deep,cond=cond))

fig = px.treemap(df, path=['parents','labels'], values='values',color='cond',
                  color_discrete_map={'A':'white','B':'royalblue','(?)':'white'})
fig.update_layout(margin = dict(t=1, l=1, r=1, b=1))
fig.update_traces(marker_line_color='grey')
fig.update_layout(uniformtext=dict(minsize=10, mode='hide'))
fig.update_traces(insidetextfont_family='Calibri')
fig.show()

nasasv = pd.read_csv(
    'na_surf_asv_sil_bulk.csv', sep=',',
    names=['Values', 'Parents', 'Labels','Deep'])

values = nasasv['Values']
labels = nasasv['Labels']
parents = nasasv['Parents']
deep = nasasv['Deep']

condition = [deep=='None',deep!='None']
choices = ['A','B']
cond = np.select(condition,choices,default='D')

df = pd.DataFrame(
    dict(labels=labels, parents=parents, values=values,deep=deep,cond=cond))

fig = px.treemap(df, path=['parents','labels'], values='values',color='cond',
                  color_discrete_map={'A':'white','B':'gold','(?)':'white'})
fig.update_layout(margin = dict(t=1, l=1, r=1, b=1))
fig.update_traces(marker_line_color='grey')
fig.update_layout(uniformtext=dict(minsize=10, mode='hide'))
fig.update_traces(insidetextfont_family='Calibri')
fig.show()

#Repeat for all replicates (NA, NP, each group, surf to deep, deep to surf, etc.)
