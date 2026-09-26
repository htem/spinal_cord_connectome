# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 13:27:15 2025

@author: Adrian Holtrup
"""

import numpy as np
import pickle
import cloudvolume as cv
import trimesh

#%%

SEG_ID = 720575940873438782
VISUALIZE = False
ROOT_PATH = '/Users/wangchu/connectivity_analysis/meshes'

vol = cv.CloudVolume(
    'graphene://https://cave.fanc-fly.com/segmentation/table/wclee_mouse_spinalcord_cltmr',
    secrets = 'REPLACE_WITH_YOUR_CAVE_TOKEN', 
    use_https=True, 
    progress=False
    )
mesh_raw = vol.mesh.get(SEG_ID)[SEG_ID]

if VISUALIZE: 
    
    import navis
    import plotly.io as pio
    pio.renderers.default='browser'
    
    mesh_as_navis = navis.MeshNeuron((mesh_raw.vertices, mesh_raw.faces))
    fig = navis.plot3d(mesh_as_navis, backend='plotly')
    fig.show()  

mesh = trimesh.base.Trimesh(vertices=mesh_raw.vertices, faces=mesh_raw.faces)

exp_ply = trimesh.exchange.ply.export_ply(mesh)

path = (
    ROOT_PATH
    + '/'
    + str(SEG_ID)
    + '.ply'
    )

with open(path, 'wb') as f:
    f.write(exp_ply)

