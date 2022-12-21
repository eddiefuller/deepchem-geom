#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 09:21:11 2022

@author: efuller
"""

import numpy as np
import matplotlib as mpl
import networkx as nx

from rdkit import Chem
import rdkit.Chem.rdmolfiles  as rdmol
from rdkit.Chem import AllChem


def get_geom( mol,
                 verb=False,
                 hydro=False,
                 align=False):
    
    if mol.GetNumAtoms() > 1:
        
        print("Computing atom curvatures and normals.")

        src=''
        kappa = []
        kappam = []
        theta = []
        dT = []
        n = []
        eta = []
        ad = []
        
        if hydro:
            mol = Chem.AddHs(mol)
            
        AllChem.EmbedMolecule(mol, useRandomCoords=True , randomSeed=0xf00d)
    
        if hydro:
            mol = Chem.RemoveHs(mol)
        
        try:
            
            conf = mol.GetConformer()
            if align == 'True':
                Chem.rdMolTransforms.CanonicalizeConformer(conf)
            #dists = Chem.rdmolops.GetDistanceMatrix(mol)
            
            
            #parse through list of nodes and compute 'curvature' at each node defined to be rotation angle 
            for i,atm in enumerate( mol.GetAtoms() ):
                ### check to see if it's pendant node - if so, handle differently
                if (i,) in  mol.GetSubstructMatches( Chem.MolFromSmarts('[D1]') ):
                    ### pendant set curvature = 0
                    if verb:
                        print("Pendant", atm.GetSymbol())
                    kappa.append([i,0])
                    n.append( [i,[0,0,0]] )
                    ad.append( [i,0] )
                else:
                    ### not pendant so compute
                    positions = conf.GetAtomPosition(i)
                    if verb:
                        print("Evaluating Atom: ", atm.GetIdx(), atm.GetSymbol(), positions.x, positions.y, positions.z)
                    #nbhd = enumerate(atm.GetNeighbors())
                    nbrs = atm.GetNeighbors()
                    
                    loop = len(nbrs)
                    localk = []
                    localdT = []
                    localN = []
                    alphas = []
                    #localIV = []
            
                    #print("Neighborhood size is:", loop)
                    if (loop>0):
                        for j in range(0,loop-1):
                            #posprev = conf.GetAtomPosition(atoms.GetIdx())
                            #posnext = conf.GetAtomPosition(atoms.GetIdx())
                           
                            #print("hi")
                            if verb:
                                print("-Prior Atom:", nbrs[j].GetIdx(), 
                                  conf.GetAtomPosition(nbrs[j].GetIdx()).x,
                                  conf.GetAtomPosition(nbrs[j].GetIdx()).y,
                                  conf.GetAtomPosition(nbrs[j].GetIdx()).z)
                            for l in range(j+1,loop):
                                if verb:
                                    print("---Post Atom:", nbrs[l].GetIdx(), 
                                      conf.GetAtomPosition(nbrs[l].GetIdx()).x,
                                      conf.GetAtomPosition(nbrs[l].GetIdx()).y,
                                      conf.GetAtomPosition(nbrs[l].GetIdx()).z)
                                tempdTpre = [positions.x - conf.GetAtomPosition(nbrs[j].GetIdx()).x,
                                          positions.y - conf.GetAtomPosition(nbrs[j].GetIdx()).y,
                                          positions.z - conf.GetAtomPosition(nbrs[j].GetIdx()).z]
                                tempdTpost = [conf.GetAtomPosition(nbrs[l].GetIdx()).x - positions.x,
                                          conf.GetAtomPosition(nbrs[l].GetIdx()).y - positions.y,
                                          conf.GetAtomPosition(nbrs[l].GetIdx()).z - positions.z]
                                
                                localdT.append( np.array(np.array(tempdTpost)-np.array(tempdTpre)) )
                                
                                localN.append( np.cross(tempdTpre, tempdTpost) )
                                
                                alphas.append( 
                                    np.arccos( 
                                        np.dot(tempdTpost,tempdTpre)/
                                        (np.linalg.norm(tempdTpre)*np.linalg.norm(tempdTpost)) 
                                        )
                                    )
                                
                                print( 'Angle:', 
                                      np.arccos( 
                                          np.dot(tempdTpost,tempdTpre)/
                                          ( np.linalg.norm(tempdTpre)*np.linalg.norm(tempdTpost) )
                                          ) 
                                      )
                                #print('Tpre:',)
                                #localIV.append( )
                
                                localk.append( np.linalg.norm(np.array(tempdTpost)-np.array(tempdTpre)))
                        kappa.append([i,np.prod(localk)])
                        kappam.append(np.mean(localk))
                        dT.append( np.array(np.sum(localdT,axis=0)))
                        n.append( [i,np.sum( localN, axis=0)] )
                        ad.append( [i,2*np.pi - np.sum(alphas)] )
                    else:
                        if verb:
                            print("Skipped", atm.GetIdx(), ": no nbhd")
                        kappa.append([i,0])
                        n.append( [i,[0,0,0]] )
                        ad.append( [i,0] )
                            
                    #dT.append(conf.GetAtomPosition(i)-conf.GetAtomPosition(i-1))
                    #kappa.append(1)
            
            #print("dT Vectors:", dT)
            #print("Total Curvature:", np.sum(kappa))
            #print("Average: ", np.mean(kappa))
            
            #print("Total Mean:", np.sum(kappam))
            #print("Average Mean:", np.mean(kappam))
            
            d = dict(); 
            d['tcurve'] = kappa #np.sum(kappa)
            #d['avgcurve']   = np.mean(kappa[1])
            #d['avgdT']   = dT
            d['N'] = n
            d['ad'] = ad
            return d
        
        except:
            return 'error'
        
    else:
        print("Single Atom Molecule - Default Values Returned.")
        d = dict(); 
        d['tcurve'] = 0
        #d['avgcurve'] = 0
        #d['avgdT']   = 0
        d['N'] = np.array([0,0,0])
        d['ad'] = 0
        return d
        
