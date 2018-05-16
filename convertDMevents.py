#!/usr/bin/env python
#############################
# Author: Martina Ferrillo  #
# Creation date: April 2018 #
#############################
import os
import sys
import getopt
import time
import lhe_parser
import ROOT
import logging
from array import array
import argparse

DMpdgcode = 9000007
# DMpdgcode = 5000521
# DMpdgcode = 22008   # Depending on the provided DM decay model
INITIAL_STATE = -1
FINAL_STATE = 1


def convertFile(input_dm_file, output_root_file='rootfile_dm.root',
                input_hadrons_file=None):
    '''Convert outputs of MadDump in LHE files format to
    ROOT files, hadronizing on the way. The input file should
    conatain DIS in name is case of DIS scattering.'''

    ESevent = True if "ES" in input_dm_file else False
    DISevent = True if "DIS" in input_dm_file else False
    assert ESevent != DISevent, "Both flags are {}".format(ESevent)
    assert DISevent == (input_hadrons_file is not None),\
        "DIS is {0} but hadron_file is {1}".format(DISevent, input_hadrons_file)

    lhe = lhe_parser.EventFile(input_dm_file)

    rootfile_dm = ROOT.TFile(output_root_file, "recreate")
    tree_dm = ROOT.TTree("dmtree", "dmtree")

    Edm = array('d', [0])
    pxdm = array('d', [0])
    pydm = array('d', [0])
    pzdm = array('d', [0])
    dmpdg = array('i', [0])
    dis = array('i', [0])
    el = array('i', [0])
    E_2ry = array('d', [0]*500)
    px_2ry = array('d', [0]*500)
    py_2ry = array('d', [0]*500)
    pz_2ry = array('d', [0]*500)
    pdg_2ry = array('i', [0]*500)
    n_2ry = array('i', [0])

    tree_dm.Branch('Edm', Edm, 'Edm,/D')
    tree_dm.Branch('pxdm', pxdm, 'pxdm,/D')
    tree_dm.Branch('pydm', pydm, 'pydm,/D')
    tree_dm.Branch('pzdm', pzdm, 'pzdm,/D')
    tree_dm.Branch('dmpdg', dmpdg, 'dmpdg,/I')
    tree_dm.Branch('dis', dis, 'dis,/I')
    tree_dm.Branch('el', el, 'el,/I')
    tree_dm.Branch('E_2ry', E_2ry, 'E_2ry,/D')
    tree_dm.Branch('px_2ry', px_2ry, 'px_2ry,/D')
    tree_dm.Branch('py_2ry', py_2ry, 'py_2ry,/D')
    tree_dm.Branch('pz_2ry', pz_2ry, 'pz_2ry,/D')
    tree_dm.Branch('pdg_2ry', pdg_2ry, 'pdg_2ry,/I')
    tree_dm.Branch('n_2ry', n_2ry, 'n_2ry,/I')

    if DISevent:
        print 'input file pythia: ', input_hadrons_file
        P8gen = ROOT.TPythia8()
        pythiagen = P8gen.Pythia8()
        pythiagen.readString("Beams:frameType = 4")
        pythiagen.readString("Beams:LHEF = {}".format(input_hadrons_file))
        pythiagen.init()

    # Cycle over events and particles within each event
    for event in lhe:
        index = 0
        for particle in event:
            el[0] = ESevent
            dis[0] = DISevent

            if particle.status == INITIAL_STATE and particle.pdg == DMpdgcode:
                Edm[0] = particle.E
                pxdm[0] = particle.px
                pydm[0] = particle.py
                pzdm[0] = particle.pz
                dmpdg[0] = particle.pdg

            if ESevent:
                if particle.status == FINAL_STATE and particle.pdg != DMpdgcode:
                    E_2ry[index] = particle.E
                    px_2ry[index] = particle.px
                    py_2ry[index] = particle.py
                    pz_2ry[index] = particle.pz
                    pdg_2ry[index] = particle.pdg
                    index += 1

        if DISevent:
            if pythiagen.next():
                for i in range(pythiagen.event.size()):
                    if pythiagen.event[i].isFinal() and\
                      pythiagen.event[i].id() != DMpdgcode:
                        E_2ry[index] = pythiagen.event[i].e()
                        px_2ry[index] = pythiagen.event[i].px()
                        py_2ry[index] = pythiagen.event[i].py()
                        pz_2ry[index] = pythiagen.event[i].pz()
                        pdg_2ry[index] = pythiagen.event[i].id()
                        index += 1
        n_2ry[0] = index
        tree_dm.Fill()

    # Write the tree into the output file and close the file
    rootfile_dm.Write()
    rootfile_dm.Close()

    print "----------------------------------------------------------------"
    print "------------------------------------RootFile successfully saved!"
    print "----------------------------------------------------------------"
    return(rootfile_dm)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select events in neutrino files",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", dest='dm_file', type=str, help="Input DM file")
    parser.add_argument("-d", dest='h_hfile', type=str, help="Input hadron file",
                        default=None)
    parser.add_argument("-o", dest='o_hfile', type=str, help="Output file",
                        default='rootfile_dm.root')
    args = parser.parse_args()
    print(args.h_hfile)
    convertFile(args.dm_file, args.o_hfile, args.h_hfile)
