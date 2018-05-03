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
import ROOT as r
import logging
from array import array
r.gROOT.ProcessLine('#include "Pythia8/LesHouches.h"')

DMpdgcode = 9000007
#DMpdgcode = 5000521
#DMpdgcode = 22008   # Depending on the provided DM decay model
pythia_out = "DMpythiaevents.lhe"


def convertFile(inputFile):    # Convert LHE file's infos into a ROOT TTree
    ESevent = False
    DISevent = False
    if inputFile.find('ES') != -1:
        ESevent = True
        events_lhefile = inputFile
    elif inputFile.find('DIS') != -1:
        DISevent = True
        events_lhefile = runPythia(inputFile)
    lhe = lhe_parser.EventFile(events_lhefile)    # Load LHE file

    rootfile_dm = r.TFile("rootfile_dm.root", "recreate")
    tree_dm = r.TTree("dmtree", "dmtree")
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

    # Cycle over events and particles within each event
    for event in lhe:
        index = 0
        for particle in event:
            if ESevent:
                el[0] = ESevent
            if DISevent:
                dis[0] = DISevent
            if particle.status == -1 and particle.pdg == DMpdgcode:
                Edm[0] = particle.E
                pxdm[0] = particle.px
                pydm[0] = particle.py
                pzdm[0] = particle.pz
                dmpdg[0] = particle.pdg
            elif particle.status == 1 and particle.pdg != DMpdgcode:
                E_2ry[index] = particle.E
                #print "energy of particle # ",index, " is", E_2ry[index]
                px_2ry[index] = particle.px
                #print "px of particle # ",index, " is ", px_2ry[index]
                py_2ry[index] = particle.py
                #print "py of particle # ",index, " is ", py_2ry[index]
                pz_2ry[index] = particle.pz
                #print "pz of particle # ",index, " is ", pz_2ry[index]
                pdg_2ry[index] = particle.pdg
                #print "pdg of particle # ",index, " is ", pdg_2ry[index]
                index += 1
                #print '--------------------------------------------------'
        n_2ry[0] = index
        tree_dm.Fill()

    # Write the tree into the output file and close the file
    rootfile_dm.Write()
    rootfile_dm.Close()

    print "----------------------------------------------------------------"
    print "------------------------------------RootFile successfully saved!"
    print "----------------------------------------------------------------"
    return(rootfile_dm)


# -----------------------------------------------------------
# Provide hadronization of free quarks in DIS events with Pythia8
def runPythia(inputFile):

    P8gen = r.TPythia8()
    pythia = P8gen.Pythia8()
    pythia.readString("LesHouches:matchInOut = off")
    pythia.readString("Beams:frameType = 4")
    pythia.readString("Beams:LHEF = "+str(inputFile))

    newLHA = r.Pythia8.LHAupFromPYTHIA8(pythia.process, pythia.info)
    newLHA.openLHEF(pythia_out)
#   pythia.readString("SoftQCD:nonDiffractive = on")
#   pythia.readString("HadronLevel:Hadronize = on")
#   pythia.readString("ProcessLevel:all = off")
#   pythia.readString(str(DMpdgcode)+":mayDecay = false")
#   print "----Pythia8 configuration: Made DM particle stable for Pythia"
    pythia.init()
    newLHA.setInit()    # Store initialization info in the LHAup object.
    newLHA.initLHEF()    # Write out this initialization info on the file.

    while not pythia.info.atEndOfFile():
        pythia.next()    # Generate events, and check whether generation failed
        newLHA.setEvent()    # Store event info in the LHAup object
        newLHA.eventLHEF()    # Write out this event info on the file

    pythia.stat()
    newLHA.updateSigma()    # Update cross section info based on MC integration during run
    newLHA.closeLHEF(True)

    print '---------------------------end pythia'
    return(pythia_out)

# ------------------------------------------------------------
