#! /usr/bin/env python

import ROOT
import sys
from DataFormats.FWLite import Events, Handle

files = ["ttbsm_414_data.root"]
events = Events (files)
handle1  = Handle ("std::vector<pat::Jet>")
handle2  = Handle ("std::vector<pat::Jet>")

# for now, label is just a tuple of strings that is initialized just
# like and edm::InputTag
label1 = ("goodPatJetsCA8PrunedPF")
label2 = ("goodPatJetsCATopTagPF")

f = ROOT.TFile("ttbsm_unittest_fwlite.root", "RECREATE")
f.cd()


# loop over events
i = 0
for event in events:
    i = i + 1
    print  '--------- Processing Event ' + str(i)


    print '---- ' + label1
    # use getByLabel, just like in cmsRun
    event.getByLabel (label1, handle1)
    # get the product
    jets1 = handle1.product()

    ijet = 0
    for jet in jets1 :
        print 'Jet {0:4.0f}, pt = {1:6.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, nda = {5:3.0f}, ptda1 = {6:6.2f}, ptda1 = {7:6.2f}'.format(
            ijet, jet.pt(), jet.eta(), jet.phi(), jet.mass(), jet.numberOfDaughters(), jet.daughter(0).pt(), jet.daughter(1).pt()
            )
        ijet += 1


    print '---- ' + label2
    # use getByLabel, just like in cmsRun
    event.getByLabel (label2, handle2)
    # get the product
    jets2 = handle2.product()

    ijet = 0
    for jet in jets2 :
        print 'Jet {0:4.0f}, pt = {1:6.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, nda = {5:3.0f}, topmass = {6:6.2f}, minmass = {7:6.2f}'.format(
            ijet, jet.pt(), jet.eta(), jet.phi(), jet.mass(), jet.numberOfDaughters(), jet.tagInfo('CATop').properties().topMass, jet.tagInfo('CATop').properties().minMass
            )
        ijet += 1


f.cd()

f.Close()
