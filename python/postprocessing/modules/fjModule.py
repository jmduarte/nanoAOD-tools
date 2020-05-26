import ROOT
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True
from array import array
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class pfcandProducer(Module):
    def __init__(self, jetSelection):
        self.jetSel = jetSelection
        self.Nparts = 40
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("fj_idx", "I", 1);
        self.out.branch("fj_msd", "F", 1);
        self.out.branch("fj_m", "F", 1);
	self.out.branch("fj_pt", "F", 1);
	self.out.branch("fj_eta", "F", 1);
 	self.out.branch("fj_phi", "F", 1);
 	self.out.branch("fj_n2b1", "F", 1);
 	self.out.branch("fj_tau21", "F", 1);
 	self.out.branch("fj_deepTagZqq", "F", 1);
 	self.out.branch("fj_deepTagWqq", "F", 1);
	self.out.branch("fj_nProngs","I",1);
        self.out.branch("PF_pt",  "F", self.Nparts);
        self.out.branch("PF_eta",  "F", self.Nparts);
        self.out.branch("PF_phi",  "F", self.Nparts);
        self.out.branch("PF_pup",  "F", self.Nparts);
        self.out.branch("PF_q",  "F", self.Nparts);
        self.out.branch("PF_trk",  "F", self.Nparts);
        self.out.branch("PF_dz",  "F", self.Nparts);

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        pfcands	 = Collection(event, "FatJetPFCands")
        jets 	 = Collection(event, "FatJet")
	genparts = Collection(event, "GenPart")

        for ij, jet in enumerate(jets):
 

            if jet.pt < 400 or jet.msoftdrop < 40 : continue

            ##Fill basic jet properties
            jpt    = jet.pt
            jeta   = jet.eta
            jphi   = jet.phi
            jmsd   = jet.msoftdrop 
            jm   = jet.mass
	    jn2b1  = jet.n2b1
	    jdeepTagZqq  = jet.deepTagZqq
	    jdeepTagWqq  = jet.deepTagWqq
	    jn2b1  = jet.n2b1
            try: jtau21 = float(jet.tau2)/float(jet.tau1)
            except: jtau21 = 0.



            ##Calculate # prongs by looping over daughters of Z'
            status = 0
            motheridx = None
            jetv  = ROOT.TLorentzVector()
            jetv.SetPtEtaPhiM(jpt, jeta, jphi, jmsd) 
            dau1 = ROOT.TLorentzVector()
            dau2 = ROOT.TLorentzVector()
	    for ig, gpart in enumerate(genparts):
                if gpart.pdgId == 55 and gpart.status > status:
                    motheridx = ig
                    status = gpart.status

	    for ig, gpart in enumerate(genparts):
                if type(motheridx) == None: break
                if gpart.genPartIdxMother == motheridx:
                    if   dau1.Pt() == 0.: dau1.SetPtEtaPhiM(gpart.pt, gpart.eta, gpart.phi, gpart.mass)
                    elif dau2.Pt() == 0.: dau2.SetPtEtaPhiM(gpart.pt, gpart.eta, gpart.phi, gpart.mass)
            
            nProngs = 0
            if dau1.Pt() > 0. and dau2.Pt() > 0.:
                if jetv.DeltaR(dau1) < 0.8: nProngs += 1
                if jetv.DeltaR(dau2) < 0.8: nProngs += 1
                #if nProngs == 2: break


            ##Fill PF candidates
            if ij == 0:  candrange = range(0, jet.nPFConstituents)
            elif ij > 0: candrange = range(jets[ij-1].nPFConstituents, jets[ij].nPFConstituents)
            pfpt  = np.zeros(self.Nparts, dtype = np.float16)
            pfeta = np.zeros(self.Nparts, dtype = np.float16)
            pfphi = np.zeros(self.Nparts, dtype = np.float16)
            pftrk = np.zeros(self.Nparts, dtype = np.float16)
            pfpup = np.zeros(self.Nparts, dtype = np.float16)
            pfq   = np.zeros(self.Nparts, dtype = np.float16)
            pfdz  = np.zeros(self.Nparts, dtype = np.float16)
            arrIdx = 0
            for ip, part in enumerate(pfcands):
                if ip not in candrange: continue
                if arrIdx == self.Nparts: break
                pfpt[arrIdx]  = part.pt*part.puppiWeight/jpt
                pfeta[arrIdx] = part.eta - jeta
                pfphi[arrIdx] = signedDeltaPhi(part.phi, jphi)
                pfpup[arrIdx] = part.puppiWeight
                pfq[arrIdx]   = part.charge
                pfdz[arrIdx]  = part.dz
                pftrk[arrIdx] = part.trkChi2 
                arrIdx += 1

	    self.out.fillBranch("fj_idx",ij)
	    self.out.fillBranch("fj_pt",jpt)
	    self.out.fillBranch("fj_eta",jeta)
            self.out.fillBranch("fj_phi",jphi)
	    self.out.fillBranch("fj_n2b1",jn2b1)
	    self.out.fillBranch("fj_tau21",jtau21)
	    self.out.fillBranch("fj_nProngs",nProngs)
	    self.out.fillBranch("fj_deepTagZqq",jdeepTagZqq)
	    self.out.fillBranch("fj_deepTagWqq",jdeepTagWqq)
	    self.out.fillBranch("fj_msd",jmsd)
	    self.out.fillBranch("fj_m",jm)
            self.out.fillBranch("PF_pt",pfpt)
            self.out.fillBranch("PF_eta",pfeta)
            self.out.fillBranch("PF_phi",pfphi)
            self.out.fillBranch("PF_pup",pfpup)
            self.out.fillBranch("PF_q",pfq)
            self.out.fillBranch("PF_dz",pfdz)
            self.out.fillBranch("PF_trk",pftrk)
            return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def signedDeltaPhi(phi1,phi2):
    dPhi = phi1 - phi2
    if(dPhi < -np.pi):
       dPhi = 2*np.pi+dPhi
    elif(dPhi > np.pi):
       dPhi = -2*np.pi+dPhi
    return dPhi

pfModule = lambda : pfcandProducer(jetSelection= lambda j : j.pt >= 400. and j.msoftdrop >= 40.) 
 
