#!/usr/bin/env python
import abc
import collections
import re
import argparse, os

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("outputfile")
  parser.add_argument("inputfile", nargs="+")
  g = parser.add_mutually_exclusive_group(required=True)
  g.add_argument("--vbf", action="store_true")
  g.add_argument("--vbf_withdecay", action="store_true")
  g.add_argument("--zh", action="store_true")
  g.add_argument("--zh_withdecay", action="store_true")
  g.add_argument("--zh_lep", action="store_true")
  g.add_argument("--zh_lep_hawk", action="store_true")
  g.add_argument("--wh_withdecay", action="store_true")
  g.add_argument("--wh_lep", action="store_true")
  g.add_argument("--wh", action="store_true")
  g.add_argument("--ggH4l", action="store_true") # for ggH 4l JHUGen and prophecy  
  g.add_argument("--ggH4lMG", action="store_true") # for ggH4l Madgraph with weights
  parser.add_argument("--use-flavor", action="store_true")
  parser.add_argument("--merge_photon", action="store_true") # for ggH 4l JHUGen and prophecy
  parser.add_argument("--calc_prodprob", action="store_true")
  parser.add_argument("--calc_decayprob", action="store_true")
  parser.add_argument("--CJLST", action="store_true")
  parser.add_argument("--reweight-to", choices="fa3-0.5")
  args = parser.parse_args()

  if os.path.exists(args.outputfile): raise IOError(args.outputfile+" already exists")
  for _ in args.inputfile:
    if not os.path.exists(_) and not args.CJLST: raise IOError(_+" doesn't exist")

from array import array
import itertools

import ROOT

from lhefile import LHEFile_JHUGenVBFVH, LHEFile_Hwithdecay, LHEFile_VHHiggsdecay,LHEFile_HwithdecayOnly, LHEFile_Offshell4l,LHEFile_StableHiggs,LHEFile_StableHiggsZHHAWK
from mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar
from pythonmelautils import MultiDimensionalCppArray, SelfDParameter, SelfDCoupling


def tlvfromptetaphim(pt, eta, phi, m):
  result = ROOT.TLorentzVector()
  result.SetPtEtaPhiM(pt, eta, phi, m)
  return result

class Event(object):
  doneinit = False
  def __init__(self, mela, daughters, associated, mothers, isgen, weight=1):
    self.daughters = daughters
    self.associated = associated
    self.mothers = mothers
    self.isgen = isgen
    self.weight = weight
    self.mela = mela
    self.mela.setInputEvent(daughters, associated, mothers, isgen)
    self.doneinit = True

  def __getattr__(self, attr):
    return getattr(self.mela, attr)

  def __setattr__(self, attr, value):
    if self.doneinit:
      return setattr(self.mela, attr, value)
    return super(Event, self).__setattr__(attr, value)



class LHEEvent(object):
  __metaclass__ = abc.ABCMeta
  def __init__(self, event, isgen):
    lines = event.split("\n")

    self.weights = {}
    for line in lines:
      if "<wgt" not in line: continue
      match = re.match("<wgt id='(.*)'>([0-9+Ee.-]*)</wgt>", line)
      if match: self.weights[match.group(1)] = float(match.group(2))
    lines = [line for line in lines if not ("<" in line or ">" in line or not line.split("#")[0].strip())]
    nparticles, _, weight, _, _, _ = lines[0].split()

    nparticles = int(nparticles)
    self.weight = float(weight)
    if nparticles != len(lines)-1:
      raise ValueError("Wrong number of particles! Should be {}, have {}".format(nparticles, len(lines)-1))

    daughters, associated, mothers = (SimpleParticleCollection_t(_) for _ in self.extracteventparticles(lines[1:], isgen))
    if not list(mothers): mothers = None
    self.daughters, self.associated, self.mothers, self.isgen = self.inputevent = InputEvent(daughters, associated, mothers, isgen)

  @abc.abstractmethod
  def extracteventparticles(cls, lines, isgen): "has to be a classmethod that returns daughters, associated, mothers"

  def __iter__(self):
    return iter(self.inputevent)


class LHEEvent_Offshell4l(LHEEvent):
  @classmethod
  def extracteventparticles(cls, lines, isgen):
    daughters, mothers, associated = [], [], []

    for line in lines:
      id, status, mother1, mother2 = (int(_) for _ in line.split()[0:4])
      if (1 <= abs(id) <= 6 or abs(id) == 21) and not isgen:
        line = line.replace(str(id), "0", 1)  #replace the first instance of the jet id with 0, which means unknown jet
      if status == -1:
        mothers.append(line)

      if ( abs(id) in ( 25) or abs(id) in (25)  ) and status == 1:
        daughters.append(line)
        print "added daught"
        flav4l = flav4l*abs(id)

        #      if ( abs(id) in (11, 12, 13, 14, 15, 16) or abs(id) in (1, 2, 3, 4, 5)  ) and status == 1:
        #       daughters.append(line)
        #       flav4l = flav4l*abs(id)
      #if abs(id) in (0, 1, 2, 3, 4, 5, 21) and status == 1:
      #  associated.append(line)
      #if abs(id) in (0, 1, 2, 3, 4, 5, 21) and status == 1:
      #  associated.append(line)

    #if len(daughters) == 4:
      #print "flavour comp:",flav4l
    #if len(daughters) != 4:
    #  raise ValueError("Wrong number of daughters (expected {}, found {})\n\n".format(4, len(daughters))+"\n".join(lines))
    if cls.nassociatedparticles is not None and len(associated) != cls.nassociatedparticles:
      raise ValueError("Wrong number of associated particles (expected {}, found {})\n\n".format(cls.nassociatedparticles, len(associated))+"\n".join(lines))
    if len(mothers) != 2:
      raise ValueError("{} mothers in the event??\n\n".format(len(mothers))+"\n".join(lines))

    if not isgen: mothers = None
    return daughters, associated, mothers

  nassociatedparticles = None



bad = False


try:
  newf = ROOT.TFile(args.outputfile, "RECREATE")
  t = ROOT.TTree("tree", "tree")

  branchnames_float = "costheta1", "costheta2", "Phi1", "costhetastar", "Phi", "HJJpz","M4L","MZ1","MZ2","costheta1d","costheta2d","Phid","costhetastard","Phi1d"
  if args.calc_prodprob or args.calc_decayprob :
    branchnames_float += ("pg1", "pg4","pg2","pg1g2","pg1g4","pg2za","pg4za","pg1g2za","pg1g4za","D0minus","D0hplus", "DCP", "Dint","D0minus_za","D0hplus_za","Dint_za","DCP_za")
  if args.zh or args.wh or args.zh_withdecay or args.wh_withdecay or args.zh_lep or args.wh_lep or args.zh_lep_hawk:
    branchnames_float += ("mV", "mVstar",    "pxj1", "pyj1", "pzj1", "Ej1",
    "pxj2", "pyj2", "pzj2", "Ej2")
    
  if args.ggH4lMG:
    branchnames_float_array=("weights",)
    num_weights=30
  if args.vbf or args.vbf_withdecay:
    branchnames_float += ("q2V1", "q2V2","Dphijj")
  branchnames_float += (
    
    "ptH", "pxH",  "pyH",  "pzH",  "EH","rapH","rapHJJ","decayMode","qfl1","qfl2","qfl1mom","qfl2mom",
    "pxj1", "pyj1", "pzj1", "Ej1",
    "pxj2", "pyj2", "pzj2", "Ej2",
    "pxph1","pyph1","pzph1","Eph1",
    "weight",
    "ptdau1","pxdau1","pydau1","pzdau1","Edau1","flavdau1", 
    "ptdau2","pxdau2","pydau2","pzdau2","Edau2","flavdau2",
    "ptdau3","pxdau3","pydau3","pzdau3","Edau3","flavdau3",
    "ptdau4","pxdau4","pydau4","pzdau4","Edau4","flavdau4",
  )

  branchnames_int = ()

  branches = {name: array("f", [0]) for name in branchnames_float}
  branches_int = {name: array("i", [0]) for name in branchnames_int}
  if args.ggH4lMG:
    branches_float = {name: array("f",[0]*num_weights) for name in branchnames_float_array}
    print branches_float,branchnames_float_array,"here"
    for name in branches_float:
      t.Branch(name, branches_float[name], name+"[30]/F")
    branches.update(branches_float)
  for name in branchnames_float:
    t.Branch(name, branches[name], name+"/F")
  for name in branchnames_int:
    t.Branch(name, branches_int[names], name+"/I")

  branches.update(branches_int)
  g4 = 0
  if args.zh:
    g4 = 0.144057
  if args.wh:
    g4 = 0.1236136
  if args.vbf:
    g4 = 0.297979
  
  
  for inputfile in args.inputfile:
    print inputfile
    inputfclass = LHEFile_Hwithdecay(inputfile,isgen=args.use_flavor)
    if args.ggH4l : 
      inputfclass = LHEFile_HwithdecayOnly(inputfile,isgen=args.use_flavor)
    if args.vbf or args.zh or args.wh or args.zh_lep or args.wh_lep  :
      inputfclass = LHEFile_StableHiggs(inputfile,isgen=args.use_flavor)
    if args.zh_withdecay or args.wh_withdecay  :
      inputfclass = LHEFile_VHHiggsdecay(inputfile,isgen=args.use_flavor)
    if args.zh_lep_hawk :
      print ("Algorithm will automaticaly merge associated FSR photons to the leptons")
      inputfclass = LHEFile_StableHiggsZHHAWK(inputfile,isgen=args.use_flavor)

    #inputfclass = LHEFile_Hwithdecay(inputfile,isgen=args.use_flavor)
    
    with inputfclass  as f:
      for i, event in enumerate(f):
        print "here"
        #debugging purposes
        if i > 10000 : 
          break
        
        #if( i % 100 == 0): 
        #  print ("Processed", i, " events",end="\r")
        
        
	### Automatically detect Had or Lep associated for VH production###
        associated_flavor = 0
        for a in event.associated:
          associated_flavor=abs(a.first)
          print associated_flavor
        if associated_flavor in [1,2,3,4,5,6]:
          if args.zh or args.zh_withdecay:
            process = TVar.Had_ZH
          elif args.wh or args.wh_withdecay:
            process = TVar.Had_WH
        if associated_flavor in [11,12,13,14,15,16]:
          if args.zh or args.zh_withdecay:
            process = TVar.Lep_ZH
          elif args.wh or args.wh_withdecay:
            process = TVar.Lep_WH
        if args.vbf or args.vbf_withdecay:
          process = TVar.JJVBF
        if args.zh_lep or args.zh_lep_hawk :
          process = TVar.Lep_ZH
        if args.wh_lep :
          process = TVar.Lep_WH
        if args.ggH4l :
          process = TVar.ZZGG

        #print i,"---- 1st -----------------------------"
        #for ll in event.daughters :
        #  print ll.second.Pt(), ll.second.Eta()

        flav4l = 1;
        for d in event.daughters: 
          flav4l = flav4l*d.first
          #print "init :",d.first,d.second.Px(),d.second.Py(),d.second.Pz(),d.second.E()
        # this section merges the EW emitted photon to the closest lepton and updates the daugthers collection for Prophecy
        # and associated for HAWK ZH, before the calculation of angles and probabilities for decay by MELA. The photon is
        #still stored in  the associated collection for references for Prophecy. For HAWK it is deleted.  


        if (args.zh_lep_hawk)  :   

          for ipho, p in enumerate(event.associated):

            if not (  p.first == 22)  :
              continue
            k = p.second
            mindr = 9999
            newlep = ROOT.TLorentzVector()

            photonvector = ROOT.TLorentzVector()
            lepp = ROOT.TLorentzVector()
            photonvector.SetPtEtaPhiM(k.Pt(), k.Eta(),k.Phi(), k.M())
            #print "photon:",p.first,p.second.Px(),p.second.Py(),p.second.Pz(),p.second.E()
            i_newlep = 0        
            for ilep, d in enumerate(event.associated):
              if not( d.first == 22 )  : 
                lep = d.second
                d_lor = ROOT.TLorentzVector()            
                d_lor.SetPtEtaPhiM(lep.Pt(), lep.Eta(),lep.Phi(), lep.M())
                dr =  k.DeltaR(d_lor)
                if dr < mindr :
                  lepp.SetPtEtaPhiM(lep.Pt(), lep.Eta(),lep.Phi(), lep.M())
                  mindr = dr
                  i_newlep = ilep                  
            newlep = lepp  + photonvector 
            event.associated.pop_back()
            event.associated[i_newlep].second = newlep
            # re-initiate the event with the new associated particles  this only works for HAWK events!!!
            eventt = event
            event.mela.setInputEvent(eventt.daughters, eventt.associated, eventt.mothers)




        if (args.merge_photon and args.ggH4l)  :   

          for p in event.associated:

            if not (  p.first == 22)  :
              continue
            k = p.second
            mindr = 9999
            newlep = ROOT.TLorentzVector()

            photonvector = ROOT.TLorentzVector()
            lepp = ROOT.TLorentzVector()
            photonvector.SetPtEtaPhiM(k.Pt(), k.Eta(),k.Phi(), k.M())
            #print "photon:",p.first,p.second.Px(),p.second.Py(),p.second.Pz(),p.second.E()
            i_newlep = 0        
            for ilep, d in enumerate(event.daughters): 
              lep = d.second
              d_lor = ROOT.TLorentzVector()            
              d_lor.SetPtEtaPhiM(lep.Pt(), lep.Eta(),lep.Phi(), lep.M())
              dr =  k.DeltaR(d_lor)
              if dr < mindr :
                lepp.SetPtEtaPhiM(lep.Pt(), lep.Eta(),lep.Phi(), lep.M())
                mindr = dr
                i_newlep = ilep  

            newlep = lepp  + photonvector 
            event.daughters[i_newlep].second = newlep
            # re-initiate the event with the new daughters this only works for Prophecy events!!!
            eventt = event
            event.mela.setInputEvent(eventt.daughters, eventt.associated, eventt.mothers)
  
                
        #Probabilities
        #event.setProcess(TVar. HSMHiggs,TVar.JHUGen,process)
        if args.calc_decayprob : 
          # decayP works only for the process below 
          #everytime you call a compute Prob function all the couplings
          #are reset and have to be redefined. 

          process = TVar.ZZINDEPENDENT
          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz2 = 1
          branches["pg2"][0] = event.computeP()
                  
          event.setProcess(TVar.HSMHiggs,TVar.JHUGen,process)
          event.ghz1 = 2      
          branches["pg1"][0] = event.computeP()
        
          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz4 = 1
          branches["pg4"][0] = event.computeP()
                            
          event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
          event.ghz1 = 1
          event.ghz4 = 1
          branches["pg1g4"][0] = event.computeP() - branches["pg1"][0] - branches["pg4"][0]
        
          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz1 = 1
          event.ghz2 = 1        
          branches["pg1g2"][0] = event.computeP() - branches["pg1"][0] - branches["pg2"][0]


          

          
          c_0hplus = 1
          c_0minus = 1 
          if ( process ==  TVar.ZZINDEPENDENT  ) : 
            c_0minus = 2.55497301342
            c_0hplus = 1.66326995046


            
          branches["D0minus"][0] = branches["pg1"][0] / (branches["pg1"][0] + c_0minus*c_0minus*branches["pg4"][0])
          branches["D0hplus"][0] = branches["pg1"][0] / (branches["pg1"][0] + c_0hplus*c_0hplus*branches["pg2"][0])
          branches["DCP"][0] = branches["pg1g4"][0] / (2 * (branches["pg1"][0] * branches["pg4"][0]) ** 0.5)
          branches["Dint"][0] = branches["pg1g2"][0] / (2 * (branches["pg1"][0] * branches["pg2"][0]) ** 0.5)
          #branches["DCP_old"][0] = branches["pg1g4"][0] / (branches["pg1"][0] + branches["pg4"][0])

        if args.calc_prodprob :
          
          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz1 = 0
          event.ghz2 = 1
          branches["pg2"][0] = event.computeProdP()
          #print "g2 :",event.ghz2,branches["pg2"][0]
          
          
          event.setProcess(TVar. HSMHiggs,TVar.JHUGen,process)
          event.ghz1 = 1
          branches["pg1"][0] = event.computeProdP()
          #print "g1 :",branches["pg1"][0]

          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz4 = 1
          branches["pg4"][0] = event.computeProdP()
          #print "g4 :",branches["pg4"][0]

                    
          event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
          event.ghz1 = 1
          event.ghz4 = 1
          branches["pg1g4"][0] = event.computeProdP() - branches["pg1"][0] - branches["pg4"][0]
          #print "g1g4 :",branches["pg1g4"][0]

          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz1 = 1
          event.ghz2 = 1        
          branches["pg1g2"][0] = event.computeProdP() - branches["pg1"][0] - branches["pg2"][0]
          #print "g1g2 :",branches["pg1g2"][0]

          #added zgamma discr Prod
          '''
          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz1 = 0
          event.ghza2 = 1        
          branches["pg2za"][0] = event.computeProdP()

          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz1 = 0
          event.ghza4 = 1        
          branches["pg4za"][0] = event.computeProdP()
          

          
          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz1 = 1
          event.ghza2 = 1        
          branches["pg1g2za"][0] = event.computeProdP() - branches["pg1"][0] - branches["pg2za"][0]
          

          event.setProcess(TVar.SelfDefine_spin0,TVar.JHUGen,process)
          event.ghz1 = 1
          event.ghza4 = 1        
          branches["pg1g4za"][0] = event.computeProdP() - branches["pg1"][0] - branches["pg4za"][0]
          '''



          c_0hplus = 1
          c_0minus = 1 
          if ( process == TVar.Had_ZH ) :
            c_0hplus = 0.130395173298
            c_0minus = 0.104503154335
            c_0hplusza = 0.130395173298
            c_0minusza = 0.104503154335
          if ( process == TVar.JJVBF  ) :   
            c_0minus = 0.297979440554
            c_0hplus = 0.271880048944
          if ( process ==  TVar.ZZGG  ) : 
            c_0minus = 2.55497301342
            c_0hplus = 1.66326995046


            
          branches["D0minus"][0] = branches["pg1"][0] / (branches["pg1"][0] + c_0minus*c_0minus*branches["pg4"][0])
          branches["D0hplus"][0] = branches["pg1"][0] / (branches["pg1"][0] + c_0hplus*c_0hplus*branches["pg2"][0])
          branches["DCP"][0] = branches["pg1g4"][0] / (2 * (branches["pg1"][0] * branches["pg4"][0]) ** 0.5)
          branches["Dint"][0] = branches["pg1g2"][0] / (2 * (branches["pg1"][0] * branches["pg2"][0]) ** 0.5)

          '''
          branches["D0minus_za"][0] = branches["pg1"][0] / (branches["pg1"][0] + c_0minusza*c_0minusza*branches["pg4za"][0])
          branches["D0hplus_za"][0] = branches["pg1"][0] / (branches["pg1"][0] + c_0hplusza*c_0hplusza*branches["pg2za"][0])
          branches["DCP_za"][0] = branches["pg1g4za"][0] / (2 * (branches["pg1"][0] * branches["pg4za"][0]) ** 0.5)
          branches["Dint_za"][0] = branches["pg1g2za"][0] / (2 * (branches["pg1"][0] * branches["pg2za"][0]) ** 0.5)
          ''' 


          #branches["DCP_old"][0] = branches["pg1g4"][0] / (branches["pg1"][0] + branches["pg4"][0])
  





          
        if args.zh or args.wh or args.zh_lep or args.wh_lep or args.zh_lep_hawk:
          branches["mV"][0], branches["mVstar"][0], branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVHAngles(process)
        elif args.zh_withdecay or args.wh_withdecay :
          branches["mV"][0], branches["mVstar"][0], branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVHAngles(process)
          branches["M4L"][0], branches["MZ1"][0], branches["MZ2"][0], branches["costheta1d"][0],branches["costheta2d"][0], branches["Phid"][0], branches["costhetastard"][0], branches["Phi1d"][0]= event.computeDecayAngles()
          #branches["mV"][0] = sum((particle.second for particle in event.associated), ROOT.TLorentzVector()).M()
          #branches["mVstar"][0] = sum((particle.second for particle in itertools.chain(event.daughters, event.associated)), ROOT.TLorentzVector()).M()
        elif args.vbf:
          branches["q2V1"][0], branches["q2V2"][0], branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVBFAngles()
          branches["HJJpz"][0] = sum((particle.second for particle in itertools.chain(event.daughters, event.associated)), ROOT.TLorentzVector()).Pz()

          pj1 = event.associated[0].second
          pj2 = event.associated[1].second
          if pj1.Pt() > pj2.Pt() :
            branches["Dphijj"][0] = pj1.DeltaPhi(pj2)
          else:
            branches["Dphijj"][0] = pj2.DeltaPhi(pj1)
        elif args.vbf_withdecay:
          branches["q2V1"][0], branches["q2V2"][0], branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVBFAngles()
          branches["HJJpz"][0] = sum((particle.second for particle in itertools.chain(event.daughters, event.associated)), ROOT.TLorentzVector()).Pz()
          branches["M4L"][0], branches["MZ1"][0], branches["MZ2"][0], branches["costheta1d"][0],branches["costheta2d"][0], branches["Phid"][0], branches["costhetastard"][0], branches["Phi1d"][0]= event.computeDecayAngles()

	
	elif args.ggH4l or args.ggH4lMG:
	  branches["M4L"][0], branches["MZ1"][0], branches["MZ2"][0], branches["costheta1d"][0],branches["costheta2d"][0], branches["Phid"][0], branches["costhetastard"][0], branches["Phi1d"][0]= event.computeDecayAngles()
        #print i,branches["M4L"][0], branches["MZ1"][0], branches["MZ2"][0],len(event.associated)
        pH = sum((particle.second for particle in event.daughters), ROOT.TLorentzVector())
        if args.ggH4l:
          if( len(event.associated) > 0):
            ph1 = event.associated[0].second
            pHph = pH + ph1 
            branches["pxph1"][0] = pHph.Px()
            branches["pyph1"][0] = pHph.Py()
            branches["pzph1"][0] = pHph.Pz()
            branches["Eph1"][0]  = pHph.E()

        else:     
          branches["ptH"][0] = pH.Pt()
          branches["pxH"][0] = pH.Px()
          branches["pyH"][0] = pH.Py()
          branches["pzH"][0] = pH.Pz()
          branches["EH"][0] = pH.E()
          branches["rapH"][0] = pH.Rapidity()
                
        if not args.vbf and not args.zh and not args.wh and not args.zh_lep and not args.wh_lep and not args.zh_lep_hawk: 
          pdau1 = event.daughters[0].second
          branches["ptdau1"][0] = pdau1.Pt()
          branches["pxdau1"][0] = pdau1.Px()
          branches["pydau1"][0] = pdau1.Py()
          branches["pzdau1"][0] = pdau1.Pz()
          branches["Edau1"][0] = pdau1.E()
          branches["flavdau1"][0] = event.daughters[0].first
          
          
          pdau2 = event.daughters[1].second
          branches["ptdau2"][0] = pdau2.Pt()
          branches["pxdau2"][0] = pdau2.Px()
          branches["pydau2"][0] = pdau2.Py()
          branches["pzdau2"][0] = pdau2.Pz()
          branches["Edau2"][0] = pdau2.E()
          branches["flavdau2"][0] = event.daughters[1].first
          
          pdau3 = event.daughters[2].second
          branches["ptdau3"][0] = pdau3.Pt()
          branches["pxdau3"][0] = pdau3.Px()
          branches["pydau3"][0] = pdau3.Py()
          branches["pzdau3"][0] = pdau3.Pz()
          branches["Edau3"][0] = pdau3.E()
          branches["flavdau3"][0] = event.daughters[2].first
          
          
          pdau4 = event.daughters[3].second
          branches["ptdau4"][0] = pdau4.Pt()
          branches["pxdau4"][0] = pdau4.Px()
          branches["pydau4"][0] = pdau4.Py()
          branches["pzdau4"][0] = pdau4.Pz()
          branches["Edau4"][0] = pdau4.E()
          branches["flavdau4"][0] = event.daughters[3].first
          
        

        if args.ggH4l:
          #add FSR photon to the root file 
          if( len(event.associated) > 0):
            ph1 = event.associated[0].second
            branches["pxph1"][0] = ph1.Px()
            branches["pyph1"][0] = ph1.Py()
            branches["pzph1"][0] = ph1.Pz()
            branches["Eph1"][0] = ph1.E()

        if args.vbf or args.zh or args.wh  or args.zh_lep_hawk:     
          pj1 = event.associated[0].second
          branches["pxj1"][0] = pj1.Px()
          branches["pyj1"][0] = pj1.Py()
          branches["pzj1"][0] = pj1.Pz()
          branches["Ej1"][0] = pj1.E()          
          pj2 = event.associated[1].second
          branches["pxj2"][0] = pj2.Px()
          branches["pyj2"][0] = pj2.Py()
          branches["pzj2"][0] = pj2.Pz()
          branches["Ej2"][0] = pj2.E()
          phjj = pH + pj1 + pj2
          branches["rapHJJ"][0] = phjj.Rapidity()
        
	if args.ggH4lMG:
          key_sorted=sorted(event.weights)
          list_weights=[]
          for key in key_sorted:
            list_weights.append(event.weights.get(key,None))
          for j in range(len(list_weights)):
            branches["weights"][j] = list_weights[j]
          print branches["weights"]
	else:
          branches["weight"][0] = event.weight
              # print "FIlling!"
        t.Fill()
      print "Processed", i+1, "events"

  newf.Write()
except:
  bad = True
  raise
finally:
  if bad:
    try:
      os.remove(args.outputfile)
    except:
      pass
