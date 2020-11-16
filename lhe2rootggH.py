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
  g.add_argument("--zh", action="store_true")
  g.add_argument("--wh", action="store_true")
  g.add_argument("--ggH4l", action="store_true") # for ggH 4l JHUGen
  g.add_argument("--ggH4lMG", action="store_true") # for ggH4l Madgraph with weights
  parser.add_argument("--use-flavor", action="store_true")
  parser.add_argument("--CJLST", action="store_true")
  parser.add_argument("--reweight-to", choices="fa3-0.5")
  args = parser.parse_args()

  if os.path.exists(args.outputfile): raise IOError(args.outputfile+" already exists")
  for _ in args.inputfile:
    if not os.path.exists(_) and not args.CJLST: raise IOError(_+" doesn't exist")

from array import array
import itertools

import ROOT

from lhefile import LHEFile_JHUGenVBFVH, LHEFile_Hwithdecay, LHEFile_Offshell4l
from mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar

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
      if ( abs(id) in (11, 12, 13, 14, 15, 16) or abs(id) in (1, 2, 3, 4, 5)  ) and status == 1:
        daughters.append(line)
        flav4l = flav4l*abs(id)
      #if abs(id) in (0, 1, 2, 3, 4, 5, 21) and status == 1:
      #  associated.append(line)
      #if abs(id) in (0, 1, 2, 3, 4, 5, 21) and status == 1:
      #  associated.append(line)

    #if len(daughters) == 4:
      #print "flavour comp:",flav4l
    if len(daughters) != 4:
      raise ValueError("Wrong number of daughters (expected {}, found {})\n\n".format(4, len(daughters))+"\n".join(lines))
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
  #if args.zh or args.wh:
  #  branchnames_float += ("mV", "mVstar")
  if args.ggH4lMG:
    branchnames_float_array=("weights",)
    num_weights=30
  if args.vbf:
    branchnames_float += ("q2V1", "q2V2")
  branchnames_float += (
    "pg1", "pg4", "pg4pd", "pg1g4", "D0minus", "DCP", "DCP_old",
    "ptH", "pxH",  "pyH",  "pzH",  "EH","rapH","rapHJJ","decayMode","qfl1","qfl2","qfl1mom","qfl2mom",
    "pxj1", "pyj1", "pzj1", "Ej1",
    "pxj2", "pyj2", "pzj2", "Ej2",
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

  if args.zh:
    g4 = 0.144057
  if args.wh:
    g4 = 0.1236136
  if args.vbf:
    g4 = 0.297979
  if args.ggH4lMG:
    g4 = 0 ## PlaceHolder!!

  #fileclass = LHEEvent_Offshell4l

  for inputfile in args.inputfile:
    print inputfile
    with LHEFile_Hwithdecay(inputfile, isgen=args.use_flavor) as f:
      
      for i, event in enumerate(f):
        

        print "Processed", i, "events"
        
        if event.daughters == 0 : 
          continue

        #if i != 0 and i % 1000 == 0:
        #if i > 20000 : 
        #  break
        
        if args.zh:
          process = TVar.Had_ZH
        elif args.wh:
          process = TVar.Had_WH
        elif args.vbf:
          process = TVar.JJVBF


        #Print daughters info
        
        flav4l = 1; 
        for d in event.daughters: 
         # print d.first
          flav4l = flav4l*d.first
        '''  
        for ind,d in enumerate(event.associated,1): 
          #print d.first," ",d 
          if( ind == 1 ): branches["qfl1"][0] = d.first
          if( ind == 2 ): branches["qfl2"][0] = d.first
          
        for ind,d in enumerate(event.mothers,1): 
          #print d.first," ",d 
          if( ind == 1 ): branches["qfl1mom"][0] = d.first
          if( ind == 2 ): branches["qfl2mom"][0] = d.first
        branches["decayMode"][0] = flav4l
        
         
        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1 = 1
        branches["pg1"][0] = event.computeProdP()
        val0 = event.computeProdP()
        '''

        event.setProcess(TVar. HSMHiggs,TVar.JHUGen,TVar.ZZGG)
        event.ghz1 = 1
        #branches["pg1"][0] = event.computeProdDecP()

        #event.setProcess(TVar. HSMHiggs,TVar.MCFM,TVar.JJVBF)
        event.ghz2 = 0
        event.ghz1 = 0
        event.ghz4 = 0 #0.29
        
        #branches["pg4pd"][0] = event.computeProdDecP()
        #branches["pg4"][0] = event.computeProdP()
        #branches["pg4pd"][0] = event.computeProdDecP()
        branches["pg4"][0] = event.computeP()

        '''
        event.setProcess(TVar. HSMHiggs,TVar.MCFM,TVar.JJVBF)
        event.ghz4 = 1
        val2 =  event.computeProdDecP()

        event.setProcess(TVar. HSMHiggs,TVar.MCFM,TVar.JJVBF)
        event.ghz1 = 1
        event.ghz2 = 0
        event.ghz4 = 0
        val5 =  event.computeProdDecP()
        '''

        #print "weights :",val1," ",val5," ",val2,"  ",branches["pg4pd"][0]
        

  

        #event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        #event.ghz4 = g4
      
#        event.setProcess(TVar.SelfDefine_spin0, TVar.MCFM, TVar.JJEW)
 #       event.ghz1 = 1
 #       branches["pg1"][0] = event.computeProdDecP()

#        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1 = 1
        event.ghz4 = g4

        #branches["pg1g4"][0] = event.computeProdP() - branches["pg1"][0] - branches["pg4"][0]

        if branches["pg1"][0] == branches["pg4"][0] == 0:
        #  for _ in event.daughters: print _.first,; _.second.Print()
        #  for _ in event.associated: print _.first,; _.second.Print()
          continue

        #branches["D0minus"][0] = branches["pg1"][0] / (branches["pg1"][0] + branches["pg4"][0])
        #branches["DCP"][0] = branches["pg1g4"][0] / (2 * (branches["pg1"][0] * branches["pg4"][0]) ** 0.5)
        #branches["DCP_old"][0] = branches["pg1g4"][0] / (branches["pg1"][0] + branches["pg4"][0])

        if args.zh or args.wh:

          branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVHAngles(process)
          branches["mV"][0] = sum((particle.second for particle in event.associated), ROOT.TLorentzVector()).M()
          branches["mVstar"][0] = sum((particle.second for particle in itertools.chain(event.daughters, event.associated)), ROOT.TLorentzVector()).M()
        elif args.vbf:
          branches["q2V1"][0], branches["q2V2"][0], branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVBFAngles()
          branches["HJJpz"][0] = sum((particle.second for particle in itertools.chain(event.daughters, event.associated)), ROOT.TLorentzVector()).Pz()

          branches["M4L"][0], branches["MZ1"][0], branches["MZ2"][0], branches["costheta1d"][0],branches["costheta2d"][0], branches["Phid"][0], branches["costhetastard"][0], branches["Phi1d"][0]= event.computeDecayAngles()

	
	elif args.ggH4l or args.ggH4lMG:
	  branches["M4L"][0], branches["MZ1"][0], branches["MZ2"][0], branches["costheta1d"][0],branches["costheta2d"][0], branches["Phid"][0], branches["costhetastard"][0], branches["Phi1d"][0]= event.computeDecayAngles()

        pH = sum((particle.second for particle in event.daughters), ROOT.TLorentzVector())
        branches["ptH"][0] = pH.Pt()
        branches["pxH"][0] = pH.Px()
        branches["pyH"][0] = pH.Py()
        branches["pzH"][0] = pH.Pz()
        branches["EH"][0] = pH.E()
        branches["rapH"][0] = pH.Rapidity()
                
        
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
        



        '''
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
        '''
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
