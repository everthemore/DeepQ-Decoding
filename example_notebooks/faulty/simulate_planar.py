import sys
import os
import time
import math
import copy

import planar_lattice
import perfect_matching

home = os.environ['HOME']

class Timer():
   def __enter__(self): self.start = time.time()
   def __exit__(self, *args): print(time.time() - self.start)


def squashMatching(size, error_type, matching):

   # Convert error_type to channel index
   channel= 0 if error_type=="X" else 1

   # Initialize an array of (2*size+1) [1]'s in sets of (2*size+1)
   flip_array=[[1]*(2*size+1) for _ in range(2*size+1)]

   #
   for [(pt,p0,p1),(qt,q0,q1)] in matching:
      flips = []


      if (p0 in [-1,size*2+1] and q0 in [-1,size*2+1]) or (p1 in [-1,size*2+1] and q1 in [-1,size*2+1]):
         flips+=[]

      else:
         s0=int(math.copysign(1,q0-p0))
         s1=int(math.copysign(1,q1-p1))

         range0=range(1,abs(q0-p0),2)
         range1=range(1,abs(q1-p1),2)
         for x in range1:
            flips+=[[p0,p1+s1*x]]
         for y in range0:
            flips+=[[p0+s0*y,q1]]

      for flip in flips:
         flip_array[flip[0]][flip[1]]*=-1

   return flip_array


###
###          RUN CODE
###
def run3Drandom(size=4, tSteps=5, p=0.05, pLie=0.00, timespace=[1,1], showTextArray=False):

   # The surface code lattice
   L =planar_lattice.PlanarLattice(size)
   # Parity lattice
   PL=planar_lattice.PlanarLattice3D(size)

   # For each time-slice
   for i in range(tSteps):
      # Apply random X&Z errors to the slice
      L.applyRandomErrors(p,p)
      # Measure the stabilizers (with fauly error rate)
      L.measurePlaquettes(pLie)
      L.measureStars(pLie)
      # Add the slice to the 3D lattice
      PL.addMeasurement(L)

   # Finally, add another layer with perfect measurements
   L.measurePlaquettes(0)
   L.measureStars(0)
   PL.addMeasurement(L)

   if showTextArray==True: L.showArrayText("errors","X")

   # Find the errors
   PL.findAnyons()

   # Perfect matching on the plaquettes
   matchingX = perfect_matching.match_planar_3D(size,"plaquette",PL.anyon_positions_P,timespace)
   # Perfect matching on the starts
   matchingZ = perfect_matching.match_planar_3D(size,"star",PL.anyon_positions_S,timespace)

   # Reformat the matching
   flipsX = squashMatching(size,"X",matchingX)
   flipsZ = squashMatching(size,"Z",matchingZ)

   # Apply the correction...
   L.apply_flip_array("Z",flipsZ)
   L.apply_flip_array("X",flipsX)

   # ...and return the result of a logical error measurement
   return L.measure_logical()
