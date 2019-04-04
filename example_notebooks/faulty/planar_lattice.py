import random
import math
import copy
import itertools
import numpy as np

# uncomment to be able to use the .showArray() function
#import matplotlib.pyplot as plt


"""

INFO:
   Errors stored in an array, each entry represents either
   a qubit: [xError?,zError?] with xError,zError in {1,-1}
   a stabiliser: star or plaquette depending on location with S in {1,-1}


Data structure:
       top          top
  Q -- STAR -- Q -- STAR -- Q --
        |            |
 left   |            |
 PLAQ   Q    PLAQ    Q    PLAQ
        |            |
        |            |
  Q -- STAR -- Q -- STAR -- Q --
        |            |
 left   |            |
 PLAQ   Q    PLAQ    Q    PLAQ
        |            |
        |            |
  Q -- STAR -- Q -- STAR -- Q --
        |            |
        |            |


"""

class PlanarLattice:

    """
    Planar lattice class
    """

    def __init__(self,size):

        self.size=size

        self.N_Q=2*size*size+2*size+1
        self.N_full_P=size*(size-1)
        self.N_edge_P=size

        self.positions_anyons_S=None
        self.positions_anyons_P=None

        ## Define basic qubit and stabilizer positions

        self.positions_Q=[(x,y) for x in range(2*size+1) for y in range((x%2),2*size+1,2) ]
        self.positions_P=[(x,y) for x in range(0,2*size+1,2) for y in range(1,2*size+1,2)]
        self.positions_S=[(x,y) for x in range(1,2*size+1,2) for y in range(0,2*size+1,2)]

        ## Define position lists by stabilizer type

        self.positions_full_S=[(x,y) for x in range(2,2*size-1,2) for y in range(1,2*size+1,2)]
        self.positions_edge_S_T=[(0,y) for y in range(1,2*size,2)]
        self.positions_edge_S_B=[(2*size,y) for y in range(1,2*size,2)]
        self.positions_full_P=[(x,y) for x in range(1,2*size,2) for y in range(2,2*size-1,2)]
        self.positions_edge_P_L=[(x,0) for x in range(1,2*size,2)]
        self.positions_edge_P_R=[(x,2*size) for x in range(1,2*size,2)]

        ## Defining list of positions for two interspersed round of plaquette measurement

        self.positions_full_P_1=[(x,y) for x in range(1,2*size,2) for y in range((x+1)%4+2,2*size,4)]
        self.positions_full_P_2=[(x,y) for x in range(1,2*size,2) for y in range((x-1)%4+2,2*size,4)]
        self.positions_edge_P_L_1=[(x,0) for x in range(1,2*size,4)]
        self.positions_edge_P_L_2=[(x,0) for x in range(3,2*size,4)]
        self.positions_edge_P_R_1= [(x,2*size) for x in range(2*(size%2)+1,2*size,4)]
        self.positions_edge_P_R_2=[(x,2*size) for x in range(2*((size+1)%2)+1,2*size,4)]

        self.positions_full_S_1=[(x,y)  for x in range(2,2*size-1,2) for y in range((x+2)%4+1,2*size,4)]
        self.positions_full_S_2=[(x,y)  for x in range(2,2*size-1,2) for y in range((x%4)+1,2*size,4)]
        self.positions_edge_S_T_1=[(0,y) for y in range(1,2*size,4)]
        self.positions_edge_S_T_2=[(0,y) for y in range(3,2*size,4)]
        self.positions_edge_S_B_1=[(2*size,y) for y in range(2*(size%2)+1,2*size,4)]
        self.positions_edge_S_B_2=[(2*size,y) for y in range(2*((size+1)%2)+1,2*size,4)]

        ## Initialise array

        self.array=[[1]*(2*self.size+1) for _ in range(2*self.size+1)]

        for p0,p1 in self.positions_Q: self.array[p0][p1]=[1,1]
        for p0,p1 in self.positions_P: self.array[p0][p1]=1
        for p0,p1 in self.positions_S: self.array[p0][p1]=1




    # Methods for Displaying the state of the array
    #==============================================


    def showArray(self,arrayType,channel=0):

        c=0 if channel=="X" else 1

        if arrayType=="errors":
            print_array=[[x[c] if isinstance(x,list) else 0 for x in row]for row in self.array]

        if arrayType=="stabilizers":
            print_array=[[x if isinstance(x,int) else 0 for x in row] for row in self.array]


        plt.imshow(print_array)
        plt.show()


    def showArrayText(self,arrayType="errors",channel=0):

        if channel in ["X","x",0]: c=0
        elif channel in ["Z","z",1]: c=1
        else:
            raise ValueError('%s is not a valid channel for showArrayText(), channel must be "X" or "Z '%(channel,))

        if arrayType in ["error","errors","Errors","Error"]:

            print_array = [[str(x[c]) if isinstance(x,list) else '.' for x in row] for row in self.array]
            print_array = [[channel if x=='-1' else x for x in row] for row in print_array]


        elif arrayType in ["stabilizers","stabs","stabilisers","stabilizer","stabiliser"]:

            print_array = [[str(x) if isinstance(x,int) else '.' for x in row] for row in self.array]

        elif arrayType in ["all","both"]:

            print_array = [[str(x[c]) if isinstance(x,list) else ('.' if x==1 else '#') for x in row] for row in self.array]

        else:
            raise ValueError('%s is not a valid arrayType for showArrayText()'%(arrayType,))

        print_array = [['&' if (x=='#' and i%2==0)  else x for x in print_array[i]] for i in range(len(print_array))]


    def applyRandomErrors(self,pX,pZ):
        """ applies random X and Z errors to every qubit in the array

        Parameters:
        ----------
        pX -- probability of X error
        pZ -- probability of Z error

        """
        for q0,q1 in self.positions_Q:
            rand1=random.random()
            rand2=random.random()

            if rand1<pX:
                self.array[q0][q1][0]*=-1
            if rand2<pZ:
                self.array[q0][q1][1]*=-1


    def applyRandomErrorsXYZ(self,pX,pY,pZ):

        for p0,p1 in self.positions_Q:
            if random.random()<pX:
                self.array[p0][p1][0]*=-1
            if random.random()<pY:
                self.array[p0][p1][0]*=-1
                self.array[p0][p1][1]*=-1
            if random.random()<pZ:
                self.array[p0][p1][1]*=-1



    def findAnyons(self):

        anyon_positions_x=[]
        anyon_positions_z=[]

        for i in range(self.N_full_P):
            if self.full_S[i]==-1:
                anyon_positions_x+=[self.positions_full_S[i]]

        for i in range(self.N_edge_P):
            if self.edge_S_T[i]==-1:
                anyon_positions_x+=[self.positions_edge_S_T[i]]
            if self.edge_S_B[i]==-1:
                anyon_positions_x+=[self.positions_edge_S_B[i]]

        for i in range(self.N_full_P):
            if self.full_P[i]==-1:
                anyon_positions_z+=[self.positions_full_P[i]]

        for i in range(self.N_edge_P):
            if self.edge_P_L[i]==-1:
                anyon_positions_z+=[self.positions_edge_P_L[i]]
            if self.edge_P_R[i]==-1:
                anyon_positions_z+=[self.positions_edge_P_R[i]]


        self.positions_anyons_S=anyon_positions_x
        self.positions_anyons_P=anyon_positions_z



    # MEASUREMENT
    #============

    def updateStabilizer(self,p0,p1,stabQubits,channel,pLie):
        stab=1
        for s0,s1 in stabQubits:
            stab*=self.array[s0][s1][channel]

        rand = random.random()
        if rand<pLie: stab*=-1
        self.array[p0][p1]=stab


    def measurePlaquettes(self,pLie=0):


        for p0,p1 in self.positions_full_P:
            stabQubits=((p0,p1-1),(p0,p1+1),(p0-1,p1),(p0+1,p1))
            self.updateStabilizer(p0,p1,stabQubits,0,pLie)

        for p0,p1 in self.positions_edge_P_L:

             stabQubits = ((p0,p1+1),(p0-1,p1),(p0+1,p1))
             self.updateStabilizer(p0,p1,stabQubits,0,pLie)

        for p0,p1 in self.positions_edge_P_R:
             stabQubits = ((p0,p1-1),(p0-1,p1),(p0+1,p1))
             self.updateStabilizer(p0,p1,stabQubits,0,pLie)


    def measureStars(self,pLie=0):

        for p0,p1 in self.positions_full_S:
            stabQubits=((p0,p1-1),(p0,p1+1),(p0-1,p1),(p0+1,p1))
            self.updateStabilizer(p0,p1,stabQubits,1,pLie)

        for p0,p1 in self.positions_edge_S_T:
            stabQubits = ((p0,p1-1),(p0,p1+1),(p0+1,p1))
            self.updateStabilizer(p0,p1,stabQubits,1,pLie)

        for p0,p1 in self.positions_edge_S_B:
            stabQubits = ((p0,p1-1),(p0,p1+1),(p0-1,p1))
            self.updateStabilizer(p0,p1,stabQubits,1,pLie)


    # MEASUREMENT ACCORDING TO AN ERROR VECTOR
    #=========================================

    def stabilizer(self,channel,pos,error,sType,round12,stabilizersNotComplete=0):

        if round12 not in [1,2]: raise ValueError('round12 must be either 1 or 2')

        if channel in ["plaquette","P"]: c=0
        elif channel in ["star","Star","S"]:c=1
        else: raise ValueError('%s is not a valid channel for the stabilizer'%(channel))


        p0,p1=pos

        [up,down,left,right] = [(p0-1,p1),(p0+1,p1),(p0,p1-1),(p0,p1+1)]

        if sType =="F": stabQubits = (left,right,up,down)
        elif sType =="L": stabQubits = (right,up,down)
        elif sType =="R": stabQubits = (left,up,down)
        elif sType =="T": stabQubits = (left,right,down)
        elif sType =="B": stabQubits = (left,right,up)
        else:
            raise ValueError(' pType must be a valid plaquette type: F,R,L,T or B ')

        order = range(len(stabQubits))


        # Add option that a stabilizer doesn't get measured. The stabilizer value stays the same as the previous round and no errors are applied

        if random.random()>stabilizersNotComplete:

            lie,err=error

            if round12==2:
                stab=lie
                for q0,q1 in stabQubits:
                    stab*=self.array[q0][q1][c]
                    self.array[p0][p1]=stab

            random.shuffle(order)
            errQubits=[stabQubits[order[0]],stabQubits[order[1]]]


            if err!=[[1,1],[1,1]]:
                for i in range(len(errQubits)):
                    q=errQubits[i]
                    self.array[q[0]][q[1]][0]*=err[i][0]
                    self.array[q[0]][q[1]][1]*=err[i][1]

            if round12==1:
                stab=lie
                for q0,q1 in stabQubits:
                    stab*=self.array[q0][q1][c]
                self.array[p0][p1]=stab


    def apply_matching(self,error_type,matching):
        """ applies appropriate flips to the array to bring it back to the codespace using the given matching """

        channel=0 if error_type=="X" else 1

        flips=[]

        for pair in matching:

            [p0,p1]=pair[0]
            [q0,q1]=pair[1]

            if channel==1 and (p1==-1 or p1==self.size*2+1)and(q1==-1 or q1==self.size*2+1):
                flips+=[]
            elif channel==0 and (p0==-1 or p0==self.size*2+1)and(q0==-1 or q0==self.size*2+1):
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
            self.array[flip[0]][flip[1]][channel]*=-1


    def apply_flip_array(self,channel,flip_array):

        # Convert channel to index
        c=0 if channel=="X" else 1

        # For each qubit, apply a flip to the qubit error channel if needed
        for (x0,x1) in self.positions_Q:
            self.array[x0][x1][c]*=flip_array[x0][x1]


    def measure_logical(self):

        logical_x=1
        logical_z=1
        positions_x=[[0,x] for x in range(0,2*self.size+1,2)]
        positions_z=[[y,0] for y in range(0,2*self.size+1,2)]

        for pos in positions_x:
            logical_x*=self.array[pos[0]][pos[1]][0]
        for pos in positions_z:
            logical_z*=self.array[pos[0]][pos[1]][1]

        return [logical_x,logical_z]
















######################################################################






class PlanarLattice3D:

    def __init__(self,size):

        self.size=size
        self.N=size*(size+1) # number of stabilizers

        self.syndrome_P=[1]*self.N
        self.syndrome_S=[1]*self.N

        self.parity_array_P=[]
        self.parity_array_S=[]

        self.time=len(self.parity_array_P)

        self.definite_array_P=[[1]*self.N]

        self.positions_S=[(x,y) for x in range(0,2*size+1,2) for y in range(1,2*size+1,2)];
        self.positions_P=[(x,y) for x in range(1,2*size+1,2) for y in range(0,2*size+1,2)]

    def getTime(self):
        self.time=len(self.parity_array_P)


    def showTopLayer(self):

        print_array=[[0 for x in range(2*self.size+1)] for y in range(2*self.size+1)]


        for i in range(len(self.positions_S)):
            Spos=self.positions_S[i]
            Ppos=self.positions_P[i]
            print_array[Spos[0]][Spos[1]]=self.parity_array_S[-1][i]
            print_array[Ppos[0]][Ppos[1]]=self.parity_array_P[-1][i]


        plt.imshow(print_array)
        plt.show()


    def addMeasurement(self,lat):

        plaquette_layer=[1 for x in range(self.N)]
        star_layer=[1 for x in range(self.N)]

        for i in range(self.N):

            Ppos=self.positions_P[i]
            Spos=self.positions_S[i]

            plaquette_layer[i]=lat.array[Ppos[0]][Ppos[1]]
            star_layer[i]=lat.array[Spos[0]][Spos[1]]


        new_syndrome_P=copy.copy(plaquette_layer)
        new_syndrome_S=copy.copy(star_layer)


        for i in range(self.N):
            plaquette_layer[i]*=self.syndrome_P[i]
            star_layer[i]*=self.syndrome_S[i]

        self.parity_array_S+=[star_layer]
        self.parity_array_P+=[plaquette_layer]

        self.syndrome_P=new_syndrome_P
        self.syndrome_S=new_syndrome_S

    def findAnyons(self):

        self.getTime()

        anyon_positions_x=()
        anyon_positions_z=()

        for t in range(self.time):

            anyon_positions_x_t=()
            anyon_positions_z_t=()

            for i in range(self.N):


                if self.parity_array_P[t][i]==-1:
                    anyon_positions_x_t+=((t,)+(self.positions_P[i]),)

                if self.parity_array_S[t][i]==-1:
                    anyon_positions_z_t+=((t,)+(self.positions_S[i]),)

            anyon_positions_x+=(anyon_positions_x_t,)
            anyon_positions_z+=(anyon_positions_z_t,)

        self.anyon_positions_x=anyon_positions_x
        self.anyon_positions_z=anyon_positions_z

        self.anyon_positions_P = anyon_positions_x
        self.anyon_positions_S = anyon_positions_z
