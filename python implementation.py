## here are the necessary functions to calculate the persistent cup length diagram associated to a choice of cocycle representatives.

#import libraries

import matplotlib.pyplot as plt
import oineus as oin
import numpy as np
from copy import deepcopy

# functions: cup_Product (performs cup product on cocycle representatives, Algorithm 1 in "Persistent Cup-Length")
#            get_Pivots (pulls out the pivot rows of the reduced coboundary matrix R, necessary for determining whether a given cochain is a coboundary)
#            matrix_Multiplication_for_Sparse_Arrays (Used in evaluating U.y in Algorithm 2 of "Persistent Cup-Length")
#            get_Row_Operation_Matrix (Algorithm 3 from "Persistent Cup-Length" with correction to the last line of the algorithm)
#            get_Cocycle_Reps (Takes the barcodes from the oineus implementation and assigns them their cocycle representative by constructing a Cocycle object)
#            check_Coboundary (Gives the condition for the while loop in line 14 of Algorithm 2 in "Persistent Cup-Length")
#            Persistent_Cup_Diagram (Algorithm 2 in "Persistent Cup-Length" which produces the Persistent Cup-Length Diagram)



def cup_Product(sigma_1,sigma_2, cells):
    #cells is a list of simplices, sigma_i is an element of the cohomology ring given by the ids of the simplices.
    #initialize output
    sigma = []
    #Check for 0 input
    if (sigma_1 == [] or sigma_2 == []):
        return sigma
    #Efficiencies to implement later: Dimension bound to skip a loop
    for i in range(0,len(sigma_1)):
            a = cells[sigma_1[i]].vertices
            for j in range(0,len(sigma_2)):
                b = cells[sigma_2[j]].vertices
                if(a[-1] == b[0]):
                    c = a+b[1:len(b)]
                    #Check if the result is in cells. 
                    #Efficiency to implement later: Reduce range to just those simplices of the right dimension. maybe using dim first dim last in filtration
                    for k in range(0,len(cells)):
                        if(c==cells[k].vertices):
                            sigma.append(k)
    sigma.sort()
    return sigma

def get_Pivots(Rmatrix):
    #Get the pivot columns from the reduced coboundary matrix.
    #Initalise output
    Pivots_R = []
    for i in range(0,len(Rmatrix)):
        if(Rmatrix[i] != []):
            Pivots_R.append(Rmatrix[i][-1])
    return Pivots_R

def matrix_Multiplication_for_Sparse_Arrays(A,B):
    #Given the sparse matrix notation of oineus which lists columns with a nonzero entry assuming Z_2 coefficients, compute their matrix multiplication A.B
    #I am assuming size A = size B
    AB = []
    m = len(B)
    for j in range(m):
        C = set([])
        for i in B[j]:
            C = C^set(A[i])
        AB.append(list(C))
    for i in range(m):
        AB[i].sort()
    return AB

def get_Row_Operation_Matrix( R):
    #Given reduced matrix R get the Row Operation Matrix U needed to check for coboundaries
    #Initialise output
    U = []
    Rm = deepcopy(R)
    m = len(Rm)
    for j in range(m):
        #make U the identity matrix
        U.append([j])
    for j in range(m):
        #Go over each column to reduce
        Rm[j].sort()
        verts = Rm[j]
        #Check for empty column
        if (verts != []):
            E = []
            for n in range(m):
                #make U the identity matrix
                E.append([n])
            #select row of Rmatrix
            for i in verts[0:len(verts)-1]:
                #Update each column of Rmatrix
                for k in range(m):
                    if (verts[-1] in Rm[k]):
                        if (i in Rm[k]):
                            #remove i from Rmatrix[k]
                            Rm[k].remove(i)
                        else:
                            Rm[k].append(i)
                #Update U to account for the above row operation            
                E[verts[-1]].append(i)
            U = matrix_Multiplication_for_Sparse_Arrays(E,U)
    return U


class Cocycle:
    def __init__(self, simplices, birth, death):
        self.simplices = simplices
        self.birth = birth
        self.death = death
        self.life = [birth, death]
        self.birth_index = 0
        self.death_index = 0
        
def get_Cocycle_Reps(dcmp,fil,V, dim_bound):
    #ensure dcmp has been reduced before calling this function
    Cocycle_Reps = []
    for i in range(1,dim_bound+1):
        dgms = dcmp.diagram(fil).in_dimension(i)
        indx_dgms = dcmp.diagram(fil).index_diagram_in_dimension(i)
        for j in range(len(dgms)):
            sigma_idx = fil.id_by_sorted_id(indx_dgms[j][0])
            new_Cocycle = Cocycle(V[len(V)-1-sigma_idx],dgms[j][0],dgms[j][1])
            Cocycle_Reps.append(new_Cocycle)
    return Cocycle_Reps

def check_Coboundary(Pivots, U, sigma, length):
    vector = [[]]
    for k in sigma:
        if k >= length:
            vector[0].append(k)
    image = matrix_Multiplication_for_Sparse_Arrays(U,vector)
    result = []
    for k in image[0]:
        if k>= length:
            result.append(k)
    sigma_is_coboundary = set(result).issubset(Pivots)
    return sigma_is_coboundary        

def Persistent_Cup_Diagram(dim_bound, Cosimplices, Rmatrix, Representative_Cocycles ):
    #representative cocycles must have birth time and death time 
    #they should be ordered first by increasing death time and then in increasing birth time
    B = [[]]
    #Note that this puts the index one below the index listed in Algorithm 2
    for barcode in Representative_Cocycles:
        if len(Cosimplices[barcode.simplices[0]].vertices) != 1:
            B[0].append(barcode)
    b_time = []
    d_time = []
    for barcode in B[0]:
        if (barcode.birth in b_time) == False:
            barcode.birth_index = len(b_time)
            b_time.append(barcode.birth)
        else:
            for k in range(len(b_time)):
                if barcode.birth == b_time[k]:
                    barcode.birth_index = k
        if (barcode.death in d_time) == False:
            barcode.death_index = len(d_time)
            d_time.append(barcode.death)
        else: 
            for k in range(len(d_time)):
                if barcode.death == d_time[k]:
                    barcode.death_index = k
    #The index stuff is to avoid having float == float statements later on though there are some in determining the indices
    m = len(Cosimplices)
    Pivots = get_Pivots(Rmatrix)
    l = 1
    U = get_Row_Operation_Matrix(Rmatrix)
    A_0 = np.zeros((len(b_time),len(d_time)))
    A_1 = np.zeros((len(b_time),len(d_time)))
    A = []
    A.append(A_0)
    A.append(A_1)
    for barcode in B[0]:
        i = barcode.birth_index
        j = barcode.death_index
        A[1][i,j] = 1
    while ((A[l] != A[l-1]).any()) and (l <= dim_bound-1):
        A.append(deepcopy(A[l]))
        B.append([])
        for barcode_1 in B[0]:
            for barcode_2 in B[l-1]:
                sigma = cup_Product(barcode_1.simplices,barcode_2.simplices, Cosimplices)
                if sigma == []:
                    continue
                death_index = min(barcode_1.death_index, barcode_2.death_index)
                index_set = []
                for k in range(len(b_time)):
                    if b_time[k] <= d_time[death_index]:
                        index_set.append(k)
                birth_index = max(index_set)
                s = 0
                for tau in Cosimplices:
                    if tau.value <= b_time[birth_index]:
                        s+= 1
                while (check_Coboundary(Pivots, U, sigma, m+1-s) == True) and (birth_index >= 0):
                    birth_index = birth_index-1
                    s = 0
                    for tau in Cosimplices:
                        if tau.value <= b_time[birth_index]:
                            s+= 1
                if b_time[birth_index] < d_time[death_index]:
                    new_Cocycle = Cocycle(sigma,b_time[birth_index], d_time[death_index])
                    new_Cocycle.birth_index = birth_index
                    new_Cocycle.death_index = death_index
                    B[l].append(deepcopy(new_Cocycle))
                    A[l+1][birth_index,death_index] = l+1
        l+=1
    return [A[l],b_time,d_time]

