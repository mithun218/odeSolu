
# cython: language_level=3, boundscheck=False

import numpy as np

cpdef long int factorial(long int n):
    return np.prod(range(1,n+1))


cdef class adomianMat:

    cdef:
        long int j,k,l,m,N,DEtoArrayLLen,DEtoArrayLtermLen,DEtoArrayNLLen,iniCondLen,lastSeriestermPos,orderDegz,lastCalTerm,lastTermInr
        bint isBefrAssn
        
        #Using memoryviews#
        double[:] emptyNDM
        double[:] emptyNAM1
        double[:] emptyNAM2
        double[:] zeroMNAddclt
        double[:] zeroMNMulclt
        double[:] zeroSoluM
        double[:] oneMN
        double[:] soluMroll

        double[:] derivM
        double[:] integM
        
        double[:] soluM
        double[:,:] DEtoArrayL, DEtoArrayNL
        
    def __cinit__(self, _N,_DEtoArrayLLen,_DEtoArrayLtermLen,_DEtoArrayNLLen,_iniCondLen,_lastSeriestermPos,_orderDegz,
                _soluM,_DEtoArrayL,_DEtoArrayNL):
        self.N = _N
        self.DEtoArrayLLen = _DEtoArrayLLen
        self.DEtoArrayLtermLen = _DEtoArrayLtermLen
        self.DEtoArrayNLLen = _DEtoArrayNLLen
        self.iniCondLen = _iniCondLen
        self.lastSeriestermPos = _lastSeriestermPos
        self.orderDegz = _orderDegz
                
        self.soluM = _soluM
        self.DEtoArrayL = _DEtoArrayL
        self.DEtoArrayNL = _DEtoArrayNL
        
        self.emptyNDM = np.empty(_N,dtype=np.double) 
        self.emptyNAM1 = np.empty(_N,dtype=np.double)
        self.emptyNAM2 = np.empty(_N,dtype=np.double)
        self.zeroMNAddclt = np.empty(_N,dtype=np.double) # addition collector
        self.zeroMNMulclt = np.empty(_N,dtype=np.double) # production collector
        self.zeroSoluM = np.empty(_N,dtype=np.double) # to copy solution matrix
        self.oneMN = np.ones(_N,dtype=np.double)

        self.derivM = np.array(range(_N),dtype=np.double)
        self.integM = 1./(np.array(range(1,(_N+1)),dtype=np.double))
        self.integM[-1] = 0 # last elements are assigned to 0 for avoiding large number multipication        


    cdef double[:] doDerivM(self,double[:] matIn, double _Nt): # _Nt times differentiation
        for self.l in range(int(_Nt)):
            matIn = np.roll(np.multiply(self.derivM, matIn), -1)
        return matIn
    
    cdef double[:] doIntegM(self,double[:] matIn,long int _Nt): # _Nt times integration
        for self.l in range(1, _Nt+1): 
            matIn = np.roll(np.multiply(self.integM, matIn), 1)
        return matIn
    
    cdef double[:] adomianMuN(self,double[:] _soluM, double _Nt): # this will create adomian matrix of _soluM^(_Nt-1) for the matrix elements
        self.emptyNAM1[:] = _soluM     # after the element position lastSeriestermPos
        for self.m in range(int(_Nt)-2):
            for self.l in range(self.lastSeriestermPos-self.orderDegz+1,-1,-1):
                self.emptyNAM1[self.l] = np.multiply(_soluM[:(self.l+1)],np.array([np.flip(self.emptyNAM1[:(self.l+1)])])).sum()
        return self.emptyNAM1

    cdef double[:] adomianU1U2(self,double[:] U1,double[:] U2): # this will create adomian matrix between the matrix U1 and U2
        np.asarray(self.emptyNAM2).fill(0)
        for self.l in range(self.lastSeriestermPos-self.orderDegz+1,-1,-1):
            self.emptyNAM2[self.l] = np.multiply(U1[:(self.l+1)],np.array([np.flip(U2[:(self.l+1)])])).sum()
        return self.emptyNAM2 
    
    cpdef seriesTerms(self):
        '''initializetion (here we use self.orderDegz in place of seriesTermCalALoop in first term calculation in series)'''
        #N = 40
        #self.soluM = np.array([ 1.,  1 ,0.5 ,0, 0,0,0])
        # Calculating series terms from linear part (self.DEtoArrayL) in r.h.s. (We need not to calculate Adomian Matrix)
        np.asarray(self.zeroMNAddclt).fill(0) 
        np.asarray(self.zeroMNMulclt).fill(0)

        self.isBefrAssn = False # it indicates whether self.zeroMNMulclt has been assigned before



        self.lastTermInr = 1 #It increases self.lastSeriestermPos by 1
        self.lastSeriestermPos = self.orderDegz-1 #last series term position in self.soluM. At first it is equal to order of diff. equ.
        while self.lastTermInr<self.N:        
            for self.j in range(self.DEtoArrayLLen): # cal. are performed in summation term
                for self.k in range(self.DEtoArrayLtermLen): # cal. are performed in a product term
                    if(self.k==0 and self.DEtoArrayL[self.j][self.k]!=0): # constant in product term 
                        if(np.asarray(self.DEtoArrayL[self.j][2:]).sum()!=0): # from 2 as dpndt var must be present
                            self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = self.DEtoArrayL[self.j][self.k]
                            self.isBefrAssn = True
                        elif(self.lastSeriestermPos == self.orderDegz-1
                            and np.asarray(self.DEtoArrayL[self.j][1:]).sum()==0):# one time derivation during running of the loop
                            self.zeroMNMulclt[0] = self.DEtoArrayL[self.j][self.k]
                            self.isBefrAssn = True
                        else:
                            continue

                    elif(self.k==1 and self.DEtoArrayL[self.j][self.k]!=0): # only indpnd^n in product term
                        if( np.asarray(self.DEtoArrayL[self.j][2:]).sum()==0 and 
                           self.lastSeriestermPos == self.DEtoArrayL[self.j][self.k]+self.orderDegz-1):# one time calculation for polynomial terms  already 
                                                                                # present in input diff. equ. during running of the loop
                            if(self.DEtoArrayL[self.j][0]==0): # checself.king nonzero number at self.DEtoArrayL[self.j][0]
                                self.zeroMNMulclt[int(self.DEtoArrayL[self.j][self.k])] = float(1) # coefficient is 1
                            else:
                                self.zeroMNMulclt[int(self.DEtoArrayL[self.j][self.k])] = self.DEtoArrayL[self.j][0] # coefficient is of self.DEtoArrayL[self.j][0]
                            self.isBefrAssn = True

                    elif(self.k==2 and self.DEtoArrayL[self.j][self.k]!=0): # dpnd in product term
                        if(not self.isBefrAssn):
                            if(self.DEtoArrayL[self.j][1]==0):
                                self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = self.soluM[self.lastSeriestermPos-self.orderDegz+1]
                            else:
                                self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = np.roll(self.soluM, int(self.DEtoArrayL[self.j][1]))[self.lastSeriestermPos-self.orderDegz+1]
                            self.isBefrAssn = True
                        else:
                            if(self.DEtoArrayL[self.j][1]==0): # dpndnd and indpnd^n in product term
                                self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = (self.soluM[self.lastSeriestermPos-self.orderDegz+1]
                                                                                *self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1])
                            else:
                                self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = (np.roll(self.soluM, int(self.DEtoArrayL[self.j][1]))[self.lastSeriestermPos-self.orderDegz+1]
                                                                                *self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1])

                    elif(self.k>2 and self.DEtoArrayL[self.j][self.k]!=0): # higher order derivative terms
                        if(not self.isBefrAssn):
                            self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = self.doDerivM(self.soluM, self.k-2)[self.lastSeriestermPos-self.orderDegz+1]
                            self.isBefrAssn = True
                        else:
                            if(self.DEtoArrayL[self.j][1]==0): # dpndnd and indpnd^n in product term
                                self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = (self.doDerivM(self.soluM, self.k-2)[self.lastSeriestermPos-self.orderDegz+1]
                                                                                *self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1])
                            else:
                                self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = (np.roll(self.doDerivM(self.soluM, self.k-2), int(self.DEtoArrayL[self.j][1]))[self.lastSeriestermPos-self.orderDegz+1]
                                                                                *self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1])                                  

                    else:
                        continue
      
                if(self.DEtoArrayLLen):
                    self.zeroMNAddclt = np.add(self.zeroMNMulclt, self.zeroMNAddclt) # adding production terms in addition collector
                    np.asarray(self.zeroMNMulclt).fill(0) # clearing production collector
                    self.isBefrAssn = False


            # Calculating series terms from nonlinear part (self.DEtoArrayNL) in r.h.s. (We need to calculate Adomian Matrix)
            self.lastCalTerm = self.lastSeriestermPos-self.orderDegz+1
            for self.j in range(self.DEtoArrayNLLen): # cal. are performed in summation term
                for self.k in range(self.DEtoArrayLtermLen-1, -1, -1): # cal. are performed in a product term
                    if(self.k==0 and self.DEtoArrayNL[self.j][self.k]!=0): # constant in product term 
                        if(not self.isBefrAssn):
                            self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = self.DEtoArrayNL[self.j][self.k]
                            self.isBefrAssn = True
                        else:
                            self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1] = self.DEtoArrayNL[self.j][self.k]*self.zeroMNMulclt[self.lastSeriestermPos-self.orderDegz+1]


                    elif(self.k==2 and self.DEtoArrayNL[self.j][self.k]!=0): # dpnd^n or dpnd in product term. Here adomian cal ends.
                        if(not self.isBefrAssn):
                            self.zeroMNMulclt[:] = self.adomianMuN(self.soluM, self.DEtoArrayNL[self.j][self.k])
                            if(self.DEtoArrayNL[self.j][1]==0):
                                self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.soluM[:self.lastCalTerm+1], np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()
                            else:
                                if(self.lastCalTerm<self.DEtoArrayNL[self.j][1]): #self.lastCalTerm is the power of last term. If it is smaller than
                                    self.zeroMNMulclt[self.lastCalTerm] = 0  # self.DEtoArrayNL[self.j][1], it is assined to zero.
                                else:
                                    self.soluMroll = np.roll(self.soluM, int(self.DEtoArrayNL[self.j][1]))
                                    self.soluMroll[:int(self.DEtoArrayNL[self.j][1])] = 0
                                    self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.soluMroll[:self.lastCalTerm+1], 
                                                            np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()
                            self.isBefrAssn = True
                            if(self.lastCalTerm!=0):
                                self.zeroMNMulclt[:self.lastCalTerm] = 0
                        else:
                            if(self.DEtoArrayNL[self.j][self.k]==1):                                     
                                if(self.DEtoArrayNL[self.j][1]==0):
                                    self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.soluM[:self.lastCalTerm+1], np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()
                                else:
                                    if(self.lastCalTerm<self.DEtoArrayNL[self.j][1]): #self.lastCalTerm is the power of last term. If it is smaller than
                                        self.zeroMNMulclt[self.lastCalTerm] = 0  # self.DEtoArrayNL[self.j][1], it is assined to zero.
                                    else:
                                        self.soluMroll = np.roll(self.soluM, int(self.DEtoArrayNL[self.j][1]))
                                        self.soluMroll[:int(self.DEtoArrayNL[self.j][1])] = 0
                                        self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.soluMroll[:self.lastCalTerm+1], 
                                                                np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()                                                  

                                if(self.lastCalTerm!=0):
                                    self.zeroMNMulclt[:self.lastCalTerm] = 0
                            else:
                                self.zeroMNMulclt[:] = self.adomianU1U2(self.adomianMuN(self.soluM, self.DEtoArrayNL[self.j][self.k]), self.zeroMNMulclt)
                                if(self.DEtoArrayNL[self.j][1]==0):
                                    self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.soluM[:self.lastCalTerm+1], np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()
                                else:
                                    if(self.lastCalTerm<self.DEtoArrayNL[self.j][1]): #self.lastCalTerm is the power of last term. If it is smaller than
                                        self.zeroMNMulclt[self.lastCalTerm] = 0  # self.DEtoArrayNL[self.j][1], it is assined to zero.
                                    else:
                                        self.soluMroll = np.roll(self.soluM, int(self.DEtoArrayNL[self.j][1]))
                                        self.soluMroll[:int(self.DEtoArrayNL[self.j][1])] = 0
                                        self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.soluMroll[:self.lastCalTerm+1], 
                                                                np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()

                                if(self.lastCalTerm!=0):
                                    self.zeroMNMulclt[:self.lastCalTerm] = 0


                    elif(self.k>2 and self.DEtoArrayNL[self.j][self.k]!=0): # higher order derivative terms 
                        np.asarray(self.emptyNDM).fill(0)
                        self.emptyNDM[:self.lastCalTerm+1] = (self.doDerivM(self.soluM, self.k-2)[:self.lastCalTerm+1])
                        if(not self.isBefrAssn and np.asarray(self.DEtoArrayNL[self.j][2:self.k]).sum()!=0):
                            if(self.DEtoArrayNL[self.j][self.k]==1):
                                self.zeroMNMulclt[:] = self.emptyNDM
                            else:
                                self.zeroMNMulclt[:] = self.adomianMuN(self.emptyNDM, self.DEtoArrayNL[self.j][self.k]+1)
                            self.isBefrAssn = True

                        elif(not self.isBefrAssn and np.asarray(self.DEtoArrayNL[self.j][2:self.k]).sum()==0):
                            self.zeroMNMulclt[:] = self.adomianMuN(self.emptyNDM, self.DEtoArrayNL[self.j][self.k])
                            if(self.DEtoArrayNL[self.j][1]==0):
                                self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.emptyNDM[:self.lastCalTerm+1], np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()
                            else:
                                if(self.lastCalTerm<self.DEtoArrayNL[self.j][1]): #self.lastCalTerm is the power of last term. If it is smaller than
                                    self.zeroMNMulclt[self.lastCalTerm] = 0  # self.DEtoArrayNL[self.j][1], it is assined to zero.
                                else:
                                    self.soluMroll = np.roll(self.emptyNDM, int(self.DEtoArrayNL[self.j][1]))
                                    self.soluMroll[:int(self.DEtoArrayNL[self.j][1])] = 0
                                    self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.soluMroll[:self.lastCalTerm+1], 
                                                            np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()                                      

                            if(self.lastCalTerm!=0):
                                self.zeroMNMulclt[:self.lastCalTerm] = 0
                            self.isBefrAssn = True

                        elif(self.isBefrAssn and np.asarray(self.DEtoArrayNL[self.j][2:self.k]).sum()!=0):
                            if(self.DEtoArrayNL[self.j][self.k]==1):
                                self.zeroMNMulclt[:] = self.adomianU1U2(self.emptyNDM, self.zeroMNMulclt)
                            else:
                                self.zeroMNMulclt[:] = self.adomianU1U2(self.adomianMuN(self.emptyNDM, self.DEtoArrayNL[self.j][self.k]+1), self.zeroMNMulclt)

                        elif(self.isBefrAssn and np.asarray(self.DEtoArrayNL[self.j][2:self.k]).sum()==0):
                            if(self.DEtoArrayNL[self.j][self.k]!=1):
                                self.zeroMNMulclt[:] = self.adomianU1U2(self.adomianMuN(self.emptyNDM, self.DEtoArrayNL[self.j][self.k]), self.zeroMNMulclt)      
                            if(self.DEtoArrayNL[self.j][1]==0):
                                self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.emptyNDM[:self.lastCalTerm+1], np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()
                            else:
                                if(self.lastCalTerm<self.DEtoArrayNL[self.j][1]): #self.lastCalTerm is the power of last term. If it is smaller than
                                    self.zeroMNMulclt[self.lastCalTerm] = 0  # self.DEtoArrayNL[self.j][1], it is assined to zero.
                                else: 
                                    self.soluMroll = np.roll(self.emptyNDM, int(self.DEtoArrayNL[self.j][1]))
                                    self.soluMroll[:int(self.DEtoArrayNL[self.j][1])] = 0
                                    self.zeroMNMulclt[self.lastCalTerm] = np.multiply(self.soluMroll[:self.lastCalTerm+1], 
                                                            np.flip(self.zeroMNMulclt[:self.lastCalTerm+1])).sum()                                          

                            if(self.lastCalTerm!=0):
                                self.zeroMNMulclt[:self.lastCalTerm] = 0

                    else:
                        continue

                if(self.DEtoArrayNLLen):
                    self.zeroMNMulclt[self.lastCalTerm+1:] = 0
                    self.zeroMNAddclt = np.add(self.zeroMNMulclt, self.zeroMNAddclt) # adding production terms in addition collector
                    np.asarray(self.zeroMNMulclt).fill(0) # clearing production collector    
                    self.isBefrAssn = False


            self.soluM = np.add(self.soluM, self.doIntegM(self.zeroMNAddclt,self.iniCondLen)) # final result is collected in self.soluM
            np.asarray(self.zeroMNAddclt).fill(0) # clearing summation collector
            self.lastSeriestermPos = self.orderDegz-1 + self.lastTermInr# after one term cal. last pos. in self.soluM is increased by 1
            self.lastTermInr += 1 # to inrease last pos. in self.soluM 
            
    def solution(self):
        return np.asarray(self.soluM)