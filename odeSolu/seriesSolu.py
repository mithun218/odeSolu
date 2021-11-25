import numpy as np
#import numpy.polynomial.polynomial as nppoly
from sympy import *
from odeSolu.adomianMat import adomianMat

def seriesSolu(_inputDE,_idpdDpdVar,_iniCond,_N):
    inputDE = expand(_inputDE)
    idpdDpdVar = _idpdDpdVar
    iniCond = _iniCond # this is the initial conditions
    N = _N
    inputDEz_ = symbols('inputDEz_')
    w1_ = Wild('w1_')
    w2_ = Wild('w2_')
    inputDEz = inputDE.subs(diff(idpdDpdVar[1], (idpdDpdVar[0], ode_order(inputDE,idpdDpdVar[1]))), inputDEz_)

    if(inputDEz.is_polynomial(inputDEz_)):
        orderDeg = (ode_order(inputDE,idpdDpdVar[1]), degree(inputDEz, inputDEz_))   # and independent variables
        inputDEzargs = inputDEz.args
        inputDEzargsDict = {}

        #isDerivInRhs = False # at least one Derivative must be present at r.h.s.
        for j in inputDEzargs:

            if((j).has(inputDEz_)): # the term containing inputDEz_ with its degree is collected in inputDEzargsDict
                inputDEzargsDict[j] = degree(j, inputDEz_)
            #if((j).has(Derivative)): # checking at least one Derivative term other than highest derivative term
                #isDerivInRhs = True       
        if(len(inputDEzargsDict)!=0):
            if(len(inputDEzargsDict)>1):
                inputDEzargsDict = dict(sorted(inputDEzargsDict.items(), key=lambda x: x[1], reverse=False)) # sorting in desending order
            inputDEzargsDictKeys = inputDEzargsDict.keys()
            if(inputDEzargsDict[next(iter(inputDEzargsDictKeys))]==1): # isolating higher deriv. term with degree 1 and r.h.s term 
                inputDErhs = expand((-(inputDE -  next(iter(inputDEzargsDictKeys)).subs(inputDEz_, diff(idpdDpdVar[1], (idpdDpdVar[0], orderDeg[0])))))/ # is assigned to inputDErhs
                                  (next(iter(inputDEzargsDictKeys))).coeff(inputDEz_))

                ### Checking polynomials w.r.t indp. var., diff. terms in diff. eq. ###########
                temInputDErhs = inputDErhs
                diffSymDict = {orderDeg[0]+1:idpdDpdVar[0]}

                for j in range(orderDeg[0], -1, -1):
                    diffSymDict[j] = sympify("_inputDEz"+str(j))
                    temInputDErhs = temInputDErhs.subs( diff(idpdDpdVar[1], (idpdDpdVar[0], j)), diffSymDict[j])
                isPoly = True

                for j in range(orderDeg[0]+2):
                    if(not temInputDErhs.is_polynomial(diffSymDict[j])):
                        isPoly = False
                        break

                if(not isPoly):
                    print("ODE cannot be solved for non-polynomial terms.")
                    return
            else:
                print("ODE cannot be solved for the power of the highest derivative term.")
                return
        else:
            print("ODE cannot be solved.")
            return
    else:
        print("ODE cannot be solved for non-polynomial terms.")
        return

    ###################################
    DEtoList = []
    if(isinstance(inputDErhs,Add)):
        
        coeffAdd = inputDErhs.as_coeff_add()
        
        L1 = len(coeffAdd[1])
        for j in range(L1):
            if((not isinstance(coeffAdd[1][j],Symbol)) 
               and coeffAdd[1][j]!=idpdDpdVar[1] 
               and (not (coeffAdd[1][j]).match(idpdDpdVar[1]**w1_))
               and (not (coeffAdd[1][j]).match(idpdDpdVar[0]**w1_))
               and (not (coeffAdd[1][j]).is_number)): # added1
                
                if((coeffAdd[1][j]).has(Derivative)): 
                    if(not (coeffAdd[1][j]).match(diff(idpdDpdVar[1],idpdDpdVar[0],ode_order(coeffAdd[1][j],idpdDpdVar[1]))**w2_)):
                        
                        tem = list(coeffAdd[1][j].args)
                        if(not (tem[0]).is_number):
                            tem.append(sympify('1'))
                            DEtoList.append(tuple(tem))
                        else:
                            DEtoList.append(tuple(tem))
                    else:
                        DEtoList.append((sympify('1'),coeffAdd[1][j]))

                else:
                    tem = list(coeffAdd[1][j].args)
                    if(not (tem[0]).is_number):
                        tem.append(sympify('1'))
                        DEtoList.append(tuple(tem))
                    else:
                        DEtoList.append(tuple(tem)) 

            else:
                if((coeffAdd[1][j]).is_number):
                    DEtoList.append((coeffAdd[1][j],))
                else:
                    DEtoList.append((sympify('1'),coeffAdd[1][j]))
                    
        DEtoList.append((coeffAdd[0],))
        
    else:
        if((inputDErhs.is_number)):
            DEtoList.append((inputDErhs,))   
        elif((not isinstance(inputDErhs,Symbol)) 
             and inputDErhs!=idpdDpdVar[1]
             and (not (inputDErhs).match(idpdDpdVar[1]**w1_))
             and (not (inputDErhs).match(idpdDpdVar[0]**w1_))
             and (not (inputDErhs).is_number)): # added1
            
            if((inputDErhs).has(Derivative)): 
                if(not (inputDErhs).match(diff(idpdDpdVar[1],idpdDpdVar[0],ode_order(inputDErhs,idpdDpdVar[1]))**w2_)):
                    tem = list(inputDErhs.args)
                    if(not (tem[0]).is_number):
                        tem.append(sympify('1'))
                        DEtoList.append(tuple(tem))
                    else:
                        DEtoList.append(tuple(tem))

                else:
                    DEtoList.append((sympify('1'),inputDErhs))

            else:
                tem = list(inputDErhs.args)
                if(not (tem[0]).is_number):
                    tem.append(sympify('1'))
                    DEtoList.append(tuple(tem))
                else:
                    DEtoList.append(tuple(tem))
                    
        else:
            if((inputDErhs).is_number):
                DEtoList.append((inputDErhs,))
            else:
                DEtoList.append((sympify('1'),inputDErhs))
            
            DEtoList.append((inputDErhs,))

    ############################################

    L1 = len(DEtoList)
    L2 = 3+ode_order(inputDE,idpdDpdVar[1])
    DEtoArrayL = np.zeros((L1,L2))
    DEtoArrayNL = np.zeros((L1,L2))
    for j in range(L1):
        L3 = DEtoList[j]
        NLchk = 0

        for k in range(len(L3)): # nonlinearities are checked in this loop

            if(m1:=(DEtoList[j][k]).match(idpdDpdVar[1]**w1_)): # nonlinearity for dep^w1_
                if(m1[w1_]>1):
                    NLchk = 2
                    break  

            elif((DEtoList[j][k]).has(Derivative)):
                # nonlinearity for diff(dep,ind,w1_)^w2
                if(m1:=(DEtoList[j][k]).match(diff(idpdDpdVar[1],idpdDpdVar[0],ode_order(DEtoList[j][k],idpdDpdVar[1]))**w2_)):
                    if(m1[w2_]>1):
                        NLchk = 2
                        break
            else:
                pass
            
            if(DEtoList[j][k].has(idpdDpdVar[1])): # counting dependent variables for checking nonlinearity
                NLchk += 1


        for k in range(len(L3)):

            if(NLchk>1):

                if((DEtoList[j][k]).is_number): # number
                    DEtoArrayNL[j][0] = lambdify((),DEtoList[j][k],modules='numpy')()

                elif(m1:=(DEtoList[j][k]).match(idpdDpdVar[0]**w1_)): # independent^w1_.
                    DEtoArrayNL[j][1] = lambdify((),m1[w1_],modules='numpy')()

                elif(m1:=(DEtoList[j][k]).match(idpdDpdVar[1]**w1_)): # dependent^w1_
                    DEtoArrayNL[j][2] = lambdify((),m1[w1_],modules='numpy')()

                elif((DEtoList[j][k]).has(Derivative)):
                    #sympy match funtion does not work with wild symbol in Derivative arguments#
                    if(m1:=(DEtoList[j][k]).match(diff(idpdDpdVar[1],idpdDpdVar[0],
                                                       ode_order(DEtoList[j][k],idpdDpdVar[1]))**w2_)): # diff(dep,ind,w1_)^w2_
                        DEtoArrayNL[j][2+ode_order(DEtoList[j][k],idpdDpdVar[1])] = lambdify((),m1[w2_],modules='numpy')()

                else:
                    print("ODE cannot be solved for the parameter(s) in the ODE.")
                    return
            else:
                if((DEtoList[j][k]).is_number): # number
                    DEtoArrayL[j][0] = lambdify((),DEtoList[j][k],modules='numpy')()

                elif(m1:=(DEtoList[j][k]).match(idpdDpdVar[0]**w1_)): # independent^w1_.
                    DEtoArrayL[j][1] = lambdify((),m1[w1_],modules='numpy')()

                elif(m1:=(DEtoList[j][k]==idpdDpdVar[1])): # dependent (sympy match function does not work without wild symbol, so we use ==)
                    DEtoArrayL[j][2] = lambdify((),1,modules='numpy')()

                elif((DEtoList[j][k]).has(Derivative)):
                    if((DEtoList[j][k]==diff(idpdDpdVar[1],idpdDpdVar[0],ode_order(DEtoList[j][k],idpdDpdVar[1])))): # diff(dep,ind,w1_)
                        DEtoArrayL[j][2+ode_order(DEtoList[j][k],idpdDpdVar[1])] = lambdify((),1,modules='numpy')()
                else:
                    print("ODE cannot be solved for the parameter(s) in the ODE.")
                    return

    DEtoArrayL = DEtoArrayL[~np.all(DEtoArrayL == 0, axis=1)] # deleting row containing only zero elements
    DEtoArrayNL = DEtoArrayNL[~np.all(DEtoArrayNL == 0, axis=1)] # deleting row containing only zero elements

    DEtoArrayL = np.double(DEtoArrayL) # converting double number
    DEtoArrayNL = np.double(DEtoArrayNL) # converting double number



    DEtoArrayLLen = len(DEtoArrayL)
    DEtoArrayNLLen = len(DEtoArrayNL)

    DEtoArrayLtermLen = 0 
    if(DEtoArrayLLen):
        DEtoArrayLtermLen = len(DEtoArrayL[0])# array of length for each tuple present in DEtoArrayL and DEtoArrayNL
    elif(DEtoArrayNLLen):
        DEtoArrayLtermLen = len(DEtoArrayNL[0])

    

    max_xPow = np.max(np.append(DEtoArrayL[:, 1], DEtoArrayNL[:, 1])) # maximum power of independent
    #seriesTermCalALoop = max_xPow + orderDeg[0] # number of series terms calculation in a loop

    lastSeriestermPos = orderDeg[0] # last series term position in soluM. At first it is equal to order of diff. equ.. 

          

    #########################################
    if(len(iniCond)<orderDeg[0]):
        print('Please provide sufficient initial conditions.')
        return
    else:
        soluM = np.zeros(N)
        iniCondLen = len(iniCond)
        iniCondM = np.empty(iniCondLen)
        for j in range(iniCondLen):
            orderj = ode_order(iniCond[j][0], idpdDpdVar[1])
            if(orderj == 0):
                soluM[0] = iniCond[j][1]
            else: 
                soluM[orderj] = iniCond[j][1]/factorial(orderj)
                
                
    mat = adomianMat(N,DEtoArrayLLen,DEtoArrayLtermLen,DEtoArrayNLLen,iniCondLen,lastSeriestermPos,orderDeg[0],
                soluM,DEtoArrayL,DEtoArrayNL)

    mat.seriesTerms()
    print('Successfully solved')   
    return mat.solution()
    