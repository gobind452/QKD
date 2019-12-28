import numpy as np 
from cvxopt import glpk,matrix,solvers,spmatrix,sparse
from base import Sequence,getErrors
import gc
from scipy.sparse import random
import sys
import time

# Setting parameters
keyLength = 1000
errorFraction = 0.6
syndromeFraction = 0.50
checkMatrixDensity = 0.01 # Important 
timeLimit = 50000
similarityEve = 0.3 # Expected errors in Eve will (errorFraction*similarityEve + (1-errorFraction)*(1-similarityEve))*keyLength

def getBobAndEveErrors(): # Returns a tuple (BobErros,EveErrors) (If timeout then it is equal to the -1)
	errorLength = int(errorFraction*keyLength)
	syndromeLength = int(syndromeFraction*keyLength)
	checkMatrix = random(syndromeLength,keyLength,density=checkMatrixDensity,data_rvs=lambda length: np.random.randint(low=1,high=2,size=(length)))
	aliceKey = Sequence() 
	aliceKey.randomiseSequence(keyLength) # Random Alice Key
	syndrome = np.reshape(checkMatrix.dot(aliceKey),(-1,1)) # Calculate Syndrome
	bobKey = Sequence(aliceKey)
	bobKey.induceErrors(errorLength) # Induce errors in Bob
	print("Errors in Bob",errorLength)
	bobKey = np.reshape(np.array(bobKey),(-1,1))
	c = np.subtract(np.equal(bobKey,0).astype(int),np.equal(bobKey,1).astype(int)).astype(float)
	c = matrix(c)
	A = spmatrix(V=checkMatrix.data,I=checkMatrix.row,J=checkMatrix.col,size=(syndromeLength,keyLength))
	b = matrix(syndrome)
	gc.collect()
	start = time.time()
	status,x = glpk.ilp(c=c,G=matrix(np.zeros([1,keyLength])),h=matrix(np.array([[float(0)]])),A=A,b=b,B=set([i for i in range(keyLength)]),options={'tm_lim': timeLimit,'show_progress':False,'msg_lev':'GLP_MSG_OFF'})
	end = time.time()
	sol = list(x)
	sol = [round(e) for e in sol]
	if status == "undefined":
		bobErrors = -1
	else:
		bobErrors = getErrors(aliceKey,sol[:keyLength])
	similarIndexes = list(np.random.choice(a=range(len(bobKey)),size=(int(len(bobKey)*similarityEve)),replace=False))
	eveKey = Sequence(list([bobKey[i]^1 for i in range(keyLength)]))
	print("Errors in Eve",keyLength-len(similarIndexes))
	for i in similarIndexes:
		eveKey[i] = bobKey[i][0]
	gc.collect()
	eveKey = np.reshape(np.array(eveKey),(-1,1))
	c = np.subtract(np.equal(eveKey,0).astype(int),np.equal(eveKey,1).astype(int)).astype(float)
	A = spmatrix(V=checkMatrix.data,I=checkMatrix.row,J=checkMatrix.col,size=(syndromeLength,keyLength))
	#status,x = glpk.ilp(c=c,G=matrix(np.zeros([1,keyLength])),h=matrix(np.array([[float(0)]])),A=A,b=b,B=set([i for i in range(keyLength)]),options={'tm_lim': timeLimit,'show_progress':False,'msg_lev':'GLP_MSG_OFF'})
	gc.collect()
	sol = list(x)
	sol = [round(e) for e in sol]
	if status == "undefined":
		eveErrors = -1
	else:
		eveErrors = getErrors(aliceKey,sol[:keyLength])
	return bobErrors,eveErrors,end-start

print(getBobAndEveErrors())
