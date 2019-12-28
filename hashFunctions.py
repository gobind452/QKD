import numpy as np

primesFile = "primes.txt"
primes = {}

def loadPrimes():
	global primes
	with open(primesFile,'r') as f:
		for line in f:
			x = line.rstrip().split(',')
			primes[int(x[0])] = int(x[1])
	return

def polyHash(message,key,outputLength,upperBound): # Optimise modulo large primes
	l = int((upperBound+1)/outputLength)
	appendLength = l*outputLength- len(message)
	if appendLength >= 1:
		message = message + "1" + (appendLength-1)*"0"
	prime = pow(2,outputLength)+ primes[outputLength]
	hashResult = 0
	for i in range(l):
		num = message[i*outputLength:(i+1)*outputLength]
		num = (int(num,2)*int(key[i]))% prime
		hashResult = hashResult + num
	return np.binary_repr(hashResult)

def toeplitzHash(message,key):
	alpha = len(message)
	beta = len(key)-len(message)+1
	key = [int(key[i]) for i in range(len(key))]
	matrix = []
	for i in range(beta):
		matrix.append(key[beta-i-1:beta-i+alpha-1])
	matrix = np.array(matrix)
	message = np.reshape(np.array([int(message[i]) for i in range(len(message))]),(-1,1))
	return np.remainder(np.matmul(matrix,message),2)

def oneTimePad(message,key):
	if len(message) != len(key):
		raise Exception("Incorrect size for one time pad")
	return np.binary_repr(int(message,2)^int(key,2))