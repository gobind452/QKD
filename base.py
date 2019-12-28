import numpy as np

class Sequence(list): # List subclass for data sequences
	def randomiseSequence(self,length): # Generate a random sequence
		self[:] = Sequence(np.random.randint(low=0,high=2,size=(length)))
	
	def induceErrors(self,numErrors,errorIndexes=-1): # Induce errors at specified indexes or random if no index specified
		if errorIndexes == -1:
			errorIndexes = np.random.choice(a=range(len(self)),size=(numErrors),replace=False)
		for errorIndex in errorIndexes:
			self[errorIndex] = self[errorIndex]^1
		return
 
class Matrix(np.ndarray):
	def __new__(cls,inputArray):
		obj = np.asarray(inputArray).view(cls)
		return obj

	def __array_finalize__(self, obj):
		return 

	def add(self,x):
		x = np.array(x)
		return Matrix(self^x)

	def matmul(self,x):
		x = np.array(x)
		return Matrix(np.matmul(self,x)%2)

class Binary(): # Utility class for binary-decimal operations
	def getBinary(number):
		return bin(number)[2:]

	def getOnesPosition(number):
		binary = Binary.getBinary(number)
		binary = binary[::-1]
		answer = []
		for i in range(len(binary)):
			if binary[i] == '1':
				answer.append(i)
		return answer

	def generateDecimal(binary):
		if type(binary) == Sequence:
			binary = [str(e) for e in binary]
			binary = "".join(binary)
		return int(binary,2)

def getErrors(a,b):
	return sum([int(a[i])^int(b[i]) for i in range(len(a))])