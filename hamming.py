from base import Sequence
from base import Binary
import numpy as np

class Hamming():
	def lengthHammingBits(length): # Returns the number of parity bits for data sequence for length
		hammingBits = 0
		while 1:
			if np.power(2,hammingBits)-hammingBits-1>=length:
				return hammingBits
			hammingBits = hammingBits+1

	def encodeHamming(seq): # Encodes a data sequence and returns the new stream
		hammingBits = Hamming.lengthHammingBits(len(seq))
		encoded = Sequence([0 for i in range(hammingBits+len(seq))])
		pos = 0
		for i in range(len(encoded)):
			ones = Binary.getOnesPosition(i+1)
			if len(ones) == 1:
				encoded[i] = 0
			else:
				encoded[i] = seq[pos]
				pos = pos + 1
				for one in ones:
					encoded[np.power(2,one)-1] = encoded[np.power(2,one)-1]^encoded[i]
		return encoded

	def correctData(seq): # Corrects a data sequence for 1 bit error
		hammingBits = len(Binary.getBinary(len(seq)))
		actualParityBits = Sequence([0 for i in range(hammingBits)])
		receivedParityBits = Sequence()
		for i in range(len(seq)):
			ones = Binary.getOnesPosition(i+1)
			if len(ones) == 1:
				receivedParityBits.append(seq[i])
			else:
				for one in ones:
					actualParityBits[one] = actualParityBits[one]^seq[i]
		errorPosition = Sequence()
		for i in range(hammingBits):
			errorPosition.append(actualParityBits[i]^receivedParityBits[i])
		errorPosition.reverse()
		errorPosition = Binary.generateDecimal(errorPosition)
		if errorPosition > 0:
			seq[errorPosition-1] = seq[errorPosition-1]^1
		return seq

	def extractData(seq): # Extracts data from a encoded sequence
		data = Sequence()
		for i in range(len(seq)):
			ones = Binary.getOnesPosition(i+1)
			if len(ones) > 1:
				data.append(seq[i])
		return data

if __name__ == "__main__":
	x = Sequence()
	x.randomiseSequence(8) # Random data stream of 5 bits
	print(x) # Original data
	encoded = Hamming.encodeHamming(x) # Encoded stream (Alice)
	print(encoded)
	encoded.induceErrors(1) # Random 1 bit error
	print(encoded)
	corrected = Hamming.correctData(encoded) # Bob corrects data
	print(Hamming.extractData(corrected)) # Extracted original data (Should be equal to x)