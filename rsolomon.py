from utils import *
import reedsolomon as rdsolomonlib
import random

__all__ = ['Rsolomon']

"""
  rs.Codec(n, k)
  n is the total number of symbols per codeword, including data
    and parity.  
  k is the number of data symbols per codeword.
  symsize is max 8, so n should be less than 2 ** symsize - 1 < 256,

  if m=8
  then according to n = 3*k - 2*m
  k could be at most 90. then total number of parities = 255-90 
  so lets divide original to chunks of 90
  that means a list length of 90 minutias, if we have 8 number of correct matches then we are able to 
  regenerate the whole set 
  
  new:
  This means that the encoder takes k data symbols of s bits each and adds parity symbols to make an n symbol codeword. 
  There are n-k parity symbols of s bits each. A Reed-Solomon decoder can correct up to t symbols that contain errors in a codeword, where 2t = n-k.
  t = (n-k)/2 max numbers of errors / n will be adjusted accordingly
  rm = k-t=> rm = k-(n-k)/2
  n = 3k - 2m => n = 3k - 2rm 
"""
_rsymsize=8 #default, cannot be greater than 8
class Rsolomon:
    def __init__(self, _rm = 127):
        self.rsymsize = _rsymsize
        self.rm = _rm  # m correct matches will fix all,
        self.max_rn = 2 ** self.rsymsize - 1  # max 255
        self.min_rn = 2
        # n = 3*k - 2*m , taken from thesis paper
        self.kmax = int((self.max_rn + 2 * self.rm) / 3)  # kmax is 169, so number of data part is 169 (k), parity part is 86
        #print "max_rn: ", self.max_rn, " symsize: ", self.rsymsize, " kmax: ", self.kmax  # n is max 255

    """
        Encode reed solomon for a given list of minutia
    """
    def encodeSolomon(self, l, prune_extra = False):
        l2 = l
        k = len(l2)  # must be <= kmax
        if prune_extra and k > self.kmax:
            print "In encodeSolomon number of list: "+str(k)+ " must be less than or equal to kmax: " + str(self.kmax)
            print "That is why we are pruning extra points beginning from last. need a splitting solution!!"
            l2 = l[0:self.kmax]
            k = len(l2)
        if k > self.kmax:
            raise "In encodeSolomon number of list: "+str(k)+ " must be less than or equal to kmax: " + str(self.kmax)
        n = (3 * k) - (2 * self.rm)
        if n <= k:
            n = k*(self.rm)
        #print "k: ", k, " n: ", n, " rm: ", self.rm
        rcodec = rdsolomonlib.Codec(n, k, self.rsymsize)
        encoded = rcodec.encodechunks(l2)
        #print "encoded: ", encoded
        #checksums = list(encoded[k:])  # this is the parity list
        checksums = encoded[k:]  # this is the parity list
        return k, checksums, n # returns the data part length, and parity

    """
            Decode reed solomon for a given list of minutia / minutia list may not be correct
    """
    def decodeSolomon(self, data_part, k, checksum, n):
        encoded = tuple(data_part)+checksum
        assert k <= self.kmax, "In encodeSolomon length of list must be less than or equal to kmax: " + self.kmax
        #print "k: ", k, " n: ", n, " rm: ", self.rm
        rcodec = rdsolomonlib.Codec(n, k, self.rsymsize)
        decoded, corrections = rcodec.decodechunks(encoded)
        return (decoded, corrections)

    """
    Test encode reed solomon for a given list of minutia
    """
    def encodeSolomonTest(self, l):
        k = len(l)  # must be <= kmax
        n = (3 * k) - (2 * self.rm)
        assert k <= self.kmax, "In encodeSolomon length of list must be less than or equal to kmax: "+self.kmax
        print "k: ", k, " n: ", n, " rm: ", self.rm
        rcodec = rdsolomonlib.Codec(n, k, self.rsymsize)
        rcodec2 = rdsolomonlib.Codec(n, k, self.rsymsize)
        #for decoding either rcodec or rcodec2 can be used.

        # exit()
        encoded = rcodec.encodechunks(l)
        print "encoded: ", encoded
        checksums = encoded[k:]  # this is the parity list
        print "checksum-parities: ", checksums
        decoded, corrections = rcodec.decodechunks(encoded)
        print "decoded: ", decoded
        assert cmpT(decoded, l) == True, " Cannot decode!"
        brokenl = list(encoded)
        brokenl = brokenl[0:k]
        random.shuffle(brokenl)
        brokenl = brokenl + list(checksums)
        for i in range(0, len(l) - self.rm):
        #for i in range(0, len(l)):
            print "i: ", i
            brokenl[i] = 'x' * 8
            broken = tuple(brokenl)
            print "broken: ", broken
            #decoded, corrections = rcodec.decodechunks(broken)
            ereasures = range(0,i+1)
            decoded, corrections = rcodec2.decodechunks(broken, ereasures)
            print "to-be-decoded: ", broken
            print "corrections: ", corrections
            assert cmpT(decoded, l)== True, "Cannot decode!"
            print "decoded[0:k]: ", decoded[0:k]
            print ""


#n=k*2 + 8

"""
If P is the number of parity symbols, R is the number of errors (corrupted bytes), 
and S is the number of erasures, RS coding guarantees recoverability if the following equation is true::

2R + S <= P
 Codec c = Codec(7, 5)
 can recover if 2 bytes are missing
 P=2

"""


