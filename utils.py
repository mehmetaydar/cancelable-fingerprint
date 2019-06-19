import random
import hashlib
import base64

__all__ = ['matrixMultiply', 'cmpT','cmpLT', 'chunks', 'genTransformationMatrix',
           'getHashOfListOfTuples', 'splitStringEqualLength','getSplittedHashes',
           'testSets', 'getHash_base64_grabFirst8', 'printList']

def matrixMultiply(X, Y):
    result = [[sum(a * b for a, b in zip(X_row, Y_col)) for Y_col in zip(*Y)] for X_row in X]
    return result

def cmpT(t1, t2):#compre two tuples
  return sorted(t1) == sorted(t2)
def cmpLT(l1, l2):#compre list of tuples
  return sorted(l1) == sorted(l2)
def chunks(l, lchunk):
    """Yield successive lchunk-sized chunks from l."""
    for i in range(0, len(l), lchunk):
        yield l[i:i + lchunk]

"""generate lenCxlenC transformation matrix values randomly selected either 0 or 1  """
def genTransformationMatrix(lenC):
    """
    if
    C=[
    [1,2,3,4]
    ]
    M=[
        [0,0,0,0],
        [0,1,0,0],
        [1,0,0,1],
        [0,0,1,0]
    ]
    so generate lenCxlenC matrix, in which we have a 1 in only of the columns. the rest is 0
    """
    M = []
    for i in range(1, lenC+1): #rows
        M.append( [0] * lenC)
    for j in range(0, lenC):
        index1_i = random.randint(0, lenC-1)
        M[index1_i][j] = 1
    return M

def getHash_base64_original(o):
    h = hashlib.md5(o).digest()
    h = base64.b64encode(h)
    return h

"""
since in rsolomon we are restricted in 8 char we are doing this stripping.
splitting does not always work in rsolomon because of the restriction of number of data items. k must be less than kmax issue.
this might cause collusion. neeeds to be tested
"""
def getHash_base64_grabFirst8(o):
    h = hashlib.md5(o).digest()
    h = base64.b64encode(h)
    h = h[0:8]
    return h

def getHashOfListOfTuples(l, pad = True):
    if pad:
        paddedminutialisttest = [','.join(map(str, ptuple)).rjust(8, '0') for ptuple in l]
    else:
        paddedminutialisttest = l
    #minutiae_set_hash = hashlib.md5(','.join( sorted(paddedminutialisttest) )).hexdigest()

    #minutiae_set_hash = hashlib.md5(','.join(sorted(paddedminutialisttest))).digest()
    #minutiae_set_hash = base64.b64encode(minutiae_set_hash)
    #return minutiae_set_hash
    return getHash_base64_original(','.join(sorted(paddedminutialisttest)))
def getHashOfListOfTuples_grabFirst8(l, pad = True):
    if pad:
        paddedminutialisttest = [','.join(map(str, ptuple)).rjust(8, '0') for ptuple in l]
    else:
        paddedminutialisttest = l
    #minutiae_set_hash = hashlib.md5(','.join( sorted(paddedminutialisttest) )).hexdigest()

    #minutiae_set_hash = hashlib.md5(','.join(sorted(paddedminutialisttest))).digest()
    #minutiae_set_hash = base64.b64encode(minutiae_set_hash)
    #return minutiae_set_hash
    return getHash_base64_grabFirst8(','.join(sorted(paddedminutialisttest)))

def splitStringEqualLength(x, size):
    #lenx = 24, size = 8, num = 3
    num = len(x)/size
    chunks, chunk_size = len(x), len(x) / num
    return [x[i:i + chunk_size] for i in range(0, chunks, chunk_size)]

def getSplittedHashes(hashes0, hash_size):
    hashes = []
    for hash in hashes0:
        for hash_part in splitStringEqualLength(hash, hash_size):
            hashes.append(hash_part)
    return hashes

def testSets(s1, s2):
    assert len(s1) == len(s2), "List length are not equal! " + len(s1)+ ' != '+ len(s2)
    print "Len list: ", len(s1)
    s1 = set(s1)
    s2 = set(s2)
    intersect = s1.intersection(s2)
    print "Len intersection: ", len(intersect)

#def alignListsForReedSolomonDecode(l1, l2):


def printList(ls, desc):
    print desc + ":"
    for l in ls:
        print l