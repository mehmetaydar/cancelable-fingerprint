from utils import *
import reedsolomon as rs
import random
from rsolomon2 import *
from cartesian2 import *
import os
import math
import json
import cv2
import os
import numpy
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import distance

class Fvs:
    def __init__(self):
        self.num_cart = 0.0
        self.tot_cart_sim = 0.0
        self.avg_sim = 0.0
        self.num_cart_hash_match = 0.0

        #self.testWithoutRecoveringOriginalMinutiaPoints()
        #self.testWithRecoveringOriginalMinutiaPoints()
        #self.testWithRecovering_HashOfOriginalMinutiaPoints(rm_divide=4)
        pass

    def testWithRecovering_HashOfOriginalMinutiaPoints(self, des_original, des_candidate, rm_divide1 = 2, rm_divide2 = 2, prune_extra_rsolomon = False):
        cart = Cartesian(des_original)
        #Reed solomon per cartesian blocks
        original_rd_conf = {} #keep
        #cart.Cprime_list_minutias#keep
        rd_kmax1 = int(255/ (3-2/rm_divide1))
        rd_kmax2 = int(255 / (3 - 2 / rm_divide2))
        for c, minutia_list in cart.C_list_minutias.iteritems():
            paddedminutialist = [','.join( map(str, ptuple) ).rjust(8,'0') for ptuple in minutia_list]

            """a temp fix, need to have a better reed solomon library"""
            if len(paddedminutialist) > rd_kmax1:
                paddedminutialist = sorted(paddedminutialist)[0:rd_kmax1]

            """reed solomon per hashes of minutia points per cartesian blocks"""
            hashesof_paddedminutialist = [getHash_base64_grabFirst8(padded_minutia) for padded_minutia in paddedminutialist]
            hash_order = {}
            order = 0
            for hash0 in hashesof_paddedminutialist:
                hash_order[hash0] = order
                order = order + 1
            splitted_hashesof_paddedminutialist = getSplittedHashes(hashesof_paddedminutialist, 8) #need to do this because of max symsize of the reed solomon is 8

            total_hash_of_splitted_hashesof_paddedminutialist = getHashOfListOfTuples(splitted_hashesof_paddedminutialist, pad=False)
            rm = max(int(len(splitted_hashesof_paddedminutialist) / rm_divide1), 2)
            rs = Rsolomon(_rm=rm)
            if len(splitted_hashesof_paddedminutialist) >rs.kmax:
                print "In encodeSolomon number of list: " + str(len(splitted_hashesof_paddedminutialist)) + " must be less than or equal to kmax: " + str(
                    rs.kmax)
                pass
            k, checksums, n = rs.encodeSolomon(splitted_hashesof_paddedminutialist, prune_extra_rsolomon)
            original_rd_conf[c] = {'checksums': checksums, 'rm': rm, 'k': k, 'n': n, 'hash_order': hash_order,
                                   #'number_of_original_hashes': len(hash_order),
                                   'number_of_original_hashes': len(hashesof_paddedminutialist),
                                   'hash': total_hash_of_splitted_hashesof_paddedminutialist}

        """reed solomon per overall hashes of cartesian blocks"""
        hash_order = {}
        hashes0 = [v['hash'] for c, v in original_rd_conf.iteritems()]
        order = 0
        for hash0 in hashes0:
            hash_order[hash0] = order
            order = order + 1
        hashes = getSplittedHashes(hashes0, 8)
        total_hash = getHashOfListOfTuples(hashes, pad=False)
        rm = max(int(len(hashes) / rm_divide2), 2)
        rs = Rsolomon(_rm=rm)
        k, checksums, n = rs.encodeSolomon(hashes, prune_extra_rsolomon)
        hash_rd_conf = {'checksums': checksums, 'rm': rm, 'k': k, 'n': n,
                        #'number_of_original_hashes': len(hash_order),
                        'number_of_original_hashes': len(hashes0),
                        'hash': total_hash, 'hash_order': hash_order}

        Cgenerated = {}
        #assert des_candidate.bounds == cart.bounds, "Original and candiate bounds must be equal"
        cart2 = Cartesian(des_candidate, W=cart.W, H=cart.H, bounds=cart.bounds, M=cart.M, keep_transforms=True) #notice that keep_transforms=False is true

        self.num_cart += len(cart.Cprime_list_minutias)
        """for c, minutia_list2 in cart2.Cprime_list_minutias.iteritems():
            print "cartesian block: ", c
            if c in cart.Cprime_list_minutias:  # candidate
                minutia_list1 = cart.Cprime_list_minutias[c]  # original

                print "minutialist1: ", minutia_list1
                print "minutialist2: ", minutia_list2
                matches = set([])
                for p in sorted(minutia_list1):
                    for q in sorted(minutia_list2):
                        if q in matches:
                            continue
                        if self.checkIfMinutiaMatch(p, q):
                            # print(cart2.Creverse)
                            c_original = cart2.Creverse[hash(tuple(q))]
                            if c_original not in Cgenerated:
                                Cgenerated[c_original] = set([])
                            original_minutia = cart2.transformMinutiaToNewBlock(p, c, c_original, back_to_original=True)
                            Cgenerated[c_original].add(original_minutia)
                            matches.add(q)
                            break"""


        for c, minutia_list2 in cart2.Cprime_list_minutias.iteritems():
            print "cartesian block: ", c
            if c in cart.Cprime_list_minutias: #candidate
                minutia_list1 = cart.Cprime_list_minutias[c] #original


                """s1 = np.array( [(x[0], x[1]) for x in minutia_list1])
                s2 = np.array([(x[0], x[1]) for x in minutia_list2])
                print np.min(cdist(s1,s2))
                print distance.cdist(s1, s2,  'euclidean').min(axis=1)"""

                #print "minutialist1: ", sorted(minutia_list1)
                #print "minutialist2: ", sorted(minutia_list2)
                matchesq = set([])
                matchesp = set([])
                for p in minutia_list1:
                    if tuple(p) in matchesp:
                        continue
                    dist = 99999999999999999
                    bestq = None
                    for q in minutia_list2:
                        if tuple(q) in matchesq:
                            continue
                        if self.checkIfMinutiaMatch(p, q) and self.minutiaDistance(p, q) < dist:
                            bestq = q
                    if bestq is not None:
                        c_original = cart2.Creverse[hash(tuple(bestq))]
                        if c_original not in Cgenerated:
                            Cgenerated[c_original] = set([])
                        original_minutia = cart2.transformMinutiaToNewBlock(p, c, c_original, back_to_original=True )
                        Cgenerated[c_original].add( tuple(original_minutia) )
                        matchesq.add(tuple(bestq))
                        matchesp.add(tuple(p))
                        #print "Match: ",p, bestq

        genhashes = []
        for c, minutia_list in Cgenerated.iteritems():
            #print "generated:"
            #print c, sorted(list(minutia_list))

            paddedminutialist = [','.join(map(str, ptuple)).rjust(8, '0') for ptuple in minutia_list]
            hashesof_paddedminutialist = [getHash_base64_grabFirst8(padded_minutia) for padded_minutia in paddedminutialist]
            splitted_hashesof_paddedminutialist = getSplittedHashes(hashesof_paddedminutialist, 8)
            genhash = getHashOfListOfTuples(splitted_hashesof_paddedminutialist, pad=False)

            print "generated-hash: ", genhash

            if c not in cart.C_list_minutias:
                print "Cartesian number: ", c, " does not exists in original fingerprint template!!"
                continue

            #print "original:"
            #print c, sorted(cart.C_list_minutias[c])
            #orghash = getHashOfListOfTuples(cart.C_list_minutias[c]) #later you need to make sure that you pre-store the hash
            orghash = original_rd_conf[c]['hash']
            print "original-hash: ", orghash

            if genhash == orghash:
                print "c generated and original minutias are equal \n"
                self.tot_cart_sim += 1.0
                self.num_cart_hash_match += 1.0
            else:
                print "generated and original minutias are NOT equal"
                print "trying reed solomon error correction codes"

                rm = original_rd_conf[c]['rm']
                checksums = original_rd_conf[c]['checksums']
                k = original_rd_conf[c]['k']

                #k_number_of_original_hashes = k / 3  # this is because we previously splitted 24 length of hashes to 8 length
                k_number_of_original_hashes = original_rd_conf[c]['number_of_original_hashes']

                n = original_rd_conf[c]['n']
                hash_order = original_rd_conf[c]['hash_order']
                new_hashesof_paddedminutialist = []
                for index in range(0, k_number_of_original_hashes):  # if available minutia is less than data part then add dummy data
                    #new_hashesof_paddedminutialist.append('x' * 24)  # 24 is the length of the base64 hash
                    new_hashesof_paddedminutialist.append('x' * 8)  # 8 is the length of the base64 hash / grabbed the first 8
                num_match_hash = 0.0
                for hash0 in hashesof_paddedminutialist:
                    if hash0 in hash_order:
                        order = hash_order[hash0]
                        if order < len(new_hashesof_paddedminutialist):
                            new_hashesof_paddedminutialist[order] = hash0
                            num_match_hash +=1.0
                        else:
                            print "Possible error. order < len(new_hashesof_paddedminutialist), ", order, len(new_hashesof_paddedminutialist)
                hashesof_paddedminutialist = new_hashesof_paddedminutialist
                splitted_hashesof_paddedminutialist = getSplittedHashes(hashesof_paddedminutialist, 8)

                self.avg_sim = self.tot_cart_sim/self.num_cart
                print "Avg-sim: ", self.avg_sim
                rs = Rsolomon(_rm=rm)
                print "to-be-decoded: ", splitted_hashesof_paddedminutialist
                try:
                    decoded, corrections = rs.decodeSolomon(splitted_hashesof_paddedminutialist, k, checksums, n)
                    print "corrections: ", corrections
                    print "decoded ", decoded
                    genhash = getHashOfListOfTuples(list(decoded), pad = False)

                    assert genhash==orghash, "Generated and original hashes are not equal after decode!"

                    print "c hashes are equal \n"
                    self.tot_cart_sim += 1.0
                    self.num_cart_hash_match += 1.0
                except Exception as e:
                    print "c hashes are not equal ! reed-solomon-cannot-decode:", str(e), " \n"
                    self.tot_cart_sim += num_match_hash/ len(new_hashesof_paddedminutialist)
            genhashes.append(genhash)

        k = hash_rd_conf['k']

        #k_number_of_original_hashes = k/3 # this is because we previously splitted 24 length of hashes to 8 length
        k_number_of_original_hashes = hash_rd_conf['number_of_original_hashes']

        n = hash_rd_conf['n']
        hash_order = hash_rd_conf['hash_order']
        new_genhashes = []
        for index in range(0, k_number_of_original_hashes):  # if available minutia is less than data part then add dummy data
            new_genhashes.append('x' * 24) # 24 is the length of the base64 hash
            # new_genhashes.append('x' * 8)  # 8 is the length of the base64 hash / grabbed the first 8
        for hash0 in genhashes:
            if hash0 in hash_order:
                new_genhashes[hash_order[hash0]] = hash0
        genhashes = new_genhashes
        genhashes = getSplittedHashes(genhashes, 8)

        try:
            rs = Rsolomon(_rm=hash_rd_conf['rm'])
            print "to-be-decoded: ", genhashes
            decoded, corrections = rs.decodeSolomon(genhashes, k, hash_rd_conf['checksums'], n)
            print "corrections: ", corrections
            print "decoded ", decoded

            total_hash_decoded = getHashOfListOfTuples(list(decoded), pad=False)
            assert hash_rd_conf['hash'] == total_hash_decoded, "Decoded hash is not equal to the original hash"
            print "Final Original hash and decoded hash are equal: ", hash_rd_conf['hash'], " equals ", total_hash_decoded, "\n"
            return True
        except Exception as e:
            print "Final Original hash and decoded hash are NOT equal:  ! reed-solomon-cannot-decode:", str(e), "\n"
        return False

    def minutiaDistance(self, p, q, match = 'eucludian'):
        if p[2] != q[2]: return 9999999
        if match == 'scalar':
            dist = abs(p[0]- q[0])+abs(p[1] - q[1])
        elif match == 'eucludian':
            dist = math.sqrt(sum([(a - b) ** 2 for a, b in zip(p, q)]))
        return dist
    def checkIfMinutiaMatch(self, p, q):
        rang1 = math.atan2(p[1], p[0])
        rang2 = math.atan2(q[1], q[0])
        #check minutia types, and distances
        if p[2] == q[2] and \
                abs(p[0]- q[0]) <= xmatchthreshold and \
                abs(p[1] - q[1]) <= ymatchthreshold \
                and abs(rang1-rang2) <= anglethreshold:
        #if        abs(rang1 - rang2) <= anglethreshold:
            return True
        return False

_prune_extra_rsolomon = False
xmatchthreshold = 12.0
ymatchthreshold = 12.0
anglethreshold = 1
#_rm_divide1 = 4
#_rm_divide2 = 5
_rm_divide1 = 4
_rm_divide2 = 4
ffiles = []
for i in range(1, 4):
    for j in range(1, 4):
        if i < 10:
            finputminutia = '/Users/Admin/fvs/samples/DB3_B/10' + str(i) + '_' + str(j) + '.tif'
        else:
            finputminutia = '/Users/Admin/fvs/samples/DB3_B/1' + str(i) + '_' + str(j) + '.tif'
        ffiles.append({
            'i': i,
            'j': j,
            'image': finputminutia
        })
errors = 0
totaltry = 0
tp = tn = fp = fn = 0.0
actual = False
tps =[]
fps =[]
fns =[]
tns =[]
tot_avg_sim = 0.0
avg_avg_sim = 0.0
for index1 in range(0, len(ffiles)):
    f1 = ffiles[index1]
    for index2 in range(0, len(ffiles)):
        if index1 == index2: continue
        f2 = ffiles[index2]
        finputminutia = f1['image']
        finputminutia2 = f2['image']
        finputminutia_dump = "/Users/Admin/fvs/samples/dumps2/" + os.path.basename(finputminutia) + ".yaml"
        finputminutia2_dump = "/Users/Admin/fvs/samples/dumps2/" + os.path.basename(finputminutia2) + ".yaml"
        with open(finputminutia_dump) as f:
            des_original = json.load(f)
        with open(finputminutia2_dump) as f:
            des_candidate = json.load(f)
        actual = False
        if f1['i'] == f2['i']:
            actual = True
        print "Original: ", finputminutia
        print "Candidate: ", finputminutia2
        if not os.path.isfile(finputminutia) or not os.path.isfile(finputminutia2):
            print "Problem on files. Either them does not exists."
            continue
        match = False
        """try:
            fvs = Fvs()
            match = fvs.testWithRecovering_HashOfOriginalMinutiaPoints(des_original, des_candidate, rm_divide=_rm_divide, prune_extra_rsolomon = _prune_extra_rsolomon)
            #exit()
        except Exception as ex:
            print "Exception: ", ex.message
            errors = errors+1"""
        fvs = Fvs()
        match = fvs.testWithRecovering_HashOfOriginalMinutiaPoints(des_original, des_candidate
                                                                   , rm_divide1=_rm_divide1, rm_divide2=_rm_divide2,
                                                                   prune_extra_rsolomon=_prune_extra_rsolomon)
        print "num-cart: ", fvs.num_cart, "num-cart-match: ", fvs.num_cart_hash_match
        print "processed: ", (f1['i'], f1['j']), (f2['i'], f2['j'])
        print "--next-pair--\n"
        tot_avg_sim+= fvs.avg_sim
        totaltry = totaltry + 1
        if match and actual:
            tp = tp+1.0
            tps.append( (finputminutia, finputminutia2) )
        elif match and not actual:
            fp = fp+1.0
            fps.append( (finputminutia, finputminutia2) )
        elif not match and actual:
            fn = fn+1.0
            fns.append( (finputminutia, finputminutia2) )
        elif not match and not actual:
            tn = tn+1.0
            tns.append( (finputminutia, finputminutia2) )
if tp+fp > 0:
    precision = tp / (tp+fp)
else: precision = None

if tp+fn > 0:
    recall = tp / (tp+fn)
else: recall = None

if precision is not None and recall is not None and (precision+recall) > 0:
    fmeasure = 2*((precision*recall) / (precision+recall))
else: fmeasure = None

printList(tps, "True Positives")
printList(fps, "False Positives")
printList(tns, "True Negatives")
printList(fns, "False Negatives")

avg_avg_sim = tot_avg_sim/totaltry

print "\nNum-try: ", totaltry, "\tNum errors: ", errors,
print "xmatchthreshold: ", xmatchthreshold, "\tymatchthreshold: ", xmatchthreshold, "rmdivide: ", _rm_divide1, _rm_divide2
print "Precision: ", precision, "\tRecall: ", recall, "\tFmeasure: ", fmeasure, "\n"
print "avg_avg_sim: ", avg_avg_sim

""" 
1. investigate if shamir secret can work here
since each cartesian block will have a different hash
so if I know some hashes, can i recover all hashes.
or
2. re-do reed-solomon for cartesian blocks hash values
"""
"""
TODO:
do this for the candidate fingerpring image too
from candidate fingerprint image do matching and reversing
and then do reed-solomon decoding for each of the cartesian blocks
"""

"""
minutiaelist=[]
with open(finputminutia) as f:
    lines = f.readlines()
for line in lines:
    minutiaelist.append(line.strip('\n').rjust(8,'0'))
    l = line.strip('\n').split(',')
rs = Rsolomon()
for l in chunks(minutiaelist, rs.kmax):
    rs.encodeSolomon(l)
"""


