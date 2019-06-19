from utils import *
import reedsolomon as rs
import random
from rsolomon import *
from cartesian import *
import os
import math

class Fvs:
    def __init__(self):
        #self.testWithoutRecoveringOriginalMinutiaPoints()
        #self.testWithRecoveringOriginalMinutiaPoints()
        #self.testWithRecovering_HashOfOriginalMinutiaPoints(rm_divide=4)
        pass

    """
          performs reed solomon encode on the transformed cartesian blocks (does it for each block)
          does the matching on transformed cartesian blocks.
          does NOT keep the original c(cartesian rect numbers) transformation on the candidate image
          if matches, using the list of matches for each cartesian block,
            it does the reed solomon encodings/decodings
          checks the hashes on the transformed cartesian blocks
          Use this for cancellable fingerprint templates
    """
    def testWithoutRecoveringOriginalMinutiaPoints(self):
        cart = Cartesian(finputminutia)
        #Reed solomon per transformed cartesian blocks
        transformed_rd_conf = {} #keep
        hash_rd_conf = {}  # keep reed solomon conf for the hash values of the transformed cartesian blocks
        #cart.Cprime_list_minutias#keep
        for c, minutia_list in cart.Cprime_list_minutias.iteritems():
            paddedminutialist = [','.join( map(str, ptuple) ).rjust(8,'0') for ptuple in minutia_list]
            hash = getHashOfListOfTuples(paddedminutialist, pad=False)
            #if there is only half of matches then we can still recover them
            rm = int(len(paddedminutialist)/2)
            rs = Rsolomon(_rm = rm )

            rs.encodeSolomonTest(paddedminutialist)

            k, checksums, n = rs.encodeSolomon(paddedminutialist)
            transformed_rd_conf[c] = {'checksums': checksums, 'rm': rm, 'k': k, 'n': n, 'hash': hash}

        hashes0 = [v['hash'] for c,v in transformed_rd_conf.iteritems()]
        hashes = getSplittedHashes(hashes0, 8)
        total_hash = getHashOfListOfTuples(hashes, pad=False)
        rm = int(len(hashes) / 2)
        rs = Rsolomon(_rm=rm)
        k, checksums, n = rs.encodeSolomon(hashes)
        hash_rd_conf = {'checksums': checksums, 'rm': rm, 'k': k, 'n': n, 'hash': total_hash}


        Cgenerated = {}
        cart2 = Cartesian(finputminutia2, W=cart.W, H=cart.H, bounds=cart.bounds, M=cart.M, keep_transforms=False) #notice that keep_transforms=False is false
        for c, minutia_list2 in cart2.Cprime_list_minutias.iteritems():
            print "cartesian block: ", c
            if c in cart.Cprime_list_minutias: #transformed cartesian block in the candidate template
                minutia_list1 = cart.Cprime_list_minutias[c] ##transformed cartesian block in the original template

                print "minutialist1: ", minutia_list1
                print "minutialist2: ", minutia_list2
                matches = set([])
                for p in sorted(minutia_list1):
                    for q in sorted(minutia_list2):
                        if q in matches:
                            continue
                        if self.checkIfMinutiaMatch(p, q):
                            if c not in Cgenerated:
                                Cgenerated[c] = set([])
                            Cgenerated[c].add(p)
                            matches.add(q)
                            break

        genhashes = []
        for c, minutia_list in Cgenerated.iteritems():
            print "generated:"
            print c, sorted(list(minutia_list))
            genhash = getHashOfListOfTuples(minutia_list)
            print "generated-hash: ", genhash
            if c in cart.Cprime_list_minutias:
                print "original-transformed:"
                print c, sorted(cart.Cprime_list_minutias[c])
                #orghash = getHashOfListOfTuples(cart.Cprime_list_minutias[c])
                orghash = transformed_rd_conf[c]['hash']
                print "original-transformed-hash: ", orghash

                if genhash == orghash:
                    print "generated and original minutias are equal"
                else:
                    print "generated and original minutias are NOT equal"
                    print "trying reed solomon error correction codes"
                    paddedminutialist = [','.join(map(str, ptuple)).rjust(8, '0') for ptuple in minutia_list]

                    rm = transformed_rd_conf[c]['rm']
                    checksums = transformed_rd_conf[c]['checksums']
                    k = transformed_rd_conf[c]['k']
                    n = transformed_rd_conf[c]['n']

                    leng = len(paddedminutialist)
                    if leng < k:
                        for index in range(0, k-leng): #if available minutia is less than data part then add dummy data
                            paddedminutialist.append('x' * 8)
                    elif leng > k: #if available minutia is more than data part then remove elements from the last
                        for index in range(0, leng-k):
                            del paddedminutialist[-1]

                    paddedminutialisttest = [','.join(map(str, ptuple)).rjust(8, '0') for ptuple in
                                             cart.Cprime_list_minutias[c]]

                    testSets(paddedminutialist, paddedminutialisttest)

                    rs = Rsolomon(_rm=rm)
                    print "to-be-decoded: ", paddedminutialist
                    decoded, corrections = rs.decodeSolomon(paddedminutialist, k, checksums, n)
                    print "corrections: ", corrections
                    print "decoded ", decoded
                    genhash = getHashOfListOfTuples(list(decoded), pad = False)

                    assert cmpT(decoded, paddedminutialisttest) == True, "Cannot decode!"
                    assert genhash==orghash, "Generated and original hashes are not equal after decode!"

                    print "hashes are equal"
            else:
                print "Cartesian number: ", c, " does not exists in the transformed form of the original fingerprint template!!"
            genhashes.append(genhash)

        genhashes = getSplittedHashes(genhashes, 8)
        k = hash_rd_conf['k']
        n = hash_rd_conf['n']
        leng = len(genhashes)
        if leng < k:
            for index in range(0, k - leng):  # if available minutia is less than data part then add dummy data
                genhashes.append('x' * 8)
        elif leng > k:  # if available minutia is more than data part then remove elements from the last
            for index in range(0, leng - k):
                del genhashes[-1]
        rs = Rsolomon(_rm=hash_rd_conf['rm'])
        print "to-be-decoded: ", genhashes
        decoded, corrections = rs.decodeSolomon(genhashes, k, hash_rd_conf['checksums'], n)
        print "corrections: ", corrections
        print "decoded ", decoded
        total_hash_decoded = getHashOfListOfTuples(list(decoded), pad=False)
        assert hash_rd_conf['hash'] == total_hash_decoded, "Decoded hash is not equal to the original hash"
        print "Original hash and decoded hash are equal: ", hash_rd_conf['hash'], " equals ", total_hash_decoded

    """
    WORKS
      does the matching on transformed cartesian blocks.
      performs reed solomon encode on the original cartesian blocks (does it for each block)
      keeps the original c(cartesian rect numbers) transformation on the candidate image
      if matches it recovers original minutia point using the original c(cartesian rect numbers) transformation on the candidate image from the candidate template
      it recovers all the original minutia points using reed solomon encoding/decoding
      checks hashes on the original cartesian blocks
    """
    def testWithRecoveringOriginalMinutiaPoints(self):
        """
        per each xprime,yprime in Cprime_list_minutia - also keep the index-order of the original x,y in the list when it
        went to the reed solomon encoding
        now rd works.
        calculateS the final hash
        """

        cart = Cartesian(finputminutia)
        #Reed solomon per cartesian blocks
        original_rd_conf = {} #keep
        #cart.Cprime_list_minutias#keep
        for c, minutia_list in cart.C_list_minutias.iteritems():
            paddedminutialist = [','.join( map(str, ptuple) ).rjust(8,'0') for ptuple in minutia_list]
            hash0 = getHashOfListOfTuples(paddedminutialist, pad=False)
            minutia_order = {}#for decodin g we need exact orders
            #if there is only half of matches then we can still recover them
            rm = int(len(paddedminutialist)/2)
            rs = Rsolomon(_rm = rm )
            k, checksums,n = rs.encodeSolomon(paddedminutialist)

            order = 0
            for padded_minutia in paddedminutialist:
                minutia_order[hash(padded_minutia)] = order
                order = order + 1
            original_rd_conf[c] = {'checksums': checksums, 'rm': rm, 'k': k, 'n': n, 'minutia_order': minutia_order, 'hash': hash0}

        hash_order = {}
        hashes0 = [v['hash'] for c, v in original_rd_conf.iteritems()]
        order = 0
        for hash0 in hashes0:
            hash_order[hash0] = order
            order = order + 1
        hashes = getSplittedHashes(hashes0, 8)
        total_hash = getHashOfListOfTuples(hashes, pad=False)
        rm = int(len(hashes) / 2)
        rs = Rsolomon(_rm=rm)
        k, checksums, n = rs.encodeSolomon(hashes)
        hash_rd_conf = {'checksums': checksums, 'rm': rm, 'k': k, 'n': n, 'hash': total_hash, 'hash_order': hash_order}

        Cgenerated = {}
        cart2 = Cartesian(finputminutia2, W=cart.W, H=cart.H, bounds=cart.bounds, M=cart.M, keep_transforms=True) #notice that keep_transforms=False is true

        for c, minutia_list2 in cart2.Cprime_list_minutias.iteritems():
            print "cartesian block: ", c
            if c in cart.Cprime_list_minutias: #candidate
                minutia_list1 = cart.Cprime_list_minutias[c] #original

                print "minutialist1: ", minutia_list1
                print "minutialist2: ", minutia_list2
                matches = set([])
                for p in sorted(minutia_list1):
                    for q in sorted(minutia_list2):
                        if q in matches:
                            continue
                        if self.checkIfMinutiaMatch(p, q):
                            #print(cart2.Creverse)
                            c_original = cart2.Creverse[q]
                            if c_original not in Cgenerated:
                                Cgenerated[c_original] = set([])
                            original_minutia = cart2.transformMinutiaToNewBlock(p, c, c_original )
                            Cgenerated[c_original].add( original_minutia )
                            matches.add(q)
                            break

        genhashes = []
        for c, minutia_list in Cgenerated.iteritems():
            print "generated:"
            print c, sorted(list(minutia_list))
            genhash = getHashOfListOfTuples(minutia_list)
            print "generated-hash: ", genhash

            if c not in cart.C_list_minutias:
                print "Cartesian number: ", c, " does not exists in original fingerprint template!!"
                continue

            print "original:"
            print c, sorted(cart.C_list_minutias[c])
            orghash = getHashOfListOfTuples(cart.C_list_minutias[c])
            print "original-hash: ", orghash

            if genhash == orghash:
                print "generated and original minutias are equal"
            else:
                print "generated and original minutias are NOT equal"
                print "trying reed solomon error correction codes"
                paddedminutialist = [','.join(map(str, ptuple)).rjust(8, '0') for ptuple in minutia_list]

                rm = original_rd_conf[c]['rm']
                checksums = original_rd_conf[c]['checksums']
                k = original_rd_conf[c]['k']
                n = original_rd_conf[c]['n']
                minutia_order = original_rd_conf[c]['minutia_order']

                new_paddedminutialist = []
                for index in range(0, k):  # if available minutia is less than data part then add dummy data
                    new_paddedminutialist.append('x' * 8)
                for padded_minutia in paddedminutialist:
                    hash_of = hash(padded_minutia)
                    if hash_of in minutia_order:
                        new_paddedminutialist[ minutia_order[hash_of] ] = padded_minutia
                paddedminutialist = new_paddedminutialist

                rs = Rsolomon(_rm=rm)
                print "to-be-decoded: ", paddedminutialist
                try:
                    decoded, corrections = rs.decodeSolomon(paddedminutialist, k, checksums, n)
                    print "corrections: ", corrections
                    print "decoded ", decoded
                    genhash = getHashOfListOfTuples(list(decoded), pad = False)

                    paddedminutialisttest = [','.join(map(str, ptuple)).rjust(8, '0') for ptuple in cart.C_list_minutias[c]]
                    assert cmpT(decoded, paddedminutialisttest) == True, "Cannot decode!"
                    assert genhash==orghash, "Generated and original hashes are not equal after decode!"

                    print "hashes are equal"
                except Exception as e:
                    print "hashes are not equal ! reed-solomon-cannot-decode:", str(e)


            genhashes.append(genhash)

        k = hash_rd_conf['k']
        k_number_of_original_hashes = k/3 # this is because we previously splitted 24 length of hashes to 8 length
        n = hash_rd_conf['n']
        new_genhashes = []
        for index in range(0, k_number_of_original_hashes):  # if available minutia is less than data part then add dummy data
            new_genhashes.append('x' * 24) # 24 is the length of the base64 hash
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
            print "Original hash and decoded hash are equal: ", hash_rd_conf['hash'], " equals ", total_hash_decoded
        except Exception as e:
            print "Original hash and decoded hash are NOT equal:  ! reed-solomon-cannot-decode:", str(e)

    """
     same as testWithRecoveringOriginalMinutiaPoints except we don't recover original minutia points
     we reed solomon with hashes of the original minutia points and we just recover the hashes
     this should be used. because it does not reveal minutia points from the original template
     also - the reed solomon encoding part needs to happens once before registering.
     the rest is matching part and re-generating the hash value
    """
    def testWithRecovering_HashOfOriginalMinutiaPoints(self, rm_divide = 2, prune_extra_rsolomon = False):
        cart = Cartesian(finputminutia)
        #Reed solomon per cartesian blocks
        original_rd_conf = {} #keep
        #cart.Cprime_list_minutias#keep
        for c, minutia_list in cart.C_list_minutias.iteritems():
            paddedminutialist = [','.join( map(str, ptuple) ).rjust(8,'0') for ptuple in minutia_list]

            """reed solomon per hashes of minutia points per cartesian blocks"""
            hashesof_paddedminutialist = [getHash_base64_grabFirst8(padded_minutia) for padded_minutia in paddedminutialist]
            hash_order = {}
            order = 0
            for hash0 in hashesof_paddedminutialist:
                hash_order[hash0] = order
                order = order + 1
            splitted_hashesof_paddedminutialist = getSplittedHashes(hashesof_paddedminutialist, 8) #need to do this because of max symsize of the reed solomon is 8
            total_hash_of_splitted_hashesof_paddedminutialist = getHashOfListOfTuples(splitted_hashesof_paddedminutialist, pad=False)
            rm = max(int(len(splitted_hashesof_paddedminutialist) / rm_divide), 2)
            rs = Rsolomon(_rm=rm)
            if len(splitted_hashesof_paddedminutialist) >rs.kmax:
                print "In encodeSolomon number of list: " + str(len(splitted_hashesof_paddedminutialist)) + " must be less than or equal to kmax: " + str(
                    rs.kmax)
                pass
            k, checksums, n = rs.encodeSolomon(splitted_hashesof_paddedminutialist, prune_extra_rsolomon)
            original_rd_conf[c] = {'checksums': checksums, 'rm': rm, 'k': k, 'n': n, 'hash_order': hash_order,
                                   'number_of_original_hashes': len(hash_order),
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
        rm = rm = max(int(len(hashes) / rm_divide), 2)
        rs = Rsolomon(_rm=rm)
        k, checksums, n = rs.encodeSolomon(hashes, prune_extra_rsolomon)
        hash_rd_conf = {'checksums': checksums, 'rm': rm, 'k': k, 'n': n,
                        'number_of_original_hashes': len(hash_order),
                        'hash': total_hash, 'hash_order': hash_order}

        Cgenerated = {}
        cart2 = Cartesian(finputminutia2, W=cart.W, H=cart.H, bounds=cart.bounds, M=cart.M, keep_transforms=True) #notice that keep_transforms=False is true

        for c, minutia_list2 in cart2.Cprime_list_minutias.iteritems():
            print "cartesian block: ", c
            if c in cart.Cprime_list_minutias: #candidate
                minutia_list1 = cart.Cprime_list_minutias[c] #original

                print "minutialist1: ", minutia_list1
                print "minutialist2: ", minutia_list2
                matches = set([])
                for p in sorted(minutia_list1):
                    for q in sorted(minutia_list2):
                        if q in matches:
                            continue
                        if self.checkIfMinutiaMatch(p, q):
                            #print(cart2.Creverse)
                            c_original = cart2.Creverse[q]
                            if c_original not in Cgenerated:
                                Cgenerated[c_original] = set([])
                            original_minutia = cart2.transformMinutiaToNewBlock(p, c, c_original, back_to_original=True )
                            Cgenerated[c_original].add( original_minutia )
                            matches.add(q)
                            break

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
                print "generated and original minutias are equal"
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
                for hash0 in hashesof_paddedminutialist:
                    if hash0 in hash_order:
                        new_hashesof_paddedminutialist[hash_order[hash0]] = hash0
                hashesof_paddedminutialist = new_hashesof_paddedminutialist
                splitted_hashesof_paddedminutialist = getSplittedHashes(hashesof_paddedminutialist, 8)

                rs = Rsolomon(_rm=rm)
                print "to-be-decoded: ", splitted_hashesof_paddedminutialist
                try:
                    decoded, corrections = rs.decodeSolomon(splitted_hashesof_paddedminutialist, k, checksums, n)
                    print "corrections: ", corrections
                    print "decoded ", decoded
                    genhash = getHashOfListOfTuples(list(decoded), pad = False)

                    assert genhash==orghash, "Generated and original hashes are not equal after decode!"

                    print "hashes are equal"
                except Exception as e:
                    print "hashes are not equal ! reed-solomon-cannot-decode:", str(e)
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
            print "Original hash and decoded hash are equal: ", hash_rd_conf['hash'], " equals ", total_hash_decoded
            return True
        except Exception as e:
            print "Original hash and decoded hash are NOT equal:  ! reed-solomon-cannot-decode:", str(e)
        return False


    def checkIfMinutiaMatch(self, p, q):
        rang1 = math.atan2(p[1], p[0])
        rang2 = math.atan2(q[1], q[0])
        if abs(p[0]- q[0]) <= xmatchthreshold and \
                abs(p[1] - q[1]) <= ymatchthreshold:
        #if        abs(rang1 - rang2) <= anglethreshold:
            return True
        return False

_prune_extra_rsolomon = False
xmatchthreshold = 12
ymatchthreshold = 12
anglethreshold = 0.1
_rm_divide = 2
ffiles = []
for i in range(1, 5):
    for j in range(1, 5):
        if i< 10:
            finputminutia='/Users/Admin/fvs/samples/DB2_B/out/10'+ str(i)+ '_'+str(j)+ '_out.txt'
            finputminutia = '/Users/Admin/fvs/samples/DB2_B/out/10' + str(i) + '_' + str(j) + '_out.txt'
        else:
            finputminutia = '/Users/Admin/fvs/samples/DB2_B/out/1' + str(i) + '_' + str(j) + '_out.txt'
            finputminutia = '/Users/Admin/fvs/samples/DB2_B/out/1' + str(i) + '_' + str(j) + '_out.txt'
        ffiles.append( {
            'i': i,
            'j': j,
            'finputminutia': finputminutia
        } )
errors = 0
totaltry = 0
tp = tn = fp = fn = 0.0
actual = False
tps =[]
fps =[]
fns =[]
tns =[]
for index1 in range(0, len(ffiles)):
    f1 = ffiles[index1]
    for index2 in range(0, len(ffiles)):
        if index1 == index2: continue
        f2 = ffiles[index2]
        finputminutia = f1['finputminutia']
        finputminutia2 = f2['finputminutia']
        actual = False
        if f1['i'] == f2['i']:
            actual = True
        print "Original: ", finputminutia
        print "Candidate: ", finputminutia2
        if not os.path.isfile(finputminutia) or not os.path.isfile(finputminutia2):
            print "Problem on files. Either them does not exists."
            continue
        match = False
        try:
            fvs = Fvs()
            match = fvs.testWithRecovering_HashOfOriginalMinutiaPoints(rm_divide=_rm_divide, prune_extra_rsolomon = _prune_extra_rsolomon)
            #exit()
        except Exception as ex:
            print "Exception: ", ex.message
            errors = errors+1
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

print "\nNum-try: ", totaltry, "\tNum errors: ", errors,
print "xmatchthreshold: ", xmatchthreshold, "\tymatchthreshold: ", xmatchthreshold, "rmdivide: ", _rm_divide
print "Precision: ", precision, "\tRecall: ", recall, "\tFmeasure: ", fmeasure, "\n"

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


