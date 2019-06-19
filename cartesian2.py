from utils import *
import math

__all__ = ['Cartesian']

_H = 5.0
_W = 5.0
class Cartesian:
    def xy_to_Cxy(self, x, y):
        cx_ = ((x - self.xmin) / self.dx) + 1.0
        cy_ = ((y - self.ymin) / self.dy) + 1.0
        # print("cx_, cy_", cx_, cy_)
        cx = min(int(self.W), int(math.floor(cx_)))
        cy = min(int(self.H), int(math.floor(cy_)))
        cx = min(self.W, max(1, cx))
        cy = min(self.H, max(1, cy))
        if cx <=0 or cy <= 0 or cx > self.W or cy > self.H:
            print "cx cy values are incorrect: ", cx, cy
        return cx, cy
    def Cxy_to_c(self, cx, cy):
        #c = int(((self.H-cy) * self.H) + cx )
        c = int(((cy-1) * self.W) + cx)
        cx2, cy2 = self.C_to_Cxy(c)
        if (int(cx), int(cy)) != (int(cx2), int(cy2)):
            print "not equal"
            raise Exception("(cx, cy) != (cx2, cy2) ", (cx, cy) != (cx2, cy2))
        return c
    def C_to_Cxy(self, c):
        cx = (c % self.W) if (c % self.W) > 0 else self.W
        cy = math.ceil(c / self.W)
        return cx, cy
    def __init__(self, des, W=_W, H=_H, bounds = None, M = None, keep_transforms = False):
        nminutiaelist=[] # integer minutia list
        minutiaexpositions=[]
        minutiaeypositions=[]

        # Find extremities of bounding box
        if bounds is not None:
            xmin = bounds[0]
            xmax = bounds[1]
            ymin = bounds[2]
            ymax = bounds[3]

        nminutiaelist = des['minutiaes']
        for minutia in des['minutiaes']:
            if bounds is not None:
                x = minutia[0]
                y = minutia[1]
                type = minutia[2]
                some_not_in_bount = False
                if x < xmin or y < ymin or x > xmax or y > ymax:
                    #print "minutia: ", minutia, " is out of the template bounding box!! excluding the minutia."
                    some_not_in_bount = True
                    continue
                if some_not_in_bount:
                    print "Some minutias are out of the template bounding box!! excluding them"
            minutiaexpositions.append(minutia[0])
            minutiaeypositions.append(minutia[1])

        # Find extremities of bounding box
        boundsxadd = 0 #to be secure
        boundsyadd = 0
        if bounds is None:
            xmax = max(minutiaexpositions)+boundsxadd
            ymax = max(minutiaeypositions)+boundsyadd
            xmin = min(minutiaexpositions)-boundsxadd
            ymin = min(minutiaeypositions)-boundsyadd
            bounds = (xmin, xmax, ymin, ymax)
        # Find the step lengths for the grid
        dx = (xmax-xmin)/W
        dy = (ymax-ymin)/H

        self.W = W
        self.H = H
        self.bounds = bounds
        self.xmax = xmax
        self.ymax = ymax
        self.xmin = xmin
        self.ymin = ymin
        # Find the step lengths for the grid
        self.dx = dx
        self.dy = dy

        R = {} #per a minutia point(x,y), its representative cartesian rectangle point
        C = {} #per a minutia point(x,y), its representative cartesian rectangle number
        C_list_minutias = {} # per each cartesian rectangle number the list of minutia points
        lenC = int(H*W)
        LC = range(1, lenC+1) #list of cartesian block numbers from 1 to HxW
        if M is None:
            M = genTransformationMatrix(lenC)
        Cprime = {} #transformed cartesian blocks, per a minutia point(x,y), its representative transformed cartesian rectangle number
        Cprime_list_minutias = {} # per each transformed cartesian rectangle number the list of minutia points, kind of inverse of C
        T = {} #cartesian transformation
        Treverse = {} #this will only be used when matching to get the original blocks with points
        Tnminutiaelist=[] # transformed minutia list according to the cartesian transformation
        Creverse = {}  # per a transformed minutia point(x,y), its origin cartesian rectangle number

        #generate cartesian blocks and minutia lists for each cartesian blocks
        for minutia in nminutiaelist:
            x = minutia[0]
            y = minutia[1]
            cx, cy = self.xy_to_Cxy(x, y)
            R[x, y] = [cx, cy]  # mapping it to rectangular points
            c = self.Cxy_to_c(cx, cy)
            if c > lenC:
                c = self.Cxy_to_c(cx, cy)
                raise Exception("Rectangle number: ", c, " cannot be bigger than max rectangle number: ", lenC)

            C[x, y] = c
            if C[x, y] not in C_list_minutias:
                C_list_minutias[C[x, y]] = []
            C_list_minutias[C[x, y]].append(minutia)

        #sort minutia list
        for c in C_list_minutias.keys():
            C_list_minutias[c] = sorted(C_list_minutias[c])

        #generate transformed cartesian blocks
        Cprime = matrixMultiply([LC], M)[0]
        #print "Cprime:"
        print Cprime
        for index in range(0, len(Cprime)):
            c = index+1
            cprime = Cprime[index]
            #print "cell: ",c , " => ", cprime
            T[c] = cprime
            Treverse[cprime] = c
        #print "Cartesian transformation: "
        #print (T)
        #print "Cartesian reverse transformation: "
        #print (Treverse)

        #generate Cprime_list_minutias
        #list of minutia points per transformed catesian block numbers
        #minutia points are transformed accordingly

        #print "C values:"
        #print set(C.values())

        for minutia in nminutiaelist:
            c = C[minutia[0], minutia[1]] #original cartesian block
            cprime = T[c] #transformed cartesian block
            if cprime not in Cprime_list_minutias:
                Cprime_list_minutias[cprime] = []

            #now transform x and y location for the minutia
            transformed_minutia = self.transformMinutiaToNewBlock(minutia, c, cprime)
            minutia_back_from_transformed = self.transformMinutiaToNewBlock(transformed_minutia, cprime, c,
                                                                            back_to_original=True)
            if self.checkMinutiasEqual(minutia, minutia_back_from_transformed) == False:
                print minutia, transformed_minutia, minutia_back_from_transformed
                transformed_minutia = self.transformMinutiaToNewBlock(minutia, c, cprime)
                minutia_back_from_transformed = self.transformMinutiaToNewBlock(transformed_minutia, cprime, c,
                                                                                back_to_original=True)
                raise Exception("minutial cartesian transformation - vice-versa, not equal")

            Cprime_list_minutias[cprime].append(transformed_minutia)
            try:
                #print hash(tuple(transformed_minutia))
                Creverse[ hash(tuple(transformed_minutia)) ] = c
            except Exception as ex:
                raise ex

            # sort minutia prime list
            for cprime in Cprime_list_minutias.keys():
                Cprime_list_minutias[c] = sorted(Cprime_list_minutias[cprime])

        #print "Cprime values:"
        #print set(Cprime)

        #print "C_list_minutias lengths"
        #for c in C_list_minutias:
        #    print len(C_list_minutias[c])
        #print "\nCprime_list_minutias lengths"
        #for cprime in Cprime_list_minutias:
        #    print cprime, ":", len(Cprime_list_minutias[cprime])

        self.R = R
        self.C = C
        self.C_list_minutias = C_list_minutias
        self.M = M
        self.Cprime = Cprime
        self.Cprime_list_minutias = Cprime_list_minutias
        if keep_transforms:
            self.T = T
            self.Treverse = Treverse
            self.Creverse = Creverse
        else:
            self.T = None
            self.Treverse = None
            self.Creverse = None

        """a=[
            [1,2,3,4]
        ]
        b=[
            [0,0,0,0],
            [0,1,0,0],
            [1,0,0,1],
            [0,0,1,0]
        ]
        matrixMultiply(a,b)
        """

    """
       transform x and y location for the minutia
       c is the old cartesian block
       cprime is the new one. transform pre_minutia(x, y) accordingly
    """

    def transformMinutiaToNewBlock(self, pre_minutia, c, cprime, back_to_original=False):
        if c == cprime:
            return (int(pre_minutia[0]), int(pre_minutia[1]), pre_minutia[2])
        c = float(c)
        cprime = float(cprime)
        x = float(pre_minutia[0])
        y = float(pre_minutia[1])
        type = pre_minutia[2]

        cx, cy = self.C_to_Cxy(c)
        cxprime, cyprime = self.C_to_Cxy(cprime)

        xprime0 = x + (cxprime - cx) * self.dx
        yprime0 = y + (cyprime - cy) * self.dy

        """xprime = int(xprime0)
        yprime = int(yprime0)
        xprime2 = math.ceil(xprime0)
        yprime2 = math.ceil(yprime0)
        xprime3 = math.floor(xprime0)
        yprime3 = math.floor(yprime0)"""

        # xprime4 = round(xprime0)
        # yprime4 = round(yprime0)
        # xprime = int(xprime4)
        # yprime = int(yprime4)
        # return (xprime, yprime)
        if not back_to_original:
            xprime = xprime0
            yprime = yprime0
        else:
            xprime = int(round(xprime0))
            yprime = int(round(yprime0))
        ret = [xprime, yprime, type]
        return tuple(ret)

    def checkMinutiasEqual(self, min1, min2):
        if min1[0] == min2[0] and \
                min1[1] == min2[1] and \
                min1[2] == min2[2]:
            return True
        return False