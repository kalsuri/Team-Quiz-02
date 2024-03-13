import numpy
import itertools
import math
import sympy.ntheory

########################################################################################
# The ASCIIPad takes a string (plaintext) as an input and gives a number as an output. #
# The number is represented in a ASCII like format that is padded.                     #
########################################################################################

def ASCIIPad(Message):
    newList = []
    messLen = len(Message)
    tempIter = messLen - 1
    while tempIter >= 0:
        newList.append(Message[tempIter])
        tempIter -= 1
    x = [100+ord(newList[i]) for i in range(messLen)];
    x = ZZ(x,1000);
    return(x);

def ASCIIPad(Mess,Mod):
    K = []
    for letter in Mess:
        K.append(ord(letter))
    L = Mod.ndigits()
    l = len(K)
    y = ZZ(math.floor(L/3))
    count = 0
    padded = []
    buffer = ""
    for numChar in K:
        numChar+=100
        buffer+=str(numChar)
        count+=1
        if count==y:
            padded.append(ZZ(buffer))
            buffer = ""
            count = 0
    if len(buffer)>0:
        padded.append(ZZ(buffer))
    return padded

#################################################################################################
# The input is a number and the output is the original message. If the input is not padded ASCII#
# version of a message it returns the value: "This is not a padded ASCII string"                #
#################################################################################################

def ASCIIDepad(Number):
    n = Number.ndigits() % 3;
    if (n > 0):
        print("This is not a padded ASCII string\n");
    else:
        L = [((Number - (Number % (1000^i)))/1000^i)%1000 - 100 for i in range(Number.ndigits()/3)];
        N = "";
        for i in range(Number.ndigits()/3):
            N = chr(L[i]) + N;
        return(N)

def isASCIIPadded(Number):
    N = ""
    n = Number.ndigits() % 3;
    if (n > 0):
        return False;
    L = [((Number - (Number % (1000^i)))/1000^i)%1000 - 100 for i in range(Number.ndigits()/3)];
    for i in range(Number.ndigits()/3):
        if L[i] < 0:
            return False
        if L[i] > 255:
            return False
    return True

def ECSearch (lowerbound, upperbound, Group):
    p = Group[2]
    a = Group[0]
    b = Group[1]
    if (4 * a ** 3 + 27 * b ** 2) % p == 0:
        print ("This is not an elliptic curve.")
    else:
        for j in itertools.count(lowerbound):
            if j==lowerbound-1 or j>upperbound or j>upperbound:
                return "not found"
            j=j%p
            if kronecker(j ** 3 + a * j + b, p) == 1:
                x = (j ** 3 + a * j + b) % p
                var('z')
                L=TonSh(x,p)
                y = L
                return([j,y])
             
def ECEmbed (Message, gp, tolerance):
    p = ZZ(math.floor(gp[2] / (tolerance + 1)))
    M = ASCIIPad(Message, p)
    packets = len(M)
    pointlist = [0]*packets
    for j in range(0, packets ):
        N = M[j]
        pointlist[j] = ECSearch(tolerance * N, tolerance * (N + 1) - 1, gp)
    return(pointlist)

def ECUnembed (pointlist, tolerance):
    k = len(pointlist)
    paddedasciilist=[0]*k
    for j in range(0, k ):
        pointlist[j][0]=ZZ(pointlist[j][0])
        toType = ZZ(QQ((pointlist[j][0])/tolerance).floor())
        paddedasciilist[j] = ((pointlist[j][0])/tolerance).floor()
    returnStr = ""
    for paddedItem in paddedasciilist:
        buffer = ASCIIDepad(paddedItem)
        returnStr+= buffer
    return returnStr

def TonSh (a, Prime):
    if kronecker(a, Prime) == -1:
        pass
        print ("{0} does not have a square root mod {1}".format(a, Prime))
        return None
    elif Prime % 4 == 3:
        x = power_mod(ZZ(a), ZZ((Prime+1)/4),ZZ(Prime))
        return(x)
    else:
        #################################################################################
        # Determine e so that Prime-1 = (2^e)*q for some odd number q                   # 
        #################################################################################

        e = ZZ(0)
        q = ZZ(Prime - 1)
        for j in itertools.count(1):
            if q % 2 == 0:
                e = ZZ(e + 1)
                q = ZZ(q / 2)
            else:
                break
        for i in range(1, 101):
            n = i
            if kronecker(n, Prime) == -1:
                break
        z = power_mod(ZZ(n),ZZ(q),ZZ(Prime))
        y = ZZ(z)
        r = ZZ(e)
        x = power_mod(ZZ(a),ZZ( (q-1)/2),ZZ( Prime))
        b = (a * x ** 2) % Prime
        x = (a * x) % Prime
        #for i in range(1, e + 1):
        for i in itertools.count(1):
            if b == 1:
                break
            else:
                for m in itertools.count(1):
                    t = power_mod(ZZ(b), ZZ(2^m) ,  ZZ(Prime))
                    if t == 1:
                        mm = m
                        break
                t = power_mod(ZZ(y), ZZ(2^(r - mm - 1)),ZZ(Prime))
                y = power_mod(ZZ(t), ZZ(2),ZZ(Prime))
                r = mm
                x = (x * t) % Prime
                b = (b * y) % Prime
        return(x)

def findptOrder(point,group):
    E = EllipticCurve(GF(group[2]),[group[0],group[1]])
    Ep = E(point[0],point[1])
    return Ep.additive_order()

def ECDouble (Point, Group):
    a = Group[0]
    b = Group[1]
    p = Group[2]
    if Point != []:
        x1 = Point[0]
        y1 = Point[1]
        if (4 * a ** 3 + 27 * b ** 2) % p == 0:
            print ("This is not an elliptic curve. ")
        elif Point != [] and y1 ** 2 % p != (x1 ** 3 + a * x1 + b) % p:
            print ("The point to double is not on the elliptic curve. ")
        elif y1 == 0:
            R = []
        else:
            s = mod((3 * x1 ** 2 + a) / y1 / 2, p)
            x = (s ** 2 - 2 * x1) % p
            y = (s * (x1 - x) - y1) % p
            R = [x,y]
    else:
        R = []
    return(R)

def ECAdd (Point1, Point2, Group):
    a = Group[0]
    b = Group[1]
    p = Group[2]
    if(Point1!=[]):
        x1 = Point1[0]
        y1 = Point1[1]
    if(Point2!=[]):
        x2 = Point2[0]
        y2 = Point2[1]
    if (4 * a ** 3 + 27 * b ** 2) % p == 0:
        print ("This is not an elliptic curve. ")
    elif Point1 != [] and y1 ** 2 % p != (x1 ** 3 + a * x1 + b) % p:
        print ("Point 1 is not on the elliptic curve. \n")
    elif Point2 != [] and y2 ** 2 % p != (x2 ** 3 + a * x2 + b) % p:
        print ("Point 2 is not on the elliptic curve. \n")
    else:
        if Point1 == []:
            R = Point2
        elif Point2 == []:
            R = Point1
        else:
            if x1 == x2 and 0 == (y1 + y2) % p:
                R = []
            elif x1 == x2 and y1 == y2:
                R = ECDouble(Point1, Group)
            else:
                s = ((y1 - y2) / (x1 - x2)) % p
                x = (s ** 2 - x1 - x2) % p
                y = (s * (x1 - x) - y1) % p
                R = [x,y]
    return(R)

def ECTimes (Point, scalar, Group):
    ECIdentity=[]
    if Point == ECIdentity or scalar == 0:
        return(ECIdentity)
    else:
        E=EllipticCurve(GF(Group[2]),[Group[0],Group[1]])
        EPoint = E(Point[0],Point[1])
        cgret = scalar*EPoint
        if(cgret[0]==0 and cgret[1]==1 and cgret[2]==0):
            return ECIdentity
        return([cgret[0],cgret[1]])
    
def ECInverse (Point, Group):
    if Point == []:
        return(Point)
    else:
        p = Group[2]
        x = Point[0]
        y = Point[1]
        return([x,(p - y) % p])