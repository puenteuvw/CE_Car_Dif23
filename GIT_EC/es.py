#Esta programa busca crear una clase para trabajar con superficies cubicas generadas por una curva eliptica no singular sobre un campo discreto
from EM import mod
from random import randint
from functools import partial, reduce

#Ley de composicion

def df(f,x,h=0.001):
    return (f(x+h)-f(x))/h

def ddf(f,x,h=0.001):
    return (df(f,x+h)-df(f,x))/h

def nr(f):
    x0, x1 = 0, 0
    while True:
        x1 += -f(x0)/df(f,x0)
        if abs(x1-x0) < 0.01:
            break
    return x1

Bin = lambda k : bin(k)[2:][::-1]

"""
def dr(x,ls):
    "determina si x es generado por ls y su representacion"
    s = x
    d = []
    
    for l in ls:
            while s%l == 0 :
                    s = s/l
                    d.append(l)
    return (d,True if s == 1 else False)

def ds(X):
    "determina una representacion de X en terminos de elementos primos"
    p = 2
    P = [p]

    while True:
        p += 1
        s = 0                       #Flag off
        for q in P:
                if p%q == 0:
                        s = 1       #Flag on
                        break

        if s == 1:          #elemento no primo
            continue

        else:               #elemento primo
            P.append(p)

            l, z = dr(X,P)

            if z:           #comprobacion de factorizacion total
                return(l)

#funcion de rapida exponenciacion 
def exp(N):
    def f(g,k):
        v = 1
        rb = Bin(k)
        vp = g%N

        for i in rb:
            v = v*vp if i=="1" else v
            vp = (vp**2)%N
        return v
    return f
"""

trns = lambda x : [list(row) for row in zip(*x)]

def v(X):
    #Impresion de la matriz X
    for x in X:
        print("")
        for i in x:
            print("\t",i,end = "")

#Funcion signo, du
sgn = lambda x : -1 if x<0 else (1 if x >0 else 0)

#Descomprime
uz = lambda D : [l for p in D for l in p ]

###### Producto de matrices ######
PR = lambda A, B: [[a*b for b in trns(B) ] for a in A ]

class vec:
    def __init__(self,X):
        self.cmp = list(X)

    def __repr__(self):
        return f"{self.cmp}"

    def __str__(self):
        return f"{self.cmp}"

    def __iter__(self):
        return iter(self.cmp)
    
    def __getitem__(self,index):
        if isinstance(index, int) or isinstance(index,slice):
            return self.cmp[index]
        
        else:
            raise ValueError("El indice no es valido")

    def __add__(self,X):
        return vec([x+y for x,y in zip(self,X)])

    def __radd__(self,K):
        if hasattr(K,"__int__"):
            return vec([K+x for x in self])

    def __iadd__(self,X):
        return self + X

    def __neg__(self):
        return vec([-x for x in self])

    def __sub__(self,X):
        return self+(-X)

    def __isub__(self,X):
        return self-X

    def __mul__(self,X):
        return sum([x*y for x,y in zip(self,X)])

    def __rmul__(self,K):
        if hasattr(K,"__int__"):
            return vec([K*x for x in self])
        else:
            raise ValueError("valor no admisible")

class es:
    def __init__(self, r = None, h = None, s = None, b = None, c = None, p = None, locus = None, orbs = None, json = False):
        if json:
            self.r, self.h, self.s, self.b, self.c = r, h, s, b, c
            self.car, self.locus, self.orbs = p, locus, orbs
        else:
            if any([x == None for x in [r, h, s, b, c, p]]):
                raise ValueError("Error 1")
            
            elif ((pbc := 4*(b**3) + 27*(c**2)) == 0) or ((QD := h**2-4*r*s) == 0) or (mod(QD,p).es_residuo_cuadratico()) or not(bool(r) and bool(s)):
                if pbc == 0:
                    raise ValueError("Error 2")
                elif QD == 0:
                    raise ValueError("Error 3")
                elif mod(QD,p).es_residuo_cuadratico():
                    raise ValueError("Error 4")
                elif not(bool(r) and bool(s)):
                    raise ValueError("Error 5")
            
            else:
                self.r = mod(r,p) 
                self.h = mod(h,p)
                self.s = mod(s,p)
                self.b = mod(b,p)
                self.c = mod(c,p)
                self.car = p

                self.locus = [(x,y,z) for z in range(p) for y in range(p) for x in range(p) if self.behoren((x,y,z))]+["inf"]
                self.orbs = []

    def __repr__(self):
        return f"{self.r.k}z^2 + {self.h.k}zy + {self.s.k}y^2 = x^3 + {self.b.k}x + {self.c.k} mod {self.car}"

    def __str__(self):
        return f"{self.r.k}z^2 + {self.h.k}zy + {self.s.k}y^2 = x^3 + {self.b.k}x + {self.c.k} mod {self.car}"

    def oneness(self):
        #identidad para el formato json
        return f"{self.r.k}_{self.h.k}_{self.s.k}_{self.b.k}_{self.c.k}_{self.car}"

    def behoren(self,X):
        if X!= "inf":
            x,y,z = X
            return True if self.r*(z**2) + self.h*z*y + self.s*(y**2) == x**3 + self.b*x + self.c else False
        else:
            return True
    
    def rndm(self):
        n = randint(0,len(self.locus)-1)
        return self.locus[n]
    
    def __len__(self):
        return len(self.locus)

    def sm(self, P, Q, bhrn = True):
        if bhrn:
            if(P == "inf" or Q == "inf"):
                return P if Q == "inf" else Q
            else:
                if P!= Q:
                    if P[0] != Q[0]:
                        P, Q = list(map((f := partial(mod,n = self.car)),P)), list(map(f,Q))
                        x1, y1, z1 = P
                        M = [a2-a1 for a2,a1 in zip(Q,P)]
                        m1, m2, m3 = M
                        x3 = (self.r*(m3**2) + self.h*m2*m3 + self.s*(m2**2) )/(m1**2) - x1 - Q[0]
                        y3 = (m2/m1)*(x1-x3) - y1
                        z3 = (m3/m1)*(x1-x3) - z1
                        return ([x3,y3,z3])
                    else:
                        return "inf"
                else:
                    x1, y1, z1 = P
                    if (y1 == 0) and (z1 == 0):
                        return "inf"
                    else:
                        x1, y1, z1 = mod(x1,self.car), mod(y1,self.car), mod(z1,self.car)
                        x3 = ((m := (3*(x1**2))+self.b)**2)/(4*(d := self.r*(z1**2)+self.s*(y1**2)+self.h*z1*y1)) - 2*x1
                        y3 = y1*((m/(2*d))*(x1-x3)-1)
                        z3 = z1*((m/(2*d))*(x1-x3)-1)
                        return([x3,y3,z3])

        else:
            if all(map(self.behoren,[P,Q])):
                if(P == "inf" or Q == "inf"):
                    return P if Q == "inf" else Q
                else:
                    if P!= Q:
                        if P[0] != Q[0]:
                            P, Q = list(map((f := partial(mod,n = self.car)),P)), list(map(f,Q))
                            x1, y1, z1 = P
                            M = [a2-a1 for a2,a1 in zip(Q,P)]
                            m1, m2, m3 = M
                            x3 = (self.r*(m3**2) + self.h*m2*m3 + self.s*(m2**2) )/(m1**2) - x1 - Q[0]
                            y3 = (m2/m1)*(x1-x3) - y1
                            z3 = (m3/m1)*(x1-x3) - z1
                            return ([x3,y3,z3])
                        else:
                            return "inf"
                        
                    else:

                        x1, y1, z1 = P
                        if (y1 == 0) and (z1 == 0):
                            return "inf"
                        else:
                            x1, y1, z1 = mod(x1,self.car), mod(y1,self.car), mod(z1,self.car)
                            x3 = ((m := (3*(x1**2))+self.b)**2)/(4*(d := self.r*(z1**2)+self.s*(y1**2)+self.h*z1*y1)) - 2*x1
                            y3 = y1*((m/(2*d))*(x1-x3)-1)
                            z3 = z1*((m/(2*d))*(x1-x3)-1)
                            return([x3,y3,z3])

    def pn(self,n,P):
        if n < 0:
            n *= -1
            P = P[0]+[-x for x in P[1:]]
        if n == 1:
            return P
        elif n == 2 :
            return self.sm(P,P)
        elif n > 2 :
            eb = bin(n)[:1:-1]
            bc = []
            m = P
            for j in eb[1:]:
                m = self.sm(m,m)
                if j!= "0":
                    bc.append(m)
            d = ([P] if eb[0] != "0" else []) + bc
            return reduce(self.sm,d)
        else:
            raise ValueError("el parametro n no es permisible")
    
    def orb(self,P,R = False):
        if P != "inf":
                if R:
                    Q = self.sm(P,P)
                    O = [P,Q]
                    k = 2
                    while Q != "inf":
                        Q = self.sm(Q,P)
                        O.append(Q)
                        k += 1
                    return k,O
                else:
                    Q = self.sm(P,P)
                    k = 2
                    while Q != "inf":
                        Q = self.sm(Q,P)
                        k += 1
                    return k
        else:
            raise ValueError("no existe la orbita del infinito")