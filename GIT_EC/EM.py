#Campo de enteros modulo n
from math import ceil
import json
import os

base_dir = os.path.dirname(os.path.abspath(__file__))
cache_f = os.path.join(base_dir, "qr")

os.makedirs(cache_f, exist_ok=True)

def qr_file(p):
    return os.path.join(cache_f, f"qr_mod_{p}.json")


def qrld(p):
    pth = qr_file(p)
    if os.path.exists(pth):
        with open(pth, "r") as f:
            data = json.load(f)
            return set(data)  
    return None


def qrsv(p, residuos):
    pth = qr_file(p)
    tmp = pth + ".tmp"

    with open(tmp, "w") as f:
        json.dump(list(residuos), f)  

    os.replace(tmp, pth)


def get_qr(p):
    qr = qrld(p)

    if qr is None:
        # calcular solo una vez
        qr = {
            a for a in range(p)
            if a == 0 or pow(a, (p - 1) // 2, p) == 1
        }
        qrsv(p, qr)

    return qr

def Bin(k):
    #Funcion que determina la representacion binaria de un numero entero k
    return [ 1 if i == "1" else 0 for i in bin(k)[2:]][::-1]

def d(x,y):
    #division generalizada
    q = 0
    while x >= y :
        x -= y
        q += 1
    return([q,x])

def gcd(x,y,d):
    #maximo comun divisor dada una funcion de division
    if y == 0:
        return x

    else:
        _, r = d(x,y)
        return gcd(y,r,d)

def inv(a,p):
    #inversa de un elemento sobre un campo de car p
    u, v = a, p
    x, y = 1, 0

    while u != 1:
        q = int(v/u)
        r, z = v-q*u, y-x*q
        v, u, y, x = u, r, x, z
    return x

def lnd(x,a):
    #logaritmo discreto de x en base a modulo n
    m = ceil((a.n)**(1/2))
    J = [a**j for j in range(m)]
    ai = a**(-m)
    y = x
    for i in range(m):
        try:
            j = J.index(y)
            return (i*m+j)
        except:
            y = y*ai
    return None

def ec(x):
    #criterio de euler para residuos cuadraticos
    if x != 0:
        return True if x**((x.n-1)//2) == 1 else False
    else:
        return True

class mod :
    #Esta clase define el conjunto de enteros modular con sus respectivas operaciones
    def __init__(self, k, n, json = False):
        if json:
            self.k = k
            self.n = n
        else:
            if isinstance(k,mod):
                self.k = k.k
                self.n = k.n
            else:
                self.k = k%n
                self.n = n
                
    def es_residuo_cuadratico(self):
        return self.k in get_qr(self.n)
    ######## METODOS ########
    
    #REPRESENTACION

    def __repr__(self):
        return f"{self.k} mod {self.n}"    

    def __str__(self):
        return f"{self.k} mod {self.n}"
    
    def __int__(self):
        return self.k
    
    #BOOLEANOS
    
    def __eq__(self,other):
        if isinstance(other,mod):
            return self.k == other.k if other.n == self.n else False
        
        else:
            try:
                return self.k == other % self.n
            except:
                return False
    
    def __bool__(self):
        return bool(self.k)
    
    ######## OPERACIONES NUMERICAS ########

    #SUMA
    
    def __add__(self,other):
        if isinstance(other,mod):
            if self.n == other.n :
                return mod(self.k+other.k, self.n)
            
            else:
                raise ValueError(f"{self.k} es modulo {self.n} y {other.k} es modulo {other.k}")
            
        elif isinstance(other,int) :
            return mod(self.k + other, self.n)
        
        else:
            raise ValueError("solo se aceptan operaciones con enteros u enteros modulares")
        
    def __radd__(self,other):
        return self + other
    
    def __iadd__(self,other):
        self = self + other
        return self
    
    #RESTA
    def __neg__(self):
        return mod((-self.k),self.n)
    
    def __sub__(self,other):
        if isinstance(other,mod):
            if self.n == other.n :
                return self + (-other)
            
            else:
                raise ValueError(f"{self.k} es modulo {self.n} y {other.k} es modulo {other.k}")
            
        elif isinstance(other,int) :
            return self + (-other)
        
        else:
            raise ValueError("solo se aceptan operaciones con enteros u enteros modulares")
    
    def __rsub__(self,other):
        return self - other
        
    def __isub__(self,other):
        self = self - other
        return self
    
    #MULTIPLICACION
    
    def __mul__(self,other):
        if isinstance(other,mod):
            
            if self.n == other.n :
                return mod(self.k*other.k, self.n)
            
            else:
                raise ValueError(f"{self.k} es modulo {self.n} y {other.k} es modulo {other.k}")
            
        elif isinstance(other,int) :
            return mod(self.k * other, self.n)
        
        else:
            raise ValueError("solo se aceptan operaciones con enteros u enteros modulares")
    
    def __rmul__(self,other):
        return self * other
    
    def __imul__(self,other):
        self = self * other
        return self
        
    #DIVISION
    def __truediv__(self,other):
        if isinstance(other,mod) :

            if self.n == other.n:

                if gcd(other.k,other.n,d) == 1:
                    return self * mod(inv(other.k,other.n),other.n)

                else:
                    raise ValueError(f"No hay inverso para {other}")

            else:
                raise ValueError(f"{self.k} es modulo {self.n} y {other.k} es modulo {other.k}")
            
        elif isinstance(other,int):
            if gcd(other,self.n,d) == 1:
                return self * mod(inv(other,self.n),self.n)
            
            else:
                raise ValueError("No hay inverso")
            
        else:
            raise ValueError("no es posible hacer esta operacion con elementos diferentes de enteros modulares")

    #EXPONENCIACION

    def __pow__(self,other):
        try:
            v = mod(1,self.n)
            rb = Bin(other)

            if other >= 0 :
                vp = self

            elif gcd(self.k,self.n,d) == 1:
                vp = mod(inv(self.k,self.n),self.n)

            else:
                raise TypeError

            if any(rb):
                for i in rb:
                    v = v*vp if i == 1 else v
                    vp *= vp
                return v

            else:
                return v

        except:
            raise ValueError if not isinstance(other,int) else TypeError

    #Funciones

    def sqrt(self):
        if self:
            if ec(self):
                q = self.n-1
                z = mod(2,self.n)
                S = 0
                Q = q

                while True:
                    if Q%2 == 0:
                        S += 1
                        Q = Q//2
                    else:
                        break
                while ec(z):
                    z += 1

                M = S
                c = z**Q
                t = self**Q
                R = self**((Q+1)//2)

                while t!= 1:
                    for i in range(1,M):
                        if t**(2**i) == 1:
                            b = c**(2**(M-i-1))
                            M = i
                            c = b**2
                            t = t*c
                            R = R*b                        
                            break
                return R

            else:
                raise ValueError("no existe raiz")
        else:
            return 0

