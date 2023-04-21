# to exit
import sys

# to parse the input
import argparse

# to estimate the RLWE hardness
from latticeEstimator.estimator import *

# to write the parameters into a csv file
import csv

# to check if an integer is prime
import miller_rabin

# to compute binary logarithms of large integers
from sage.functions.log import logb

def log2(x):
    return numerical_approx(logb(x,2))

# to silently run functions
class silent(object):
    def __init__(self,stdout = None, stderr = None):
        self.devnull = open(os.devnull,'w')
        self._stdout = stdout or self.devnull or sys.stdout
        self._stderr = stderr or self.devnull or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
        self.devnull.close()


# PARSING INPUT PARAMETERS
parser = argparse.ArgumentParser(description='Finds a secure set of parameters.')

parser.add_argument('lamb', type=int, help="Security parameter.")
parser.add_argument('n', type=int, help="Dimension of the lattice.")
parser.add_argument('q', type=int, help="Prime modulus.")
parser.add_argument('d', type=int, help="Number of splitting factors of x^n+1 mod q.")
parser.add_argument('-o','--output', type=str, help="Output filename.")
args = parser.parse_args()

lamb, n, q, d = [args.lamb, args.n, args.q, args.d]

# INPUT PARAMETERS TESTS
if not (n and (not(n & (n - 1))) ):
    print("n = {0:d} is not a power of 2".format(n))
    sys.exit(1)

if not miller_rabin.miller_rabin(int(q)):
    print("q = {0:d} is not a prime number".format(q))
    sys.exit(1)

if not ((d and (not(d & (d - 1)))) and n>=d and d>1):
    print("d = {0:d} is not a power of 2 such that n ≥ d > 1".format(d))
    sys.exit(1)

if not (q >= 37):
    print("q = {0:d} is not greater or equal than 37".format(q))
    sys.exit(1)

if (q - (2*d+1))%(4*d):
    print("q = {0:d} is not equivalent to 2d+1 mod 4d for d={1:d}".format(q,d))
    sys.exit(1)

# AUXILIARY FUNCTIONS
def vecBoundedPrToBoundedPr(b,d):
    # outputs a such that some event happens at least once in d samples
    # with a probability of at most 2^-b if the probability of it
    # in a single sample is at most 2^-a
    return b + log2(d)

def sigmaFromB(a,B):
    # outputs sigma such that a sample from D_sigma
    # has absolute value greater than B
    # with a probability lower or equal to 2^-a
    return B*math.sqrt(log2(math.e)/(2*(a+1)))

def bitsec(params,estimateLevel=1):
    # estimates the hardness of a RLWE problem
    if estimateLevel == 0:
        try: x = LWE.primal_usvp(params)
        except: x = {'rop': float('inf')}
        try: y = LWE.dual(params)
        except: y =  {'rop': float('inf')}
        if math.isinf(min(x['rop'],y['rop'])): return 0
        bits = log2(min(x['rop'],y['rop']))
    elif estimateLevel == 1:
        try: x = LWE.primal_usvp(params)
        except: x = {'rop': float('inf')}
        try: y = LWE.dual_hybrid(params, mitm_optimization=True)
        except: y = {'rop': float('inf')}
        try: z = LWE.primal_bdd(params, red_shape_model="gsa")
        except: z = {'rop': float('inf')}
        if math.isinf(min(x['rop'],y['rop'],z['rop'])): return 0
        bits = log2(min(x['rop'],y['rop'],z['rop']))
    else:
        with silent():
            try: x = LWE.estimate(params)
            except: x = {'Error':{'rop': float('inf')}}
            if math.isinf(min([x[y]['rop'] for y in x])): return 0
            bits = log2(min([x[y]['rop'] for y in x]))
    return bits

# FUNCTIONS COMPUTING OPTIMAL PARAMETERS
def bestDeltaOL(lamb,q):
    delta = ceil(lamb/(logb(2*q/(q+1),2)))
    while not ((2*lamb - 1)*log2((2*lamb - 1)/(delta*(1-1/q))) + (delta - 2*lamb + 1)*log2(q*(delta-2*lamb+1)/delta) >= 2*lamb):
        delta += 1
    return delta

def bestDeltaM(lamb,q):
    delta = ceil(lamb/(logb(2*q*q/(q*q+3*q-2),2)))
    while not ((2*lamb - 1)*log2((2*lamb - 1)/(delta*(1-1/q))) + (delta - 2*lamb + 1)*log2(q*(delta-2*lamb+1)/delta) >= 2*lamb):
        delta += 1
    return delta

def bestK(lamb,n,q,d,B):
    return ceil((lamb+2*n*log2(q))/(n*(log2(q)/d - log2(4*B-1))))

def bestSigma(lamb,n,q,d,B,k,estimateLevel=1):
    sigma = sigmaFromB(a=vecBoundedPrToBoundedPr(b=lamb,d=k*n),B=(B-1))
    # we truncate it to two decimals
    sigma = float(floor(sigma*100)/100)
    params = LWE.Parameters(n=n, q=q, Xs=ND.UniformMod(q),Xe=ND.DiscreteGaussian(sigma),m=k*n)
    bits = bitsec(params,estimateLevel)
    if bits >= lamb: return sigma
    else: return None

def bestB(lamb,n,q,d,estimateLevel=1,minB=2):
    B = minB
    sigma = None
    while (sigma is None) and (d < log2(q)/(log2(4*B-1))):
        print("λ={0:d} | n={1:<5d} | q~2^{2:<3d} | d={3:<5d} | B=2^{4:<3d}"
                              .format(lamb,n,math.ceil(log2(q)),d,int(log2(B))),
                              end="\r")
        k = bestK(lamb,n,q,d,B)
        sigma = bestSigma(lamb,n,q,d,B,k,estimateLevel)
        if sigma is None:
            B *= 2
    if not (d < log2(q)/(log2(4*B-1))):
        return None
    else:
        return B

# COMPUTING THE PARAMETERS
B = 2
for level in range(3):
    B = bestB(lamb,n,q,d,estimateLevel=level,minB=B)
    if B is None:
        print("\nNot found")
        break

if B is not None:
    k = bestK(lamb,n,q,d,B)
    sigma = bestSigma(lamb,n,q,d,B,k,estimateLevel=2)
    deltaOL = bestDeltaOL(lamb,q)
    deltaM = bestDeltaM(lamb,q)
    print([lamb,n,q,d,k,sigma,B,deltaOL,deltaM])
    with open(args.output, 'a') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow([lamb,n,q,d,k,sigma,B,deltaOL,deltaM])
else:
    sys.exit(1)
