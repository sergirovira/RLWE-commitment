import math

def commitment(param):
    lamb,n,q,d,k,sigma,B,deltaOL,deltaM = param
    return n*math.ceil(math.log2(q))*k

def opening(param):
    lamb,n,q,d,k,sigma,B,deltaOL,deltaM = param
    com_size = 256
    com_op_size = 8*math.ceil(lamb/8.)
    seed_size = 8*math.ceil(lamb/8.)
    # ch_size = math.ceil(math.log2(q))
    res = 0.0
    res += 2*com_size # commitments c1 and c2
    # res += ch_size # alpha challenge
    res += (math.log2(B)+1)*2*n*k*math.ceil(math.log2(q)) # g
    # res += 1 # second challenge bit
    # if b challenge is 0
    res += 0.5*seed_size # pi
    res += 0.5*(k*n*math.ceil(math.log2(q))) # y
    res += 0.5*(n*math.ceil(math.log2(q))) # s
    res += 0.5*com_op_size # d1
    # if b challenge is 1
    res += 0.5*((math.log2(B)+1)*2*n*k) # e'
    res += 0.5*com_op_size # d2
    return res*deltaOL

def linear(param):
    lamb,n,q,d,k,sigma,B,deltaOL,deltaM = param
    com_size = 256
    com_op_size = 8*math.ceil(lamb/8.)
    seed_size = 8*math.ceil(lamb/8.)
    # ch_size = math.ceil(math.log2(q))
    res = 0.0
    res += 2*com_size # commitments c1 and c2
    # res += ch_size # alpha challenge
    res += 3*(math.log2(B)+1)*2*n*k*math.ceil(math.log2(q)) # g
    # res += 1 # second challenge bit
    # if b challenge is 0
    res += 0.5*3*seed_size # pi
    res += 0.5*3*(k*n*math.ceil(math.log2(q))) # y
    res += 0.5*3*(n*math.ceil(math.log2(q))) # s
    res += 0.5*com_op_size # d1
    # if b challenge is 1
    res += 0.5*(3*(math.log2(B)+1)*2*n*k) # e'
    res += 0.5*com_op_size # d2
    return res*deltaOL

def multiplicative(param):
    lamb,n,q,d,k,sigma,B,deltaOL,deltaM = param
    com_size = 256
    com_op_size = 8*math.ceil(lamb/8.)
    seed_size = 8*math.ceil(lamb/8.)
    # ch_size = math.ceil(math.log2(q))
    res = 0.0
    res += 4*com_size # commitments c1, c2, c3 and c4
    # res += 2*ch_size # alpha and beta challenges
    res += 3*(math.log2(B)+1)*2*n*k*math.ceil(math.log2(q)) # g
    res += com_size # commitment c5
    # res += 1 # second challenge bit
    # if b challenge is 0
    res += 0.5*3*seed_size # pi
    res += 0.5*3*(k*n*math.ceil(math.log2(q))) # y
    res += 0.5*5*(n*math.ceil(math.log2(q))) # t
    res += 0.5*3*com_op_size # d1, d4 and d5
    # if b challenge is 1
    res += 0.5*(3*(math.log2(B)+1)*2*n*k) # e'
    res += 0.5*3*(n*math.ceil(math.log2(q))) # mu
    res += 0.5*3*com_op_size # d2, d3 and d5
    return res*deltaM
