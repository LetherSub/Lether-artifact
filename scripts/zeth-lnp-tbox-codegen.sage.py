

# This file was *autogenerated* from the file zeth-lnp-tbox-codegen.sage
from sage.all_cmdline import *   # import sage library

_sage_const_512 = Integer(512); _sage_const_8 = Integer(8); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_4 = Integer(4); _sage_const_128 = Integer(128); _sage_const_1p0044 = RealNumber('1.0044'); _sage_const_1p64 = RealNumber('1.64'); _sage_const_14 = Integer(14); _sage_const_5 = Integer(5); _sage_const_64 = Integer(64); _sage_const_6p36 = RealNumber('6.36'); _sage_const_1p6 = RealNumber('1.6'); _sage_const_12193p83 = RealNumber('12193.83'); _sage_const_140 = Integer(140); _sage_const_129 = Integer(129); _sage_const_59 = Integer(59); _sage_const_147 = Integer(147); _sage_const_337 = Integer(337); _sage_const_256 = Integer(256); _sage_const_3 = Integer(3); _sage_const_2p5 = RealNumber('2.5'); _sage_const_2p25 = RealNumber('2.25'); _sage_const_127 = Integer(127); _sage_const_26 = Integer(26); _sage_const_377 = Integer(377); _sage_const_41 = Integer(41); _sage_const_1p55 = RealNumber('1.55'); _sage_const_13 = Integer(13)
import mpmath as mp
from mpmath import mpf, nstr
import sys

# XXX update to new LWE estimator
# sys.path.append('../third_party/lattice-estimator/')
# from estimator import *  # nopep8

# cannot print more than 822
mp.mp.prec = _sage_const_512 
prec = _sage_const_8   # precision for nstr



if not 'verbose' in globals():
    verbose = _sage_const_1 
if not 'code' in globals():
    code = _sage_const_1 
if not 'loaded' in globals():
    loaded = _sage_const_0 


if not loaded:
    assert len(sys.argv) == _sage_const_2 
    params_file = sys.argv[_sage_const_1 ]

load("codegen.sage")


# bound B
def Bound_():
    global stdev2
    global m2
    global d
    global eta
    global D
    global kmsis
    global gamma
    return stdev2 * mp.sqrt(_sage_const_2  * m2 * d) + mpf(eta) * _sage_const_2  ** (D-_sage_const_1 ) * mp.sqrt(kmsis*d) + (gamma * mp.sqrt(kmsis * d))/mpf(_sage_const_2 )


# bound B1 on zbar1
def Bound1():
    global stdev1
    global m1
    global Z
    global d
    # XXX compare B1 page 150 [1]
    return mpf(_sage_const_2 ) * stdev1 * mp.sqrt(_sage_const_2  * (m1 + Z) * d)


# bound B2
def Bound2():
    return mpf(_sage_const_2 ) * Bound_()


# bound B on the extracted MSIS solution
def Bound():
    global eta
    return _sage_const_4  * mpf(eta) * mp.sqrt(Bound1() ** _sage_const_2  + Bound2() ** _sage_const_2 )


# constants, dont change
KAPPA = _sage_const_128     # security param, bit security
DELTA128 = _sage_const_1p0044   # root hermite factor for 128-bit security
# number of irreducible factors of X^d + 1 modulo each q_i,  q_i = 2l+1 (mod 4l)
L = _sage_const_2 
T = _sage_const_1p64   # for KAPPA=128
NADDS = _sage_const_128   # chose P big enough for this many additions


# default values rejection sampling
rejs1 = _sage_const_1 
rejs2 = _sage_const_2 
rejs3 = _sage_const_2 
rejs4 = _sage_const_2 
rejs5 = _sage_const_2 

gamma1 = _sage_const_14 
gamma2 = _sage_const_1 
gamma3 = _sage_const_5 
gamma4 = _sage_const_5 
gamma5 = _sage_const_1 

if not loaded:
    load(params_file)


# default values if no norms are proven
if not 'nbin' in globals():
    nbin = _sage_const_0 
if not 'B' in globals():
    B = []
if not 'n' in globals():
    n = []
if not 'Bprime' in globals():
    Bprime = _sage_const_0 
if not 'nprime' in globals():
    nprime = _sage_const_0 


assert rejs1 == _sage_const_1 
assert rejs2 == _sage_const_1  or rejs2 == _sage_const_2 
assert rejs3 == _sage_const_2 
assert rejs4 == _sage_const_2 
gamma1 = mpf(gamma1)
gamma2 = mpf(gamma2)
gamma3 = mpf(gamma3)
gamma4 = mpf(gamma4)
gamma5 = mpf(gamma5)

# X^d + 1 mod q1,q2 must split into L=2 irreductible factors.
# q1, q2 odd primes, q1, q2 = 2L+1 mod 4L, q1 < q2

alpha = mpf(alpha)
B = [mpf(x) for x in B]

if d not in [_sage_const_64 , _sage_const_128 ]:
    err("d not in [64,128]")
log2d = log(d, _sage_const_2 )

n_div = _sage_const_1   # number of divisors of q
if 'log2q1' in globals():
    if log2q1 >= log2q:
        err("log(q1) > log(q) = log (q1 * q2)")
    n_div = _sage_const_2 
else:
    log2q1_ = log2q

Z = len(B)  # number of exact l2 norm proofs

nex = sum(n) + nbin + Z
alpha3 = mp.sqrt(sum(([x ** _sage_const_2  for x in B])) + (nbin + Z)*d)

if alpha3 == _sage_const_0 :
    alpha3 = _sage_const_1     

approx_proof = _sage_const_1 
alpha4 = mp.sqrt((n_cols + m_rows) * d * (_sage_const_6p36 *_sage_const_1p6 )**_sage_const_2  + n_cols * d)

alpha5 = mp.sqrt(num_rec * d * (_sage_const_12193p83 *_sage_const_1p6 )**_sage_const_2 )
#if nprime > 0 and Bprime > 0:
#    approx_proof = 1
#    Bprime = mpf(Bprime)
#    alpha4 = Bprime
#elif nprime != 0 or Bprime != 0:  # either both > 0, or both == 0
#    err("Invalid approximate proof params")

# number of repetitions for boosting soundness, we assume lambda is even
lmbda = _sage_const_2  * ceil(KAPPA/(_sage_const_2 *log2q1_))

D = _sage_const_0        # dropping low-order bits of t_A
gamma = _sage_const_0    # dropping low-order bits of w

# challenge space
if d == _sage_const_64  and L == _sage_const_2  and log2q1_ >= _sage_const_4 :
    omega = _sage_const_8 
    eta = _sage_const_140 
    Csize = _sage_const_2  ** _sage_const_129 
elif d == _sage_const_128  and L == _sage_const_2  and log2q1_ >= _sage_const_4 :
    omega = _sage_const_2 
    eta = _sage_const_59 
    Csize = _sage_const_2  ** _sage_const_147 
else:
    err("challenge space undefined")

# sample from [-omega,omega] <=> sample from [0,2*omega] - omega
omega_bits = ceil(log(_sage_const_2 *omega+_sage_const_1 , _sage_const_2 ))

# standard deviations for standard rejection sampling
stdev1 = gamma1 * mpf(eta) * mp.sqrt(alpha ** _sage_const_2  + Z*d)
stdev2 = mpf(_sage_const_0 )  # set later (depends on length of randomness s2)

# standard deviations for bimodal rejection sampling
stdev3 = gamma3 * mp.sqrt(_sage_const_337 ) * alpha3
stdev4 = gamma4 * mp.sqrt(_sage_const_337 ) * alpha4
stdev5 = gamma5 * mp.sqrt(_sage_const_337 ) * alpha5

# XXX
stdev1 = round_stdev(stdev1)
gamma1 = stdev1 / (mpf(eta) * mp.sqrt(alpha ** _sage_const_2  + Z*d))

stdev3 = round_stdev(stdev3)
gamma3 = stdev3 / (mp.sqrt(_sage_const_337 ) * alpha3)

stdev4 = round_stdev(stdev4)
gamma4 = stdev4 / (mp.sqrt(_sage_const_337 ) * alpha4)

stdev5 = round_stdev(stdev5)
gamma5 = stdev5 / (mp.sqrt(_sage_const_337 ) * alpha5)
# XXX

if gamma1 <= _sage_const_0 :
    err("gamma1 is negative")
if gamma2 <= _sage_const_0 :
    err("gamma2 is negative")
if gamma3 <= _sage_const_0 :
    err("gamma3 is negative")
if gamma4 <= _sage_const_0 :
    err("gamma4 is negative")
if gamma5 <= _sage_const_0 :
    err("gamma5 is negative")


nu = _sage_const_1       # randomness vector s2 with coefficients between -nu and nu
kmlwe = _sage_const_0            # MLWE dim, to be determined
easy_mlwe_dim = _sage_const_0    # lower bound for MLWE dim
hard_mlwe_dim = _sage_const_64   # guess for upper bound for MLWE dim
# find upper actual bound (and possibly improve lower bound)
while True:
    delta_mlwe = get_delta_mlwe(nu, hard_mlwe_dim, d, _sage_const_2  ** log2q)
    if delta_mlwe <= DELTA128:
        # print(f"MLWE dim {hard_mlwe_dim}: hard")
        break
    # print(f"MLWE dim {hard_mlwe_dim}: easy")
    easy_mlwe_dim = hard_mlwe_dim
    hard_mlwe_dim *= _sage_const_2 
# binary search for smallest MLWE dimension that is still hard
while True:
    kmlwe = (easy_mlwe_dim + hard_mlwe_dim) / _sage_const_2 
    delta_mlwe = get_delta_mlwe(nu, kmlwe, d, _sage_const_2  ** log2q)
    if delta_mlwe <= DELTA128:
        # print(f"MLWE dim {kmlwe} : hard")
        hard_mlwe_dim = kmlwe
    else:
        # print(f"MLWE dim {kmlwe} : easy")
        easy_mlwe_dim = kmlwe
    if hard_mlwe_dim == easy_mlwe_dim + _sage_const_1 :
        kmlwe = hard_mlwe_dim
        # print(f"found MLWE dim : {kmlwe}")
        break


# Find an appropriate Module-SIS dimension
kmsis = _sage_const_0    # dimension of the MSIS problem
while True:
    kmsis += _sage_const_1 
    # we use the packing optimisation from Section 5.3 XXX (paper)
    m2 = kmlwe + kmsis + l + lmbda/_sage_const_2  + _sage_const_256 /d + _sage_const_1  + approx_proof * _sage_const_256 /d + _sage_const_1  + add_ARP * (_sage_const_256 /d + _sage_const_1 )
    stdev2 = gamma2 * mpf(eta) * mpf(nu) * mp.sqrt(m2 * d)
    stdev2 = round_stdev(stdev2)  # XXX
    gamma2 = stdev2 / (mpf(eta) * mpf(nu) * mp.sqrt(m2 * d))  # XXX
    # print(f"d {d}")
    # print(f"2^log2q {2^log2q}")
    # print(f"kmsis {kmsis}")
    # print(f"Bound {Bound()}")
    # print(f"delta {get_delta_msis(Bound(), kmsis, d, 2 ** log2q)}")
    if get_delta_msis(Bound(), kmsis, d, _sage_const_2  ** log2q) < DELTA128 and Bound() < _sage_const_2  ** log2q:
        break

# Find the largest possible gamma which makes the MSIS solution still small.
gamma = _sage_const_2  ** log2q
while True:       # searching for right gamma
    gamma /= _sage_const_2 
    if get_delta_msis(Bound(), kmsis, d, _sage_const_2  ** log2q) < DELTA128 and Bound() < _sage_const_2  ** log2q:
        break

# Find exact values for q, q1 and gamma:
done = false
# q1 = 2L+1 mod 4L
q1 = ceil((_sage_const_2  ** log2q1_)/(_sage_const_4 *L)) * _sage_const_4 *L + _sage_const_2 *L+_sage_const_1 
q1 -= _sage_const_4 *L
while done == False:
    q1 += _sage_const_4 *L
    while is_prime(q1) == False:
        q1 += _sage_const_4 *L
    if n_div == _sage_const_1 :
        q = q1
    elif n_div == _sage_const_2 :
        # q2 = 2L+1 mod 4L
        q2 = ceil((_sage_const_2  ** log2q)/(_sage_const_4 *L*q1)) * _sage_const_4 *L + _sage_const_2 *L + _sage_const_1 
        while is_prime(q2) == False:
            q2 += _sage_const_4 *L
        q = q1 * q2
    else:
        assert n_div == _sage_const_1  or n_div == _sage_const_2 
    Div_q = divisors(q-_sage_const_1 )
    for i in Div_q:
        # find a divisor which is close to gamma
        if gamma*_sage_const_4 /_sage_const_5  < i and i <= gamma and is_even(i):
            gamma = i
            done = True
m = (q-_sage_const_1 ) / gamma

# print("************************ q1 = ", q1)

# Check q,q1,q2
if n_div == _sage_const_2 :
    if q1.divides(q) == False:
        err("q1 is not a divisor of q")
    q2 = q / q1
    if q1 <= _sage_const_3  or is_prime(q1) == False:
        err("q1 is not an odd prime")
    if q2 <= _sage_const_3  or is_prime(q2) == False:
        err("q2 is not an odd prime")
    if not q1 < q2:
        err("q1 is not less than q2")
    if q1 % (_sage_const_4 *L) != _sage_const_2 *L+_sage_const_1 :
        err(f"q1 != {_sage_const_2 *L+_sage_const_1 } mod {_sage_const_4 *L}")
    if q2 % (_sage_const_4 *L) != _sage_const_2 *L+_sage_const_1 :
        err(f"q2 != {_sage_const_2 *L+_sage_const_1 } mod {_sage_const_4 *L}")
elif n_div == _sage_const_1 :
    if q <= _sage_const_3  or is_prime(q) == False:
        err("q is not an odd prime")
    if q % (_sage_const_4 *L) != _sage_const_2 *L+_sage_const_1 :
        err(f"q != {_sage_const_2 *L+_sage_const_1 } mod {_sage_const_4 *L}")
else:
    assert n_div == _sage_const_1  or n_div == _sage_const_2 


# Find the largest possible D which makes the MSIS solution small
D = log2q
while True:
    D -= _sage_const_1 
    if get_delta_msis(Bound(), kmsis, d, q) < DELTA128 and Bound() < _sage_const_2  ** log2q and _sage_const_2  ** (D-_sage_const_1 )*omega*d < gamma:
        break


# update MLWE root hermite factor with exact q
delta_mlwe = get_delta_mlwe(nu, kmlwe, d, q)

M1 = std_gamma2M(gamma1)
if rejs2 == _sage_const_2 :
    M2 = bim_gamma2M(gamma2)  # one-time commitments ([1], chapter 7)
elif rejs2 == _sage_const_1 :
    M2 = std_gamma2M(gamma2)

M3 = bim_gamma2M(gamma3)
M4 = bim_gamma2M(gamma4)
M5 = bim_gamma2M(gamma5)
rate = M1 * M2 * M3 * ((_sage_const_1 -approx_proof) + approx_proof*M4) * M5

# compute proof size
nonshort_bits = kmsis * d *     (ceil(log(q, _sage_const_2 )) - D) + (l + _sage_const_256 /d + _sage_const_1  +
                             approx_proof * _sage_const_256 /d + lmbda + _sage_const_1 ) * d * ceil(log(q, _sage_const_2 ))
challenge_bits = ceil(log(_sage_const_2 *omega+_sage_const_1 , _sage_const_2 )) * d
short_bits1 = (m1 + Z) * d * (ceil(log(stdev1, _sage_const_2 ) + _sage_const_2p5 )) +     (m2 - kmsis) * d * (ceil(log(stdev2, _sage_const_2 ) + _sage_const_2p5 ))
short_bits2 = _sage_const_256  * (ceil(log(stdev3, _sage_const_2 ) + _sage_const_2p5 )) +     approx_proof * _sage_const_256  * (ceil(log(stdev4, _sage_const_2 ) + _sage_const_2p5 )) + add_ARP * _sage_const_256  * (ceil(log(stdev5, _sage_const_2 ) + _sage_const_2p5 ))
hint_bits = _sage_const_2p25  * kmsis * d
proof_bits = nonshort_bits + challenge_bits +     short_bits1 + short_bits2 + hint_bits

printv(f"auto-generated by lnp-tbox.sage.")
printv(f"")

# check completeness conditions ([1], theorem 6.4.1)
if not ((m1 + Z) * d >= _sage_const_5 *KAPPA and m2 * d >= _sage_const_5 *KAPPA):
    err("protocol not complete")
ecorr = _sage_const_1  - _sage_const_1 /(M1*M2*M3*M4) + _sage_const_2  ** (-_sage_const_127 )
printv(
    f"protocol is statistically complete with correctness error >= 1 - 2^({floor(log(_sage_const_1 -ecorr, _sage_const_2 ))})")

# check simulatability conditions ([1], theorem 6.4.2)
if not (kmlwe >= _sage_const_0  and kmlwe == m2 - kmsis - l - lmbda/_sage_const_2  - _sage_const_256 /d - approx_proof * _sage_const_256 /d - _sage_const_2  - add_ARP * (_sage_const_256 /d + _sage_const_1 )):
    err("protocol not simulatable")
printv(
    f"protocol is simulatable under MLWE({kmlwe},{kmsis+l+lmbda/_sage_const_2 +_sage_const_256 /d + approx_proof * _sage_const_256 /d + _sage_const_2 },[-{nu},{nu}])")  # XXX extended MLWE - yes because stdev 2 is bimodal?

# check knowledge-soundness conditions ([1], theorem 6.4.3)
t = mp.sqrt(_sage_const_1  - mp.log(_sage_const_2  ** (-KAPPA)) / mpf(_sage_const_128 ))  # Figure 6.3 [1]
Barp = mpf(_sage_const_2 ) * mp.sqrt(_sage_const_256 /_sage_const_26 ) * t * stdev3
psi = mpf(_sage_const_2 ) * gamma4 * mp.sqrt(_sage_const_377 *_sage_const_2 *KAPPA)  # ~ gamma4 * 621.33
if not q >= _sage_const_41  * nex * d * Barp:
    err("protocol not knowledge-sound: cannot use lemma 3.2.5")
if not q > Barp ** _sage_const_2  + Barp * sqrt(nbin * d):
    err("protocol not knowledge-sound: cannot prove Ps1 * s1 + Pm * m + f has binary coefficients")
if not q > Barp ** _sage_const_2  + Barp * sqrt(d):
    err(
        "protocol not knowledge-sound: cannot prove theta[1], ..., theta[Z] have binary coefficients")
for i in range(len(B)):
    if not q > _sage_const_3  * B[i] ** _sage_const_2  + Barp ** _sage_const_2 :
        err(
            f"protocol not knowledge-sound: cannot prove l2(Es[{i}] * s1 + Em[{i}] * m +  v[{i}]) <= B[{i}]")
eknow = mpf(_sage_const_1 )/mpf(_sage_const_2 *Csize) + mpf(q1) ** (mpf(-d)/mpf(L)) + mpf(
    q1) ** (mpf(-lmbda)) + mpf(_sage_const_2 ) ** (mpf(-_sage_const_128 )) + mpf(_sage_const_2 ) ** (mpf(-_sage_const_256 ))
printv(
    f"protocol is knowledge-sound with knowledge error <= 2^({nstr(mp.ceil(mp.log(eknow,_sage_const_2 )),prec)}) under MSIS({kmsis},{m1+m2},2^({nstr(mp.log(Bound(),_sage_const_2 ), prec)}))")


#################################
# print params
# print comments

printv(f"")
printv(f"Ring")
printv(f"degree d = {d}")
printv(f"modulus q = {q}, log(q) ~ {nstr(mp.log(q,_sage_const_2 ),prec)}")
if n_div == _sage_const_1 :
    printv(f"factors q = q1")
if n_div == _sage_const_2 :
    printv(f"modulus factors q = q1 * q2")
    printv(f"q1 = {q1}, log(q1) ~ {nstr(mp.log(q1,_sage_const_2 ),prec)}")
    printv(f"q2 = {q1}, log(q2) ~ {nstr(mplog(q2,_sage_const_2 ),prec)}")
else:
    assert n_div == _sage_const_1  or n_div == _sage_const_2 
printv(f"")
printv(f"Compression")
printv(f"D = {D}")
printv(f"gamma = {gamma}, log(gamma) ~ {nstr(mp.log(gamma,_sage_const_2 ),prec)}")
printv(f"")
printv(f"Dimensions of secrets")
printv(f"s1: m1 = {m1}")
printv(f"m: l = {l}")
printv(f"s2: m2 = {m2}")
printv(f"")
printv(f"Size of secrets")
printv(f"l2(s1) <= alpha = {nstr(alpha,prec)}")
printv(f"m unbounded")
printv(f"s2 uniform in [-nu,nu] = [{-nu},{nu}]")
printv(f"")
printv(f"Norm proofs")
if nbin == _sage_const_0 :
    printv(f"binary: no")
elif nbin > _sage_const_0 :
    printv(f"binary: yes (dimension: {nbin})")
else:
    assert nbin >= _sage_const_0 
if Z == _sage_const_0 :
    printv(f"exact euclidean: no")
elif Z > _sage_const_0 :
    printv(
        f"exact euclidean: yes (dimensions: {n}, bounds: [{', '.join([nstr(x,prec) for x in B])}])")
else:
    assert Z >= _sage_const_0 
if approx_proof == _sage_const_0 :
    printv(f"approximate infinity: no")
elif approx_proof == _sage_const_1 :
    printv(
        f"approximate infinity: yes (psi: {nstr(psi,prec)}, dimension: {nprime}, bound: {nstr(Bprime,prec)})")
else:
    assert approx_proof == _sage_const_0  or approx_proof == _sage_const_1 
printv(f"")
printv(f"Challenge space")
printv(
    f"c uniform in [-omega,omega] = [{-omega},{omega}], o(c)=c, sqrt(l1(o(c)*c)) <= eta = {eta}")  # XXX square root or 2*k-th root ?
printv(f"")
printv(f"Standard deviations")
printv(
    f"stdev1 = {stdev1}, log(stdev1/1.55) = {mp.log(stdev1/mpf(_sage_const_1p55 ),_sage_const_2 )}")
printv(
    f"stdev2 = {stdev2}, log(stdev2/1.55) = {mp.log(stdev2/mpf(_sage_const_1p55 ),_sage_const_2 )}")
printv(
    f"stdev3 = {stdev3}, log(stdev3/1.55) = {mp.log(stdev3/mpf(_sage_const_1p55 ),_sage_const_2 )}")
if approx_proof == _sage_const_1 :
    printv(
        f"stdev4 = {stdev4}, log(stdev4/1.55) = {mp.log(stdev4/mpf(_sage_const_1p55 ),_sage_const_2 )}")
printv(
    f"stdev5 = {stdev5}, log(stdev4/1.55) = {mp.log(stdev5/mpf(_sage_const_1p55 ),_sage_const_2 )}")
printv(f"")
printv(f"Repetition rate")
printv(f"M1 = {nstr(M1,prec)}")
printv(f"M2 = {nstr(M2,prec)}")
printv(f"M3 = {nstr(M3,prec)}")
if approx_proof == _sage_const_1 :
    printv(f"M4 = {nstr(M4,prec)}")
printv(f"M5 = {nstr(M5,prec)}")
printv(f"total = {nstr(rate, prec)}")
printv(f"")
printv(f"Security")
printv(f"MSIS dimension: {kmsis}")
printv(
    f"MSIS root hermite factor: {nstr(get_delta_msis(Bound(), kmsis, d, q), prec)}")
printv(f"MLWE dimension: {kmlwe}")
printv(f"MLWE root hermite factor: {nstr(mpf(delta_mlwe), prec)}")
printv(f"")
printv(f"Proof size")
printv(f"~ {nstr(proof_bits / (_sage_const_2 **_sage_const_13 ), prec)} KiB")
printv(f"")

log2stdev1 = int(mp.log(stdev1/mpf(_sage_const_1p55 ), _sage_const_2 ))
log2stdev2 = int(mp.log(stdev2/mpf(_sage_const_1p55 ), _sage_const_2 ))
log2stdev3 = int(mp.log(stdev3/mpf(_sage_const_1p55 ), _sage_const_2 ))
log2stdev4 = int(mp.log(stdev4/mpf(_sage_const_1p55 ), _sage_const_2 ))
log2stdev5 = int(mp.log(stdev5/mpf(_sage_const_1p55 ), _sage_const_2 ))

q_nlimbs = int2limbs(q, -_sage_const_1 )[_sage_const_1 ]
if m % _sage_const_2  == _sage_const_0 :
    mby2 = m / _sage_const_2 
else:
    mby2 = _sage_const_0 

minP = min_P(d, q, NADDS)
moduli = moduli_list(nbit, d, minP)[_sage_const_0 ]
P = prod(moduli)
assert P >= minP
nmoduli = len(moduli)
Pmodq = redc(Mod(P, q), q)
Ppmodq = []
Ppmodq_str = []
for i in range(len(moduli)):
    Ppmodq_str += [f"{name}_Ppmodq_{i}"]
    Ppmodq += [redc(Mod(P/moduli[i], q), q)]
Ppmodq_array = strlist2ptrarray(Ppmodq_str)

Bz3sqr = floor((T*stdev3*sqrt(_sage_const_256 )) ** _sage_const_2 )
Bz4 = floor((T*stdev4*sqrt(_sage_const_256 )) ** _sage_const_2 )
Bz5 = floor((T*stdev5*sqrt(_sage_const_256 )) ** _sage_const_2 )

#############################################
# print code


out = ""
out += f"""
#include "lazer.h"
{int_t(f"{name}_q", q)}
{int_t(f"{name}_qminus1", q - _sage_const_1 )}
{int_t(f"{name}_m", (q - _sage_const_1 ) / gamma, q_nlimbs)}
{int_t(f"{name}_mby2", (q - _sage_const_1 ) / gamma / _sage_const_2 , q_nlimbs)}
{int_t(f"{name}_gamma", gamma, q_nlimbs)}
{int_t(f"{name}_gammaby2", gamma / _sage_const_2 , q_nlimbs)}
{int_t(f"{name}_pow2D", _sage_const_2 **D, q_nlimbs)}
{int_t(f"{name}_pow2Dby2", _sage_const_2 **D / _sage_const_2 , q_nlimbs)}
{int_t(f"{name}_Bsq", floor(Bound_()**_sage_const_2 ), _sage_const_2 *q_nlimbs)}
{int_t(f"{name}_scM1", int(mp.nint(mpf(_sage_const_2 **_sage_const_128 ) * M1)))}
{int_t(f"{name}_scM2", int(mp.nint(mpf(_sage_const_2 **_sage_const_128 ) * M2)))}
{int_t(f"{name}_scM3", int(mp.nint(mpf(_sage_const_2 **_sage_const_128 ) * M3)))}
{int_t(f"{name}_scM4", int(mp.nint(mpf(_sage_const_2 **_sage_const_128 ) * M4)))}
{int_t(f"{name}_scM5", int(mp.nint(mpf(_sage_const_2 **_sage_const_128 ) * M5)))}
{int_t(f"{name}_stdev1sq", int(mp.nint(stdev1**_sage_const_2 )), _sage_const_2 *q_nlimbs)}
{int_t(f"{name}_stdev2sq", int(mp.nint(stdev2**_sage_const_2 )), _sage_const_2 *q_nlimbs)}
{int_t(f"{name}_stdev3sq", int(mp.nint(stdev3**_sage_const_2 )), _sage_const_2 *q_nlimbs)}
{int_t(f"{name}_stdev4sq", int(mp.nint(stdev4**_sage_const_2 )), _sage_const_2 *q_nlimbs)}
{int_t(f"{name}_stdev5sq", int(mp.nint(stdev5**_sage_const_2 )), _sage_const_2 *q_nlimbs)}
{int_t(f"{name}_inv2", redc(_sage_const_1 /_sage_const_2  % q, q))}
{int_t(f"{name}_inv4", redc(_sage_const_1 /_sage_const_4  % q, q))}
static const unsigned int {name}_n[{Z}] = {intlist2intarray(n)};
{int_t(f"{name}_Bz3sqr", Bz3sqr, _sage_const_2 *q_nlimbs)}
{int_t(f"{name}_Bz4", Bz4, _sage_const_2 *q_nlimbs)}
{int_t(f"{name}_Bz5", Bz5, _sage_const_2 *q_nlimbs)}
{int_t(f"{name}_Pmodq", Pmodq, q_nlimbs)}
"""
# {int_t(f"{name}_alphasq", int(mp.nint(alpha^2)))}
# {int_t(f"{name}_alpha3sq", int(mp.nint(alpha3^2)))}
# {int_t(f"{name}_Bprimesq", int(mp.nint(Bprime^2)))}
l2Bsq_str = []
for i in range(len(B)):
    l2Bsq_str += [f"{name}_l2Bsq{i}"]
    out += int_t(f"{name}_l2Bsq{i}", int(mp.nint(B[i] ** _sage_const_2 )), q_nlimbs) + "\n"
l2Bsq = strlist2ptrarray(l2Bsq_str)
for i in range(len(moduli)):
    out += int_t(f"{name}_Ppmodq_{i}", Ppmodq[i], q_nlimbs) + f"\n"
out += f"""static const int_srcptr {name}_l2Bsq[] = {l2Bsq};
static const int_srcptr {name}_Ppmodq[] = {Ppmodq_array};
static const polyring_t {name}_ring = {{{{{name}_q, {d}, {ceil(log(q-_sage_const_1 ,_sage_const_2 ))}, {log2d}, moduli_d{d}, {nmoduli}, {name}_Pmodq, {name}_Ppmodq, {name}_inv2}}}};
static const dcompress_params_t {name}_dcomp = {{{{ {name}_q, {name}_qminus1, {name}_m, {name}_mby2, {name}_gamma, {name}_gammaby2, {name}_pow2D, {name}_pow2Dby2, {D}, {m % _sage_const_2 }, {ceil(log(m,_sage_const_2 ))} }}}};
static const abdlop_params_t {name}_tbox = {{{{ {name}_ring, {name}_dcomp, {m1 + Z}, {m2}, {l}, {_sage_const_256 /d * num_norm + _sage_const_1  + add_ARP * (_sage_const_256 /d + _sage_const_1 ) + lmbda/_sage_const_2  + _sage_const_1 }, {kmsis}, {name}_Bsq, {nu}, {omega}, {omega_bits}, {eta}, {rejs1}, {log2stdev1}, {name}_scM1, {name}_stdev1sq, {rejs2}, {log2stdev2}, {name}_scM2, {name}_stdev2sq}}}};
static const abdlop_params_t {name}_quad_eval_ = {{{{ {name}_ring, {name}_dcomp, {m1 + Z}, {m2}, {l + _sage_const_256 /d * num_norm + _sage_const_1  + add_ARP * (_sage_const_256 /d + _sage_const_1 )}, {lmbda/_sage_const_2  + _sage_const_1 }, {kmsis}, {name}_Bsq, {nu}, {omega}, {omega_bits}, {eta}, {rejs1}, {log2stdev1}, {name}_scM1, {name}_stdev1sq, {rejs2}, {log2stdev2}, {name}_scM2, {name}_stdev2sq}}}};
static const abdlop_params_t {name}_quad_many_ = {{{{ {name}_ring, {name}_dcomp, {m1 + Z}, {m2}, {l + _sage_const_256 /d * num_norm + _sage_const_1  + lmbda/_sage_const_2  + add_ARP * (_sage_const_256 /d + _sage_const_1 )}, {_sage_const_1 }, {kmsis}, {name}_Bsq, {nu}, {omega}, {omega_bits}, {eta}, {rejs1}, {log2stdev1}, {name}_scM1, {name}_stdev1sq, {rejs2}, {log2stdev2}, {name}_scM2, {name}_stdev2sq}}}};
static const lnp_quad_eval_params_t {name}_quad_eval = {{{{ {name}_quad_eval_, {name}_quad_many_, {lmbda}}}}};
static const lnp_tbox_params_t {name} = {{{{ {name}_tbox, {name}_quad_eval, {nbin}, {name}_n, {nprime}, {Z}, {nex}, {rejs3}, {log2stdev3}, {name}_scM3, {name}_stdev3sq, {rejs4}, {log2stdev4}, {name}_scM4, {name}_stdev4sq, {name}_Bz3sqr, {name}_Bz4, &{name}_l2Bsq[0], {name}_inv4, {int(ceil(proof_bits / _sage_const_8 ))}UL }}}};
static const add_lnp_tbox_params_t add_tbox_para = {{{{ {add_ARP}, {bit_proof}, {num_rec}, {bit_index}, {bit_length}, {rejs5}, {log2stdev5}, {name}_scM5, {name}_stdev5sq,  {name}_Bz5 }}}};
"""

printc(out)

# sys.exit(int(0))

# [1] Lattice-Based Zero-Knowledge Proofs Under a Few Dozen Kilobytes
# https://doi.org/10.3929/ethz-b-000574844

