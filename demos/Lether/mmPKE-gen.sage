import re
import random
#from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianSampler

from random import gauss

from math import floor, ceil

from copy import deepcopy


def to_poly_list(flat_list, num_polys, deg):
    return [Rq(flat_list[i * deg:(i + 1) * deg]) for i in range(num_polys)]

def R2_to_poly_list(flat_list, num_polys, deg):
    return [Rq2(flat_list[i * deg:(i + 1) * deg]) for i in range(num_polys)]

def centered_mod(c, q):
    return ((int(c) + q // 2) % q) - q // 2

def centered_poly(poly, q):
    return [centered_mod(c, q) for c in poly.list()]

def discrete_gaussian_sample(sigma=1.0):
    return int(round(gauss(0, sigma)))

def round_away(x):
    return int(x + 0.5) if x >= 0 else int(x - 0.5)


def int_to_bin_list(x, bit_len=None):
    if x < 0:
        raise ValueError("error")
    bin_str = bin(x)[2:]  
    bin_list = [int(b) for b in bin_str[::-1]]  
    if bit_len is not None:
        while len(bin_list) < bit_len:
            bin_list.append(0)
    return bin_list

def decrypt(c, ct ,sk):
    acc = Rq(0)
    global m
    global d
    global eval_t
    for i in range(m):
        acc += c[i] * sk[i]
    mn = ct - acc
    mn_list = mn.list()
    m_list = []
    de_list = []
    sum_m = 0
    for i in range(d):
        m = (round_away(float(mn_list[i]) / float(q / (2 * eval_t + 2)))) % (2 * eval_t + 2)
        m_list.append(m)
        de = x - m * (2 * eval_t + 2)
        de_list.append(de)
    for i in range(32):
        sum_m += 2**i * (m_list[i] - eval_t)
    print(m_list)
    print(sum_m)


def decompose(elem, bound):
    length = ceil(log(bound, 2))
    #bins = [0] * length
    beta = bound - 2**(length - 1)
    if elem >= beta:
        tmp = elem - beta
        flag = 1
    else: 
        tmp = elem
    bin_list = int_to_bin_list(tmp, length - 1)
    bin_list += [flag]
    return bin_list, beta, length



q = 70368744177829
global d
d = 64 
global m 
global n   
m = 26
n = 26

num_rec = 16

sigma = 1.0
center = 0

global eval_t
eval_t = 61


R.<x> = PolynomialRing(Integers(q))
Rq = R.quotient(x^d + 1, 'xbar')
xbar = Rq.gen()

A_list = [randint(-q//2, q//2) for _ in range(d * m * n)]

B_list = [randint(-q//2, q//2) for _ in range(d * num_rec * n)]

s_list = [discrete_gaussian_sample(sigma=6.36) for _ in range(d * (m + n))]
s2_list = [discrete_gaussian_sample(sigma=6.36) for _ in range(d * (m + n))]

y_list = [discrete_gaussian_sample(sigma=12191.83) for _ in range(d * num_rec)]
y2_list = [discrete_gaussian_sample(sigma=12191.83) for _ in range(d )]

c2_list = [randint(-q//2, q//2) for _ in range(d * m * num_rec)]

ct2_list = [randint(-q//2, q//2) for _ in range(d * num_rec)]

while(1):
    mo_list = [randint(0, 1) for _ in range(30)] + [0]*(d-30)
    mp_list = [randint(0, 1) for _ in range(30)] + [0]*(d-30)
    sum_mo = 0
    sum_mp = 0
    for i in range(32):
        sum_mo += 2**i * mo_list[i]
        sum_mp += 2**i * mp_list[i]
    if sum_mo >= sum_mp:
        mo_list = [x + eval_t for x in mo_list]
        break

print("old message:")
print(mo_list)
print(sum_mo)

print("update message:")
print(mp_list)
print(sum_mp)

index_s = randint(0, num_rec - 1)

while(1):
    index_r = randint(0, num_rec - 1)
    if index_r != index_s:
        break
    if num_rec == 1:
        break

sk_s_list = [randint(0, 1) for _ in range(d * m)]
sk_e_list = [randint(0, 1) for _ in range(d * n)]

IdR_list = [0] * (d * num_rec**2)

set_random_seed(hash("Lether0"))
hat_p = 887
hat_q = hat_p * 2
n_pri = 2

AH_list = [randint(-1 * hat_q/2, hat_q/2-1) for _ in range(d * m * n_pri)]

S.<x> = PolynomialRing(Integers(hat_q))
Rq2 = S.quotient(x^d + 1, 'xbar2')
xbar2 = Rq2.gen()

AH_polys = R2_to_poly_list(AH_list, m * n_pri, d)
AH_poly = [[AH_polys[i * m + j] for j in range(m)] for i in range(n_pri)]
sh_poly = R2_to_poly_list(sk_s_list, m, d)[:m]

cmod = centered_mod(round_away(float (q) / float (2 * eval_t + 2)), q)
for i in range(num_rec):
    for j in range(num_rec):
        if i == j:
            IdR_list[(i * num_rec + j) * d] = cmod

A_polys = to_poly_list(A_list, m * n, d)
A_poly = [[A_polys[i * n + j] for j in range(n)] for i in range(m)]

A_T_poly = [[A_poly[i][j] for i in range(m)] for j in range(n)]


B_polys = to_poly_list(B_list, num_rec * n, d)
B_poly = [[B_polys[i * n + j] for j in range(n)] for i in range(num_rec)]

c2_polys = to_poly_list(c2_list, num_rec * m, d)
c2_poly = [[c2_polys[i * m + j] for j in range(m)] for i in range(num_rec)]
ct2_poly = to_poly_list(ct2_list, num_rec, d)[0:num_rec]


s_poly = to_poly_list(s_list[:d * n], n, d)[:n]
s2_poly = to_poly_list(s2_list[:d * n], n, d)[:n]
e_poly = to_poly_list(s_list[d * n:], m, d)[:m]
e2_poly = to_poly_list(s2_list[d * n:], m, d)[:m]

y_poly = to_poly_list(y_list, num_rec, d)[:num_rec]
y2_poly = Rq(y2_list)


mp = Rq(mp_list)
mo = Rq(mo_list)

m_poly = [ ((i == index_r) - (i == index_s)) * mp  for i in range(num_rec)]

sk_s_poly = to_poly_list(sk_s_list, m, d)[:m]
sk_e_poly = to_poly_list(sk_e_list, n, d)[:n]

t_poly = []
for i in range(m):
    acc = Rq(0)
    for j in range(n):
        acc += A_poly[i][j] * s_poly[j]
    acc += e_poly[i]
    t_poly.append(acc)

t2_poly = []
for i in range(m):
    acc = Rq(0)
    for j in range(n):
        acc += A_poly[i][j] * s2_poly[j]
    acc += e2_poly[i]
    t2_poly.append(acc)


pk_poly = []
for i in range(n):
    acc = Rq(0)
    for j in range(m):
        acc += A_T_poly[i][j] * sk_s_poly[j]
    acc += sk_e_poly[i]
    pk_poly.append(acc)

B_poly[index_s] = pk_poly

ct_poly = []
for i in range(num_rec):
    acc = Rq(0)
    for j in range(n):
        acc += B_poly[i][j] * s_poly[j]
    acc += y_poly[i] + cmod * m_poly[i]
    ct_poly.append(acc)

for i in range(m):
    c2_poly[index_s][i] = t2_poly[i]

ct2 = Rq(0)
i = index_s
for j in range(n):
    ct2 += B_poly[i][j] * s2_poly[j]
ct2 += y2_poly + cmod * mo

ct2_poly[index_s] = ct2


for i in range(num_rec):
    for j in range(m):
        c2_poly[i][j] += t_poly[j]
    ct2_poly[i] += ct_poly[i]


acc = Rq(0)
for i in range(m):
    acc += c2_poly[index_s][i] * sk_s_poly[i]
mn = ct2_poly[index_s] - acc
mn_list = mn.list()
dec_m_list = []
de_list = []
sum_m = 0
for i in range(d):
    x = (round_away(float(mn_list[i]) / float(q / (2 * eval_t + 2)))) % (2 * eval_t + 2)
    dec_m_list.append(x)
    de = mn_list[i] - x * round_away( float (q) / float (2 * eval_t + 2))
    cval = centered_mod(de, q) + round_away( float (q) / float (4 * eval_t + 4))
    de_list.append(cval)
for i in range(32):
    sum_m += 2**i * (dec_m_list[i] - eval_t)
print("dec_m_list = ", dec_m_list)
print("rest message = ", sum_m)

bal_bin_list = int_to_bin_list(sum_m, 64)
print("bal_bin_list", bal_bin_list)

tmp_sum = 0
for i in range(32):
    tmp_sum += 2**i * bal_bin_list[i]
print("tmp_sum = ", tmp_sum)

bin_list, beta, length = decompose(dec_m_list[0], 2*eval_t + 2 - 1)
mbin_list = [0] * (d * length)
for i in range(d):
    bin_list, beta, length = decompose(dec_m_list[i], 2*eval_t + 2 - 1)
    for j in range(length):
        mbin_list[d*j + i] = bin_list[j]

print("beta = ", beta)
print("length = ", length)

delta_m_list = [0] * (d * length)
for i in range(length):
    if i < length - 1:
        delta_m_list[i * d] = 2**i * round_away( float (q) / float (2 * eval_t + 2))
    else:
        delta_m_list[i * d] = beta * round_away( float (q) / float (2 * eval_t + 2))
 
    

bin_list, beta2, length2 = decompose(de_list[0], round_away( float (q) / float (2 * eval_t + 2)))
debin_list = [0] * (d * length2)
for i in range(d):
    bin_list, beta2, length2 = decompose(de_list[i], round_away( float (q) / float (2 * eval_t + 2)))
    for j in range(length2):
        debin_list[d*j + i] = bin_list[j]

print("beta2 = ", beta2)
print("length2 = ", length2)

shift_factor = -1 * round_away( float (q) / float (4 * eval_t + 4))

ct2_poly = [-1 * x for x in ct2_poly]

delta_h_list = [0]*(length2 * d)
for i in range(length2):
    if i < length2 - 1:
        delta_h_list[i*d] = 2**i
    else:
        delta_h_list[i*d] = beta2

mes_factor_list = [0]*d
mes_factor_list[0] = round_away( float (q) / float (2 * eval_t + 2))


shift_factor_list = [shift_factor]*d


bar_C_list = sum([centered_poly(p, q) for row in c2_poly for p in row], [])

minus_bar_ct_list = sum([centered_poly(ct, q) for ct in ct2_poly], [])

B_T_poly = [[B_poly[i][j] for i in range(num_rec)] for j in range(n)]

tag_poly =[]
for i in range(n_pri):
    acc = Rq2(0)
    for j in range(m):
        acc += AH_poly[i][j] * sh_poly[j]
    tag_poly.append(acc)

tag_list = sum([centered_poly(t, hat_q) for t in tag_poly], []) 


tag_list = [floor(float(x) / 2.0) for x in tag_list]

tag_p_poly = R2_to_poly_list(tag_list, n_pri, d)[:n_pri]

r_e_poly = []
for i in range(n_pri):
    acc = tag_poly[i] - 2 * tag_p_poly[i]
    r_e_poly.append(acc)

r_e_list = sum([r_e.list() for r_e in r_e_poly], [])

inv_hat_q = centered_mod(1/hat_q % q, q)

AH_list = [centered_mod( inv_hat_q * x, q) for x in AH_list]

Id_inv_list = [0] * (n_pri * n_pri * d)
for i in range(n_pri):
    for j in range(n_pri):
        if i == j:
            Id_inv_list[(i * n_pri + j) * d] = -1 * inv_hat_q

tag_list = [centered_mod(2 * inv_hat_q * x, q) for x in tag_list]

t_list = sum([centered_poly(t, q) for t in t_poly], [])  
ct_list = sum([centered_poly(ct, q) for ct in ct_poly], [])  
pk_list = sum([centered_poly(pk, q) for pk in pk_poly], [])  

AT_list = sum([centered_poly(p, q) for row in A_T_poly for p in row], [])

B_list = []
B_list = sum([centered_poly(p, q) for row in B_poly for p in row], [])

BT_list = sum([centered_poly(p, q) for row in B_T_poly for p in row], [])

m_list = sum([centered_poly(m, q) for m in m_poly], [])

bin_s_list = [0] * (d * num_rec)
bin_s_list[index_s * d] = 1

bin_r_list = [0] * (d * num_rec)
bin_r_list[index_r * d] = 1


# XXX negative BT_list for further proof
BT_list = [-x for x in BT_list]


s_list_half = s_list[:len(s_list) // 2]
e_list_half = s_list[len(s_list) // 2:]

output_file = "data.h"

def format_c_array(name, data):
    lines = [f"static int64_t {name}[] = {{"]
    for i in range(0, len(data), 8):  
        line = "  " + ", ".join(str(x) for x in data[i:i+8])
        if i + 8 < len(data):
            line += ","
        lines.append(line)
    lines.append("};\n")
    return "\n".join(lines)

with open(output_file, "w") as f:
    f.write("#include <stdint.h>\n\n")
    f.write(format_c_array("A_", A_list))
    f.write(format_c_array("s_", s_list_half))
    f.write(format_c_array("t_", t_list))
    f.write(format_c_array("B_", B_list))
    f.write(format_c_array("nBT_", BT_list))
    f.write(format_c_array("ct_", ct_list))
    f.write(format_c_array("e_", e_list_half))
    f.write(format_c_array("y_", y_list))
    f.write(format_c_array("mes_", m_list))
    f.write(format_c_array("IdR_", IdR_list))
    f.write(format_c_array("AT_", AT_list))
    f.write(format_c_array("pk_", pk_list))
    f.write(format_c_array("sks_", sk_s_list))
    f.write(format_c_array("ske_", sk_e_list))
    f.write(format_c_array("bin_s_", bin_s_list))
    f.write(format_c_array("bin_r_", bin_r_list))
    f.write(format_c_array("mp_", mp_list))
    f.write(format_c_array("AH_", AH_list))
    f.write(format_c_array("IdQ_", Id_inv_list))
    f.write(format_c_array("tag_", tag_list))
    f.write(format_c_array("re_", r_e_list))
    f.write(format_c_array("barC_", bar_C_list))
    f.write(format_c_array("bar_minus_ct_", minus_bar_ct_list))
    f.write(format_c_array("sf_", shift_factor_list))
    f.write(format_c_array("mf_", mes_factor_list))
    f.write(format_c_array("delta_h_", delta_h_list))
    f.write(format_c_array("dec_m_", dec_m_list))
    f.write(format_c_array("bin_h_", debin_list))
    f.write(format_c_array("bin_m_", mbin_list))
    f.write(format_c_array("delta_factor_m_", delta_m_list))
    f.write(format_c_array("bin_bal_", bal_bin_list))
    

print(f" output loaded to {output_file}")

