import sys

d = 64  
PSI = 3358 


vname = "paramsTP"
name = f"_{vname}"  

deg = 64

length2 = 40

length = 7

n_prime = 2

num_rec = 16

m_rows = 26

n_cols = 26 

eval_t = 61

beta_m = 59
nrows = m_rows + n_cols +num_rec + n_cols + n_prime

ncols = n_cols + num_rec + m_rows + num_rec * 2 + 1 + n_prime + length2 + length + 1
 

# XXX XXX indicate regular proof, 1 for only ARP, 2 for 1 ARP and 1 bin/exact norm proof
num_norm = 2

# indicate additional ARP
add_ARP = 1

# indicate integer proof
bit_proof = 2
bit_index = n_cols + num_rec + m_rows
bit_length = num_rec * 2

mes_index = n_cols
mp_index = n_cols + num_rec + m_rows + num_rec * 2

bin_h_index = n_cols + num_rec + m_rows + num_rec * 2 + 1 + n_prime 
bin_m_index = n_cols + num_rec + m_rows + num_rec * 2 + 1 + n_prime + length2
bin_bal_index = n_cols + num_rec + m_rows + num_rec * 2 + 1 + n_prime + length2 + length
balance_range = 32

nbin = m_rows + num_rec * 2 + 1 + n_prime + length2 + length


wdim = ncols   # dimension of witness vector

# XXX 
wpart = [list(range(0, n_cols + num_rec)), list(range(n_cols + num_rec, wdim))]

wnsub = len(wpart)      

wl2 = [0, 0]
wbin = [0, 1]


wlinf = 20


mod = 70368744178037

k = deg/d
P = ( mod -1)/2
S = wlinf  # bound on linf(s)
E = 0   # bound on linf(e), n/a here


log2q = 46

# find dimension of binary vector
nbin_ = 0
for i in range(wnsub):
    if wbin[i] == 1:
        nbin_ += len(wpart[i])
nbin = k * nbin_       # set length of vector with binary coefficients

# find dimensions of vectors bounded in l2 norm
n_ = []
B_ = []
for i in range(wnsub):
    if wl2[i] > 0:
        n_ += [len(wpart[i])]
        B_ += [wl2[i]]
n = [i * k for i in n_]
B = B_.copy()
n = []
B = []

# find dimension of vector bounded in linf norm
nprime = k * nrows
Bprime = 15.9*5.3     # set linf norm bound


# copy l2 bound list and insert naive l2 norm bounds for binary vectors
wl2_ = wl2.copy()
for i in range(wnsub):
    if wl2_[i] == 0:
        wl2_[i] = sqrt((1 ** 2) * deg * len(wpart[i]))

# XXX
wl2_ = [0, sqrt(deg * nbin)]


# Search for smallest proof size:
# 1. Put the whole witness in the Ajtai part.
# 2. Find the subvector v of largest l2-norm in the Ajtai part and move v to the BDLOP part.
# 3. If the proof size decreased, go back to 2. Otherwise, move v back to the Ajtai part and stop.
ajtai = [i for i in range(wnsub)]  # list of subvectors in Ajtai part
bdlop = []  # list of subvectors in BDLOP part

#***** give below value
m1 = wdim     # set length of vector s1 (committed in Ajtai part)
l = 0         # set length of vector m (committed in BDLOP part)
# set l2-norm bound on vector s1

#***** give this alpha value
alpha = ceil(sqrt( n_cols * d * (6.36*1.6)**2) +  (m_rows + 3 + length + length2 + 1) * d + 2)


##########################
## compute parameters
##########################

verbose = 1
code = 1
loaded = 1
codegen_err = 0
load("zeth-lnp-tbox-codegen.sage")




##########################
## output parameters 
##########################

# indices of elements of w that go to s1
s1_indices = []
Ps_indices = []
Es_indices = []
for i in ajtai:
    if wbin[i] == 1:
        Ps_indices += list(range(len(s1_indices),
                                 len(s1_indices)+len(wpart[i])))
    if wl2[i] > 0:
        Es_indices += [list(range(len(s1_indices),
                                  len(s1_indices)+len(wpart[i])))]
    s1_indices += wpart[i]
# indices of elements of w that go to m
m_indices = []
Em_indices = []
for i in bdlop:
    # cant put binary subvecs in ajtai (Pm not implemented)
    assert wbin[i] == 0
    if wl2[i] > 0:
        Em_indices += [list(range(len(s1_indices),
                                  len(s1_indices)+len(wpart[i])))]
    m_indices += wpart[i]


out = ""

if Ps_indices == []:
    matPs_nrows = 0
    vname_Ps = f"NULL"
else:
    out += f"static const unsigned int {vname}_Ps[{len(Ps_indices)}] = {intlist2intarray(Ps_indices)};\n"
    matPs_nrows = len(Ps_indices)
    vname_Ps = f"{vname}_Ps"

vname_Es = []
for i in range(len(Es_indices)):
    _Es_indices = Es_indices[i]
    if _Es_indices == []:
        vname_Es += [f"NULL"]
    else:
        out += f"static const unsigned int {vname}_Es{i}[{len(_Es_indices)}] = {intlist2intarray(_Es_indices)};\n"
        vname_Es += [f"{vname}_Es{i}"]
if Es_indices == []:
    strEs = "NULL"
    strEs_nrows = "NULL"
else:
    strEs = f"{vname}_Es"
    strEs_nrows = f"{vname}_Es_nrows"
    out += f"static const unsigned int *{vname}_Es[{len(Es_indices)}] = {{ "
    for i in range(len(Es_indices)):
        out += f"{vname_Es[i]}, "
    out += f"}};\n"
    out += f"static const unsigned int {vname}_Es_nrows[{len(Es_indices)}] = {intlist2intarray([len(i) for i in Es_indices])};\n"

vname_Em = []
for i in range(len(Em_indices)):
    _Em_indices = Em_indices[i]
    if _Em_indices == 0:
        vname_Em += [f"NULL"]
    else:
        out += f"static const unsigned int {vname}_Em{i}[{len(_Em_indices)}] = {intlist2intarray(_Em_indices)};\n"
        vname_Em += [f"{vname}_Em{i}"]
if Em_indices == []:
    strEm = "NULL"
    strEm_nrows = "NULL"
else:
    strEm = f"{vname}_Es"
    strEm_nrows = f"{vname}_Es_nrows"
    out += f"static const unsigned int *{vname}_Em[{len(Em_indices)}] = {{ "
    for i in range(len(Em_indices)):
        out += f"{vname_Em[i]}, "
    out += f"}};\n"
    out += f"static const unsigned int {vname}_Em_nrows[{len(Em_indices)}] = {intlist2intarray([len(i) for i in Em_indices])};\n"

if l > 0:
    out += f"static const unsigned int {vname}_m_indices[{len(m_indices)}] = {intlist2intarray(m_indices)};\n"
    vname_m_indices = f"{vname}_m_indices"
else:
    vname_m_indices = f"NULL"

mod = 887*2
# note: here to avoid output error in pinv, we temporarily set mod = 70368744176221 which is different with read q
out += f"""
{int_t(f"{vname}_p", mod)}
{int_t(f"{vname}_pinv", redc(1/mod % q, q))}
static const unsigned int {vname}_s1_indices[{len(s1_indices)}] = {intlist2intarray(s1_indices)};
static const lin_params_t {vname} = {{
{{ {name}, {deg}, {vname}_p, {vname}_pinv, {k}, 
{vname}_s1_indices, {len(s1_indices)},
{vname_m_indices}, {len(m_indices)}, 
{vname_Ps}, {matPs_nrows}, 
{strEs}, {strEs_nrows}, 
{strEm}, {strEm_nrows} }}
}};
"""
printc(out)

##########################
## output parameters finished
##########################
