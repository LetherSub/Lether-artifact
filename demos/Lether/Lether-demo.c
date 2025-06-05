#include "data.h"
#include "lazer.h"
#include "params.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

const int num_rec = 16;
const int num_m = 26;
const int num_n = 26;
const int64_t mod = 70368744177829;
const int n_prime = 2;

const int length2 = 40;
const int length = 7;

polymat_t barC;
polyvec_t bar_ct, delat_h, delat_m;
poly_t sf_elem;

size_t 
prover (uint8_t *proof, polyvec_t s, polymat_t A, polyvec_t t,
        const uint8_t pp[32])
{
  size_t len = 0;
  lin_prover_state_t prover;
  lin_prover_init (prover, pp, paramsTP);
  lin_prover_set_statement (prover, A, t);
  lin_prover_set_witness (prover, s);
  lin_prover_prove (prover, proof, &len, NULL);
  lin_prover_clear (prover);

  return len;
}

int
verifier (const uint8_t *proof, polymat_t A, polyvec_t t, const uint8_t pp[32])
{
  lin_verifier_state_t verifier;
  lin_verifier_init (verifier, pp, paramsTP);
  lin_verifier_set_statement (verifier, A, t);
  int accept = lin_verifier_verify (verifier, proof, NULL);
  lin_verifier_clear (verifier);
  return accept;
}

static double get_sec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return ((double) tv.tv_sec) + 1E-6*((double) tv.tv_usec);
}

int
main (void)
{
  lazer_init ();

  size_t proof_sz = 0;  
  int deg = 64;

  INT_T (p, 1);

  int_set_i64 (p, mod);           
  POLYRING_T (Rp, p, deg);      
  const uint8_t pp[32] = { 0 }; 


  double tim, proof_tim = 0.0, verify_tim = 0.0;

  
  polymat_t A, A1, B1, Id, IdR, M_z, A1T, nB1T, AH, IdQ, barC; 
  polyvec_t s, s_r, s_m, s_sk, s_bs, s_br, s_mp, s_re, s_barm, s_binh, s_binm, s_bin_bal;         
  polyvec_t t, tag;        
  polyvec_t t1;         
  polyvec_t ct1;         
  polyvec_t pk;
  polyvec_t v_zero;         

  polymat_alloc (A, Rp, num_m + num_n + num_rec + num_n + n_prime, num_n + num_rec + num_m + num_rec * 2 + 1 + n_prime + length2 + length + 1);
  polyvec_alloc (s, Rp, num_n + num_rec + num_m + num_rec * 2 + 1 + n_prime + length2 + length + 1);
  polyvec_alloc (t, Rp, num_m + num_n + num_rec + num_n + n_prime);
  
  polymat_alloc (barC, Rp, num_rec, num_m);
  polyvec_alloc (bar_ct, Rp, num_rec);
  polyvec_alloc (delat_h, Rp, length2);
  polyvec_alloc (delat_m, Rp, length);
  poly_alloc(sf_elem, Rp);

  polymat_set_i64(barC, barC_);

  polyvec_set_coeffvec_i64(bar_ct, bar_minus_ct_);
  polyvec_set_coeffvec_i64(delat_h, delta_h_);
  polyvec_set_coeffvec_i64(delat_m, delta_factor_m_);
  poly_set_coeffvec_i64(sf_elem, sf_);

  printf ("load statement As=t, s small ... ");

  polymat_set_zero(A);

  polymat_get_submat (Id, A, 0, 0, num_n, num_n, 1, 1);
  polymat_set_one (Id);
  polymat_get_submat (A1, A, num_n, 0, num_m, num_n, 1, 1);
  polymat_set_i64 (A1, A_);
  polymat_get_submat (A1T, A, num_n + num_m, num_n + num_rec, num_n, num_m, 1, 1);
  polymat_set_i64 (A1T, AT_);
  polymat_get_submat (nB1T, A, num_n + num_m, num_n + num_rec + num_m, num_n, num_rec, 1, 1);
  polymat_set_i64 (nB1T, nBT_);
  polymat_get_submat (AH, A, num_n + num_m + num_n, num_n + num_rec, n_prime, num_m, 1, 1);
  polymat_set_i64 (AH, AH_);
  polymat_get_submat (IdQ, A, num_n + num_m + num_n, num_n + num_rec + num_m + num_rec * 2 + 1, n_prime, n_prime, 1, 1);
  polymat_set_i64 (IdQ, IdQ_);
  polymat_get_submat (B1, A, num_n + num_m + num_n + n_prime, 0, num_rec, num_n, 1, 1);
  polymat_set_i64 (B1, B_);
  polymat_get_submat (IdR, A, num_n + num_m + num_n + n_prime, num_n, num_rec, num_rec, 1, 1);
  polymat_set_i64 (IdR, IdR_);

  polyvec_get_subvec(s_r, s, 0, num_n, 1);
  polyvec_set_coeffvec_i64 (s_r, s_);
  polyvec_get_subvec(s_m, s, num_n, num_rec, 1);
  polyvec_set_coeffvec_i64 (s_m, mes_);
  polyvec_get_subvec(s_sk, s, num_n + num_rec, num_m, 1);
  polyvec_set_coeffvec_i64 (s_sk, sks_);
  polyvec_get_subvec(s_bs, s, num_n + num_rec + num_m, num_rec, 1);
  polyvec_set_coeffvec_i64 (s_bs, bin_s_);
  polyvec_get_subvec(s_br, s, num_n + num_rec + num_m + num_rec, num_rec, 1);
  polyvec_set_coeffvec_i64 (s_br, bin_r_);
  polyvec_get_subvec(s_mp, s, num_n + num_rec + num_m + num_rec * 2, 1, 1);
  polyvec_set_coeffvec_i64 (s_mp, mp_);
  polyvec_get_subvec(s_re, s, num_n + num_rec + num_m + num_rec * 2 + 1, n_prime, 1);
  polyvec_set_coeffvec_i64 (s_re, re_);
  polyvec_get_subvec(s_binh, s, num_n + num_rec + num_m + num_rec * 2 + 1 + n_prime, length2, 1);
  polyvec_set_coeffvec_i64 (s_binh, bin_h_);
  polyvec_get_subvec(s_binm, s, num_n + num_rec + num_m + num_rec * 2 + 1 + n_prime + length2, length, 1);
  polyvec_set_coeffvec_i64 (s_binm, bin_m_);
  polyvec_get_subvec(s_bin_bal, s, num_n + num_rec + num_m + num_rec * 2 + 1 + n_prime + length2 + length, 1, 1);
  polyvec_set_coeffvec_i64 (s_bin_bal, bin_bal_);

  polyvec_set_zero(t);
  polyvec_get_subvec(t1, t, num_n, num_m, 1);
  polyvec_set_coeffvec_i64 (t1, t_);
  polyvec_get_subvec(tag, t, num_n + num_m + num_n, n_prime, 1);
  polyvec_set_coeffvec_i64 (tag, tag_);
  polyvec_get_subvec(ct1, t, num_n + num_m + num_n + n_prime, num_rec, 1);
  polyvec_set_coeffvec_i64 (ct1, ct_);


  polyvec_neg_self (t);

  printf ("[OK]\n");

  printf ("prover generates PoK of s: As-t=0, s small ... ");
  
  uint8_t proof[500000]; // some upper bound on proof size
  
  tim = get_sec();
  proof_sz = prover (proof, s, A, t, pp);
  tim = get_sec() - tim;
  proof_tim += tim;
  
  printf("[OK] proof: %zu bytes. Time: %.6f sec.\n", proof_sz, tim);
  
  printf ("verifier verifies proof ... ");
  
  tim = get_sec();
  int accept = verifier (proof, A, t, pp);
  tim = get_sec() - tim;
  verify_tim += tim;

  printf ("%s Time: %.6f sec.\n", accept ? "[OK]" : "[FAILED]", tim);
  
  polymat_free (A);
  polyvec_free (s);
  polyvec_free (t);
  return 0;
}
