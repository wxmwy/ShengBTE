#include <cuda_runtime.h>

// Ind
void run_cuda_ind_plus_wrapper_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, int* _N_plus, double* energy, double* velocity, double2* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* Indof2ndPhonon_plus, int* Indof3rdPhonon_plus, double* Gamma_plus, double* WP3_plus, int* types, double* masses, int* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
#ifdef DEBUG
    printf("==== Begin running Ind_plus CUDA acceleration ====\n");
#endif
    run_cuda_ind_plus_(_rank, _mm,  _nband,  _scalebroad,  _nptk,  _T,  rlattvec,  _N_plus, energy,  velocity,  eigenvect,  _Nlist,  list, _Ntri, Phi, R_j, R_k, Index_i,  Index_j,  Index_k,  IJK, Indof2ndPhonon_plus, Indof3rdPhonon_plus, Gamma_plus,  WP3_plus,  types,  masses,  _onlyharmonic, Ngrid,  _nelements,  _natoms);
#ifdef DEBUG
    printf("==== Finish running Ind_plus CUDA acceleration ====\n");
#endif
}

void run_cuda_ind_minus_wrapper_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, int* N_minus, double* energy, double* velocity, double2* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* Indof2ndPhonon_minus, int* Indof3rdPhonon_minus, double* Gamma_minus, double* WP3_minus, int* types, double* masses, int* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
#ifdef DEBUG
    printf("==== Begin running Ind_minus CUDA acceleration ====\n");
#endif
    run_cuda_ind_minus_(_rank, _mm,  _nband,  _scalebroad,  _nptk,  _T,  rlattvec,  N_minus, energy,  velocity,  eigenvect,  _Nlist,  list, _Ntri, Phi, R_j, R_k, Index_i,  Index_j,  Index_k,  IJK, Indof2ndPhonon_minus, Indof3rdPhonon_minus, Gamma_minus,  WP3_minus,  types,  masses,  _onlyharmonic, Ngrid,  _nelements,  _natoms);
#ifdef DEBUG
    printf("==== Finish running Ind_minus CUDA acceleration ====\n");
#endif
}

void init_cuda_ind_wrapper_(int* _rank, int* _nband, int* _nptk, double* rlattvec, int* _N, double* energy, double* velocity, double2* eigenvect, int* _Nlist, int *list, int* _Ntri, double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* types, double* masses, int* Ngrid, int* _nelements, int _natoms)
{
    init_cuda_ind_(_rank, _nband, _nptk, rlattvec, _N, energy, velocity, eigenvect, _Nlist, list, _Ntri, Phi, R_j, R_k, Index_i, Index_j, Index_k, IJK, types, masses, Ngrid, _nelements, _natoms);
}
void finalize_cuda_ind_wrapper_() { finalize_cuda_ind_(); }



//NP

/*void run_cuda_np_plus_wrapper_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, double* energy, double* velocity, int* _Nlist, int* list, int* IJK, int* N_plus, double* P_plus, int* types, double* masses, int* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
#ifdef DEBUG
    printf("==== Begin running NP_plus CUDA acceleration ====\n");
#endif
    run_cuda_np_plus_(_rank, _mm,  _nband,  _scalebroad,  _nptk,  _T,  rlattvec,  energy,  velocity,  _Nlist,  list,  IJK, N_plus,  P_plus,  types,  masses,  _onlyharmonic, Ngrid,  _nelements,  _natoms);
#ifdef DEBUG
    printf("==== Finish running NP_plus CUDA acceleration ====\n");
#endif
}

void run_cuda_np_minus_wrapper_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, double* energy, double* velocity, int* _Nlist, int* list, int* IJK, int* N_minus, double* P_minus, int* types, double* masses, int* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
#ifdef DEBUG
    printf("==== Begin running NP_minus CUDA acceleration ====\n");
#endif
    run_cuda_np_minus_(_rank, _mm,  _nband,  _scalebroad,  _nptk,  _T,  rlattvec,  energy,  velocity,  _Nlist,  list, IJK, N_minus,  P_minus,  types,  masses,  _onlyharmonic, Ngrid,  _nelements,  _natoms);
#ifdef DEBUG
    printf("==== Finish running NP_minus CUDA acceleration ====\n");
#endif
}

void init_cuda_np_wrapper_(int* _rank, int* _nband, int* _nptk, double* rlattvec, double* energy, double* velocity, int* _Nlist, int *list, int* IJK, int* types, double* masses, int* Ngrid, int* _nelements, int _natoms)
{
    init_cuda_np_(_rank, _nband, _nptk, rlattvec, energy, velocity, _Nlist, list, IJK, types, masses, Ngrid, _nelements, _natoms);
}
void finalize_cuda_np_wrapper_() { finalize_cuda_np_(); }*/




// RTA
void run_cuda_rta_plus_wrapper_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, double* energy, double* velocity, double2* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, double* Gamma_plus, double* WP3_plus, int* types, double* masses, int* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
#ifdef DEBUG
    printf("==== Begin running RTA_plus CUDA acceleration ====\n");
#endif
    run_cuda_rta_plus_(_rank, _mm,  _nband,  _scalebroad,  _nptk,  _T,  rlattvec,  energy,  velocity,  eigenvect,  _Nlist,  list, _Ntri, Phi, R_j, R_k, Index_i,  Index_j,  Index_k,  IJK, Gamma_plus,  WP3_plus,  types,  masses,  _onlyharmonic, Ngrid,  _nelements,  _natoms);
#ifdef DEBUG
    printf("==== Finish running RTA_plus CUDA acceleration ====\n");
#endif
}

void run_cuda_rta_minus_wrapper_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, double* energy, double* velocity, double2* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, double* Gamma_minus, double* WP3_minus, int* types, double* masses, int* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
#ifdef DEBUG
    printf("==== Begin running RTA_minus CUDA acceleration ====\n");
#endif
    run_cuda_rta_minus_(_rank, _mm,  _nband,  _scalebroad,  _nptk,  _T,  rlattvec,  energy,  velocity,  eigenvect,  _Nlist,  list, _Ntri, Phi, R_j, R_k, Index_i,  Index_j,  Index_k,  IJK, Gamma_minus,  WP3_minus,  types,  masses,  _onlyharmonic, Ngrid,  _nelements,  _natoms);
#ifdef DEBUG
    printf("==== Finish running RTA_minus CUDA acceleration ====\n");
#endif
}

void init_cuda_rta_wrapper_(int* _rank, int* _nband, int* _nptk, double* rlattvec, double* energy, double* velocity, double2* eigenvect, int* _Nlist, int *list, int* _Ntri, double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* types, double* masses, int* Ngrid, int* _nelements, int _natoms)
{
    init_cuda_rta_(_rank, _nband, _nptk, rlattvec, energy, velocity, eigenvect, _Nlist, list, _Ntri, Phi, R_j, R_k, Index_i, Index_j, Index_k, IJK, types, masses, Ngrid, _nelements, _natoms);
}
void finalize_cuda_rta_wrapper_() { finalize_cuda_rta_(); }
