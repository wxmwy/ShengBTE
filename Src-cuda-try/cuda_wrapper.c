#include <cuda_runtime.h>
void run_cuda_rta_plus_wrapper_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, double* energy, double* velocity, double2* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, double2* Gamma_plus, double* WP3_plus, int* types, double* masses, int* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
#ifdef DEBUG
    printf("==== Begin running CUDA acceleration ====\n");
#endif
    run_cuda_rta_plus_(_rank, _mm,  _nband,  _scalebroad,  _nptk,  _T,  rlattvec,  energy,  velocity,  eigenvect,  _Nlist,  list, _Ntri, Phi, R_j, R_k, Index_i,  Index_j,  Index_k,  IJK, Gamma_plus,  WP3_plus,  types,  masses,  _onlyharmonic, Ngrid,  _nelements,  _natoms);
#ifdef DEBUG
    printf("==== Finish running CUDA acceleration ====\n");
#endif
}

void run_cuda_rta_minus_wrapper_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, double* energy, double* velocity, double2* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, double2* Gamma_minus, double* WP3_minus, int* types, double* masses, int* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
#ifdef DEBUG
    printf("==== Begin running CUDA acceleration ====\n");
#endif
    run_cuda_rta_minus_(_rank, _mm,  _nband,  _scalebroad,  _nptk,  _T,  rlattvec,  energy,  velocity,  eigenvect,  _Nlist,  list, _Ntri, Phi, R_j, R_k, Index_i,  Index_j,  Index_k,  IJK, Gamma_minus,  WP3_minus,  types,  masses,  _onlyharmonic, Ngrid,  _nelements,  _natoms);
#ifdef DEBUG
    printf("==== Finish running CUDA acceleration ====\n");
#endif
}

void init_cuda_rta_wrapper_(int* _rank, int* _nband, int* _nptk, double* rlattvec, double* energy, double* velocity, double2* eigenvect, int* _Nlist, int *list, int* _Ntri, double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* types, double* masses, int* Ngrid, int* _nelements, int _natoms)
{
    init_cuda_rta_(_rank, _nband, _nptk, rlattvec, energy, velocity, eigenvect, _Nlist, list, _Ntri, Phi, R_j, R_k, Index_i, Index_j, Index_k, IJK, types, masses, Ngrid, _nelements, _natoms);
}
void finalize_cuda_rta_wrapper_() { finalize_cuda_rta_(); }
