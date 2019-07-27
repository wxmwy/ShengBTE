#include <math.h>
#include "cuda_settings.h"
#include "cuda_debug.h"
extern "C" {
const double hbar=1.05457172647e-22;
const double hbarp=hbar*1e22;
const double Kb=1.380648813e-23;
const double pi=3.141592653589793238;
//#define DEBUG 1
#define EPS 1e-8
//double2* cmplx_debug;
//#define NO_CUDA
double* d_debug;
double* debug;

struct double_complex { double x,y;};

__device__ inline double cuda_atomic_add(double* address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do
	{
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
	}
	while (assumed != old); // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	return __longlong_as_double(old);
}

#ifndef NO_CUDA
__device__
#endif
// Return the base broadening (without prefactor) for a mode.
double base_sigma(double v0, double v1, double v2, int ngrid0, int ngrid1, int ngrid2, double* rlattvec) {
  double r1 = (rlattvec[0]*v0+rlattvec[1]*v1+rlattvec[2]*v2)/ngrid0;
  r1 = r1*r1;
  double r2 = (rlattvec[3]*v0+rlattvec[4]*v1+rlattvec[5]*v2)/ngrid1;
  r2 = r2*r2;
  double r3 = (rlattvec[6]*v0+rlattvec[7]*v1+rlattvec[8]*v2)/ngrid2;
  r3 = r3*r3;
  return sqrt((r1+r2+r3)/6.);
}

#ifndef NO_CUDA
__device__
#endif
// Compute one of the matrix elements involved in the calculation of Ind_plus.
double_complex Vp_plus(int i,int j,int k,int q,int qprime,int qdprime,
        double realqprime0, double realqprime1, double realqprime2,
        double realqdprime0, double realqdprime1, double realqdprime2,
        double_complex* eigenvect, int Ntri, double *Phi, double* R_j, double* R_k,
        int* Index_i,int* Index_j, int* Index_k,
        double* masses, int* types, int nptk, int Nbands) 
{
  double_complex Vp_plus_res;

  int ll;
  int rr;
  int ss;
  int tt;
  double temp;
  double_complex prefactor;
  double_complex Vp0;

  Vp_plus_res.x=0.;
  Vp_plus_res.y=0.;
  //do ll=1,Ntri
  for(ll=0;ll<Ntri;++ll) {
      double val1 = realqprime0*R_j[3*ll] + realqprime1*R_j[1+3*ll] + realqprime2*R_j[2+3*ll];
      double val2 = realqdprime0*R_k[3*ll] + realqdprime1*R_k[1+3*ll] + realqdprime2*R_k[2+3*ll];
      double tmp1_x = cos(val1);
      double tmp1_y = sin(val1);
      double tmp2_x = cos(-val2);
      double tmp2_y = sin(-val2);
      double tmp3_r = 1./sqrt(masses[types[Index_i[ll]-1]-1]*
                        masses[types[Index_j[ll]-1]-1]*masses[types[Index_k[ll]-1]-1]);
      prefactor.x = tmp3_r*tmp1_x;
      //prefactor.x = realqprime0*R_j[3*ll];
      //prefactor.x = R_j[3*ll];
      //prefactor.x = R_j[3*ll]+R_j[3*ll+1]+R_j[3*ll+2];
      prefactor.y = tmp3_r*tmp1_y;
      temp = prefactor.x;
      prefactor.x = prefactor.x*tmp2_x - prefactor.y*tmp2_y;
      prefactor.y = temp*tmp2_y + prefactor.y*tmp2_x;

     Vp0.x=0.;
     Vp0.y=0.;
     for(rr=0;rr<3;++rr) {
         for(ss=0;ss<3;++ss) {
             for(tt=0;tt<3;++tt) {
                int idx_1 = q+nptk*(i+Nbands*(tt+3*(Index_i[ll]-1)));
                int idx_2 = qprime+nptk*(j+Nbands*(ss+3*(Index_j[ll]-1)));
                int idx_3 = qdprime+nptk*(k+Nbands*(rr+3*(Index_k[ll]-1)));
                //double _tmp_x = Phi[tt+ss*3+rr*9+ll*27];
                //double _tmp_y = 0;
                double _tmp_x = Phi[tt+ss*3+rr*9+ll*27] * eigenvect[idx_1].x;
                double _tmp_y = Phi[tt+ss*3+rr*9+ll*27] * eigenvect[idx_1].y;
                temp = _tmp_x;
				_tmp_x = _tmp_x*eigenvect[idx_2].x - _tmp_y*eigenvect[idx_2].y;
				_tmp_y = _tmp_y*eigenvect[idx_2].x + temp*eigenvect[idx_2].y;
                // conjg
				temp = _tmp_x;
                _tmp_x = _tmp_x*eigenvect[idx_3].x + _tmp_y*eigenvect[idx_3].y;
                _tmp_y = _tmp_y*eigenvect[idx_3].x - temp*eigenvect[idx_3].y;
                Vp0.x += _tmp_x;
                Vp0.y += _tmp_y;
             }
         }
     }
     //Vp_plus_res.x=Vp_plus_res.x+(prefactor.x);
     //Vp_plus_res.y=Vp_plus_res.y+(prefactor.y);
     Vp_plus_res.x=Vp_plus_res.x+(prefactor.x*Vp0.x-prefactor.y*Vp0.y);
     Vp_plus_res.y=Vp_plus_res.y+(prefactor.x*Vp0.y+prefactor.y*Vp0.x);
  }
  return Vp_plus_res;
}

#ifndef NO_CUDA
__device__
#endif

double_complex Vp_minus(int i,int j,int k,int q,int qprime,int qdprime,
        double realqprime0, double realqprime1, double realqprime2,
        double realqdprime0, double realqdprime1, double realqdprime2,
        double_complex* eigenvect, int Ntri, double *Phi, double* R_j, double* R_k,
        int* Index_i,int* Index_j, int* Index_k,
        double* masses, int* types, int nptk, int Nbands) 
{
  double_complex Vp_minus_res;

  int ll;
  int rr;
  int ss;
  int tt;
  double temp;
  double_complex prefactor;
  double_complex Vp0;

  Vp_minus_res.x=0.;
  Vp_minus_res.y=0.;
  //do ll=1,Ntri
  for(ll=0;ll<Ntri;++ll) {
      double val1 = realqprime0*R_j[3*ll] + realqprime1*R_j[1+3*ll] + realqprime2*R_j[2+3*ll];
      double val2 = realqdprime0*R_k[3*ll] + realqdprime1*R_k[1+3*ll] + realqdprime2*R_k[2+3*ll];
      double tmp1_x = cos(-val1);
      double tmp1_y = sin(-val1);
      double tmp2_x = cos(-val2);
      double tmp2_y = sin(-val2);
      double tmp3_r = 1./sqrt(masses[types[Index_i[ll]-1]-1]*
                        masses[types[Index_j[ll]-1]-1]*masses[types[Index_k[ll]-1]-1]);
      prefactor.x = tmp3_r*tmp1_x;
      prefactor.y = tmp3_r*tmp1_y;
      temp = prefactor.x;
      prefactor.x = prefactor.x*tmp2_x - prefactor.y*tmp2_y;
      prefactor.y = temp*tmp2_y + prefactor.y*tmp2_x;

     Vp0.x=0.;
     Vp0.y=0.;
     for(rr=0;rr<3;++rr) {
         for(ss=0;ss<3;++ss) {
             for(tt=0;tt<3;++tt) {
                int idx_1 = q+nptk*(i+Nbands*(tt+3*(Index_i[ll]-1)));
                int idx_2 = qprime+nptk*(j+Nbands*(ss+3*(Index_j[ll]-1)));
                int idx_3 = qdprime+nptk*(k+Nbands*(rr+3*(Index_k[ll]-1)));
                //double _tmp_x = Phi[tt+ss*3+rr*9+ll*27];
                //double _tmp_y = 0;
                double _tmp_x = Phi[tt+ss*3+rr*9+ll*27] * eigenvect[idx_1].x;
                double _tmp_y = Phi[tt+ss*3+rr*9+ll*27] * eigenvect[idx_1].y;
                temp = _tmp_x;
				        _tmp_x = _tmp_x*eigenvect[idx_2].x + _tmp_y*eigenvect[idx_2].y;
				        _tmp_y = _tmp_y*eigenvect[idx_2].x - temp*eigenvect[idx_2].y;
                // conjg
				        temp = _tmp_x;
                _tmp_x = _tmp_x*eigenvect[idx_3].x + _tmp_y*eigenvect[idx_3].y;
                _tmp_y = _tmp_y*eigenvect[idx_3].x - temp*eigenvect[idx_3].y;
                Vp0.x += _tmp_x;
                Vp0.y += _tmp_y;
             }
         }
     }
     Vp_minus_res.x=Vp_minus_res.x+(prefactor.x*Vp0.x-prefactor.y*Vp0.y);
     Vp_minus_res.y=Vp_minus_res.y+(prefactor.x*Vp0.y+prefactor.y*Vp0.x);
  }
  return Vp_minus_res;
}


#ifndef NO_CUDA

/*__global__ void np_plus_kernel(int Nbands, double scalebroad, int nptk, double T, double* rlattvec, double* energy, double* velocity, double_complex* eigenvect, int Nlist, int* list, int Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, double* Gamma_plus, double* WP3_plus, int* types, double* masses, bool onlyharmonic, int Ngrid0,int Ngrid1,int Ngrid2, int q1, int q2, int q3, double omega,int i, int ll, double* d_debug) {
    int j = threadIdx.x+blockIdx.x*blockDim.x;
    int ii = threadIdx.y+blockIdx.y*blockDim.y;
    int k = threadIdx.z+blockIdx.z*blockDim.z;
    if(j==0&&ii==0&&k==0) {
        *Gamma_plus = 0.0;
        *WP3_plus = 0.0;
    }
    __syncthreads();
    if(j>=Nbands || ii>=nptk || k>=Nbands) return;
       
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];
    //realqprime=matmul(rlattvec,qprime/dble(ngrid))
    double realqprime0 = (rlattvec[0]*(qprime0/(double)Ngrid0))+(rlattvec[3]*(qprime1/(double)Ngrid1))+(rlattvec[6]*(qprime2/(double)Ngrid2));
    //double realqprime0 = (rlattvec[0]*(qprime0/__int2double_rn(Ngrid0)))+(rlattvec[3]*(qprime1/__int2double_rn(Ngrid1)))+(rlattvec[6]*(qprime2/__int2double_rn(Ngrid2)));
    double realqprime1 = (rlattvec[1]*(qprime0/(double)Ngrid0))+(rlattvec[4]*(qprime1/(double)Ngrid1))+(rlattvec[7]*(qprime2/(double)Ngrid2));
    double realqprime2 = (rlattvec[2]*(qprime0/(double)Ngrid0))+(rlattvec[5]*(qprime1/(double)Ngrid1))+(rlattvec[8]*(qprime2/(double)Ngrid2));

    //if(j==0 && k==0)    d_debug[ii] = rlattvec[3]*(qprime1/(double)Ngrid1) + rlattvec[6]*(qprime2/(double)Ngrid2);
    //if(j==0 && k==0)  {
    //     d_debug[ii] = 0;
    //     for(int dd=0;dd<Ntri;++dd)   d_debug[ii] += (R_j[3*dd]*realqprime0);
    //}

    double omegap=energy[ii+j*nptk];
    double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
    //--------BEGIN absorption process-----------
    int qdprime0 = (q1+qprime0)%Ngrid0;
    int qdprime1 = (q2+qprime1)%Ngrid1;
    int qdprime2 = (q3+qprime2)%Ngrid2;
    // realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
    double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
    double omegadp = energy[ss+k*nptk];
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
        if(abs(omega+omegap-omegadp)<=(2.*sigma)) {
            double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            double tmp1 = omega+omegap-omegadp;
            double WP3=(fBEprime-fBEdprime)*
                        exp(-(tmp1*tmp1)/(sigma*sigma))/
                        sigma/sqrt(pi)/(omega*omegap*omegadp);
            //cuda_atomic_add(WP3_plus, WP3);
            atomicAdd(WP3_plus, WP3);
            if (!onlyharmonic) { 
                double_complex Vp=Vp_plus(i,j,k,list[ll]-1,ii,ss,
                                   realqprime0,realqprime1,realqprime2,
                                   realqdprime0,realqdprime1,realqdprime2,
                                   eigenvect,Ntri,Phi,R_j,R_k,
                                   Index_i,Index_j,Index_k,
                                   masses, types, nptk, Nbands);
		double absVp = (Vp.x*Vp.x+Vp.y*Vp.y);
		//cuda_atomic_add(Gamma_plus,hbarp*pi/4.*WP3*absVp);
		//double absVp = Vp.x+Vp.y;
		//double absVp = Vp.x;
                //cuda_atomic_add(Gamma_plus,absVp);
                atomicAdd(Gamma_plus,hbarp*pi/double(4)*WP3*absVp);
		//double absVp = sqrt(Vp.x*Vp.x+Vp.y*Vp.y);
                //cuda_atomic_add(Gamma_plus,hbarp*pi/4.*WP3*(absVp*absVp));
            }
        }
    }
}
*/

__global__ void ind_plus_kernel(int Nbands, double scalebroad, int nptk, double T, double* rlattvec, int* N_plus, double* energy, double* velocity, double_complex* eigenvect, int Nlist, int* list, int Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* Indof2ndPhonon_plus, int* Indof3rdPhonon_plus, double* Gamma_plus, double* WP3_plus,  int* types, double* masses, bool onlyharmonic, int Ngrid0,int Ngrid1,int Ngrid2, int q1, int q2, int q3, double omega,int i, int ll, double* d_debug) {
    int j = threadIdx.x+blockIdx.x*blockDim.x;
    int ii = threadIdx.y+blockIdx.y*blockDim.y;
    int k = threadIdx.z+blockIdx.z*blockDim.z;
    
    if(j==0&&ii==0&&k==0) {
        *WP3_plus = 0.0;
    }
    __syncthreads();
    if(j>=Nbands || ii>=nptk || k>=Nbands) return;
    
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];

    double realqprime0 = (rlattvec[0]*(qprime0/(double)Ngrid0))+(rlattvec[3]*(qprime1/(double)Ngrid1))+(rlattvec[6]*(qprime2/(double)Ngrid2));
    double realqprime1 = (rlattvec[1]*(qprime0/(double)Ngrid0))+(rlattvec[4]*(qprime1/(double)Ngrid1))+(rlattvec[7]*(qprime2/(double)Ngrid2));
    double realqprime2 = (rlattvec[2]*(qprime0/(double)Ngrid0))+(rlattvec[5]*(qprime1/(double)Ngrid1))+(rlattvec[8]*(qprime2/(double)Ngrid2));


    double omegap=energy[ii+j*nptk];
    double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
    //--------BEGIN absorption process-----------
    int qdprime0 = (q1+qprime0)%Ngrid0;
    int qdprime1 = (q2+qprime1)%Ngrid1;
    int qdprime2 = (q3+qprime2)%Ngrid2;

    double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
    double omegadp = energy[ss+k*nptk];
    int count = k+Nbands*(ii+j*nptk);
    Indof2ndPhonon_plus[count] = 0;
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
        if(abs(omega+omegap-omegadp)<=(2.*sigma)) {
        	
        	Indof2ndPhonon_plus[count]=ii*Nbands+j+1;
            Indof3rdPhonon_plus[count]=ss*Nbands+k+1;
            

            double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            double_complex Vp=Vp_plus(i,j,k,list[ll]-1,ii,ss,
                                   realqprime0,realqprime1,realqprime2,
                                   realqdprime0,realqdprime1,realqdprime2,
                                   eigenvect,Ntri,Phi,R_j,R_k,
                                   Index_i,Index_j,Index_k,
                                   masses, types, nptk, Nbands);
            double tmp1 = omega+omegap-omegadp;
            double WP3=(fBEprime-fBEdprime)*
                        exp(-(tmp1*tmp1)/(sigma*sigma))/
                        sigma/sqrt(pi)/(omega*omegap*omegadp);
            atomicAdd(WP3_plus, WP3);
            double absVp2 = Vp.x*Vp.x+Vp.y*Vp.y;

            Gamma_plus[count] = hbarp*pi/double(4)*WP3*absVp2*5.60626442*(1.e8)/nptk;

        }
    }
}

__global__ void ind_minus_kernel(int Nbands, double scalebroad, int nptk, double T, double* rlattvec, int* N_minus, double* energy, double* velocity, double_complex* eigenvect, int Nlist, int* list, int Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* Indof2ndPhonon_minus, int* Indof3rdPhonon_minus, double* Gamma_minus, double* WP3_minus, int* types, double* masses, bool onlyharmonic, int Ngrid0,int Ngrid1,int Ngrid2, int q1, int q2, int q3, double omega,int i, int ll, double* d_debug) {
    int j = threadIdx.x+blockIdx.x*blockDim.x;
    int ii = threadIdx.y+blockIdx.y*blockDim.y;
    int k = threadIdx.z+blockIdx.z*blockDim.z;
    if(j==0&&ii==0&&k==0) {
        
        *WP3_minus = 0.0;
    }
    __syncthreads();
    if(j>=Nbands || ii>=nptk || k>=Nbands) return;
    
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];

    double realqprime0 = (rlattvec[0]*(qprime0/(double)Ngrid0))+(rlattvec[3]*(qprime1/(double)Ngrid1))+(rlattvec[6]*(qprime2/(double)Ngrid2));
    double realqprime1 = (rlattvec[1]*(qprime0/(double)Ngrid0))+(rlattvec[4]*(qprime1/(double)Ngrid1))+(rlattvec[7]*(qprime2/(double)Ngrid2));
    double realqprime2 = (rlattvec[2]*(qprime0/(double)Ngrid0))+(rlattvec[5]*(qprime1/(double)Ngrid1))+(rlattvec[8]*(qprime2/(double)Ngrid2));


    double omegap=energy[ii+j*nptk];
    double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
    //--------BEGIN absorption process-----------
    int qdprime0 = (q1-qprime0+Ngrid0)%Ngrid0;
    int qdprime1 = (q2-qprime1+Ngrid1)%Ngrid1;
    int qdprime2 = (q3-qprime2+Ngrid2)%Ngrid2;

    double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
    double omegadp = energy[ss+k*nptk];
    int count = k+Nbands*(ii+j*nptk);
    Indof2ndPhonon_minus[count] = 0;
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
        if(abs(omega-omegap-omegadp)<=(2.*sigma)) {
        	Indof2ndPhonon_minus[count]=ii*Nbands+j+1;
            Indof3rdPhonon_minus[count]=ss*Nbands+k+1;
            

            double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            double_complex Vp=Vp_minus(i,j,k,list[ll]-1,ii,ss,
                                   realqprime0,realqprime1,realqprime2,
                                   realqdprime0,realqdprime1,realqdprime2,
                                   eigenvect,Ntri,Phi,R_j,R_k,
                                   Index_i,Index_j,Index_k,
                                   masses, types, nptk, Nbands);
            double tmp1 = omega-omegap-omegadp;
            double WP3=(fBEprime+fBEdprime+1)*
                        exp(-(tmp1*tmp1)/(sigma*sigma))/
                        sigma/sqrt(pi)/(omega*omegap*omegadp);
            atomicAdd(WP3_minus, WP3);
            double absVp2 = Vp.x*Vp.x+Vp.y*Vp.y;

            Gamma_minus[count] = hbarp*pi/double(4)*WP3*absVp2*5.60626442*(1.e8)/nptk;

        }
    }
}

__global__ void rta_plus_kernel(int Nbands, double scalebroad, int nptk, double T, double* rlattvec, double* energy, double* velocity, double_complex* eigenvect, int Nlist, int* list, int Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, double* Gamma_plus, double* WP3_plus, int* types, double* masses, bool onlyharmonic, int Ngrid0,int Ngrid1,int Ngrid2, int q1, int q2, int q3, double omega,int i, int ll, double* d_debug) {
    int j = threadIdx.x+blockIdx.x*blockDim.x;
    int ii = threadIdx.y+blockIdx.y*blockDim.y;
    int k = threadIdx.z+blockIdx.z*blockDim.z;
    if(j==0&&ii==0&&k==0) {
        *Gamma_plus = 0.0;
        *WP3_plus = 0.0;
    }
    __syncthreads();
    if(j>=Nbands || ii>=nptk || k>=Nbands) return;
       
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];
    //realqprime=matmul(rlattvec,qprime/dble(ngrid))
    double realqprime0 = (rlattvec[0]*(qprime0/(double)Ngrid0))+(rlattvec[3]*(qprime1/(double)Ngrid1))+(rlattvec[6]*(qprime2/(double)Ngrid2));
    //double realqprime0 = (rlattvec[0]*(qprime0/__int2double_rn(Ngrid0)))+(rlattvec[3]*(qprime1/__int2double_rn(Ngrid1)))+(rlattvec[6]*(qprime2/__int2double_rn(Ngrid2)));
    double realqprime1 = (rlattvec[1]*(qprime0/(double)Ngrid0))+(rlattvec[4]*(qprime1/(double)Ngrid1))+(rlattvec[7]*(qprime2/(double)Ngrid2));
    double realqprime2 = (rlattvec[2]*(qprime0/(double)Ngrid0))+(rlattvec[5]*(qprime1/(double)Ngrid1))+(rlattvec[8]*(qprime2/(double)Ngrid2));

    //if(j==0 && k==0)    d_debug[ii] = rlattvec[3]*(qprime1/(double)Ngrid1) + rlattvec[6]*(qprime2/(double)Ngrid2);
    //if(j==0 && k==0)  {
    //     d_debug[ii] = 0;
    //     for(int dd=0;dd<Ntri;++dd)   d_debug[ii] += (R_j[3*dd]*realqprime0);
    //}

    double omegap=energy[ii+j*nptk];
    double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
    //--------BEGIN absorption process-----------
    int qdprime0 = (q1+qprime0)%Ngrid0;
    int qdprime1 = (q2+qprime1)%Ngrid1;
    int qdprime2 = (q3+qprime2)%Ngrid2;
    // realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
    double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
    double omegadp = energy[ss+k*nptk];
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
        if(abs(omega+omegap-omegadp)<=(2.*sigma)) {
            double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            double tmp1 = omega+omegap-omegadp;
            double WP3=(fBEprime-fBEdprime)*
                        exp(-(tmp1*tmp1)/(sigma*sigma))/
                        sigma/sqrt(pi)/(omega*omegap*omegadp);
            //cuda_atomic_add(WP3_plus, WP3);
            atomicAdd(WP3_plus, WP3);
            if (!onlyharmonic) { 
                double_complex Vp=Vp_plus(i,j,k,list[ll]-1,ii,ss,
                                   realqprime0,realqprime1,realqprime2,
                                   realqdprime0,realqdprime1,realqdprime2,
                                   eigenvect,Ntri,Phi,R_j,R_k,
                                   Index_i,Index_j,Index_k,
                                   masses, types, nptk, Nbands);
		double absVp = (Vp.x*Vp.x+Vp.y*Vp.y);
		//cuda_atomic_add(Gamma_plus,hbarp*pi/4.*WP3*absVp);
		//double absVp = Vp.x+Vp.y;
		//double absVp = Vp.x;
                //cuda_atomic_add(Gamma_plus,absVp);
                atomicAdd(Gamma_plus,hbarp*pi/double(4)*WP3*absVp);
		//double absVp = sqrt(Vp.x*Vp.x+Vp.y*Vp.y);
                //cuda_atomic_add(Gamma_plus,hbarp*pi/4.*WP3*(absVp*absVp));
            }
        }
    }
}

__global__ void rta_minus_kernel(int Nbands, double scalebroad, int nptk, double T, double* rlattvec, double* energy, double* velocity, double_complex* eigenvect, int Nlist, int* list, int Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, double* Gamma_minus, double* WP3_minus, int* types, double* masses, bool onlyharmonic, int Ngrid0,int Ngrid1,int Ngrid2, int q1, int q2, int q3, double omega,int i, int ll, double* d_debug) {
    int j = threadIdx.x+blockIdx.x*blockDim.x;
    int ii = threadIdx.y+blockIdx.y*blockDim.y;
    int k = threadIdx.z+blockIdx.z*blockDim.z;
    if(j==0&&ii==0&&k==0) {
        *Gamma_minus = 0.0;
        *WP3_minus = 0.0;
    }
    __syncthreads();
    if(j>=Nbands || ii>=nptk || k>=Nbands) return;
       
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];
    //realqprime=matmul(rlattvec,qprime/dble(ngrid))
    double realqprime0 = (rlattvec[0]*(qprime0/(double)Ngrid0))+(rlattvec[3]*(qprime1/(double)Ngrid1))+(rlattvec[6]*(qprime2/(double)Ngrid2));
    //double realqprime0 = (rlattvec[0]*(qprime0/__int2double_rn(Ngrid0)))+(rlattvec[3]*(qprime1/__int2double_rn(Ngrid1)))+(rlattvec[6]*(qprime2/__int2double_rn(Ngrid2)));
    double realqprime1 = (rlattvec[1]*(qprime0/(double)Ngrid0))+(rlattvec[4]*(qprime1/(double)Ngrid1))+(rlattvec[7]*(qprime2/(double)Ngrid2));
    double realqprime2 = (rlattvec[2]*(qprime0/(double)Ngrid0))+(rlattvec[5]*(qprime1/(double)Ngrid1))+(rlattvec[8]*(qprime2/(double)Ngrid2));

    //if(j==0 && k==0)    d_debug[ii] = rlattvec[3]*(qprime1/(double)Ngrid1) + rlattvec[6]*(qprime2/(double)Ngrid2);
    //if(j==0 && k==0)  {
    //     d_debug[ii] = 0;
    //     for(int dd=0;dd<Ntri;++dd)   d_debug[ii] += (R_j[3*dd]*realqprime0);
    //}

    double omegap=energy[ii+j*nptk];
    double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
    //--------BEGIN absorption process-----------
    int qdprime0 = (q1-qprime0+Ngrid0)%Ngrid0;
    int qdprime1 = (q2-qprime1+Ngrid1)%Ngrid1;
    int qdprime2 = (q3-qprime2+Ngrid2)%Ngrid2;
    // realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
    double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
    double omegadp = energy[ss+k*nptk];
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
        if(abs(omega-omegap-omegadp)<=(2.*sigma)) {
            double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            double tmp1 = omega-omegap-omegadp;
            double WP3=(fBEprime+fBEdprime+1)*
                        exp(-(tmp1*tmp1)/(sigma*sigma))/
                        sigma/sqrt(pi)/(omega*omegap*omegadp);
            //cuda_atomic_add(WP3_plus, WP3);
            atomicAdd(WP3_minus, WP3);
            if (!onlyharmonic) { 
                double_complex Vp=Vp_minus(i,j,k,list[ll]-1,ii,ss,
                                   realqprime0,realqprime1,realqprime2,
                                   realqdprime0,realqdprime1,realqdprime2,
                                   eigenvect,Ntri,Phi,R_j,R_k,
                                   Index_i,Index_j,Index_k,
                                   masses, types, nptk, Nbands);
		double absVp = (Vp.x*Vp.x+Vp.y*Vp.y);

        atomicAdd(Gamma_minus,hbarp*pi/double(4)*WP3*absVp);

            }
        }
    }
}
#endif
// global device pointers
double* dGamma;
double* dWP3;
double* dvelocity, *denergy, *dPhi, *dR_j, *dR_k, *dmasses;
double_complex* deigenvect;
double* drlattvec;
int* dtypes;
int *dIndex_i,*dIndex_j,*dIndex_k,*dIJK;
int* dlist;

int* dN;
double* dP;
int *dIndof2ndPhonon, *dIndof3rdPhonon;
// =======================
// now to do
void init_cuda_ind_(int* _rank, int* _nband, int* _nptk, double* rlattvec, int* _N, double* energy, double* velocity, double_complex* eigenvect, int* _Nlist, int *list, int* _Ntri, double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* types, double* masses, int* Ngrid, int* _nelements, int* _natoms)
{
    int nelements = *_nelements;
    int natoms = *_natoms;
    int nptk = *_nptk;
    int rank = *_rank;
    int Nbands = *_nband;
    int NList = *_Nlist;
    int Ntri = *_Ntri;
    int N_plus = *_N;
    HANDLE_ERROR(cudaSetDevice(rank%GPUNUM));
#ifdef DEBUG
        printf("nptk=%d, Nbands=%d, Ntri=%d, nelements=%d, natoms=%d, NList=%d, N_plus=%d\n", nptk, Nbands, Ntri, nelements, natoms, NList, N_plus);
        printf("Assumed GPU Memory Usage : %lf GB\n", (double)(sizeof(double)*(9+nptk*Nbands*4+33*Ntri+nelements+2+nptk*Nbands*Nbands)+sizeof(double_complex)*(nptk*Nbands*Nbands)+sizeof(int)*(natoms+Ntri*3+nptk*3+NList+2+2*nptk*Nbands*Nbands))/(double)1e9);
#endif
        HANDLE_ERROR(cudaMalloc((void**)&drlattvec, sizeof(double)*9));
        HANDLE_ERROR(cudaMalloc((void**)&dN, sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&denergy,sizeof(double)*nptk*Nbands));
        HANDLE_ERROR(cudaMalloc((void**)&dvelocity,sizeof(double)*nptk*Nbands*3));
        HANDLE_ERROR(cudaMalloc((void**)&deigenvect,sizeof(double_complex)*nptk*Nbands*Nbands));
        HANDLE_ERROR(cudaMalloc((void**)&dPhi,sizeof(double)*27*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dR_j,sizeof(double)*3*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dR_k,sizeof(double)*3*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dmasses,sizeof(double)*nelements));
        HANDLE_ERROR(cudaMalloc((void**)&dtypes,sizeof(int)*natoms));
        HANDLE_ERROR(cudaMalloc((void**)&dIndex_i,sizeof(int)*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dIndex_j,sizeof(int)*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dIndex_k,sizeof(int)*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dIJK,sizeof(int)*nptk*3));
        //HANDLE_ERROR(cudaMalloc((void**)&dIndof2ndPhonon,sizeof(int)*N_plus));
        //HANDLE_ERROR(cudaMalloc((void**)&dIndof3rdPhonon,sizeof(int)*N_plus));
        //HANDLE_ERROR(cudaMalloc((void**)&dGamma,sizeof(double)*N_plus));
        HANDLE_ERROR(cudaMalloc((void**)&dIndof2ndPhonon,sizeof(int)*nptk*Nbands*Nbands));
        HANDLE_ERROR(cudaMalloc((void**)&dIndof3rdPhonon,sizeof(int)*nptk*Nbands*Nbands));
        HANDLE_ERROR(cudaMalloc((void**)&dGamma,sizeof(double)*nptk*Nbands*Nbands));
        HANDLE_ERROR(cudaMalloc((void**)&dWP3,sizeof(double)));
        HANDLE_ERROR(cudaMalloc((void**)&dlist,sizeof(int)*NList));

        HANDLE_ERROR(cudaMalloc((void**)&d_debug,sizeof(double)*nptk));
        debug = (double*)malloc(sizeof(double)*nptk);
	//HANDLE_ERROR(cudaMalloc((void**)&cmplx_debug,sizeof(double2)*nptk*Nbands*Nbands*Ntri));
#ifdef DEBUG
        printf("^^ CUDA MALLOC FINISH ^^\n");
#endif
        HANDLE_ERROR(cudaMemcpy(dlist, list, sizeof(int)*NList,cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(drlattvec, rlattvec, sizeof(double)*9, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dN, _N, sizeof(int)*1, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(denergy, energy, sizeof(double)*nptk*Nbands, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dvelocity, velocity, sizeof(double)*nptk*Nbands*3, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(deigenvect, eigenvect, sizeof(double_complex)*nptk*Nbands*Nbands, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dPhi, Phi, sizeof(double)*27*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dR_j, R_j, sizeof(double)*3*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dR_k, R_k, sizeof(double)*3*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dmasses, masses, sizeof(double)*nelements, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dtypes, types, sizeof(int)*natoms, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dIndex_i, Index_i, sizeof(int)*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dIndex_j, Index_j, sizeof(int)*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dIndex_k, Index_k, sizeof(int)*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dIJK, IJK, sizeof(int)*nptk*3, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaStreamSynchronize(0));
#ifdef DEBUG
        printf("^^ CUDA MEMCPY FINISH ^^\n");
#endif
}

/*void init_cuda_np_(int* _rank, int* _nband, int* _nptk, double* rlattvec, double* energy, double* velocity, int* _Nlist, int *list, int* IJK, int* types, double* masses, int* Ngrid, int* _nelements, int* _natoms)
{
    int nelements = *_nelements;
    int natoms = *_natoms;
    int nptk = *_nptk;
    int rank = *_rank;
    int Nbands = *_nband;
    int NList = *_Nlist;
    HANDLE_ERROR(cudaSetDevice(rank%GPUNUM));
#ifdef DEBUG
        printf("nptk=%d, Nbands=%d, nelements=%d, natoms=%d, NList=%d\n", nptk, Nbands, nelements, natoms, NList);
        printf("Assumed GPU Memory Usage : %lf GB\n", (double)(sizeof(double)*(9+nptk*Nbands*4+nelements+2)+sizeof(int)*(natoms+nptk*3+NList+1))/(double)1e9);
#endif
        HANDLE_ERROR(cudaMalloc((void**)&drlattvec, sizeof(double)*9));
        HANDLE_ERROR(cudaMalloc((void**)&denergy,sizeof(double)*nptk*Nbands));
        HANDLE_ERROR(cudaMalloc((void**)&dvelocity,sizeof(double)*nptk*Nbands*3));
        HANDLE_ERROR(cudaMalloc((void**)&dmasses,sizeof(double)*nelements));
        HANDLE_ERROR(cudaMalloc((void**)&dtypes,sizeof(int)*natoms));
        HANDLE_ERROR(cudaMalloc((void**)&dIJK,sizeof(int)*nptk*3));

        HANDLE_ERROR(cudaMalloc((void**)&dN,sizeof(int)));
        HANDLE_ERROR(cudaMalloc((void**)&dP,sizeof(double)));

        HANDLE_ERROR(cudaMalloc((void**)&dlist,sizeof(int)*NList));

        HANDLE_ERROR(cudaMalloc((void**)&d_debug,sizeof(double)*nptk));
        debug = (double*)malloc(sizeof(double)*nptk);
	//HANDLE_ERROR(cudaMalloc((void**)&cmplx_debug,sizeof(double2)*nptk*Nbands*Nbands*Ntri));
#ifdef DEBUG
        printf("^^ CUDA MALLOC FINISH ^^\n");
#endif
        HANDLE_ERROR(cudaMemcpy(dlist, list, sizeof(int)*NList,cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(drlattvec, rlattvec, sizeof(double)*9, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(denergy, energy, sizeof(double)*nptk*Nbands, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dvelocity, velocity, sizeof(double)*nptk*Nbands*3, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dmasses, masses, sizeof(double)*nelements, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dtypes, types, sizeof(int)*natoms, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dIJK, IJK, sizeof(int)*nptk*3, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaStreamSynchronize(0));
#ifdef DEBUG
        printf("^^ CUDA MEMCPY FINISH ^^\n");
#endif
}*/

void init_cuda_rta_(int* _rank, int* _nband, int* _nptk, double* rlattvec, double* energy, double* velocity, double_complex* eigenvect, int* _Nlist, int *list, int* _Ntri, double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* types, double* masses, int* Ngrid, int* _nelements, int* _natoms)
{
    int nelements = *_nelements;
    int natoms = *_natoms;
    int nptk = *_nptk;
    int rank = *_rank;
    int Nbands = *_nband;
    int NList = *_Nlist;
    int Ntri = *_Ntri;
    HANDLE_ERROR(cudaSetDevice(rank%GPUNUM));
#ifdef DEBUG
        printf("nptk=%d, Nbands=%d, Ntri=%d, nelements=%d, natoms=%d, NList=%d\n", nptk, Nbands, Ntri, nelements, natoms, NList);
        printf("Assumed GPU Memory Usage : %lf GB\n", (double)(sizeof(double)*(9+nptk*Nbands*4+33*Ntri+nelements+3)+sizeof(double_complex)*(nptk*Nbands*Nbands)+sizeof(int)*(natoms+Ntri*3+nptk*3+NList))/(double)1e9);
#endif
        HANDLE_ERROR(cudaMalloc((void**)&drlattvec, sizeof(double)*9));
        HANDLE_ERROR(cudaMalloc((void**)&denergy,sizeof(double)*nptk*Nbands));
        HANDLE_ERROR(cudaMalloc((void**)&dvelocity,sizeof(double)*nptk*Nbands*3));
        HANDLE_ERROR(cudaMalloc((void**)&deigenvect,sizeof(double_complex)*nptk*Nbands*Nbands));
        HANDLE_ERROR(cudaMalloc((void**)&dPhi,sizeof(double)*27*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dR_j,sizeof(double)*3*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dR_k,sizeof(double)*3*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dmasses,sizeof(double)*nelements));
        HANDLE_ERROR(cudaMalloc((void**)&dtypes,sizeof(int)*natoms));
        HANDLE_ERROR(cudaMalloc((void**)&dIndex_i,sizeof(int)*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dIndex_j,sizeof(int)*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dIndex_k,sizeof(int)*Ntri));
        HANDLE_ERROR(cudaMalloc((void**)&dIJK,sizeof(int)*nptk*3));
        HANDLE_ERROR(cudaMalloc((void**)&dGamma,sizeof(double)));
        HANDLE_ERROR(cudaMalloc((void**)&dWP3,sizeof(double)));
        HANDLE_ERROR(cudaMalloc((void**)&dlist,sizeof(int)*NList));

        HANDLE_ERROR(cudaMalloc((void**)&d_debug,sizeof(double)*nptk));
        debug = (double*)malloc(sizeof(double)*nptk);
	//HANDLE_ERROR(cudaMalloc((void**)&cmplx_debug,sizeof(double2)*nptk*Nbands*Nbands*Ntri));
#ifdef DEBUG
        printf("^^ CUDA MALLOC FINISH ^^\n");
#endif
        HANDLE_ERROR(cudaMemcpy(dlist, list, sizeof(int)*NList,cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(drlattvec, rlattvec, sizeof(double)*9, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(denergy, energy, sizeof(double)*nptk*Nbands, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dvelocity, velocity, sizeof(double)*nptk*Nbands*3, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(deigenvect, eigenvect, sizeof(double_complex)*nptk*Nbands*Nbands, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dPhi, Phi, sizeof(double)*27*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dR_j, R_j, sizeof(double)*3*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dR_k, R_k, sizeof(double)*3*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dmasses, masses, sizeof(double)*nelements, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dtypes, types, sizeof(int)*natoms, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dIndex_i, Index_i, sizeof(int)*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dIndex_j, Index_j, sizeof(int)*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dIndex_k, Index_k, sizeof(int)*Ntri, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(dIJK, IJK, sizeof(int)*nptk*3, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaStreamSynchronize(0));
#ifdef DEBUG
        printf("^^ CUDA MEMCPY FINISH ^^\n");
#endif
}

void run_cuda_ind_plus_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, int* _N, double* energy, double* velocity, double_complex* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* Indof2ndPhonon_plus, int* Indof3rdPhonon_plus, double* Gamma_plus, double* WP3_plus, int* types, double* masses, bool* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
    int nelements = *_nelements;
    int natoms = *_natoms;
    bool onlyharmonic = *_onlyharmonic;
    
    int nptk = *_nptk;
    double scalebroad = *_scalebroad;
    int rank = *_rank;
    int mm = *_mm;
    int Nbands = *_nband;
    int NList = *_Nlist;
    int Ntri = *_Ntri;
    double T = *_T;
    int N_plus = *_N;
    int N_plus_count=0;
    *WP3_plus=0.;

    int i=(mm-1)%Nbands;
    int ll=int((mm-1)/Nbands);
    int* q = &IJK[(list[ll]-1)*3];
    double omega=energy[list[ll]-1+i*nptk];
    if(omega!=0) {
#ifdef NO_CUDA
    for(int j=0;j<Nbands;++j) for(int ii=0;ii<nptk;++ii) {
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];
    int Ngrid0 = Ngrid[0];
    int Ngrid1 = Ngrid[1];
    int Ngrid2 = Ngrid[2];
    int q1 = q[0];
    int q2 = q[1];
    int q3 = q[2];
    //realqprime=matmul(rlattvec,qprime/dble(ngrid))
    double realqprime0 = rlattvec[0]*(qprime0/(double)Ngrid0)+rlattvec[3]*(qprime1/(double)Ngrid1)+rlattvec[6]*(qprime2/(double)Ngrid2);
    double realqprime1 = rlattvec[1]*(qprime0/(double)Ngrid0)+rlattvec[4]*(qprime1/(double)Ngrid1)+rlattvec[7]*(qprime2/(double)Ngrid2);
    double realqprime2 = rlattvec[2]*(qprime0/(double)Ngrid0)+rlattvec[5]*(qprime1/(double)Ngrid1)+rlattvec[8]*(qprime2/(double)Ngrid2);
    double omegap=energy[ii+j*nptk];
    double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
    //--------BEGIN absorption process-----------
    for(int k=0;k<Nbands;++k) {
    int qdprime0 = (q1+qprime0)%Ngrid0;
    int qdprime1 = (q2+qprime1)%Ngrid1;
    int qdprime2 = (q3+qprime2)%Ngrid2;
    // realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
    double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
    double omegadp = energy[ss+k*nptk];
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
        if(abs(omega+omegap-omegadp)<=(2.*sigma)) {

          
            Indof2ndPhonon_plus[N_plus_count]=ii*Nbands+j+1;
            Indof3rdPhonon_plus[N_plus_count]=ss*Nbands+k+1;           
//printf( "wy-test-Indof %d %d %d\n", N_plus_count,Indof2ndPhonon_plus[N_plus_count], Indof3rdPhonon_plus[N_plus_count]);

            double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            double_complex Vp=Vp_plus(i,j,k,list[ll]-1,ii,ss,
                                   realqprime0,realqprime1,realqprime2,
                                   realqdprime0,realqdprime1,realqdprime2,
                                   eigenvect,Ntri,Phi,R_j,R_k,
                                   Index_i,Index_j,Index_k,
                                   masses, types, nptk, Nbands);
            double tmp1 = omega+omegap-omegadp;
            double WP3=(fBEprime-fBEdprime)*
                        exp(-(tmp1*tmp1)/(sigma*sigma))/
                        sigma/sqrt(pi)/(omega*omegap*omegadp);
            *WP3_plus += WP3;
            double absVp2 = Vp.x*Vp.x+Vp.y*Vp.y;
            Gamma_plus[N_plus_count] = hbarp*pi/double(4)*WP3*absVp2*5.60626442*(1.e8)/nptk;
	    	N_plus_count+=1;
            
        }
    }
    }}
#else
        //HANDLE_ERROR(cudaMemset(dGamma,0,sizeof(double)));
        //HANDLE_ERROR(cudaMemset(dWP3,0,sizeof(double)));
#ifdef DEBUG
        printf("^^ CUDA MEMSET FINISH ^^\n");
#endif
       dim3 blk( (Nbands+4-1)/4, (nptk+32-1)/32, (Nbands+4-1)/4 );
       dim3 blksize( 4,32,4 );
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL ^^\n");
#endif
       ind_plus_kernel<<<blk,blksize>>>(Nbands, scalebroad, nptk, T, drlattvec, dN, denergy, dvelocity, deigenvect, NList, dlist, Ntri, dPhi, dR_j, dR_k,  dIndex_i, dIndex_j, dIndex_k, dIJK, dIndof2ndPhonon, dIndof3rdPhonon, dGamma, dWP3, dtypes, dmasses, onlyharmonic, Ngrid[0],Ngrid[1], Ngrid[2],q[0],q[1],q[2], omega, i, ll, d_debug);
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL SYNC ^^\n");
#endif
       HANDLE_ERROR(cudaGetLastError());
       HANDLE_ERROR(cudaStreamSynchronize(0));
       HANDLE_ERROR(cudaGetLastError());

#ifdef DEBUG
       printf("^^ finish CUDA KERNEL ^^\n");
       printf("^^ begin CUDA MEMCPY OUT wp3_plus ^^\n");
#endif
       HANDLE_ERROR(cudaMemcpy(WP3_plus,dWP3,sizeof(double),cudaMemcpyDeviceToHost));
#ifdef DEBUG
       printf("^^ begin CUDA MEMCPY OUT gamma_plus ^^\n");
#endif
       int temp_Indof2ndPhonon_plus[Nbands*nptk*Nbands];
       int temp_Indof3rdPhonon_plus[Nbands*nptk*Nbands];
       double temp_Gamma_plus[Nbands*nptk*Nbands];

       HANDLE_ERROR(cudaMemcpy(temp_Gamma_plus,dGamma,sizeof(double)*nptk*Nbands*Nbands,cudaMemcpyDeviceToHost));

       HANDLE_ERROR(cudaMemcpy(temp_Indof2ndPhonon_plus,dIndof2ndPhonon,sizeof(int)*nptk*Nbands*Nbands,cudaMemcpyDeviceToHost));
       HANDLE_ERROR(cudaMemcpy(temp_Indof3rdPhonon_plus,dIndof3rdPhonon,sizeof(int)*nptk*Nbands*Nbands,cudaMemcpyDeviceToHost));
       HANDLE_ERROR(cudaStreamSynchronize(0));
       int count_sum = 0;
       int count0 = Nbands*nptk*Nbands;
       for(int count1=0; count1<count0; count1++){
           	       if(temp_Indof2ndPhonon_plus[count1] != 0){
           	            Indof2ndPhonon_plus[count_sum] = temp_Indof2ndPhonon_plus[count1];
           	            Indof3rdPhonon_plus[count_sum] = temp_Indof3rdPhonon_plus[count1];
           	            Gamma_plus[count_sum] = temp_Gamma_plus[count1];
           	            count_sum++;
           	        }
       }

       //HANDLE_ERROR(cudaMemcpy(debug,cmplx_debug,sizeof(double2)*Ntri*Nbands*Nbands*nptk,cudaMemcpyDeviceToHost));
#endif
#ifdef DEBUG
       printf("^^ scale ^^\n");
#endif
       *WP3_plus=*WP3_plus/nptk;
#ifdef DEBUG
       printf("^^ CUDA ALL FINISH ^^\n");
#endif
    }
    //(*Gamma_plus)=(*Gamma_plus)*5.60626442*double(1e8)/nptk; // THz
#ifdef DEBUG
       printf("^^ ALL FINISH ^^\n");
       //printf("+++ Gamma_plus = %lf, WP3_plus = %lf\n", (*Gamma_plus), (*WP3_plus));
#endif
}


void run_cuda_ind_minus_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, int* _N, double* energy, double* velocity, double_complex* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, int* Indof2ndPhonon_minus, int* Indof3rdPhonon_minus, double* Gamma_minus, double* WP3_minus, int* types, double* masses, bool* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
    int nelements = *_nelements;
    int natoms = *_natoms;
    bool onlyharmonic = *_onlyharmonic;
    
    int nptk = *_nptk;
    double scalebroad = *_scalebroad;
    int rank = *_rank;
    int mm = *_mm;
    int Nbands = *_nband;
    int NList = *_Nlist;
    int Ntri = *_Ntri;
    double T = *_T;
    int N_minus = *_N;
    int N_minus_count=0;
    *WP3_minus=0.;

    int i=(mm-1)%Nbands;
    int ll=int((mm-1)/Nbands);
    int* q = &IJK[(list[ll]-1)*3];
    double omega=energy[list[ll]-1+i*nptk];
    if(omega!=0) {
#ifdef NO_CUDA
    for(int j=0;j<Nbands;++j) for(int ii=0;ii<nptk;++ii) {
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];
    int Ngrid0 = Ngrid[0];
    int Ngrid1 = Ngrid[1];
    int Ngrid2 = Ngrid[2];
    int q1 = q[0];
    int q2 = q[1];
    int q3 = q[2];
    //realqprime=matmul(rlattvec,qprime/dble(ngrid))
    double realqprime0 = rlattvec[0]*(qprime0/(double)Ngrid0)+rlattvec[3]*(qprime1/(double)Ngrid1)+rlattvec[6]*(qprime2/(double)Ngrid2);
    double realqprime1 = rlattvec[1]*(qprime0/(double)Ngrid0)+rlattvec[4]*(qprime1/(double)Ngrid1)+rlattvec[7]*(qprime2/(double)Ngrid2);
    double realqprime2 = rlattvec[2]*(qprime0/(double)Ngrid0)+rlattvec[5]*(qprime1/(double)Ngrid1)+rlattvec[8]*(qprime2/(double)Ngrid2);
    double omegap=energy[ii+j*nptk];
    double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
    //--------BEGIN absorption process-----------
    for(int k=0;k<Nbands;++k) {
    int qdprime0 = (q1-qprime0+Ngrid0)%Ngrid0;
    int qdprime1 = (q2-qprime1+Ngrid1)%Ngrid1;
    int qdprime2 = (q3-qprime2+Ngrid2)%Ngrid2;
    // realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
    double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
//printf("wy-test_ss %d\n",ss);
    double omegadp = energy[ss+k*nptk];
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
        if(abs(omega-omegap-omegadp)<=(2.*sigma)) {

          
            Indof2ndPhonon_minus[N_minus_count]=ii*Nbands+j+1;
            Indof3rdPhonon_minus[N_minus_count]=ss*Nbands+k+1;
            

            double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            double_complex Vp=Vp_minus(i,j,k,list[ll]-1,ii,ss,
                                   realqprime0,realqprime1,realqprime2,
                                   realqdprime0,realqdprime1,realqdprime2,
                                   eigenvect,Ntri,Phi,R_j,R_k,
                                   Index_i,Index_j,Index_k,
                                   masses, types, nptk, Nbands);
            double tmp1 = omega-omegap-omegadp;
            double WP3=(fBEprime+fBEdprime+1)*
                        exp(-(tmp1*tmp1)/(sigma*sigma))/
                        sigma/sqrt(pi)/(omega*omegap*omegadp);
            *WP3_minus += WP3;
            
            double absVp2 = Vp.x*Vp.x+Vp.y*Vp.y;

            Gamma_minus[N_minus_count] = hbarp*pi/double(4)*WP3*absVp2*5.60626442*(1.e8)/nptk;
            N_minus_count+=1;
            
        }
    }
    }}
#else
        //HANDLE_ERROR(cudaMemset(dGamma,0,sizeof(double)));
        //HANDLE_ERROR(cudaMemset(dWP3,0,sizeof(double)));
#ifdef DEBUG
        printf("^^ CUDA MEMSET FINISH ^^\n");
#endif
       dim3 blk( (Nbands+4-1)/4, (nptk+32-1)/32, (Nbands+4-1)/4 );
       dim3 blksize( 4,32,4 );
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL ^^\n");
#endif
       ind_minus_kernel<<<blk,blksize>>>(Nbands, scalebroad, nptk, T, drlattvec, dN, denergy, dvelocity, deigenvect, NList, dlist, Ntri, dPhi, dR_j, dR_k,  dIndex_i, dIndex_j, dIndex_k,  dIJK, dIndof2ndPhonon, dIndof3rdPhonon, dGamma, dWP3, dtypes, dmasses, onlyharmonic, Ngrid[0],Ngrid[1], Ngrid[2],q[0],q[1],q[2], omega, i, ll, d_debug);
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL SYNC ^^\n");
#endif
       HANDLE_ERROR(cudaGetLastError());
       HANDLE_ERROR(cudaStreamSynchronize(0));
       HANDLE_ERROR(cudaGetLastError());

#ifdef DEBUG
       printf("^^ finish CUDA KERNEL ^^\n");
       printf("^^ begin CUDA MEMCPY OUT wp3_plus ^^\n");
#endif
       HANDLE_ERROR(cudaMemcpy(WP3_minus,dWP3,sizeof(double),cudaMemcpyDeviceToHost));
#ifdef DEBUG
       printf("^^ begin CUDA MEMCPY OUT gamma_plus ^^\n");
#endif

       int temp_Indof2ndPhonon_minus[Nbands*nptk*Nbands];
       int temp_Indof3rdPhonon_minus[Nbands*nptk*Nbands];
       double temp_Gamma_minus[Nbands*nptk*Nbands];

       HANDLE_ERROR(cudaMemcpy(temp_Gamma_minus,dGamma,sizeof(double)*nptk*Nbands*Nbands,cudaMemcpyDeviceToHost));

       HANDLE_ERROR(cudaMemcpy(temp_Indof2ndPhonon_minus,dIndof2ndPhonon,sizeof(int)*nptk*Nbands*Nbands,cudaMemcpyDeviceToHost));
       HANDLE_ERROR(cudaMemcpy(temp_Indof3rdPhonon_minus,dIndof3rdPhonon,sizeof(int)*nptk*Nbands*Nbands,cudaMemcpyDeviceToHost));
       HANDLE_ERROR(cudaStreamSynchronize(0));
       int count_sum = 0;
       int count0 = Nbands*nptk*Nbands;
       for(int count1=0; count1<count0; count1++){
           	       if(temp_Indof2ndPhonon_minus[count1] != 0){
           	            Indof2ndPhonon_minus[count_sum] = temp_Indof2ndPhonon_minus[count1];
           	            Indof3rdPhonon_minus[count_sum] = temp_Indof3rdPhonon_minus[count1];
           	            Gamma_minus[count_sum] = temp_Gamma_minus[count1];
           	            count_sum++;
           	        }
       }
       HANDLE_ERROR(cudaStreamSynchronize(0));
       //HANDLE_ERROR(cudaMemcpy(debug,cmplx_debug,sizeof(double2)*Ntri*Nbands*Nbands*nptk,cudaMemcpyDeviceToHost));
#endif
#ifdef DEBUG
       printf("^^ scale ^^\n");
#endif
       *WP3_minus=*WP3_minus*(5e-1)/nptk;
#ifdef DEBUG
       printf("^^ CUDA ALL FINISH ^^\n");
#endif
    }
    //(*Gamma_plus)=(*Gamma_plus)*5.60626442*double(1e8)/nptk; // THz
#ifdef DEBUG
       printf("^^ ALL FINISH ^^\n");
       //printf("+++ Gamma_plus = %lf, WP3_plus = %lf\n", (*Gamma_plus), (*WP3_plus));
#endif
}

/*void run_cuda_np_plus_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, double* energy, double* velocity, int* _Nlist, int* list, int* IJK, int* N_plus, double* P_plus, int* types, double* masses, bool* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
    int nelements = *_nelements;
    int natoms = *_natoms;
    bool onlyharmonic = *_onlyharmonic;
    
    int nptk = *_nptk;
    double scalebroad = *_scalebroad;
    int rank = *_rank;
    int mm = *_mm;
    int Nbands = *_nband;
    int NList = *_Nlist;
    double T = *_T;
    *N_plus=0;
    *P_plus=0.;

    int i=(mm-1)%Nbands;
    int ll=int((mm-1)/Nbands);
    int* q = &IJK[(list[ll]-1)*3];
    double omega=energy[list[ll]-1+i*nptk];
    if(omega!=0) {
#ifdef NO_CUDA
    for(int j=0;j<Nbands;++j) for(int ii=0;ii<nptk;++ii) {
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];
    int Ngrid0 = Ngrid[0];
    int Ngrid1 = Ngrid[1];
    int Ngrid2 = Ngrid[2];
    int q1 = q[0];
    int q2 = q[1];
    int q3 = q[2];
    //realqprime=matmul(rlattvec,qprime/dble(ngrid))
    //double realqprime0 = rlattvec[0]*(qprime0/(double)Ngrid0)+rlattvec[3]*(qprime1/(double)Ngrid1)+rlattvec[6]*(qprime2/(double)Ngrid2);
    //double realqprime1 = rlattvec[1]*(qprime0/(double)Ngrid0)+rlattvec[4]*(qprime1/(double)Ngrid1)+rlattvec[7]*(qprime2/(double)Ngrid2);
    //double realqprime2 = rlattvec[2]*(qprime0/(double)Ngrid0)+rlattvec[5]*(qprime1/(double)Ngrid1)+rlattvec[8]*(qprime2/(double)Ngrid2);
    double omegap=energy[ii+j*nptk];
    //double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
    //--------BEGIN absorption process-----------
    for(int k=0;k<Nbands;++k) {
    int qdprime0 = (q1+qprime0)%Ngrid0;
    int qdprime1 = (q2+qprime1)%Ngrid1;
    int qdprime2 = (q3+qprime2)%Ngrid2;
    // realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
    //double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    //double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    //double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
//printf("wy-test_ss %d\n",ss);
    double omegadp = energy[ss+k*nptk];
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
        if(abs(omega+omegap-omegadp)<=(2.*sigma)) {
            //double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            *N_plus += 1;
            double tmp1 = omega+omegap-omegadp;
            *P_plus += exp(-(tmp1*tmp1)/(sigma*sigma))/
                        (sigma*sqrt(pi)*nptk*nptk*Nbands*Nbands*Nbands);
        }
    }
    }}
#else
        //HANDLE_ERROR(cudaMemset(dGamma,0,sizeof(double)));
        //HANDLE_ERROR(cudaMemset(dWP3,0,sizeof(double)));
#ifdef DEBUG
        printf("^^ CUDA MEMSET FINISH ^^\n");
#endif
       dim3 blk( (Nbands+4-1)/4, (nptk+32-1)/32, (Nbands+4-1)/4 );
       dim3 blksize( 4,32,4 );
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL ^^\n");
#endif
       //wy to do
       rta_plus_kernel<<<blk,blksize>>>(Nbands, scalebroad, nptk, T, drlattvec, denergy, dvelocity, deigenvect, NList, dlist, Ntri, dPhi, dR_j, dR_k,  dIndex_i, dIndex_j, dIndex_k, dIJK, dGamma,dWP3, dtypes, dmasses, onlyharmonic, Ngrid[0],Ngrid[1], Ngrid[2],q[0],q[1],q[2], omega, i, ll, d_debug);
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL SYNC ^^\n");
#endif
       HANDLE_ERROR(cudaGetLastError());
       HANDLE_ERROR(cudaStreamSynchronize(0));
       HANDLE_ERROR(cudaGetLastError());

#ifdef DEBUG
       printf("^^ finish CUDA KERNEL ^^\n");
       printf("^^ begin CUDA MEMCPY OUT wp3_plus ^^\n");
#endif
       //wy to do
       HANDLE_ERROR(cudaMemcpy(WP3_plus,dWP3,sizeof(double),cudaMemcpyDeviceToHost));
#ifdef DEBUG
       printf("^^ begin CUDA MEMCPY OUT gamma_plus ^^\n");
#endif
       //wy to do
       HANDLE_ERROR(cudaMemcpy(Gamma_plus,dGamma,sizeof(double),cudaMemcpyDeviceToHost));
	   HANDLE_ERROR(cudaStreamSynchronize(0));
       //HANDLE_ERROR(cudaMemcpy(debug,cmplx_debug,sizeof(double2)*Ntri*Nbands*Nbands*nptk,cudaMemcpyDeviceToHost));
#endif
#ifdef DEBUG
       printf("^^ scale ^^\n");
#endif
#ifdef DEBUG
       printf("^^ CUDA ALL FINISH ^^\n");
#endif
    }
#ifdef DEBUG
       printf("^^ ALL FINISH ^^\n");
       //printf("+++ Gamma_plus = %lf, WP3_plus = %lf\n", (*Gamma_plus), (*WP3_plus));
#endif
}*/

// RTA-only version of Ind_plus.
void run_cuda_rta_plus_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, double* energy, double* velocity, double_complex* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, double* Gamma_plus, double* WP3_plus, int* types, double* masses, bool* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
    int nelements = *_nelements;
    int natoms = *_natoms;
    bool onlyharmonic = *_onlyharmonic;
    
    int nptk = *_nptk;
    double scalebroad = *_scalebroad;
    int rank = *_rank;
    int mm = *_mm;
    int Nbands = *_nband;
    int NList = *_Nlist;
    int Ntri = *_Ntri;
    double T = *_T;
    *Gamma_plus=0.;
    *WP3_plus=0.;

    int i=(mm-1)%Nbands;
    int ll=int((mm-1)/Nbands);
    int* q = &IJK[(list[ll]-1)*3];
    double omega=energy[list[ll]-1+i*nptk];
    if(omega!=0) {
#ifdef NO_CUDA
    for(int j=0;j<Nbands;++j) for(int ii=0;ii<nptk;++ii) {
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];
    int Ngrid0 = Ngrid[0];
    int Ngrid1 = Ngrid[1];
    int Ngrid2 = Ngrid[2];
    int q1 = q[0];
    int q2 = q[1];
    int q3 = q[2];
    //realqprime=matmul(rlattvec,qprime/dble(ngrid))
    double realqprime0 = rlattvec[0]*(qprime0/(double)Ngrid0)+rlattvec[3]*(qprime1/(double)Ngrid1)+rlattvec[6]*(qprime2/(double)Ngrid2);
    double realqprime1 = rlattvec[1]*(qprime0/(double)Ngrid0)+rlattvec[4]*(qprime1/(double)Ngrid1)+rlattvec[7]*(qprime2/(double)Ngrid2);
    double realqprime2 = rlattvec[2]*(qprime0/(double)Ngrid0)+rlattvec[5]*(qprime1/(double)Ngrid1)+rlattvec[8]*(qprime2/(double)Ngrid2);
    double omegap=energy[ii+j*nptk];
    double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
    //--------BEGIN absorption process-----------
    for(int k=0;k<Nbands;++k) {
    int qdprime0 = (q1+qprime0)%Ngrid0;
    int qdprime1 = (q2+qprime1)%Ngrid1;
    int qdprime2 = (q3+qprime2)%Ngrid2;
    // realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
    double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
//printf("wy-test_ss %d\n",ss);
    double omegadp = energy[ss+k*nptk];
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
        if(abs(omega+omegap-omegadp)<=(2.*sigma)) {
            double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            double tmp1 = omega+omegap-omegadp;
            double WP3=(fBEprime-fBEdprime)*
                        exp(-(tmp1*tmp1)/(sigma*sigma))/
                        sigma/sqrt(pi)/(omega*omegap*omegadp);
	    *WP3_plus += WP3;
            if (!onlyharmonic) {
                double_complex Vp=Vp_plus(i,j,k,list[ll]-1,ii,ss,
                                   realqprime0,realqprime1,realqprime2,
                                   realqdprime0,realqdprime1,realqdprime2,
                                   eigenvect,Ntri,Phi,R_j,R_k,
                                   Index_i,Index_j,Index_k,
                                   masses, types, nptk, Nbands);
		//double absVp = (Vp.x*Vp.x+Vp.y*Vp.y);
		//cuda_atomic_add(Gamma_plus,hbarp*pi/4.*WP3*absVp);
		//double absVp = Vp.x+Vp.y;
		//double absVp = Vp.x;
                //*Gamma_plus += absVp;
		double absVp2 = Vp.x*Vp.x+Vp.y*Vp.y;
		//printf("wy_test_Vp %f\n", absVp2);
		*Gamma_plus += hbarp*pi/4*WP3*absVp2;
	
		
		//double absVp = sqrt(Vp.x*Vp.x+Vp.y*Vp.y);
                //cuda_atomic_add(Gamma_plus,hbarp*pi/4.*WP3*(absVp*absVp));
            }
        }
    }
    }}
#else
        //HANDLE_ERROR(cudaMemset(dGamma,0,sizeof(double)));
        //HANDLE_ERROR(cudaMemset(dWP3,0,sizeof(double)));
#ifdef DEBUG
        printf("^^ CUDA MEMSET FINISH ^^\n");
#endif
       dim3 blk( (Nbands+4-1)/4, (nptk+32-1)/32, (Nbands+4-1)/4 );
       dim3 blksize( 4,32,4 );
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL ^^\n");
#endif
       rta_plus_kernel<<<blk,blksize>>>(Nbands, scalebroad, nptk, T, drlattvec, denergy, dvelocity, deigenvect, NList, dlist, Ntri, dPhi, dR_j, dR_k,  dIndex_i, dIndex_j, dIndex_k, dIJK, dGamma,dWP3, dtypes, dmasses, onlyharmonic, Ngrid[0],Ngrid[1], Ngrid[2],q[0],q[1],q[2], omega, i, ll, d_debug);
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL SYNC ^^\n");
#endif
       HANDLE_ERROR(cudaGetLastError());
       HANDLE_ERROR(cudaStreamSynchronize(0));
       HANDLE_ERROR(cudaGetLastError());
/*       
       HANDLE_ERROR(cudaMemcpy(debug,d_debug,sizeof(double)*nptk,cudaMemcpyDeviceToHost));
       int count = 0; int icount=0;
       double sum1=0.0, sumr=0.0;
       for(int ii=0;ii<nptk;++ii) {
           int qprime0 = IJK[ii*3  ];
           int qprime1 = IJK[ii*3+1];
           int qprime2 = IJK[ii*3+2];
           int Ngrid0 = Ngrid[0];
           int Ngrid1 = Ngrid[1];
           int Ngrid2 = Ngrid[2];
           //double ref = (rlattvec[2]*(qprime0/(double)Ngrid0))+(rlattvec[5]*(qprime1/(double)Ngrid1))+(rlattvec[8]*(qprime2/(double)Ngrid2));
           double rr = ((rlattvec[0]*(qprime0/(double)Ngrid0))+(rlattvec[3]*(qprime1/(double)Ngrid1))+(rlattvec[6]*(qprime2/(double)Ngrid2)));
           double ref=0;
           for(int dd=0;dd<Ntri;++dd) ref+=rr*R_j[dd*3];
           //double ref = rlattvec[3]*(qprime1/(double)Ngrid1)+rlattvec[6]*(qprime2/(double)Ngrid2);
           //double ref = rlattvec[3]*(qprime1/(double)Ngrid1);
           //double ref = (rlattvec[0]*(qprime0/(double)Ngrid0))+(rlattvec[3]*(qprime1/(double)Ngrid1));
           if((debug[ii]!=0 || ref!=0) && icount++<15) printf("%d : gpu %e, ref %e\n",ii, debug[ii], ref);
           if((ref==0.0&&fabs(debug[ii])>1e-8)||(ref!=0.0 && fabs(debug[ii]-ref)/ref > 1e-8)) {
                printf("Error at %d : gpu %e, ref %e, error %e\n", ii, debug[ii], ref, fabs(debug[ii]-ref));
                ++count;
                if(count>15) exit(-1);
           }
           sum1 += debug[ii];
           sumr += ref;
       }
       printf("%e, %e\n",sum1, sumr);*/
#ifdef DEBUG
       printf("^^ finish CUDA KERNEL ^^\n");
       printf("^^ begin CUDA MEMCPY OUT wp3_plus ^^\n");
#endif
       HANDLE_ERROR(cudaMemcpy(WP3_plus,dWP3,sizeof(double),cudaMemcpyDeviceToHost));
#ifdef DEBUG
       printf("^^ begin CUDA MEMCPY OUT gamma_plus ^^\n");
#endif
       HANDLE_ERROR(cudaMemcpy(Gamma_plus,dGamma,sizeof(double),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaStreamSynchronize(0));
       //HANDLE_ERROR(cudaMemcpy(debug,cmplx_debug,sizeof(double2)*Ntri*Nbands*Nbands*nptk,cudaMemcpyDeviceToHost));
#endif
#ifdef DEBUG
       printf("^^ scale ^^\n");
#endif
       *WP3_plus=*WP3_plus/nptk;
#ifdef DEBUG
       printf("^^ CUDA ALL FINISH ^^\n");
#endif
    }
    (*Gamma_plus)=(*Gamma_plus)*5.60626442*double(1e8)/nptk; // THz
#ifdef DEBUG
       printf("^^ ALL FINISH ^^\n");
       printf("+++ Gamma_plus = %lf, WP3_plus = %lf\n", (*Gamma_plus), (*WP3_plus));
#endif
}


void run_cuda_rta_minus_(int* _rank, int* _mm, int* _nband, double* _scalebroad, int* _nptk, double* _T, double* rlattvec, double* energy, double* velocity, double_complex* eigenvect, int* _Nlist, int* list, int* _Ntri,double* Phi, double* R_j, double* R_k, int* Index_i, int* Index_j, int* Index_k, int* IJK, double* Gamma_minus, double* WP3_minus, int* types, double* masses, bool* _onlyharmonic, int* Ngrid, int* _nelements, int* _natoms)
{
    int nelements = *_nelements;
    int natoms = *_natoms;
    bool onlyharmonic = *_onlyharmonic;
    
    int nptk = *_nptk;
    double scalebroad = *_scalebroad;
    int rank = *_rank;
    int mm = *_mm;
    int Nbands = *_nband;
    int NList = *_Nlist;
    int Ntri = *_Ntri;
    double T = *_T;
    *Gamma_minus=0.;
    *WP3_minus=0.;

    int i=(mm-1)%Nbands;
    int ll=int((mm-1)/Nbands);
    int* q = &IJK[(list[ll]-1)*3];
    double omega=energy[list[ll]-1+i*nptk];
    if(omega!=0) {
#ifdef NO_CUDA
    for(int j=0;j<Nbands;++j) for(int ii=0;ii<nptk;++ii) {
    int qprime0 = IJK[ii*3  ];
    int qprime1 = IJK[ii*3+1];
    int qprime2 = IJK[ii*3+2];
    int Ngrid0 = Ngrid[0];
    int Ngrid1 = Ngrid[1];
    int Ngrid2 = Ngrid[2];
    int q1 = q[0];
    int q2 = q[1];
    int q3 = q[2];
    //realqprime=matmul(rlattvec,qprime/dble(ngrid))
    double realqprime0 = rlattvec[0]*(qprime0/(double)Ngrid0)+rlattvec[3]*(qprime1/(double)Ngrid1)+rlattvec[6]*(qprime2/(double)Ngrid2);
    double realqprime1 = rlattvec[1]*(qprime0/(double)Ngrid0)+rlattvec[4]*(qprime1/(double)Ngrid1)+rlattvec[7]*(qprime2/(double)Ngrid2);
    double realqprime2 = rlattvec[2]*(qprime0/(double)Ngrid0)+rlattvec[5]*(qprime1/(double)Ngrid1)+rlattvec[8]*(qprime2/(double)Ngrid2);
    double omegap=energy[ii+j*nptk];
    double fBEprime = 1./(exp(hbar*omegap/Kb/T)-1.);
//printf("wy-test-fBEprime %f\n", fBEprime);
    //--------BEGIN absorption process-----------
    for(int k=0;k<Nbands;++k) {
    int qdprime0 = (q1-qprime0+Ngrid0)%Ngrid0;
    int qdprime1 = (q2-qprime1+Ngrid1)%Ngrid1;
    int qdprime2 = (q3-qprime2+Ngrid2)%Ngrid2;
    // realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
    double realqdprime0 = rlattvec[0]*(qdprime0/(double)Ngrid0)+rlattvec[3]*(qdprime1/(double)Ngrid1)+rlattvec[6]*(qdprime2/(double)Ngrid2);
    double realqdprime1 = rlattvec[1]*(qdprime0/(double)Ngrid0)+rlattvec[4]*(qdprime1/(double)Ngrid1)+rlattvec[7]*(qdprime2/(double)Ngrid2);
    double realqdprime2 = rlattvec[2]*(qdprime0/(double)Ngrid0)+rlattvec[5]*(qdprime1/(double)Ngrid1)+rlattvec[8]*(qdprime2/(double)Ngrid2);
    int ss = qdprime0+Ngrid0*(qdprime1+qdprime2*Ngrid1);
//printf("wy-test_ss %d\n",ss);
    double omegadp = energy[ss+k*nptk];
//printf("wy-test_omegadp %f\n",omegap);
    if (omegap!=0 && omegadp!=0) {
        double sigma=scalebroad*base_sigma( velocity[ii+j*nptk]-velocity[ss+k*nptk], velocity[ii+j*nptk+nptk*Nbands]-velocity[ss+k*nptk+nptk*Nbands], velocity[ii+j*nptk+2*nptk*Nbands]-velocity[ss+k*nptk+2*nptk*Nbands], Ngrid0, Ngrid1, Ngrid2, rlattvec);
//printf("wy-test-sigma %f\n", sigma);
        if(abs(omega-omegap-omegadp)<=(2.*sigma)) {
            double fBEdprime=1.0/(exp(hbar*omegadp/Kb/T)-1.0);
            double tmp1 = omega-omegap-omegadp;
            double WP3=(fBEprime+fBEdprime+1)*
                        exp(-(tmp1*tmp1)/(sigma*sigma))/
                        sigma/sqrt(pi)/(omega*omegap*omegadp);
	    	*WP3_minus += WP3;
//printf("wy-test-WP3 %f\n", WP3);
            if (!onlyharmonic) {
                double_complex Vp=Vp_minus(i,j,k,list[ll]-1,ii,ss,
                                   realqprime0,realqprime1,realqprime2,
                                   realqdprime0,realqdprime1,realqdprime2,
                                   eigenvect,Ntri,Phi,R_j,R_k,
                                   Index_i,Index_j,Index_k,
                                   masses, types, nptk, Nbands);
		//double absVp = (Vp.x*Vp.x+Vp.y*Vp.y);
		//cuda_atomic_add(Gamma_plus,hbarp*pi/4.*WP3*absVp);
		//double absVp = Vp.x+Vp.y;
		//double absVp = Vp.x;
                //*Gamma_plus += absVp;
		double absVp2 = Vp.x*Vp.x+Vp.y*Vp.y;
		//printf("wy_test_Vp %f\n", absVp2);
		*Gamma_minus += hbarp*pi/double(4)*WP3*absVp2;
	
		
		//double absVp = sqrt(Vp.x*Vp.x+Vp.y*Vp.y);
                //cuda_atomic_add(Gamma_plus,hbarp*pi/4.*WP3*(absVp*absVp));
            }
        }
    }
    }}
#else
        //HANDLE_ERROR(cudaMemset(dGamma,0,sizeof(double)));
        //HANDLE_ERROR(cudaMemset(dWP3,0,sizeof(double)));
#ifdef DEBUG
        printf("^^ CUDA MEMSET FINISH ^^\n");
#endif
       dim3 blk( (Nbands+4-1)/4, (nptk+32-1)/32, (Nbands+4-1)/4 );
       dim3 blksize( 4,32,4 );
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL ^^\n");
#endif
       rta_minus_kernel<<<blk,blksize>>>(Nbands, scalebroad, nptk, T, drlattvec, denergy, dvelocity, deigenvect, NList, dlist, Ntri, dPhi, dR_j, dR_k,  dIndex_i, dIndex_j, dIndex_k, dIJK, dGamma,dWP3, dtypes, dmasses, onlyharmonic, Ngrid[0],Ngrid[1], Ngrid[2],q[0],q[1],q[2], omega, i, ll, d_debug);
#ifdef DEBUG
       printf("^^ begin CUDA KERNEL SYNC ^^\n");
#endif
       HANDLE_ERROR(cudaGetLastError());
       HANDLE_ERROR(cudaStreamSynchronize(0));
       HANDLE_ERROR(cudaGetLastError());
/*       
       HANDLE_ERROR(cudaMemcpy(debug,d_debug,sizeof(double)*nptk,cudaMemcpyDeviceToHost));
       int count = 0; int icount=0;
       double sum1=0.0, sumr=0.0;
       for(int ii=0;ii<nptk;++ii) {
           int qprime0 = IJK[ii*3  ];
           int qprime1 = IJK[ii*3+1];
           int qprime2 = IJK[ii*3+2];
           int Ngrid0 = Ngrid[0];
           int Ngrid1 = Ngrid[1];
           int Ngrid2 = Ngrid[2];
           //double ref = (rlattvec[2]*(qprime0/(double)Ngrid0))+(rlattvec[5]*(qprime1/(double)Ngrid1))+(rlattvec[8]*(qprime2/(double)Ngrid2));
           double rr = ((rlattvec[0]*(qprime0/(double)Ngrid0))+(rlattvec[3]*(qprime1/(double)Ngrid1))+(rlattvec[6]*(qprime2/(double)Ngrid2)));
           double ref=0;
           for(int dd=0;dd<Ntri;++dd) ref+=rr*R_j[dd*3];
           //double ref = rlattvec[3]*(qprime1/(double)Ngrid1)+rlattvec[6]*(qprime2/(double)Ngrid2);
           //double ref = rlattvec[3]*(qprime1/(double)Ngrid1);
           //double ref = (rlattvec[0]*(qprime0/(double)Ngrid0))+(rlattvec[3]*(qprime1/(double)Ngrid1));
           if((debug[ii]!=0 || ref!=0) && icount++<15) printf("%d : gpu %e, ref %e\n",ii, debug[ii], ref);
           if((ref==0.0&&fabs(debug[ii])>1e-8)||(ref!=0.0 && fabs(debug[ii]-ref)/ref > 1e-8)) {
                printf("Error at %d : gpu %e, ref %e, error %e\n", ii, debug[ii], ref, fabs(debug[ii]-ref));
                ++count;
                if(count>15) exit(-1);
           }
           sum1 += debug[ii];
           sumr += ref;
       }
       printf("%e, %e\n",sum1, sumr);*/
#ifdef DEBUG
       printf("^^ finish CUDA KERNEL ^^\n");
       printf("^^ begin CUDA MEMCPY OUT wp3_minus ^^\n");
#endif
       HANDLE_ERROR(cudaMemcpy(WP3_minus,dWP3,sizeof(double),cudaMemcpyDeviceToHost));
#ifdef DEBUG
       printf("^^ begin CUDA MEMCPY OUT gamma_minus ^^\n");
#endif
       HANDLE_ERROR(cudaMemcpy(Gamma_minus,dGamma,sizeof(double),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaStreamSynchronize(0));
       //HANDLE_ERROR(cudaMemcpy(debug,cmplx_debug,sizeof(double2)*Ntri*Nbands*Nbands*nptk,cudaMemcpyDeviceToHost));
#endif
#ifdef DEBUG
       printf("^^ scale ^^\n");
#endif
       (*WP3_minus)=(*WP3_minus)*double(5e-1)/nptk; //wy
#ifdef DEBUG
       printf("^^ CUDA ALL FINISH ^^\n");
#endif
    }
    (*Gamma_minus)=(*Gamma_minus)*5.60626442*double(1e8)/nptk; // THz
#ifdef DEBUG
       printf("^^ ALL FINISH ^^\n");
       printf("+++ Gamma_minus = %lf, WP3_minus = %lf\n", (*Gamma_minus), (*WP3_minus));
#endif
}


void finalize_cuda_ind_() {
#ifdef DEBUG
       printf("^^ RESET... ^^\n");
#endif
        HANDLE_ERROR(cudaFree(drlattvec));
        HANDLE_ERROR(cudaFree(dN));
        HANDLE_ERROR(cudaFree(denergy));
        HANDLE_ERROR(cudaFree(dvelocity));
        HANDLE_ERROR(cudaFree(deigenvect));
        HANDLE_ERROR(cudaFree(dPhi));
        HANDLE_ERROR(cudaFree(dR_j));
        HANDLE_ERROR(cudaFree(dR_k));
        HANDLE_ERROR(cudaFree(dmasses));
        HANDLE_ERROR(cudaFree(dtypes));
        HANDLE_ERROR(cudaFree(dIndex_i));
        HANDLE_ERROR(cudaFree(dIndex_j));
        HANDLE_ERROR(cudaFree(dIndex_k));
        HANDLE_ERROR(cudaFree(dIJK));
        HANDLE_ERROR(cudaFree(dIndof2ndPhonon));
        HANDLE_ERROR(cudaFree(dIndof3rdPhonon));
        HANDLE_ERROR(cudaFree(dGamma));
        HANDLE_ERROR(cudaFree(dWP3));
        HANDLE_ERROR(cudaFree(dlist));

    HANDLE_ERROR(cudaDeviceReset());
}


/*void finalize_cuda_np_() {
#ifdef DEBUG
       printf("^^ RESET... ^^\n");
#endif
        HANDLE_ERROR(cudaFree(drlattvec));
        HANDLE_ERROR(cudaFree(denergy));
        HANDLE_ERROR(cudaFree(dvelocity));
        HANDLE_ERROR(cudaFree(dmasses));
        HANDLE_ERROR(cudaFree(dtypes));
        HANDLE_ERROR(cudaFree(dIJK));
        HANDLE_ERROR(cudaFree(dN_plus));
        HANDLE_ERROR(cudaFree(dP_plus));
        HANDLE_ERROR(cudaFree(dlist));

    HANDLE_ERROR(cudaDeviceReset());
}*/

void finalize_cuda_rta_() {
#ifdef DEBUG
       printf("^^ RESET... ^^\n");
#endif
        HANDLE_ERROR(cudaFree(drlattvec));
        HANDLE_ERROR(cudaFree(denergy));
        HANDLE_ERROR(cudaFree(dvelocity));
        HANDLE_ERROR(cudaFree(deigenvect));
        HANDLE_ERROR(cudaFree(dPhi));
        HANDLE_ERROR(cudaFree(dR_j));
        HANDLE_ERROR(cudaFree(dR_k));
        HANDLE_ERROR(cudaFree(dmasses));
        HANDLE_ERROR(cudaFree(dtypes));
        HANDLE_ERROR(cudaFree(dIndex_i));
        HANDLE_ERROR(cudaFree(dIndex_j));
        HANDLE_ERROR(cudaFree(dIndex_k));
        HANDLE_ERROR(cudaFree(dIJK));
        HANDLE_ERROR(cudaFree(dGamma));
        HANDLE_ERROR(cudaFree(dWP3));
        HANDLE_ERROR(cudaFree(dlist));

    HANDLE_ERROR(cudaDeviceReset());
}
}
