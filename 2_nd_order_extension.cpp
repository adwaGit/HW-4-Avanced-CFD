// second order extensions : 


// 
#include <iostream>
#include <math.h>
#include <cmath>
#include <conio.h>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <fstream>
using namespace std;
using std::cout;


double Psi_plus_vector_minus_1_rho = 0;
double Psi_plus_vector_minus_1_u = 0;
double Psi_plus_vector_minus_1_p = 0;


vector<double> Psi_plus_vector_rho(i_max) ; // does not need to have i_max + 1 cell but has a 0-1 value called Psi_plus_vector_minus_1 ;
vector<double> Psi_plus_vector_u(i_max) ; 
vector<double> Psi_plus_vector_p(i_max) ; 


vector<double>Psi_minus_vector_rho(i_max + 1) ;   // needs i_max + 1 to be populated;
vector<double>Psi_minus_vector_u(i_max + 1 ) ; 
vector<double>Psi_minus_vector_p(i_max + 1 ) ; 


 
double set_psi_plus(int i)
{
double r_plus_rho = 0;
double r_plus_u = 0;
double r_plus_p = 0;

double r_plus_ghost_cell_minus_half_rho  = 0;
double r_plus_ghost_cell_minus_half_u = 0;
double r_plus_ghost_cell_minus_half_p  = 0;


double r_plus_i_max_plus_half_rho = 0 ;  
double r_plus_i_max_plus_half_u = 0 ;  
double r_plus_i_max_plus_half_p = 0 ;  



// calculates from 0 to less than i_max 
if ( 0 <= i <= i_max-1 )
{  
r_plus_rho = (rho_n[i+2] - rho_n[i+1]) / (rho_n[i+1] - rho_n[i])  ; 
r_plus_u = (u_n[i+2] - u_n[i+1]) / (u_n[i+1] - u_n[i])  ; 
r_plus_p = (p_n[i+2] - p_n[i+1]) / (p_n[i+1] - p_n[i])  ; 

Psi_plus_vector_rho[i] = (r_plus_rho + abs(r_plus_rho)) / (1 + r_plus_rho) ; 
Psi_plus_vector_u[i] = (r_plus_u + abs(r_plus_u)) / (1 + r_plus_u) ; 
Psi_plus_vector_p[i] = (r_plus_p + abs(r_plus_p)) / (1 + r_plus_p) ; 
}

set_extra_ghost_cells(); 

if (i == -1)
r_plus_ghost_cell_minus_half_rho = (rho_n[0+2-1] - rho_n[0 + 1 -1]) / (rho_n[0 + 1 -1 ] - rho_minus_1)  ; 
r_plus_ghost_cell_minus_half_u = (u_n[0+2-1] - u_n[0 + 1 -1]) / (u_n[0 + 1 -1 ] - u_minus_1);
r_plus_ghost_cell_minus_half_p = (p_n[0+2-1] - p_n[0 + 1 -1]) / (p_n[0 + 1 -1 ] - p_minus_1);

Psi_plus_vector_minus_1_rho= (r_plus_ghost_cell_minus_half_rho + abs(r_plus_ghost_cell_minus_half_rho)) / (1 + r_plus_ghost_cell_minus_half_rho) ; 
Psi_plus_vector_minus_1_u = (r_plus_ghost_cell_minus_half_u + abs(r_plus_ghost_cell_minus_half_u)) / (1 + r_plus_ghost_cell_minus_half_u) ; 
Psi_plus_vector_minus_1_p = (r_plus_ghost_cell_minus_half_p + abs(r_plus_ghost_cell_minus_half_p)) / (1 + r_plus_ghost_cell_minus_half_p) ; 



if (i == i_max)
{

r_plus_i_max_plus_half_rho = (rho_plus_2 + rho_n[i_max + 1])/ (rho_n[i_max+1] - rho_n[i_max]); 
r_plus_i_max_plus_half_u  = (u_plus_2 + u_n[i_max + 1])/ (u_n[i_max+1] - u_n[i_max]); 
r_plus_i_max_plus_half_p  = (p_plus_2 + p_n[i_max + 1])/ (p_n[i_max+1] - p_n[i_max]); 

Psi_plus_vector_rho[i_max]  = (r_plus_i_max_plus_half_rho+ abs(r_plus_i_max_plus_half_rho)) / (1 + r_plus_i_max_plus_half_rho) ; 
Psi_plus_vector_u[i_max]  = (r_plus_i_max_plus_half_u+ abs(r_plus_i_max_plus_half_u)) / (1 + r_plus_i_max_plus_half_u) ; 
Psi_plus_vector_p[i_max]  = (r_plus_i_max_plus_half_p+ abs(r_plus_i_max_plus_half_p)) / (1 + r_plus_i_max_plus_half_p) ; 

}

return 0;
}

double set_psi_minus(int i)
{

double r_minus_rho = 0;
double r_minus_u = 0;
double r_minus_p = 0;


double r_minus_ghost_cell_plus_half_rho = 0;
double r_minus_ghost_cell_plus_half_u  = 0;
double r_minus_ghost_cell_plus_half_p  = 0;


double r_minus_i_max_plus_1_point_5_rho= 0; 
double r_minus_i_max_plus_1_point_5_u= 0; 
double r_minus_i_max_plus_1_point_5_p= 0; 

r_minus_rho = (rho_n[i] - rho_n[i-1]) / (rho_n[i+1] - rho_n[i])  ; 
r_minus_u = (u_n[i] - u_n[i-1]) / (u_n[i+1] - u_n[i])  ; 
r_minus_p = (p_n[i] - p_n[i-1]) / (p_n[i+1] - p_n[i])  ; 

Psi_minus_vector_rho[i] = (r_minus_rho + abs(r_minus_rho)) / (1 + r_minus_rho) ;
Psi_minus_vector_u[i] = (r_minus_u + abs(r_minus_u)) / (1 + r_minus_u) ; 
Psi_minus_vector_p[i] = (r_minus_p + abs(r_minus_p)) / (1 + r_minus_p) ;

set_extra_ghost_cells(); 

if (i == 0)
{
r_minus_ghost_cell_plus_half_rho = (rho_n[0] - rho_minus_1) / (rho_n[0+1] - rho_n[0])  ;
r_minus_ghost_cell_plus_half_u = (u_n[0] - u_minus_1) / (u_n[0+1] - u_n[0])  ;
r_minus_ghost_cell_plus_half_p = (p_n[0] - p_minus_1) / (p_n[0+1] - p_n[0])  ;

Psi_minus_vector_rho[0] = (r_minus_ghost_cell_plus_half_rho + abs(r_minus_ghost_cell_plus_half_rho)) / (1 + r_minus_ghost_cell_plus_half_rho) ; 
Psi_minus_vector_u[0] = (r_minus_ghost_cell_plus_half_u + abs(r_minus_ghost_cell_plus_half_u)) / (1 + r_minus_ghost_cell_plus_half_u) ; 
Psi_minus_vector_p[0] = (r_minus_ghost_cell_plus_half_p + abs(r_minus_ghost_cell_plus_half_p)) / (1 + r_minus_ghost_cell_plus_half_p) ; 
}

if ( i == i_max + 1)
{
r_minus_i_max_plus_1_point_5_rho =  (rho_n[i_max +1 ] - rho_n[i_max + 1 -1]) / (rho_plus_2 - rho_n[i_max + 1]);
r_minus_i_max_plus_1_point_5_u =  (u_n[i_max +1 ] - u_n[i_max + 1 -1]) / (u_plus_2 - u_n[i_max + 1]);
r_minus_i_max_plus_1_point_5_p =  (p_n[i_max +1 ] - p_n[i_max + 1 -1]) / (p_plus_2 - p_n[i_max + 1]);

Psi_minus_vector_rho[i_max + 1]  =  (r_minus_i_max_plus_1_point_5_rho + abs(r_minus_i_max_plus_1_point_5_rho)) / (1 + r_minus_i_max_plus_1_point_5_rho) ;  
Psi_minus_vector_u[i_max + 1]  =  (r_minus_i_max_plus_1_point_5_u + abs(r_minus_i_max_plus_1_point_5_u)) / (1 + r_minus_i_max_plus_1_point_5_u) ;  
Psi_minus_vector_p[i_max + 1]  =  (r_minus_i_max_plus_1_point_5_p + abs(r_minus_i_max_plus_1_point_5_p)) / (1 + r_minus_i_max_plus_1_point_5_p) ;  

}


return 0;
}





vector<double>set_left_state_2nd_order(int i)
{

vector<double>rho_a_u_ht_M_p_vector_left(6);
std::fill(rho_a_u_ht_M_p_vector_left.begin(), rho_a_u_ht_M_p_vector_left.end(), 0.0);

double T_left = 0; 
double a_left = 0; 
double total_enthalpy_left = 0;
double M_left = 0;
double p_left = 0; 
double rho_left = 0; 
double u_left = 0;
double epsilon = 0; 
double Kappa_2nd_order = 0; 


epsilon = 0 ; 
Kappa_2nd_order = 0; 

if (i >= 1)
{ 
rho_left= rho_n[i] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_plus_vector_rho[i-1]*(rho_n[i] - rho_n[i-1])\
 + (1 + Kappa_2nd_order)*Psi_minus_vector_rho[i]*(rho_n[i+1] - rho_n[i])) ; 

u_left= u_n[i] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_plus_vector_u[i-1]*(u_n[i] - u_n[i-1])\
 + (1 + Kappa_2nd_order)*Psi_minus_vector_u[i]*(u_n[i+1] - u_n[i])) ;     

p_left= p_n[i] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_plus_vector_p[i-1]*(p_n[i] - p_n[i-1])\
 + (1 + Kappa_2nd_order)*Psi_minus_vector_p[i]*(p_n[i+1] - p_n[i])) ; 
}


 if (i == 0)
 {

rho_left= rho_n[0] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_plus_vector_minus_1_rho*(rho_n[0] - rho_minus_1)\
 + (1 + Kappa_2nd_order)*Psi_minus_vector_rho[0]*(rho_n[0+1] - rho_n[0] )) ; 

u_left= u_n[0] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_plus_vector_minus_1_u*(u_n[0] - u_minus_1)\
 + (1 + Kappa_2nd_order)*Psi_minus_vector_u[0]*(u_n[0+1] - u_n[0] )) ; 

 p_left= p_n[0] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_plus_vector_minus_1_p*(p_n[0] - p_minus_1)\
 + (1 + Kappa_2nd_order)*Psi_minus_vector_p[0]*(p_n[0+1] - p_n[0] )) ; 

 }  


T_left = p_left/(rho_left*R_univ) ; 
a_left = sqrt(Gamma*R_univ*T_left);
total_enthalpy_left = pow(a_left,2)/(Gamma-1) + 0.5*(pow(u_left,2)) ; 
M_left = u_left/a_left ; 

rho_a_u_ht_M_p_vector_left[0] = rho_left ; 
rho_a_u_ht_M_p_vector_left[1] = a_left;
rho_a_u_ht_M_p_vector_left[2] = u_left ; 
rho_a_u_ht_M_p_vector_left[3] = total_enthalpy_left; 
rho_a_u_ht_M_p_vector_left[4] = M_left; 
rho_a_u_ht_M_p_vector_left[5] = p_left; 

return rho_a_u_ht_M_p_vector_left; 

}


vector<double>set_right_state_2nd_order(int i)
{

vector<double>rho_a_u_ht_M_p_vector_right(6);
std::fill(rho_a_u_ht_M_p_vector_right.begin(), rho_a_u_ht_M_p_vector_right.end(), 0.0);

double T_right= 0; 
double a_right= 0; 
double total_enthalpy_right = 0;
double M_right = 0;
double p_right = 0; 
double rho_right = 0; 
double u_right = 0;
double epsilon = 0; 
double Kappa_2nd_order = 0; 


epsilon = 0 ; 
Kappa_2nd_order = 0; 


// calculating rho : 
rho_right = rho_n[i+1] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_minus_vector_rho[i+1]*(rho_n[i+2] - rho_n[i+1])\
 + (1 + Kappa_2nd_order)*Psi_plus_vector_rho[i]*(rho_n[i+1] - rho_n[i])) ;

// calculating u : 
u_right = u_n[i+1] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_minus_vector_u[i+1]*(u_n[i+2] - u_n[i+1])\
 + (1 + Kappa_2nd_order)*Psi_plus_vector_u[i]*(rho_n[i+1] - rho_n[i])) ;
// calculating p : 

p_right = p_n[i+1] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_minus_vector_p[i+1]*(p_n[i+2] - p_n[i+1])\
 + (1 + Kappa_2nd_order)*Psi_plus_vector_p[i]*(p_n[i+1] - p_n[i])) ;

 if (i == i_max)
 {
// calculating rho : 
rho_right = rho_n[i+1] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_minus_vector_rho[i_max+1]*(rho_plus_2 - rho_n[i_max+1])\
 + (1 + Kappa_2nd_order)*Psi_plus_vector_rho[i_max]*(rho_n[i_max+1] - rho_n[i_max])) ;

// calculating u : 
u_right = u_n[i_max+1] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_minus_vector_u[i_max+1]*(u_plus_2 - u_n[i_max+1])\
 + (1 + Kappa_2nd_order)*Psi_plus_vector_u[i_max]*(u_n[i_max+1] - u_n[i_max])) ;
// calculating p : 

p_right = p_n[i_max+1] + (epsilon/4.0)*((1 - Kappa_2nd_order)*Psi_minus_vector_p[i_max+1]*(p_plus_2 - p_n[i_max+1])\
 + (1 + Kappa_2nd_order)*Psi_plus_vector_p[i_max]*(p_n[i_max+1] - p_n[i_max])) ;

 }



T_right = p_right/(rho_right*R_univ) ; 
a_right = sqrt(Gamma*R_univ*T_right);
total_enthalpy_right = pow(a_right,2)/(Gamma-1) + 0.5*(pow(u_right,2)) ; 
M_right = u_right/a_right ; 

rho_a_u_ht_M_p_vector_right[0] = rho_right ; 
rho_a_u_ht_M_p_vector_right[1] = a_right;
rho_a_u_ht_M_p_vector_right[2] = u_right ; 
rho_a_u_ht_M_p_vector_right[3] = total_enthalpy_right; 
rho_a_u_ht_M_p_vector_right[4] = M_right; 
rho_a_u_ht_M_p_vector_right[5] = p_right; 

return rho_a_u_ht_M_p_vector_left; 

}




std::fill(Psi_plus_vector_rho.begin(), Psi_plus_vector_rho.end(), 0.0);
std::fill(Psi_plus_vector_u.begin(), Psi_plus_vector_u.end(), 0.0);
std::fill(Psi_plus_vector_p.begin(), Psi_plus_vector_p.end(), 0.0);

std::fill(Psi_minus_vector_rho.begin(), Psi_minus_vector_rho.end(), 0.0); 
std::fill(Psi_minus_vector_u.begin(), Psi_minus_vector_u.end(), 0.0); 
std::fill(Psi_minus_vector_p.begin(), Psi_minus_vector_p.end(), 0.0); 

for (int i = -1 ; i <= i_max; i++) 
{ 

set_psi_plus(i) ;

}



for (int i = 0; i <= i_max + 1; i++) 

{

set_psi_minus(i); 

}



for (int i = 0; i <= i_max; i++) 
{
vector<vector<double>>rho_a_u_ht_M_p_vector_left(i_max +1 , vector<double>(6));
vector<vector<double>>rho_a_u_ht_M_p_vector_right(i_max +1 , vector<double>(6));

rho_a_u_ht_M_p_vector_left[i] =  set_left_state_2nd_order(i); 
rho_a_u_ht_M_p_vector_right[i] = set_right_state_2nd_order(i);

}


for (int i = 0; i <= i_max; i++) 
{
van_leer_output = van_leer_flux_calculator(rho_a_u_ht_M_p_vector_left[i],rho_a_u_ht_M_p_vector_right[i]);

 van_leer_flux_i_plus_half_eq_1[i] = van_leer_output[0] ; 
 van_leer_flux_i_plus_half_eq_2[i] = van_leer_output[1] ; 
 van_leer_flux_i_plus_half_eq_3[i] = van_leer_output[2] ; 
}