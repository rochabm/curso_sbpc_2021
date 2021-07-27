#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
from scipy import integrate
import numpy as np

def immune_response_v1 (P, t, pi_v, k_v1):
    """
    Simple Model
    """
    V = P[0]

    
####################################################################################
#Equações
####################################################################################
    dV_dt = pi_v*V - k_v1*V
    
    return [dV_dt]



def immune_response_v2 (P, t, pi_v, k_v3, alpha_Ap, Ap_0):
    """
    Simple Model
    """
    V,Ap = P[0],P[1]

    
####################################################################################
#Equações
####################################################################################
    dV_dt = pi_v*V - k_v3*Ap*V
    dAp_dt = alpha_Ap*(Ap_0 - Ap)
    
    return [dV_dt, dAp_dt]


def immune_response_v3 (P, t, pi_v, k_v3, alpha_Ap, Ap_0, beta_Ap, k_ap1, k_ap2, delta_Apm):
    """
    Simple Model
    """
    V,Ap,Apm = P[0],P[1],P[2]

    
####################################################################################
#Equações
####################################################################################
    dV_dt = pi_v*V - k_v3*Apm*V
    dAp_dt = alpha_Ap*(Ap_0 - Ap) - beta_Ap*Ap*(k_ap1*(V)/(k_ap2 + V)) 
    dApm_dt = beta_Ap*Ap*(k_ap1*(V)/(k_ap2 + V)) - delta_Apm*Apm
    
    return [dV_dt, dAp_dt, dApm_dt]




def immune_response_v4 (P, t, pi_v, k_v3, alpha_Ap, Ap_0, beta_Ap, k_ap1, k_ap2, delta_Apm, k_v2, beta_Tk,pi_T,k_te1, delta_te,Tkn_0):
    """
    Simple Model
    """
    V,Ap,Apm,Tkn,Tke = P[0],P[1],P[2],P[3],P[4]

    
####################################################################################
#Equações
####################################################################################
    dV_dt = pi_v*V - k_v3*Apm*V - k_v2* V *Tke
    dAp_dt = alpha_Ap*(Ap_0 - Ap) - beta_Ap*Ap*(k_ap1*(V)/(k_ap2 + V)) 
    dApm_dt = beta_Ap*Ap*(k_ap1*(V)/(k_ap2 + V)) - delta_Apm*Apm
    dTkn_dt = beta_Tk*(Tkn_0 - Tkn) - pi_T*Apm*Tkn 
    dTke_dt = pi_T*Apm*Tkn + k_te1*Apm*Tke - delta_te*Tke
    
    return [dV_dt, dAp_dt, dApm_dt, dTkn_dt, dTke_dt]


def immune_response_v5 (P, t, pi_v, k_v3, alpha_Ap, Ap_0, beta_Ap, k_ap1, k_ap2, delta_Apm, k_v2, beta_Tk, pi_T, k_te1, delta_te, Tkn_0, k_apm, k_tk, pi_c_apm, pi_c_i, pi_c_tke, delta_c):
    """
    Simple Model
    """
    V,Ap,Apm,Tkn,Tke,I,C = P[0],P[1],P[2],P[3],P[4],P[5],P[6]

    
####################################################################################
#Equações
####################################################################################
    dV_dt = pi_v*V - k_v2*V*Tke - k_v3*V*Apm 
    dAp_dt = alpha_Ap*C*(Ap_0 - Ap) - beta_Ap*Ap*(k_ap1*(V)/(k_ap2 + V))
    dApm_dt = beta_Ap*Ap*(k_ap1*(V)/(k_ap2 + V)) - delta_Apm*Apm - k_apm * Apm * V
    dI_dt = k_apm * Apm * V + k_tk * Tke * V - delta_Apm*I
    dTkn_dt = beta_Tk*(C)*(Tkn_0 - Tkn) - pi_T*(C+1)*Apm*Tkn 
    dTke_dt = pi_T*(C+1)*Apm*Tkn + k_te1*Apm*Tke - delta_te*Tke - k_tk * Tke * V
    dC_dt = pi_c_apm *V*Apm + pi_c_i*I + pi_c_tke*V*Tke - delta_c * C
    
    return [dV_dt, dAp_dt, dApm_dt, dTkn_dt, dTke_dt,dI_dt,dC_dt]



def immune_response_completo (P, t, pi_v, c_v1, c_v2, k_v1, k_v2, alpha_Ap, beta_Ap,
                     c_ap1, c_ap2, delta_Apm, alpha_Tn, beta_tk, pi_tk, delta_tk,
                     alpha_B, pi_B1, pi_B2, beta_ps, beta_pl, beta_Bm,
                     delta_S, delta_L, gamma_bm, k_bm1, k_bm2, pi_AS,
                     pi_AL, delta_ag, delta_am, alpha_th, beta_th, pi_th, delta_th, Ap0, Thn0, Tkn0, B0, 
                     pi_c_apm, pi_c_i,pi_c_tke,delta_c, beta_apm, k_v3, beta_tke):
    """
    Simple Model
    """
    V,Ap,Apm,Thn,The,Tkn,Tke,B,Ps,Pl,Bm,A_M, A_G, I ,C = P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8], P[9], P[10], P[11], P[12], P[13], P[14]

    
####################################################################################
#Equações
####################################################################################
    dV_dt = pi_v*V - k_v1*V*A_M - k_v1*V*A_G - k_v2*V*Tke - k_v3*V*Apm #- ((c_v1*V)/(c_v2+V))
    dAp_dt = C*(Ap0 - Ap) - beta_Ap*Ap*(c_ap1*(V)/(c_ap2 + V)) #alpha_Ap*(1+C)*(Ap0 - Ap) - beta_Ap*Ap*(c_ap1*(V)/(c_ap2 + V))
    dApm_dt = beta_Ap*Ap*(c_ap1*(V)/(c_ap2 + V)) - delta_Apm*Apm - beta_apm * Apm * V
    dI_dt = beta_apm * Apm * V + beta_tke * Tke * V - delta_Apm*I
    dThn_dt = alpha_th*(Thn0 - Thn) - beta_th*Apm*Thn 
    dThe_dt = beta_th*Apm*Thn + pi_th*Apm*The - delta_th*The
    dTkn_dt = (C)*(Tkn0 - Tkn) - beta_tk*(C+1)*Apm*Tkn #alpha_Tn*(1+C)*(Tkn0 - Tkn) - beta_tk*Apm*Tkn
    dTke_dt = beta_tk*(C+1)*Apm*Tkn + pi_tk*Apm*Tke - delta_tk*Tke - beta_tke * Tke * V
    dB_dt = alpha_B*(B0 - B) + pi_B1*V*B + pi_B2*The*B - beta_ps*Apm*B - beta_pl*The*B - beta_Bm*The*B
    dPs_dt = beta_ps*Apm*B - delta_S*Ps
    dPl_dt = beta_pl*The*B - delta_L*Pl + gamma_bm*Bm 
    dBm_dt = beta_Bm*The*B + k_bm1*Bm*(1 - Bm/(k_bm2)) - gamma_bm*Bm
    dA_M_dt = pi_AS*Ps - delta_am*A_M # anticorpos vida curta
    dA_G_dt = pi_AL*Pl - delta_ag*A_G #memória imune
    dC_dt = pi_c_apm *V*Apm + pi_c_i*I + pi_c_tke*V*Tke - delta_c * C
    
    return [dV_dt, dAp_dt, dApm_dt, dThn_dt, dThe_dt, dTkn_dt, dTke_dt, dB_dt, dPs_dt, dPl_dt, dBm_dt, dA_M_dt, dA_G_dt, dI_dt, dC_dt]
