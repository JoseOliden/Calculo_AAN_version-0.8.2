import streamlit as st
import pandas as pd
import numpy as np
import openpyxl
import os
import re
import io
import json
from datetime import datetime, timedelta
from scipy.optimize import root
import sympy as sp
import matplotlib.pyplot as plt
import tempfile
import base64
# ==============================================================================
#                        FUNCIONES AUXILIARES CALCULO
# ==============================================================================

# ----------------------------- Archivos Aesp ---------------------------------#

def Aesp(Cn_i, w_i,lam,tr,td,ti,tv,e):
  # Calcula la actividad específica
  C_i = lam/(1-np.exp(-lam*tr))
  D_i = np.exp(lam*td)
  H_i = tr/tv
  S_i = 1-np.exp(-lam*ti)
  return Cn_i*D_i*C_i*H_i/(S_i*w_i) # se agregó e

# ---------------------------- Calculo de alfa --------------------------------#
def crear_df_comparadores():
  # comparadores [Au Co Mo]
  #Q0_c = np.array([15.712, 2.041, 53.1])    #
  df = pd.DataFrame(
      data=[
        [1.00,      0.02272,  52346, 0.0000001358, 15.712, 0.000002977,    5.65,  306314, 1500, 10800, 1478.0],
        [1.32,      0.0113,   36082, 0.000051573,  2.041,  0.000000004167, 136.0, 306314, 1500, 10800, 1478.0],
        [0.0000846, 0.016053, 39272, 0.002088591,  50.365, 0.00000291994,  241,   299161, 900,  10800, 866.0]
      ],
      columns=["k0", "efe", "Cn", "w", "Q0", "lambda", "Er", "t_dec", "t_real", "t_irr", "t_vivo" ],
      index=["Au", "Co", "Mo"]
    )
  return df



def cal_alfa(df_comp):
  def equations(vars, *par):
    # Define el sistema de ecuaciones para hallar alfa
    alfa = vars[0]
    Aesp_1,k0_1,e_1,Er_1,Q0_1,Aesp_2,k0_2,e_2,Er_2,Q0_2, Aesp_3,k0_3,e_3,Er_3,Q0_3 = par
    eq1 = ((1-(Aesp_2/Aesp_1)*(k0_1/k0_2)*(e_1/e_2))**(-1) - (1-(Aesp_3/Aesp_1)*(k0_1/k0_3)*(e_1/e_3))**(-1))*(Q0_1 - 0.429)/(Er_1**(alfa)) - ((1-(Aesp_2/Aesp_1)*(k0_1/k0_2)*(e_1/e_2))**(-1))*(Q0_2 - 0.429)/(Er_2**(alfa)) + ((1-(Aesp_3/Aesp_1)*(k0_1/k0_3)*(e_1/e_3))**(-1))*(Q0_3 - 0.429)/(Er_3**(alfa))
    return [eq1]
  
  # calcula alfa
  #k0_c, e_c, Q0_c, Cn_c, w_c, lam_c,Er_c, td_c, tr_c, ti_c, tv_c = par_comp
  k0_c = df_comp["k0"].to_numpy(dtype="float64")
  e_c = df_comp["efe"].to_numpy(dtype="float64")
  Q0_c = df_comp["Q0"].to_numpy(dtype="float64")
  Cn_c = df_comp["Cn"].to_numpy(dtype="float64")
  w_c = df_comp["w"].to_numpy(dtype="float64")
  lam_c = df_comp["lambda"].to_numpy(dtype="float64")
  Er_c = df_comp["Er"].to_numpy(dtype="float64")
  td_c = df_comp["t_dec"].to_numpy(dtype="float64")
  tr_c = df_comp["t_real"].to_numpy(dtype="float64")
  ti_c = df_comp["t_irr"].to_numpy(dtype="float64")
  tv_c = df_comp["t_vivo"].to_numpy(dtype="float64")

  Aesp_c = np.zeros(len(k0_c))
  for i in range(len(k0_c)):
    Aesp_c[i] = Aesp(Cn_c[i], w_c[i],lam_c[i],tr_c[i],td_c[i],ti_c[i],tv_c[i],e_c[i]) # Activida especifica del elemento comparador
  initial_guesses = [0.0]
  par = (Aesp_c[0], k0_c[0], e_c[0], Er_c[0], Q0_c[0], Aesp_c[1], k0_c[1], e_c[1], Er_c[1], Q0_c[1], Aesp_c[2], k0_c[2], e_c[2], Er_c[2], Q0_c[2])
  solution = root(equations, initial_guesses, args = par)
  alfa = solution.x

  # Calcular f
  Q0_alfa_c = np.zeros(len(k0_c))
  for i in range(len(k0_c)):
    Q0_alfa_c[i] = cal_Q0_alfa_i(Q0_c[i],Er_c[i],alfa)
  f = cal_f_alfa(Q0_alfa_c,Aesp_c,e_c,k0_c)
  
  return alfa[0], f
# ---------------------------- Calculo de f --------------------------------#
def cal_f_alfa(Q0_alfa_c,Aesp_c,e_c,k0_c):
  # calcula f
  return ((k0_c[0]/k0_c[1])*(e_c[0]/e_c[1])*Q0_alfa_c[0]  - (Aesp_c[0]/Aesp_c[1])*Q0_alfa_c[1])/(Aesp_c[0]/Aesp_c[1] - (k0_c[0]/k0_c[1])*(e_c[0]/e_c[1]))
# ------------------------ Calculo de concentración ---------------------------#

def cal_Q0_alfa_i(Q0,Er,alfa):
  # calcula Q0_alfa del elemento i
  #(Q0-0.429)/(Er**alfa) + 0.429/((2*alfa+1)*0.55**alfa) literatura
  # (Q0-0.429)/(Er**alfa) + 0.429/(2*alfa+0.55**alfa) # MACROS
  return (Q0-0.429)/(Er**alfa) + 0.429/(2*alfa+0.55**alfa)

def conc(df_muestra, w,td_i,ti_i,tv_i,tr_i, df_comp_Au, w_Au,td_c_Au,ti_c_Au,tv_c_Au,tr_c_Au, alfa, f, geometria):
  #alfa = 0.226 # Forzar valor de alfa 
  #f = 34       # Forzar valor de f
  # Comparador Au
  #k0_c_Au, e_c_Au, Q0_c_Au, Cn_c_Au, w_c_Au, lam_c_Au, Er_c_Au, td_c_Au, tr_c_Au, ti_c_Au, tv_c_Au =  par_comp_Au
  st.dataframe(df_comp_Au)
  k0_c_Au = df_comp_Au["K0"].to_numpy(dtype="float64")
  if geometria == "50 mm":
    e_c_Au = df_comp_Au["EFIGAMMA50"].to_numpy(dtype="float64")*df_comp_Au["COI ROSSBACH"].to_numpy(dtype="float64")
  if geometria == "185 mm":
    e_c_Au = df_comp_Au["EFIGAMMA185"].to_numpy(dtype="float64")*df_comp_Au["COI GAMMA185"].to_numpy(dtype="float64")
  Q0_c_Au = df_comp_Au["Q0"].to_numpy(dtype="float64")
  Cn_c_Au = df_comp_Au["Net Peak Area"].to_numpy(dtype="float64")
  w_c_Au = w_Au
  lam_c_Au = np.log(2)/df_comp_Au["t(1/2) s"].to_numpy(dtype="float64")
  Er_c_Au = df_comp_Au["EREF"].to_numpy(dtype="float64")

  Aesp_c_Au = Aesp(np.float64(Cn_c_Au[0]), w_c_Au, np.float64(lam_c_Au), tr_c_Au, td_c_Au, ti_c_Au, tv_c_Au, np.float64(e_c_Au[0]))
  Q0_alfa_c_Au = cal_Q0_alfa_i(Q0_c_Au[0],Er_c_Au[0],alfa)
  
  # muestra
  #k0_i, e_i, Q0_i, Cn_i, w_i, lamb_i, Er_i, td_i, tr_i, ti_i, tv_i = par_ele
  k0_i = df_muestra["K0"].to_numpy(dtype="float64")
  if geometria == "50 mm":
    e_i = df_muestra["EFIGAMMA50"].to_numpy(dtype="float64")*df_muestra["COI ROSSBACH"].to_numpy(dtype="float64")
    #e_i = df_muestra["efe"].to_numpy(dtype="float64")
  if geometria == "185 mm":
    e_i = df_muestra["EFIGAMMA185"].to_numpy(dtype="float64")*df_muestra["COI GAMMA185"].to_numpy(dtype="float64")
  Q0_i = df_muestra["Q0"].to_numpy(dtype="float64")
  Cn_i = df_muestra["Net Peak Area"].to_numpy(dtype="float64")
  w_i = w
  lam_i = np.log(2)/df_muestra["t(1/2) s"].to_numpy(dtype="float64")
  Er_i = df_muestra["EREF"].to_numpy(dtype="float64")

  Aesp_i = np.zeros(len(k0_i))
  Q0_alfa_i = np.zeros(len(k0_i))
  Cn_corr_i = np.zeros(len(k0_i))
  for i in range(len(k0_i)):
    Cn_corr_i[i] = corr_Cn(i, df_muestra)
    Aesp_i[i] = Aesp(np.float64(Cn_corr_i[i]), w_i, np.float64(lam_i[i]), tr_i, td_i, ti_i, tv_i, np.float64(e_i[i]))
    #Aesp_i[i] = Aesp(np.float64(Cn_i[i]), w_i, np.float64(lam_i[i]), tr_i, td_i, ti_i, tv_i, np.float64(e_i[i]))
    Q0_alfa_i[i] = cal_Q0_alfa_i(Q0_i[i],Er_i[i],alfa)

  C = np.zeros(len(k0_i))
  for i in range(len(k0_i)):
    # Calcula la concentración del elemento i en la muestra
    C[i] = (Aesp_i[i]/Aesp_c_Au)*(k0_c_Au/k0_i[i])*(e_c_Au/e_i[i])*((f + Q0_alfa_c_Au)/(f + Q0_alfa_i[i]))
  return C, Cn_corr_i

# ------------------------ Calculo de Incertidumbre ---------------------------#
def parametros_cal_U(i,df_muestra,u_e,u_k0,u_w,td_i,ti_i,tr_i,tv_i,w_i,  df_comp, df_comp_Au,u_w_c_Au,td_c_Au,ti_c_Au,tr_c_Au,tv_c_Au,w_c_Au, geom,alfa ):
  # i es indice el nucleido.
  df_unico = df_muestra.iloc[[i]]
  # ------------------------------------------------------------------------
  #alfa = 0.226
  u_alfa = 0
  # ----------------------- Valores de la muestra --------------------------#
  Cn_i = np.float64(df_unico["Net Peak Area Corr"])
  Er_i = np.float64(df_unico["EREF"]) 
  Q0_i = np.float64(df_unico["Q0"]) 
  if (geom == "50 mm"):
    e_i = np.float64(df_unico["EFIGAMMA50"])*np.float64(df_unico["COI ROSSBACH"]) 
  if (geom == "185 mm"):
    e_i = np.float64(df_unico["EFIGAMMA185"])*np.float64(df_unico["COI GAMMA185"])
  k0_i = np.float64(df_unico["K0"])
  lamb_i = np.log(2)/np.float64(df_unico["t(1/2) s"])

  # ----------------------- valores de los comparadores ---------------------------#
  Cn_c = np.float64(df_comp["Cn"]) 
  Er_c = np.float64(df_comp["Er"]) 
  Q0_c = np.float64(df_comp["Q0"]) 
  e_c = np.float64(df_comp["efe"])
  k0_c = np.float64(df_comp["k0"])
  lamb_c = np.float64(df_comp["lambda"])
  td_c = np.float64(df_comp["t_dec"])
  ti_c = np.float64(df_comp["t_irr"])
  tr_c = np.float64(df_comp["t_real"])
  tv_c = np.float64(df_comp["t_vivo"])
  w_c = np.float64(df_comp["w"])
  
  Cn_1 = Cn_c[0]
  Cn_2 = Cn_c[1]
  Er_1 = Er_c[0]
  Er_2 = Er_c[1]
  Q0_1 = Q0_c[0]
  Q0_2 = Q0_c[1]
  e_1 = e_c[0]
  e_2 = e_c[1]
  k0_1 = k0_c[0]
  k0_2 = k0_c[1]
  lamb_1 = lamb_c[0]
  lamb_2 = lamb_c[1]
  td_1 = td_c[0]
  td_2 = td_c[1]
  ti_1 = ti_c[0]
  ti_2 = ti_c[1]
  tr_1 = tr_c[0]
  tr_2 = tr_c[1]
  tv_1 = tv_c[0]
  tv_2 = tv_c[1]
  w_1 = w_c[0]
  w_2 = w_c[1]
 
  # ----------------------- valores del comparador de Au ---------------------------#
  Cn_c_Au = np.float64(df_comp_Au["Net Peak Area"]) 
  Er_c_Au = np.float64(df_comp_Au["EREF"]) 
  Q0_c_Au = np.float64(df_comp_Au["Q0"]) 
  if (geom == "50 mm"):
    e_c_Au = np.float64(df_comp_Au["EFIGAMMA50"])*np.float64(df_comp_Au["COI ROSSBACH"]) 
  if (geom == "185 mm"):
    e_c_Au = np.float64(df_comp_Au["EFIGAMMA185"])*np.float64(df_comp_Au["COI GAMMA185"])
  k0_c_Au = np.float64(df_comp_Au["K0"])
  lamb_c_Au = np.log(2)/np.float64(df_comp_Au["t(1/2) s"])


  # -------------------- Incertidumbre de muestra -------------------------------#
  #u_e = 3 # se ingresa
  #u_k0 = 2.8 # se ingresa 
  #u_w = 0.01 # se ingresa
  u_Cn_v = np.float64(df_unico["Net Peak Uncert"])  # valor
  u_Cn = 100*u_Cn_v/Cn_i # %
  u_Er = 0.0
  u_Q0 = 0.0
  u_lamb = 0.0
  u_td = 0.0
  u_ti = 0.0
  u_tr = 0.0
  u_tv = 0.0
  # --------------- Incertidumbre de los comparadores ---------------------------#
  u_Cn_1 = 0.0
  u_Cn_2 = 0.0
  u_Er_1 = 0.0
  u_Er_2 = 0.0
  u_Q0_1 = 0.0
  u_Q0_2 = 0.0
  u_alfa = 0.0
  u_e_1 = 0.0
  u_e_2 = 0.0
  u_k0_1 = 0.0
  u_k0_2 = 0.0
  u_lamb_1 = 0.0
  u_lamb_2 = 0.0
  u_td_1 = 0.0
  u_td_2 = 0.0
  u_ti_1 = 0.0
  u_ti_2 = 0.0
  u_tr_1 = 0.0
  u_tr_2 = 0.0
  u_tv_1 = 0.0
  u_tv_2 = 0.0
  u_w_1 = 0.0
  u_w_2 = 0.0
# --------------- Incertidumbres del compardor de Au---------------------------#

  u_k0_c_Au = 0.0     #
  u_e_c_Au = 0.0    #
  u_Q0_c_Au = 0.0    #
  u_Cn_c_Au_v =np.float64(df_comp_Au["Net Peak Uncert"]) 
  u_Cn_c_Au = 100*u_Cn_c_Au_v/Cn_c_Au     # Area del fotopico extraer de archivos PLA, RPT
  #u_w_c_Au = 1 # se ingresa
  u_lamb_c_Au = 0.0   #
  u_Er_c_Au = 0.0   #
  u_td_c_Au = 0.0   # Tiempo de decaimiento
  u_tr_c_Au = 0.0            # Tiempo de real
  u_ti_c_Au = 0.0      # Tiempo de irradiación
  u_tv_c_Au = 0.0   # Tiempo de vivo
    

  Val_ini = (Cn_i, Cn_1, Cn_2, Cn_c_Au, Er_i, Er_1, Er_2, Er_c_Au, Q0_i, Q0_1, Q0_2,
             Q0_c_Au, alfa, e_i, e_1, e_2, e_c_Au, k0_i, k0_1, k0_2, k0_c_Au, lamb_i,
             lamb_1, lamb_2, lamb_c_Au, td_i, td_1, td_2, td_c_Au, ti_i, ti_1, ti_2,
             ti_c_Au, tr_i, tr_1, tr_2, tr_c_Au, tv_i, tv_1, tv_2, tv_c_Au, w_i, w_1,
             w_2, w_c_Au)
  u_v_ini = (Cn_i*u_Cn/100, Cn_1*u_Cn_1/100, Cn_2*u_Cn_2/100,
             Cn_c_Au*u_Cn_c_Au/100, Er_i*u_Er/100, Er_1*u_Er_1/100,
             Er_2*u_Er_2/100, Er_c_Au*u_Er_c_Au/100, Q0_i*u_Q0/100,
             Q0_1*u_Q0_1/100, Q0_2*u_Q0_2/100, Q0_c_Au*u_Q0_c_Au/100,
             alfa*u_alfa/100, e_i*u_e/100, e_1*u_e_1/100, e_2*u_e_2/100,
             e_c_Au*u_e_c_Au/100, k0_i*u_k0/100, k0_1*u_k0_1/100,
             k0_2*u_k0_2/100, k0_c_Au*u_k0_c_Au/100, lamb_i*u_lamb/100,
             lamb_1*u_lamb_1/100, lamb_2*u_lamb_2/100, lamb_c_Au*u_lamb_c_Au/100,
             td_i*u_td/100, td_1*u_td_1/100, td_2*u_td_2/100,
             td_c_Au*u_tv_c_Au/100, ti_i*u_ti/100, ti_1*u_ti_1/100,
             ti_2*u_ti_2/100, ti_c_Au*u_ti_c_Au/100, tr_i*u_tr/100,
             tr_1*u_tr_1/100, tr_2*u_tr_2/100, tr_c_Au*u_tr_c_Au/100,
             tv_i*u_tv/100, tv_1*u_tv_1/100, tv_2*u_tv_2/100,
             tv_c_Au*u_tv_c_Au/100, w_i*u_w/100, w_1*u_w_1/100, w_2*u_w_2/100,
             w_c_Au*u_w_c_Au/100)

  return Val_ini, u_v_ini



def cal_U(Val_ini,u_v_ini):
  (Cn, Cn_1, Cn_2, Cn_c_Au, Er, Er_1, Er_2, Er_c_Au, Q0, Q0_1, Q0_2, Q0_c_Au,
   alfa, e, e_1, e_2, e_c_Au, k0, k0_1, k0_2, k0_c_Au, lamb, lamb_1, lamb_2,
   lamb_c_Au, td, td_1, td_2, td_c_Au, ti, ti_1, ti_2, ti_c_Au, tr, tr_1, tr_2,
   tr_c_Au, tv, tv_1, tv_2, tv_c_Au, w, w_1, w_2, w_c_Au) = Val_ini
  (u_Cn, u_Cn_1, u_Cn_2,u_Cn_c_Au, u_Er, u_Er_1, u_Er_2, u_Er_c_Au, u_Q0, u_Q0_1,
   u_Q0_2, u_Q0_c_Au, u_alfa, u_e, u_e_1, u_e_2, u_e_c_Au, u_k0, u_k0_1, u_k0_2,
   u_k0_c_Au, u_lamb, u_lamb_1, u_lamb_2, u_lamb_c_Au, u_td, u_td_1, u_td_2,
   u_td_c_Au, u_ti, u_ti_1, u_ti_2, u_ti_c_Au, u_tr, u_tr_1, u_tr_2, u_tr_c_Au,
   u_tv, u_tv_1, u_tv_2, u_tv_c_Au, u_w, u_w_1, u_w_2, u_w_c_Au) = u_v_ini

  # Aesp
  # [Cn, lamb, td, ti, tr, tv, w]
  Val_ini_Aesp = (Cn, lamb, td, ti, tr, tv, w)
  u_v_ini_Aesp = (u_Cn, u_lamb, u_td, u_ti, u_tr, u_tv, u_w)
  u_Aesp, Aesp  = cal_U_Aesp(Val_ini_Aesp,u_v_ini_Aesp)
  
  # Aesp_1
  # [Cn, lamb, td, ti, tr, tv, w]
  Val_ini_Aesp_1 = (Cn_1, lamb_1, td_1, ti_1, tr_1, tv_1, w_1)
  u_v_ini_Aesp_1 = (u_Cn_1, u_lamb_1, u_td_1, u_ti_1, u_tr_1, u_tv_1, u_w_1)
  u_Aesp_1, Aesp_1  = cal_U_Aesp(Val_ini_Aesp_1,u_v_ini_Aesp_1)


  # Aesp_2
  # [Cn, lamb, td, ti, tr, tv, w]
  Val_ini_Aesp_2 = (Cn_2, lamb_2, td_2, ti_2, tr_2, tv_2, w_2)
  u_v_ini_Aesp_2 = (u_Cn_2, u_lamb_2, u_td_2, u_ti_2, u_tr_2, u_tv_2, u_w_2)
  u_Aesp_2, Aesp_2  = cal_U_Aesp(Val_ini_Aesp_2,u_v_ini_Aesp_2)


  # Aesp_c_Au
  # [Cn_c_Au, lamb_c_Au, td_c_Au, ti_c_Au, tr_c_Au, tv_c_Au, w_c_Au]
  Val_ini_Aesp_c_Au = (Cn_c_Au, lamb_c_Au, td_c_Au, ti_c_Au, tr_c_Au, tv_c_Au,
                       w_c_Au)
  u_v_ini_Aesp_c_Au = (u_Cn_c_Au, u_lamb_c_Au, u_td_c_Au, u_ti_c_Au, u_tr_c_Au,
                       u_tv_c_Au, u_w_c_Au)
  u_Aesp_c_Au, Aesp_c_Au  = cal_U_Aesp(Val_ini_Aesp_c_Au,u_v_ini_Aesp_c_Au)


  # [Aesp, Aesp_1, Aesp_2, Aesp_c_Au, Er, Er_1, Er_2, Er_c_Au, Q0, Q0_1, Q0_2, Q0_c_Au, alpha, e, e_1, e_2, e_c_Au, k0, k0_1, k0_2, k0_c_Au]
  Val_ini_con = (Aesp, Aesp_1, Aesp_2, Aesp_c_Au, Er, Er_1, Er_2, Er_c_Au, Q0,
                 Q0_1, Q0_2, Q0_c_Au, alfa, e, e_1, e_2, e_c_Au, k0, k0_1, k0_2,
                 k0_c_Au)
  u_v_ini_con = (u_Aesp, u_Aesp_1, u_Aesp_2, u_Aesp_c_Au, u_Er, u_Er_1, u_Er_2,
                 u_Er_c_Au, u_Q0, u_Q0_1, u_Q0_2, u_Q0_c_Au, u_alfa, u_e, u_e_1,
                 u_e_2, u_e_c_Au, u_k0, u_k0_1, u_k0_2, u_k0_c_Au)
  # calcula incertidumbre
  #[Aesp, Aesp_1, Q0_alfa, Q0_alfa_1, e, e_1, f, k0, k0_1]
  formula_str = "(Aesp/Aesp_c_Au)*(k0_c_Au/k0)*(e_c_Au/e)*(((k0_1/k0_2)*(e_1/e_2)*((Q0_1 -0.429)/((Er_1)**alfa)+0.429/((2*alfa-1)*0.55**alfa))-(Aesp_1/Aesp_2)*((Q0_2 -0.429)/((Er_2)**alfa)+0.429/((2*alfa-1)*0.55**alfa)))/((Aesp_1/Aesp_2)-(k0_1/k0_2)*(e_1/e_2))+((Q0_c_Au -0.429)/((Er_c_Au)**alfa)+0.429/((2*alfa-1)*0.55**alfa))) / (((k0_1/k0_2)*(e_1/e_2)*((Q0_1 -0.429)/((Er_1)**alfa)+0.429/((2*alfa-1)*0.55**alfa))-(Aesp_1/Aesp_2)*((Q0_2 -0.429)/((Er_2)**alfa)+0.429/((2*alfa-1)*0.55**alfa)))/((Aesp_1/Aesp_2)-(k0_1/k0_2)*(e_1/e_2))+((Q0 -0.429)/((Er)**alfa)+0.429/((2*alfa-1)*0.55**alfa)))"
  derivadas = cal_derivadas(Val_ini_con)

  #print(Val_ini)
  # Extraer variables únicas
  try:
      variables = sorted(list(sp.sympify(formula_str).free_symbols), key=lambda x: str(x))
  except Exception as e:
      st.error(f"Error al interpretar la fórmula: {e}")
  # Entrada de valores e incertidumbres
  valores = {}
  incertidumbres = {}
  i = 0
  for var in variables:
    valor =  Val_ini_con[i]
    incertidumbre = u_v_ini_con[i]
    valores[str(var)] = valor
    incertidumbres[str(var)] = incertidumbre
    i = i + 1

  # Cálculo
  try:
    # Definir símbolos
    simbolos = {str(v): sp.Symbol(str(v)) for v in variables}
    formula_sym = sp.sympify(formula_str)
    # Calcular valor central
    #y_val = (Aesp/Aesp_c_Au)*(k0_c_Au/k0)*(e_c_Au/e)*(((k0_1/k0_2)*(e_1/e_2)*((Q0_1 -0.429)/((Er_1)**alfa)+0.429/((2*alfa-1)*0.55**alfa))-(Aesp_1/Aesp_2)*((Q0_2 -0.429)/((Er_2)**alfa)+0.429/((2*alfa-1)*0.55**alfa)))/((Aesp_1/Aesp_2)-(k0_1/k0_2)*(e_1/e_2))+((Q0_c_Au -0.429)/((Er_c_Au)**alfa)+0.429/((2*alfa-1)*0.55**alfa))) / (((k0_1/k0_2)*(e_1/e_2)*((Q0_1 -0.429)/((Er_1)**alfa)+0.429/((2*alfa-1)*0.55**alfa))-(Aesp_1/Aesp_2)*((Q0_2 -0.429)/((Er_2)**alfa)+0.429/((2*alfa-1)*0.55**alfa)))/((Aesp_1/Aesp_2)-(k0_1/k0_2)*(e_1/e_2))+((Q0 -0.429)/((Er)**alfa)+0.429/((2*alfa-1)*0.55**alfa)))"
    y_val = np.float64(formula_sym.evalf(subs=valores))
    u_y_squared = 0
    contribuciones = []
    kkk = 0
    for v in variables:
        var_name = str(v)
        sensibilidad = derivadas[kkk]
        kkk = kkk + 1
        u_i = incertidumbres[var_name]
        contrib = (np.float64(sensibilidad) * u_i)**2
        u_rel_i = u_i / valores[var_name] if valores[var_name] != 0 else np.nan
        u_y_squared += contrib
    if isinstance(u_y_squared, np.ndarray):
      u_y_squared = u_y_squared.item()
    if isinstance(u_y_squared, list):
      u_y_squared = u_y_squared[0]

    u_y = np.sqrt(u_y_squared) # incertidumbre combinada
    u_rel_y = u_y / y_val if y_val != 0 else np.nan
    u_y_por = 100*u_rel_y
    # Calcular porcentaje de contribución
    for c in contribuciones:
        c["% Contribución"] = 100 * c["Contribución a u(y)²"] / u_y_squared if u_y_squared > 0 else np.nan
  except Exception as e:
      st.error(f"Ocurrió un error en el cálculo: {e}")
  return u_y, y_val, u_y_por, u_Aesp_c_Au

def cal_U_Aesp(Val_ini,u_v_ini):
  # [Cn, lamb, t_d, ti, tr, tv, w]
  # calcula incertidumbre
  formula_str = "(Cn*exp(lamb*td)*lamb*tr)/((1-exp(-lamb*ti))*(1-exp(-lamb*tr))*w*tv)"
  # Extraer variables únicas
  try:
      variables = sorted(list(sp.sympify(formula_str).free_symbols), key=lambda x: str(x))
  except Exception as e:
      st.error(f"Error al interpretar la fórmula: {e}")
  # Entrada de valores e incertidumbres
  valores = {}
  incertidumbres = {}
  i = 0
  for var in variables:
    valor = Val_ini[i]
    incertidumbre = u_v_ini[i]
    valores[str(var)] = valor
    incertidumbres[str(var)] = incertidumbre
    i = i + 1
  # Cálculo
  try:
    # Definir símbolos
    simbolos = {str(v): sp.Symbol(str(v)) for v in variables}
    formula_sym = sp.sympify(formula_str)
    # Calcular valor central
    y_val = np.float64(formula_sym.evalf(subs=valores))
    # Derivadas parciales (sensibilidades)
    u_y_squared = 0
    contribuciones = []
    i = 0
    for v in variables:
        var_name = str(v)
        derivada = formula_sym.diff(v)

        sensibilidad = np.float64(derivada.evalf(subs=valores))
        u_i = incertidumbres[var_name]
        contrib = (sensibilidad * u_i)**2
        u_rel_i = u_i / valores[var_name] if valores[var_name] != 0 else np.nan

        contribuciones.append({
            "Variable": var_name,
            "Sensibilidad ∂y/∂x": sensibilidad,
            "Incertidumbre": u_i,
            "Incertidumbre relativa": u_rel_i,
            "Contribución a u(y)²": contrib,
        })
        u_y_squared += contrib
        i = i + 1
    if isinstance(u_y_squared, np.ndarray):
      u_y_squared = u_y_squared.item()
    if isinstance(u_y_squared, list):
      u_y_squared = u_y_squared[0]
    u_y = np.sqrt(u_y_squared) # incertidumbre combinada
    u_rel_y = u_y / y_val if y_val != 0 else np.nan

    # Calcular porcentaje de contribución
    for c in contribuciones:
        c["% Contribución"] = 100 * c["Contribución a u(y)²"] / u_y_squared if u_y_squared > 0 else np.nan
  except Exception as e:
      st.error(f"Ocurrió un error en el cálculo: {e}")
  return u_y, y_val

# -------------------------- Calculo de derivadas -----------------------------#

def cal_derivadas(Val_ini_con):
  (Aesp, Aesp_1, Aesp_2, Aesp_c_Au, Er, Er_1, Er_2, Er_c_Au, Q0, Q0_1, Q0_2,
   Q0_c_Au, alfa, e, e_1, e_2, e_c_Au, k0, k0_1, k0_2, k0_c_Au) = Val_ini_con

  d_Aesp = e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_Aesp_1 = Aesp*e_c_Au*k0_c_Au*(-((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - (-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*e_c_Au*k0_c_Au*(((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + (-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_Aesp_2 = Aesp*e_c_Au*k0_c_Au*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - Aesp_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2) + Aesp*e_c_Au*k0_c_Au*(Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + Aesp_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_Aesp_c_Au = -Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au**2*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_Er = Aesp*alfa*e_c_Au*k0_c_Au*(Q0 - 0.429)*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*Er*Er**alfa*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_Er_1 = -Aesp*alfa*e_1*e_c_Au*k0_1*k0_c_Au*(Q0_1 - 0.429)/(Aesp_c_Au*Er_1*Er_1**alfa*e*e_2*k0*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*alfa*e_1*e_c_Au*k0_1*k0_c_Au*(Q0_1 - 0.429)*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*Er_1*Er_1**alfa*e*e_2*k0*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_Er_2 = Aesp*Aesp_1*alfa*e_c_Au*k0_c_Au*(Q0_2 - 0.429)/(Aesp_2*Aesp_c_Au*Er_2*Er_2**alfa*e*k0*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) - Aesp*Aesp_1*alfa*e_c_Au*k0_c_Au*(Q0_2 - 0.429)*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2*Aesp_c_Au*Er_2*Er_2**alfa*e*k0*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_Er_c_Au = -Aesp*alfa*e_c_Au*k0_c_Au*(Q0_c_Au - 0.429)/(Aesp_c_Au*Er_c_Au*Er_c_Au**alfa*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_Q0 = -Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*Er**alfa*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_Q0_1 = Aesp*e_1*e_c_Au*k0_1*k0_c_Au/(Aesp_c_Au*Er_1**alfa*e*e_2*k0*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) - Aesp*e_1*e_c_Au*k0_1*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*Er_1**alfa*e*e_2*k0*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_Q0_2 = -Aesp*Aesp_1*e_c_Au*k0_c_Au/(Aesp_2*Aesp_c_Au*Er_2**alfa*e*k0*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*Aesp_1*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2*Aesp_c_Au*Er_2**alfa*e*k0*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_Q0_c_Au = Aesp*e_c_Au*k0_c_Au/(Aesp_c_Au*Er_c_Au**alfa*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_alfa =Aesp*e_c_Au*k0_c_Au*((-Aesp_1*(-(Q0_2 - 0.429)*np.log(Er_2)/Er_2**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/Aesp_2 + e_1*k0_1*(-(Q0_1 - 0.429)*np.log(Er_1)/Er_1**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) - (Q0_c_Au - 0.429)*np.log(Er_c_Au)/Er_c_Au**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))*(-(-Aesp_1*(-(Q0_2 - 0.429)*np.log(Er_2)/Er_2**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/Aesp_2 + e_1*k0_1*(-(Q0_1 - 0.429)*np.log(Er_1)/Er_1**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)*np.log(Er)/Er**alfa - 0.256472073324161/(0.55**alfa*(2*alfa - 1)) + 0.858/(0.55**alfa*(2*alfa - 1)**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_e = -Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e**2*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_e_1 = Aesp*e_c_Au*k0_c_Au*(-k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2) + Aesp*e_c_Au*k0_c_Au*(k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_e_2 = Aesp*e_c_Au*k0_c_Au*(-e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2**2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - e_1*k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2**2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*e_c_Au*k0_c_Au*(e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2**2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + e_1*k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2**2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_e_c_Au = Aesp*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_k0 = -Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0**2*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_k0_1 = Aesp*e_c_Au*k0_c_Au*(-e_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - e_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2) + Aesp*e_c_Au*k0_c_Au*(e_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + e_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))

  d_k0_2 = Aesp*e_c_Au*k0_c_Au*(-e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - e_1*k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*e_c_Au*k0_c_Au*(e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + e_1*k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2)

  d_k0_c_Au = Aesp*e_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1))))
  (e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  Aesp*e_c_Au*k0_c_Au*(-((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - (-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*e_c_Au*k0_c_Au*(((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + (-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  Aesp*e_c_Au*k0_c_Au*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - Aesp_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2) + Aesp*e_c_Au*k0_c_Au*(Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + Aesp_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  -Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au**2*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  Aesp*alfa*e_c_Au*k0_c_Au*(Q0 - 0.429)*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*Er*Er**alfa*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  -Aesp*alfa*e_1*e_c_Au*k0_1*k0_c_Au*(Q0_1 - 0.429)/(Aesp_c_Au*Er_1*Er_1**alfa*e*e_2*k0*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*alfa*e_1*e_c_Au*k0_1*k0_c_Au*(Q0_1 - 0.429)*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*Er_1*Er_1**alfa*e*e_2*k0*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  Aesp*Aesp_1*alfa*e_c_Au*k0_c_Au*(Q0_2 - 0.429)/(Aesp_2*Aesp_c_Au*Er_2*Er_2**alfa*e*k0*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) - Aesp*Aesp_1*alfa*e_c_Au*k0_c_Au*(Q0_2 - 0.429)*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2*Aesp_c_Au*Er_2*Er_2**alfa*e*k0*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  -Aesp*alfa*e_c_Au*k0_c_Au*(Q0_c_Au - 0.429)/(Aesp_c_Au*Er_c_Au*Er_c_Au**alfa*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  -Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*Er**alfa*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  Aesp*e_1*e_c_Au*k0_1*k0_c_Au/(Aesp_c_Au*Er_1**alfa*e*e_2*k0*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) - Aesp*e_1*e_c_Au*k0_1*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*Er_1**alfa*e*e_2*k0*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  -Aesp*Aesp_1*e_c_Au*k0_c_Au/(Aesp_2*Aesp_c_Au*Er_2**alfa*e*k0*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*Aesp_1*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_2*Aesp_c_Au*Er_2**alfa*e*k0*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  Aesp*e_c_Au*k0_c_Au/(Aesp_c_Au*Er_c_Au**alfa*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  Aesp*e_c_Au*k0_c_Au*((-Aesp_1*(-(Q0_2 - 0.429)*np.log(Er_2)/Er_2**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/Aesp_2 + e_1*k0_1*(-(Q0_1 - 0.429)*np.log(Er_1)/Er_1**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) - (Q0_c_Au - 0.429)*np.log(Er_c_Au)/Er_c_Au**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))*(-(-Aesp_1*(-(Q0_2 - 0.429)*np.log(Er_2)/Er_2**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/Aesp_2 + e_1*k0_1*(-(Q0_1 - 0.429)*np.log(Er_1)/Er_1**alfa + 0.256472073324161/(0.55**alfa*(2*alfa - 1)) - 0.858/(0.55**alfa*(2*alfa - 1)**2))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)*np.log(Er)/Er**alfa - 0.256472073324161/(0.55**alfa*(2*alfa - 1)) + 0.858/(0.55**alfa*(2*alfa - 1)**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  -Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e**2*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  Aesp*e_c_Au*k0_c_Au*(-k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2) + Aesp*e_c_Au*k0_c_Au*(k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  Aesp*e_c_Au*k0_c_Au*(-e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2**2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - e_1*k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2**2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*e_c_Au*k0_c_Au*(e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2**2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + e_1*k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2**2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  Aesp*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  -Aesp*e_c_Au*k0_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0**2*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  Aesp*e_c_Au*k0_c_Au*(-e_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - e_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2) + Aesp*e_c_Au*k0_c_Au*(e_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + e_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))),
  Aesp*e_c_Au*k0_c_Au*(-e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) - e_1*k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))) + Aesp*e_c_Au*k0_c_Au*(e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))) + e_1*k0_1*(-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(e_2*k0_2**2*(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2))**2))*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))**2),
  Aesp*e_c_Au*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0_c_Au - 0.429)/Er_c_Au**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(Aesp_c_Au*e*k0*((-Aesp_1*((Q0_2 - 0.429)/Er_2**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/Aesp_2 + e_1*k0_1*((Q0_1 - 0.429)/Er_1**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))/(e_2*k0_2))/(Aesp_1/Aesp_2 - e_1*k0_1/(e_2*k0_2)) + (Q0 - 0.429)/Er**alfa + 0.429/(0.55**alfa*(2*alfa - 1)))))





  derivadas = (d_Aesp, d_Aesp_1, d_Aesp_2, d_Aesp_c_Au, d_Er, d_Er_1, d_Er_2, d_Er_c_Au, d_Q0, d_Q0_1, d_Q0_2, d_Q0_c_Au, d_alfa, d_e, d_e_1, d_e_2, d_e_c_Au, d_k0, d_k0_1, d_k0_2, d_k0_c_Au)

  return derivadas

# ---------------- Correción por picos interferentes ---------------#

def corr_Cn(i, df_muestra):
    # i ubicación
    # df_final: todos los datos
    df_muestra["Net Peak Area"] = pd.to_numeric(df_muestra["Net Peak Area"], errors="coerce")
    df_muestra["E_INTERF"] = pd.to_numeric(df_muestra["E_INTERF"], errors="coerce")
    df_muestra["FC_GAMM"] = pd.to_numeric(df_muestra["FC_GAMM"], errors="coerce")
    
    delta = 1.0
    df_unico = df_muestra.iloc[i]
    Nucl = df_unico["Identidad_Verificada_Energia"]
    Area = df_unico["Net Peak Area"]
    Interf = df_unico["INTERF"]
    E_Interf =np.float64(df_unico["E_INTERF"])
    FC = np.float64(df_unico["FC_GAMM"])

    if (Interf == "N_A"):
      return Area
    df_filtrado = df_muestra[(df_muestra["Identidad_Verificada_Energia"] == Interf) & (df_muestra["EGKEV"].between(E_Interf - delta, E_Interf + delta))]
    if df_filtrado.empty:
      return Area

    E_in_conf = df_filtrado.iloc[0]["Net Peak Area"]
    Area = Area - E_in_conf*FC

    return Area

# ---------------- Redondeo para concentracióne e incertidumbre ---------------#

def redondear_con_incert(x, u, sig_inc):
    #x: valor nominal
    #porc_u: porcentaje de incertidumbre (ej. 3 = 3%)
    #sig_inc: cifras significativas para la incertidumbre (1 o 2)
    
    u_red = float(f"{u:.{sig_inc}g}")
   
    if u_red <= 0 or np.isnan(u_red):
      st.error(f"La incertidumbre no es válida (u_red = {u_red}).")
      x_red = u_red
    else:
      orden = int(np.floor(np.log10(abs(u_red))))
      #orden = int(np.floor(np.log10(u_red)))
      decimales = -orden
      x_red = round(x, decimales)

    return x_red, u_red
