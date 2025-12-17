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
# ------------------ RPT ---------------------------------

def limpiar(valor):
    if isinstance(valor, str):
        # elimina todos los espacios vacíos al inicio y al final
        valor = valor.strip()
    return valor

def procesar_RPT(rpt_file):
    if rpt_file is None:
        return None

    # 1. Leer directo desde uploader sin guardarlo
    contenido = rpt_file.read().decode("utf-8", errors="ignore")

    # 2. Convertir a pandas como una fila por línea
    lineas = contenido.splitlines()
    df = pd.DataFrame(lineas, columns=["linea"])

    # 3. Eliminar primeras 19 filas
    df = df.iloc[19:].reset_index(drop=True)

    # 4. Crear columna auxiliar sin espacios
    df["sin_espacios"] = df["linea"].str.lstrip()

    # 5. Filtrar patrones a eliminar
    patrones_excluir = ("Peak", "M =", "m =", "F =", "Error")
    df = df[~df["sin_espacios"].str.startswith(patrones_excluir)]

    # 6. Eliminar líneas que empiezan con exactamente 12 espacios
    df = df[~df["linea"].str.startswith("            ")]

    # 7. Eliminar líneas vacías o solo con espacios
    df = df[~df["linea"].str.strip().eq("")]

    # 8. Eliminar columna auxiliar
    df = df.drop(columns=["sin_espacios"])

    # 9. Reiniciar índices
    df = df.reset_index(drop=True)

    # 10. Quitar 1 espacio inicial (si hay) y todos los finales
    df = df.applymap(limpiar)
    
    # 11. Agregar letra al inicio si la línea no empieza con una letra
    df["linea"] = df["linea"].apply(lambda x: "X" + x if not x[:1].isalpha() else x)

    # 12. Separar columnas y asignar nombres
    
    # primera columna
    df_tipo = df["linea"].str[:1].to_frame()
    df_tipo = df_tipo.applymap(limpiar)
    df_tipo.columns = ["Tipo"]
    
    # demas columnas
    df_demas = df["linea"].str[1:].to_frame()
    df_demas = df_demas.applymap(limpiar)
    df_demas_tab = df_demas["linea"].str.split(r"\s+", n=9, expand=True)
    df_demas_tab.columns = ["Peak No.", "ROI Start", "ROI End", "Peak Centroid",
        "Energy (keV)", "Net Peak Area", "Net Peak Uncert",
        "Continuum Counts", "Tentative Nuclide"
    ]
    
    # 12. Unir df
    df_tab = pd.concat([df_tipo, df_demas_tab], axis=1)
    st.success("Archivo procesado correctamente 🚀")
    return df_tab

def Selecion_Nucleidos_muestra(df_rpt_muestras,df_Nucleidos, df_database, tol=1.5):
    df_rpt_muestras["Energy (keV)"] = pd.to_numeric(df_rpt_muestras["Energy (keV)"], errors="coerce").astype("float64")
    df_Nucleidos["E (keV)"] = pd.to_numeric(df_Nucleidos["E (keV)"], errors="coerce").astype("float64")
    elementos_validos = df_Nucleidos["Elemento"].unique()
    df_filtrado1 = df_rpt_muestras[
        (df_rpt_muestras["Tentative Nuclide"].isin(elementos_validos)) 
       ]
    df_rpt_muestras.reset_index(drop=True, inplace=True)
    df_Nucleidos.reset_index(drop=True, inplace=True)
    df_filtrado1.reset_index(drop=True, inplace=True)

    filas_filtradas = []
    Nucleidos = pd.DataFrame(columns=['Identidad_Verificada_Energia'])  
    
    for _, rango in df_Nucleidos.iterrows():
        e_min = rango['E (keV)'] - tol
        e_max = rango['E (keV)'] + tol
        nucleido = rango['Elemento']
        
        # Filtrar muestras en este rango
        mascara = (df_filtrado1['Energy (keV)'] >= e_min) & (df_filtrado1['Energy (keV)'] <= e_max)
        muestras_en_rango = df_filtrado1[mascara].copy()
        
        if not muestras_en_rango.empty:    
            filas_filtradas.append(muestras_en_rango)
            #Nucleidos = pd.concat([Nucleidos, nucleido], pd.Series(nucleido) ignore_index=True)
            #Nucleidos['Identidad_Verificada_Energia'] = nucleido
            lista_nucleidos = Nucleidos['Identidad_Verificada_Energia'].tolist()
            lista_nucleidos.append(nucleido)
            Nucleidos = pd.DataFrame() 
            Nucleidos['Identidad_Verificada_Energia'] = lista_nucleidos
    
    if not filas_filtradas:
        return pd.DataFrame()
    # Combinar todos los resultados
    df_filtrado = pd.concat(filas_filtradas, ignore_index=True)
    df_filtrado = df_filtrado.join(Nucleidos)

    # agregar propiedades
    df_unido = Extra_from_database(df_filtrado, df_database,tol)
    
    return df_unido

def Selecion_Nucleidos_Au(df_rpt_Au, df_database,tol):
    # buscar en database energía de Au
   
    En_Au = np.float64(411.8)
    Nucledo_comp = "AU-198"
 
    tol_Au = np.float64(1.0)
    df_rpt_Au["Energy (keV)"] = pd.to_numeric(df_rpt_Au["Energy (keV)"], errors="coerce").astype("float64")
    df_energy_Au = df_rpt_Au[(df_rpt_Au["Tentative Nuclide"] == Nucledo_comp ) & ((df_rpt_Au["Energy (keV)"] > En_Au - tol_Au) | (df_rpt_Au["Energy (keV)"] < En_Au + tol_Au))]
    df_energy_Au["Tentative Nuclide"] = "AU-198"
    # agregar propiedades
    df_unido = Extra_from_database(df_energy_Au, df_database,tol)

    return df_unido

def Extra_from_database(df, df_database,tol=1.5):
    df["Energy (keV)"] = pd.to_numeric(df["Energy (keV)"], errors="coerce").astype("float64")
    df_database["EGKEV"] = pd.to_numeric(df_database["EGKEV"], errors="coerce").astype("float64")
    
    df.reset_index(drop=True, inplace=True)
    df_database.reset_index(drop=True, inplace=True)
    df_database_o = df_database.copy()
   
    df_prop_nucleidos=pd.DataFrame()
    for _, rango in df.iterrows():
        e_min = rango['Energy (keV)'] - tol
        e_max = rango['Energy (keV)'] + tol
        nucleido = rango['Tentative Nuclide']
        
        # Filtrar muestras en este rango
        mascara = (df_database_o['NUCLID'] == nucleido) & (df_database_o['EGKEV'] >= e_min) & (df_database_o['EGKEV'] <= e_max)
        muestras_en_rango = df_database_o[mascara].copy()
        df_prop_nucleidos = pd.concat([df_prop_nucleidos, muestras_en_rango], ignore_index=True)
        
    # Agregar propiedades
    df = df.reset_index(drop=True)
    df_prop_nucleidos = df_prop_nucleidos.reset_index(drop=True)

    df_unido = pd.concat([df, df_prop_nucleidos], axis=1)
    
    return df_unido  

# ------------------ kos ---------------------------------

def extraer_DATE_MEA_MEAS_TIM(k0s_file):
    if k0s_file is None:
        return None

    contenido = k0s_file.getvalue().decode("utf-8", errors="ignore")
    lineas = contenido.splitlines()

    fecha = hora = tiempo_real = tiempo_vivo = None

    for i, linea in enumerate(lineas):
        l = linea.strip()

        # -------------------- DATE_MEA --------------------
        if l.startswith("$DATE_MEA"):
            if i + 1 < len(lineas):
                datos = lineas[i + 1].strip().split()
                if len(datos) >= 2:
                    fecha = datos[0]
                    hora = datos[1]

        # -------------------- MEAS_TIM --------------------
        if l.startswith("$MEAS_TIM"):
            if i + 1 < len(lineas):
                datos = lineas[i + 1].strip().split()
                if len(datos) >= 2:
                    tiempo_vivo = datos[0]
                    tiempo_real = datos[1]
                    
    return fecha, hora, tiempo_vivo, tiempo_real

#  ------------------ kos ---------------------------------
# tiempo de decaimiento y tiempo de irradiación

def Proc_Irr_Dec(f_ini, h_ini, f_fin, h_fin, f_med, hora_med, f_ini_c_Au, h_ini_c_Au, f_fin_c_Au, h_fin_c_Au, f_med_c_Au, hora_med_c_Au ): 
    # jojo
    # f_ini, h_ini: fecha y hora de inicio de irradiación 
    # f_fin, h_fin: fecha y hora de fin de irradiación 
    # f_med, hora_med: de inicio de la medición muestra 
    # f_med_c_Au, hora_med_c_Au: Fecha y hora de incio de medición Comparador Au 
    fecha1 = str(f_ini+" "+h_ini) # muestra
    fecha2 = str(f_fin+" "+h_fin) # muestra
    fecha3 = str(str(f_med)+" "+str(hora_med)) # muestra
    fecha4 = str(f_ini_c_Au+" "+h_ini_c_Au) # comparador 
    fecha5 = str(f_fin_c_Au+" "+h_fin_c_Au) # comparador
    fecha6 = str(str(f_med_c_Au)+" "+str(hora_med_c_Au)) # comparador

    
    f1 = datetime.strptime(fecha1, "%m/%d/%Y %H:%M:%S")
    f2 = datetime.strptime(fecha2, "%m/%d/%Y %H:%M:%S")
    f3 = datetime.strptime(fecha3, "%m/%d/%Y %H:%M:%S")
    f4 = datetime.strptime(fecha4, "%m/%d/%Y %H:%M:%S")
    f5 = datetime.strptime(fecha5, "%m/%d/%Y %H:%M:%S")
    f6 = datetime.strptime(fecha6, "%m/%d/%Y %H:%M:%S")


    # Diferencia en segundos

    ti_i = float((f2 - f1).total_seconds())
    td_i = float((f3 - f2).total_seconds())
    
    ti_c_Au = float((f5 - f4).total_seconds())
    td_c_Au = float((f6 - f5).total_seconds())
    
    return ti_i,td_i, ti_c_Au, td_c_Au 
