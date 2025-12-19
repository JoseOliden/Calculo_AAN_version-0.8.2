# app.py
from librerias import *
from ProcArch import *
from calAAN import *
from ProcFechas import *
#https://github.com/JoseOliden/Calculo_AAN_v2/blob/main/LOGO_IPEN.jpg
# Establecer configuración de página
st.set_page_config(
    page_title="Sistema de Análisis k0 - AAN",
    page_icon="🔬",
    layout="wide"
)

# CSS personalizado
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1E3A8A;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.8rem;
        color: #1E3A8A;
        margin-top: 2rem;
        margin-bottom: 1rem;
    }
    .result-box {
        background-color: #F3F4F6;
        padding: 1.5rem;
        border-radius: 10px;
        border-left: 5px solid #3B82F6;
        margin-bottom: 1rem;
    }
    .info-box {
        background-color: #EFF6FF;
        padding: 1rem;
        border-radius: 8px;
        border: 1px solid #93C5FD;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

st.image(
    "banner.jpg",
    width=1500,
)
# Título principal
st.markdown('<h1 class="main-header">Sistema de Análisis AAN</h1>', unsafe_allow_html=True)

# Barra lateral para navegación
st.sidebar.title("🌐 Navegación")
page = st.sidebar.radio(
    "Seleccionar sección:",
    ["📁 Carga de Datos", "⚙️ Configuración", "📊 Procesamiento", "📈 Resultados", "📄 Reporte (EN DESARROLLO)"]
)

# ============================================
# SECCIÓN 1: CARGA DE DATOS
# ============================================
if page == "📁 Carga de Datos":
    st.markdown('<h2 class="section-header">📁 Carga de Archivos</h2>', unsafe_allow_html=True)
    
    # Crear columnas para la carga de archivos
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Archivos de la muestra")
        rpt_file = st.file_uploader("Subir archivo .RPT", type=['rpt', 'RPT'], key="rpt_sample")
        if rpt_file:
            st.session_state["rpt_file"] = rpt_file
            df_resultado = procesar_RPT(st.session_state["rpt_file"])
            st.session_state["df_resultado"] = df_resultado

        if "df_resultado" in st.session_state:
            st.success(f"📄 Archivo cargado: {st.session_state["rpt_file"].name}")
            st.dataframe(st.session_state["df_resultado"])
        else:
            st.warning("⚠️ No se ha cargado archivo RPT ")

        k0s_file = st.file_uploader("Subir archivo .k0s", type=['k0s', 'K0S'], key="k0s_sample")
        if k0s_file:
            st.session_state["k0s_file"] = k0s_file
            fecha, hora, t_vivo, t_real = extraer_DATE_MEA_MEAS_TIM(k0s_file)
            st.session_state["fecha"] = fecha
            st.session_state["hora"] = hora
            st.session_state["t_vivo"] = np.float64(t_vivo)
            st.session_state["t_real"] = np.float64(t_real)

        if (
            "fecha" in st.session_state and
            "hora" in st.session_state and
            "t_vivo" in st.session_state and
            "t_real" in st.session_state
            ):
            st.success(f"📄 Archivo cargado: {st.session_state['k0s_file'].name}")
            st.write("📌 Datos extraídos del archivo")
            st.write(f"**Fecha de medición:** {st.session_state["fecha"]}")
            st.write(f"**Hora de medición:** {st.session_state["hora"]}")
            st.write(f"**Tiempo vivo (s):** {st.session_state["t_vivo"]}")
            st.write(f"**Tiempo real (s):** {st.session_state["t_real"]}")
        else:
            st.warning("⚠️ No se ha cargado archivo k0s ")
    
    with col2:
        
        st.subheader("Archivos del comparador")
        rpt_au_file = st.file_uploader("Subir archivo .RPT", type=['RPT', 'RPT'], key="rpt_au")
        if rpt_au_file:
            st.session_state["rpt_au_file"] = rpt_au_file
            df_au_resultado = procesar_RPT(st.session_state["rpt_au_file"])
            st.session_state["df_au_resultado"] = df_au_resultado

        if "df_au_resultado" in st.session_state:
            st.success(f"📄 Archivo cargado: {st.session_state["rpt_au_file"].name}")
            st.dataframe(st.session_state["df_au_resultado"])
        else:
            st.warning("⚠️ No se ha cargado archivo RPT ")
            
        k0s_au_file = st.file_uploader("Subir archivo .k0s", type=['k0s', 'K0S'], key="k0s_au")
        if k0s_au_file:
            st.session_state["k0s_au_file"] = k0s_au_file
            fecha_au, hora_au, t_vivo_au, t_real_au = extraer_DATE_MEA_MEAS_TIM(k0s_au_file)
            st.session_state["fecha_au"] = fecha_au
            st.session_state["hora_au"] = hora_au
            st.session_state["t_vivo_au"] = np.float64(t_vivo_au)
            st.session_state["t_real_au"] = np.float64(t_real_au)

        if (
            "fecha_au" in st.session_state and
            "hora_au" in st.session_state and
            "t_vivo_au" in st.session_state and
            "t_real_au" in st.session_state
            ):
            st.success(f"📄 Archivo cargado: {st.session_state['k0s_au_file'].name}")
            st.write("📌 Datos extraídos del archivo")
            st.write(f"**Fecha de medición:** {st.session_state["fecha_au"]}")
            st.write(f"**Hora de medición:** {st.session_state["hora_au"]}")
            st.write(f"**Tiempo vivo (s):** {st.session_state["t_vivo_au"]}")
            st.write(f"**Tiempo real (s):** {st.session_state["t_real_au"]}")
        else:
            st.warning("⚠️ No se ha cargado archivo k0s ")
          
    col21, col22 = st.columns(2)

    with col21:
        # Base de datos de Nucléidos
        st.subheader("🗃️ Base de datos de nucléidos")
        db_file = st.file_uploader("Subir Base de Datos (.xlsx)", type=['xlsx'], key="database")
        if db_file:
            df_file = pd.read_excel(db_file)
            st.session_state["df_file"] = df_file
        if "df_file" in st.session_state:
            st.success(f"✅ Base de datos cargada")
            st.dataframe(st.session_state["df_file"])
        else:
            st.warning("⚠️ No se ha cargado archivo")

    
    with col22:
        # Librería de Nucléidos
        st.subheader("📚 Librería de Nucléidos")
        if "ref_type" not in st.session_state:
            st.session_state["ref_type"] = "Corta (C)"
        ref_type = st.radio("Seleccionar tipo de nucléidos:", ["Corta (C)", "Media (M)", "Larga (L)"],index = ["Corta (C)", "Media (M)", "Larga (L)"].index(st.session_state["ref_type"]))
        st.session_state["ref_type"] = ref_type
        ref_files = st.file_uploader(f"Subir archivo RDN_{ref_type[0]}.xlsx", type=['xlsx'], key="reference")
        if ref_files:
            ref_files = pd.read_excel(ref_files)
            st.session_state["ref_files"] = ref_files
        if "ref_files" in st.session_state:
            st.success(f"✅ Archivo cargado")
            st.dataframe(st.session_state["ref_files"])
        else:
            st.warning("⚠️ No se ha cargado archivo")  

# ============================================
# SECCIÓN 2: CONFIGURACIÓN
# ============================================
elif page == "⚙️ Configuración":
    st.markdown('<h2 class="section-header">⚙️ Configuración del Análisis</h2>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.subheader("⚖️ Masas")
        masas, ince = st.columns(2)
        with masas:
            if "masa_muestra" in st.session_state : 
                masa_muestra = st.number_input("Masa de la muestra (g):", min_value=0.0, value = st.session_state["masa_muestra"], step=0.0001, format="%.6f")
            else: 
                masa_muestra = st.number_input("Masa de la muestra (g):", min_value=0.0, value = 0.2817, step=0.0001, format="%.6f")
            st.session_state["masa_muestra"] = np.float64(masa_muestra)
            
            masa_comparador_au = st.number_input("Masa del comparador (μg):", min_value=0.0, value=16.82, step=0.01, format="%.2f")
            st.session_state["masa_comparador_au"] = np.float64(masa_comparador_au)/1000000

        with ince:
            u_w = st.number_input("Incertidumbre masa de la muestra (%):", min_value=0.0, max_value=5.0, value=0.01, step=0.001)
            st.session_state["u_w"] = np.float64(u_w)
            u_w_Au = st.number_input("Incertidumbre masa del comparador (%):", min_value=0.0, max_value=5.0, value=0.01, step=0.01)
            st.session_state["u_w_Au"] = np.float64(u_w_Au)
    
    with col2:
        st.subheader("🕐 Irradiación de la muestra")
        col_fecha1, col_hora1 = st.columns(2)
        with col_fecha1:
            #fecha_ini = st.date_input("Fecha inicio irradiación (yyyy/mm/dd):", value=datetime(2025, 9, 26))
            fecha_ini= st.text_input("Fecha inicio (muestra) (MM/DD/AAAA):", value="09/26/2025")

            st.session_state["fecha_ini"] = fecha_ini
        with col_hora1:
            #hora_ini = st.time_input("Hora inicio irradiación:", value=datetime.strptime("08:45:00", "%H:%M:%S").time(),step=timedelta(seconds=1))
            hora_ini = st.text_input("Hora inicio (muestra) (HH:MM:SS):", value="08:45:00")
            st.session_state["hora_ini"] = hora_ini
        
        col_fecha2, col_hora2 = st.columns(2)
        with col_fecha2:
            #fecha_fin = st.date_input("Fecha fin irradiación (yyyy/mm/dd):", value=datetime(2025, 9, 26))
            fecha_fin= st.text_input("Fecha fin (muestra) (MM/DD/AAAA):", value="09/26/2025")

            st.session_state["fecha_fin"] = fecha_fin
        with col_hora2:
            #hora_fin = st.time_input("Hora fin irradiación:", value=datetime.strptime("09:45:00", "%H:%M:%S").time(),step=timedelta(seconds=1))
            hora_fin= st.text_input("Hora fin (muestra) (HH:MM:SS):", value="09:45:00")
            st.session_state["hora_fin"] = hora_fin

    with col3:
        st.subheader("🕐 Irradiación del comparador")
        col_fecha1, col_hora1 = st.columns(2)
        with col_fecha1:
            #fecha_ini = st.date_input("Fecha inicio irradiación (yyyy/mm/dd):", value=datetime(2025, 9, 26))
            fecha_ini_Au= st.text_input("Fecha inicio (comparador) (MM/DD/AAAA):", value="09/26/2025")
            st.session_state["fecha_ini_Au"] = fecha_ini_Au
        with col_hora1:
            #hora_ini = st.time_input("Hora inicio irradiación:", value=datetime.strptime("08:45:00", "%H:%M:%S").time(),step=timedelta(seconds=1))
            hora_ini_Au = st.text_input("Hora inicio (comparador) (HH:MM:SS):", value="08:45:00")
            st.session_state["hora_ini_Au"] = hora_ini_Au
        
        col_fecha2, col_hora2 = st.columns(2)
        with col_fecha2:
            #fecha_fin = st.date_input("Fecha fin irradiación (yyyy/mm/dd):", value=datetime(2025, 9, 26))
            fecha_fin_Au = st.text_input("Fecha fin (comparador) (MM/DD/AAAA):", value="09/26/2025")

            st.session_state["fecha_fin_Au"] = fecha_fin_Au
        with col_hora2:
            #hora_fin = st.time_input("Hora fin irradiación:", value=datetime.strptime("09:45:00", "%H:%M:%S").time(),step=timedelta(seconds=1))
            hora_fin_Au= st.text_input("Hora fin (comparador) (HH:MM:SS):", value="09:45:00")
            st.session_state["hora_fin_Au"] = hora_fin_Au

    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.subheader("✅ Verificación tiempos ")
        
        #Tiempos de irradiación y decaimiento de la muestra
        # Irraciación: (f_fin, h_fin) - (f_ini, h_ini)
        # Decaimiento: (f_ini, h_ini) -  (f_med, h_med) 

        #Tiempos de irradiación y decaimiento del comparador Au 
        # Se el comparador fue irradiado en un tiempo diferente el cálculo
        # Irraciación: (f_fin_Au, h_fin_Au) - (f_ini_Au, h_ini_Au)
        # Decaimiento: (f_ini_Au, h_ini_Au) -  (f_med_c_Au, hora_med_c_Au)
            
        f_ini = st.session_state["fecha_ini"]
        h_ini = st.session_state["hora_ini"]
        f_fin = st.session_state["fecha_fin"]
        h_fin = st.session_state["hora_fin"]
        f_med = st.session_state["fecha"]
        hora_med = st.session_state["hora"]
        f_ini_c_Au = st.session_state["fecha_ini_Au"]
        h_ini_c_Au = st.session_state["hora_ini_Au"]
        f_fin_c_Au = st.session_state["fecha_fin_Au"]
        h_fin_c_Au = st.session_state["hora_fin_Au"]
        f_med_c_Au = st.session_state["fecha_au"] 
        hora_med_c_Au = st.session_state["hora_au"]
             
        t_irr, t_dec, t_irr_Au, t_dec_Au = Proc_Irr_Dec(f_ini, h_ini, f_fin, h_fin, f_med, hora_med, f_ini_c_Au, h_ini_c_Au, f_fin_c_Au, h_fin_c_Au, f_med_c_Au, hora_med_c_Au)
        st.session_state["t_irr"] = t_irr
        st.session_state["t_dec"] = t_dec
        st.session_state["t_irr_Au"] = t_irr_Au
        st.session_state["t_dec_Au"] = t_dec_Au
        
        if np.float64(t_irr) > 0 :
            st.markdown( f"<span style='color:{"green"}'><b>Tiempo irradiación de la muestra (s):</b> {t_irr}</span>",unsafe_allow_html=True)
        else: 
            st.markdown(f"<span style='color:{"red"}'><b>Tiempo irradiación de la muestra (s): ERROR INGRESO DE DATOS</b></span>", unsafe_allow_html=True)

        if np.float64(t_dec) > 0 :
            st.markdown( f"<span style='color:{"green"}'><b>Tiempo decaimiento de la muestra (s):</b> {t_dec}</span>",unsafe_allow_html=True)
        else: 
            st.markdown(f"<span style='color:{"red"}'><b>Tiempo decaimiento de la muestra (s): ERROR INGRESO DE DATOS</b></span>", unsafe_allow_html=True)

        if np.float64(t_irr_Au) > 0 :
            st.markdown( f"<span style='color:{"green"}'><b>Tiempo irradiación del comparador Au (s):</b> {t_irr_Au}</span>",unsafe_allow_html=True)
        else: 
            st.markdown(f"<span style='color:{"red"}'><b>Tiempo irradiación del comparador Au (s): ERROR INGRESO DE DATOS</b></span>", unsafe_allow_html=True)      

        if np.float64(t_dec_Au) > 0 :
            st.markdown( f"<span style='color:{"green"}'><b>Tiempo decaimiento del comparador Au (s):</b> {t_dec_Au}</span>",unsafe_allow_html=True)
        else: 
            st.markdown(f"<span style='color:{"red"}'><b>Tiempo decaimiento del comparador Au (s): ERROR INGRESO DE DATOS</b></span>", unsafe_allow_html=True)      
        
    
    with col2:
        st.subheader("📐 Geometría")
        geometria = st.radio("Geometría de detección:", ["50 mm", "185 mm"])
        geometria_val = "50" if geometria == "50 mm" else "185"
        st.session_state["geometria"] = geometria

        st.subheader("⏰ Tolerancia de Energía")
        tolerancia = st.slider("Tolerancia de energía (keV):", min_value=0.1, max_value=5.0, value=1.5, step=0.1)
        st.session_state["tolerancia"] = np.float64(tolerancia)

    with col3:
        st.subheader("📊 Incertidumbres")
        u_k0 = 2.8
        u_e = 3.0
        st.write(f"**Incertidumbre k0 de la muestra (%):** {u_k0}")
        st.write(f"**Incertidumbre eficiencia de la muestra (%):** {u_e}")
        st.write(f"**Incertidumbre masa de la muestra (%):** {st.session_state["u_w"]}")
        st.write(f"**Incertidumbre masa de la comparador (%):** {st.session_state["u_w_Au"]}")
       
        st.session_state["u_k0"] = np.float64(u_k0)
        st.session_state["u_e"] = np.float64(u_e)
       
        


# ============================================
# SECCIÓN 3: PROCESAMIENTO
# ============================================
elif page == "📊 Procesamiento":
    st.markdown('<h2 class="section-header">📊 Procesamiento de Datos</h2>', unsafe_allow_html=True)
    
    if st.button("🚀 Iniciar Procesamiento", type="primary", use_container_width=True):
        with st.spinner("Procesando datos..."):
            # Aquí iría la lógica de procesamiento
            # Por ahora mostramos un ejemplo simulado
            
            # Simulación de progreso
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            steps = [
                "Procesando archivos...",
                "Calculando concentraciones...",
                "Calculando incertidumbre...",
                "Generando resultados..."
            ]
        
            for i, step in enumerate(steps):
                progress_bar.progress((i + 1) / len(steps))
                status_text.text(f"📋 {step}")
                if (step == "Procesando archivos..."):
                    st.write("Procesando archivos...")
                    
                    # Comparadores para cálculo de alfa
            
                    df_comparadores_alfa_f = crear_df_comparadores()
                    if "df_comparadores_alfa_f" not in st.session_state:
                        st.session_state["df_comparadores_alfa_f"] = crear_df_comparadores()
                    else:
                        st.session_state["df_comparadores_alfa_f"]
            
                    # Procesa comparador de Au y sus datos
                    
                    df_Au = Selecion_Nucleidos_Au(st.session_state["df_au_resultado"], st.session_state["df_file"],st.session_state["tolerancia"])
                    # Hallar los nucleidos y sus datos
                    df_filtrado_Nuclidos = Selecion_Nucleidos_muestra(st.session_state["df_resultado"],st.session_state["ref_files"], st.session_state["df_file"], st.session_state["tolerancia"])

                    #Tiempos de irradiación y decaimiento de la muestra
                    # Irraciación: (f_fin, h_fin) - (f_ini, h_ini)
                    # Decaimiento: (f_ini, h_ini) -  (f_med, h_med) 

                    #Tiempos de irradiación y decaimiento del comparador Au 
                    # Se el comparador fue irradiado en un tiempo diferente el cálculo
                    # Irraciación: (f_fin_Au, h_fin_Au) - (f_ini_Au, h_ini_Au)
                    # Decaimiento: (f_ini_Au, h_ini_Au) -  (f_med_c_Au, hora_med_c_Au)
            
                    t_irr = st.session_state["t_irr"] 
                    t_dec = st.session_state["t_dec"] 
                    t_irr_Au = st.session_state["t_irr_Au"] 
                    t_dec_Au = st.session_state["t_dec_Au"]  
                    
                    # Cálculo de f y alfa
                    alfa, f = cal_alfa(st.session_state["df_comparadores_alfa_f"])
                    # ---------forzar valores -------
                    alfa = 0.226
                    f = 34
                    st.write(f"**alfa:** {alfa}")
                    st.write(f"**f:** {f}")
                    time.sleep(1.0)
            
                if (step == "Calculando concentraciones..."):
                    st.write("Calculando concentraciones...")
                    # Calculo de la concentración
                    df_muestra = df_filtrado_Nuclidos.copy() 
                    w = st.session_state["masa_muestra"]
                    td_i = t_dec
                    ti_i = t_irr
                    tv_i = st.session_state["t_vivo"]
                    tr_i = st.session_state["t_real"]
                    df_comp_Au = df_Au.copy()
                    w_Au = st.session_state["masa_comparador_au"]
    
                    td_c_Au = t_dec_Au 
                    ti_c_Au = t_irr
                    tv_c_Au = st.session_state["t_vivo_au"]
                    tr_c_Au = st.session_state["t_real_au"]
                    geom = st.session_state["geometria"]
                    C, Cn_corr_i = conc(df_muestra, w,td_i,ti_i,tv_i,tr_i, df_comp_Au, w_Au,td_c_Au,ti_c_Au,tv_c_Au,tr_c_Au, alfa, f, geom)
                    df_muestra["Net Peak Area Corr"] = Cn_corr_i
                    df_muestra["Concentracion (PPM)"] = C*1000000
                    time.sleep(1.0)
            
                if (step == "Calculando incertidumbre..."):
                    st.write("Calculando incertidumbre...")

                    # calculo de incertidumbre
                    u_e = st.session_state["u_e"]
                    u_k0 = st.session_state["u_k0"]
                    u_w = st.session_state["u_w"]
                    u_w_c_Au = st.session_state["u_w_Au"]
                    df_comp = st.session_state["df_comparadores_alfa_f"]

                    Inc_valor = np.zeros(len(df_muestra))
                    Inc_por = np.zeros(len(df_muestra))
                    Inc_valor_red = np.zeros(len(df_muestra))
                    C_red = np.zeros(len(df_muestra))
                    for i in range(len(df_muestra)):
                        Val_ini,u_v_ini = parametros_cal_U(i,df_muestra,u_e,u_k0,u_w,td_i,ti_i,tr_i,tv_i,w,  df_comp, df_comp_Au,u_w_c_Au,td_c_Au,ti_c_Au,tr_c_Au,tv_c_Au,w_Au, geom,alfa )
                        u_y, y_val, u_y_por, simbolos = cal_U(Val_ini,u_v_ini)
                        Inc_valor[i] = 1000000*u_y
                        Inc_por[i] = round(u_y_por,2)
                        x_red, u_red = redondear_con_incert(1000000*C[i], 1000000*u_y, sig_inc = 3)
                        C_red[i] = x_red
                        Inc_valor_red[i] = u_red

                if (step == "Generando resultados..."):
                    st.write("Generando resultados...")
                    df_ejemplo = pd.DataFrame()
                    df_ejemplo["Nucleido"] =  df_muestra["NUCLID"]
                    df_ejemplo["Energía (keV)"] = df_muestra["EGKEV"]
                    df_ejemplo["Área Neto"] = df_muestra["Net Peak Area"]
                    df_ejemplo["Concentración (ppm)"] = C_red
                    df_ejemplo["Incertidumbre (ppm)"] = Inc_valor_red
                    df_ejemplo["% Incertidumbre"] = Inc_por 
                    df_ejemplo = df_ejemplo.dropna()
                    time.sleep(1.0)

            st.success("✅ Procesamiento completado!")
            status_text.text("✅ Procesamiento finalizado")
            
            # Mostrar resultados
            st.subheader("📋 Resultados del Procesamiento")
            st.dataframe(df_ejemplo, use_container_width=True)
            
            # Guardar sesión
            st.session_state['resultados'] = df_ejemplo
            st.session_state['procesado'] = True


# ============================================
# SECCIÓN 4: RESULTADOS
# ============================================
elif page == "📈 Resultados":
    st.markdown('<h2 class="section-header">📈 Visualización de Resultados</h2>', unsafe_allow_html=True)
    
    if 'resultados' in st.session_state:
        df_resultados = st.session_state['resultados']
        
        # Mostrar tabla de resultados
        st.subheader("📊 Tabla de Resultados")
        st.dataframe(df_resultados, use_container_width=True)
        
        # Gráficos
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("📈 Concentraciones por Elemento")
            fig1, ax1 = plt.subplots(figsize=(8, 5))
            bars = ax1.bar(df_resultados['Nucleido'], df_resultados['Concentración (ppm)'])
            ax1.set_ylabel('Concentración (ppm)')
            ax1.set_xlabel('Nucleido')
            ax1.set_title('Concentraciones Calculadas')
            ax1.tick_params(axis='x', rotation=45)
            
            # Añadir etiquetas de valor
            for bar in bars:
                height = bar.get_height()
                ax1.text(bar.get_x() + bar.get_width()/2., height,
                        f'{height:.2f}', ha='center', va='bottom', fontsize=9)
            
            st.pyplot(fig1)
        
        with col2:
            st.subheader("📊 Incertidumbre Relativa")
            fig2, ax2 = plt.subplots(figsize=(8, 5))
            colors = ['#FF6B6B' if x > 10 else '#4ECDC4' for x in df_resultados['% Incertidumbre']]
            bars = ax2.bar(df_resultados['Nucleido'], df_resultados['% Incertidumbre'], color=colors)
            ax2.axhline(y=10, color='r', linestyle='--', alpha=0.5, label='Límite 10%')
            ax2.set_ylabel('Incertidumbre Relativa (%)')
            ax2.set_xlabel('Nucleido')
            ax2.set_title('Incertidumbre por Elemento')
            ax2.tick_params(axis='x', rotation=45)
            ax2.legend()
            
            st.pyplot(fig2)
        
        # Estadísticas resumidas
        st.subheader("📋 Resumen Estadístico")
        col_stat1, col_stat2, col_stat3, col_stat4 = st.columns(4)
        
        with col_stat1:
            st.metric("Número de nucléidos", len(df_resultados))
        with col_stat2:
            min_conc = df_resultados['Concentración (ppm)'].min()
            st.metric("Concentración menor", f"{min_conc:.2f} ppm")
        with col_stat3:
            max_conc = df_resultados['Concentración (ppm)'].max()
            st.metric("Concentración mayor", f"{max_conc:.2f} ppm")
        with col_stat4:
            max_uncert = df_resultados['% Incertidumbre'].max()
            st.metric("Incertidumbre mayor", f"{max_uncert:.2f}%")
        
        # Botón para exportar
        st.download_button(
            label="📥 Descargar Resultados (Excel)",
            data=df_resultados.to_csv(index=False).encode('utf-8'),
            file_name="resultados_k0_analisis.csv",
            mime="text/csv",
            use_container_width=True
        )
    else:
        st.warning("⚠️ No hay resultados disponibles. Por favor, ejecute el procesamiento primero.")

# ============================================
# SECCIÓN 5: REPORTE
# ============================================
elif page == "📄 Reporte (EN DESARROLLO)":
    st.markdown('<h2 class="section-header">📄 Generación de Reporte</h2>', unsafe_allow_html=True)
    
    # Información del reporte
    col_info1, col_info2 = st.columns(2)
    
    with col_info1:
        proyecto = st.text_input("Nombre del Proyecto:", value="Evaluación Elemental por k0-INAA")
        operador = st.text_input("Nombre del Operador:", value="José Oliden")
        laboratorio = st.text_input("Laboratorio:", value="Laboratorio de Análisis por Activación Neutrónica")
    
    with col_info2:
        muestra_id = st.text_input("ID de Muestra:", value="6824a2131025G50")
        fecha_analisis = st.date_input("Fecha de Análisis:", value=datetime.now())
        metodo = st.selectbox("Método:", ["k0-INAA", "k0-EDXRF", "k0-PIXE"])
    
    # Parámetros del reporte
    st.subheader("⚙️ Configuración del Reporte")
    incluir_graficos = st.checkbox("Incluir gráficos", value=True)
    incluir_datos_crudos = st.checkbox("Incluir datos crudos", value=False)
    formato = st.radio("Formato del reporte:", ["PDF", "HTML", "Word"], horizontal=True)
    
    # Vista previa
    st.subheader("👁️ Vista Previa del Reporte")
    if st.button("🔄 Generar Vista Previa", type="secondary"):
        with st.expander("📋 Contenido del Reporte", expanded=True):
            st.markdown(f"""
            ## Reporte de Análisis k0
            
            ### Información General
            - **Proyecto:** {proyecto}
            - **Operador:** {operador}
            - **Laboratorio:** {laboratorio}
            - **ID Muestra:** {muestra_id}
            - **Fecha de Análisis:** {fecha_analisis.strftime('%d/%m/%Y')}
            - **Método:** {metodo}
            
            ### Parámetros de Análisis
            - **Geometría:** 50 mm
            - **Comparadores:** Au, Co, Mo
            - **Fecha Irradiación:** 26/09/2025 08:45:00 - 26/09/2025 09:45:00
            - **Masa muestra:** 0.2817 g
            - **Masa comparador Au:** 16.82 μg
            
            ### Resumen de Resultados
            - **Número de elementos detectados:** 17
            - **Concentración promedio:** 514.2 ppm
            - **Incertidumbre promedio:** 6.3%
            
            ### Próximos pasos
            1. Verificar resultados
            2. Validar con estándares
            3. Archivar reporte
            """)
    
    # Botón para generar reporte completo
    if st.button("🖨️ Generar Reporte Completo", type="primary", use_container_width=True):
        st.success("✅ Reporte generado exitosamente!")
        st.info("📄 El reporte se ha generado y está listo para descargar")
        
        # Crear un archivo de ejemplo (en realidad sería un PDF generado)
        reporte_texto = f"""
        REPORTE DE ANÁLISIS k0-INAA
        ============================
        
        Proyecto: {proyecto}
        Operador: {operador}
        Laboratorio: {laboratorio}
        Muestra ID: {muestra_id}
        Fecha: {fecha_analisis.strftime('%d/%m/%Y')}
        
        RESULTADOS:
        -----------
        
        Este es un reporte de ejemplo generado por el sistema.
        
        Para generar el reporte PDF completo, se necesitaría implementar
        la biblioteca ReportLab o similar.
        """
        
        st.download_button(
            label="📥 Descargar Reporte (.txt)",
            data=reporte_texto.encode('utf-8'),
            file_name=f"reporte_{muestra_id}.txt",
            mime="text/plain",
            use_container_width=True
        )

# Pie de página
st.markdown("---")
st.markdown(
    """
    <div style='text-align: center; color: #6B7280;'>
        <p>Sistema de Análisis k0 - AAN Versión 0.8.2 - beta | Desarrollado para análisis por activación neutrónica</p>
        <p>© 2025 Laboratorio de Técnicas Analíticas - Instituto Peruano de Energía Nuclear</p>
    </div>
    """,
    unsafe_allow_html=True
)
