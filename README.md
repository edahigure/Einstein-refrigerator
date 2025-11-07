# Einstein-refrigerator
Sistema de Refrigeración por Absorción - Simulación Termodinámica
📋 Descripción General
Este programa simula un sistema de refrigeración por absorción que utiliza una mezcla ternaria de agua (H₂O), amoníaco (NH₃) y butano (C₄H₁₀). El modelo incluye:

Cálculos termodinámicos usando la ecuación de estado Patel-Teja

Flash isotérmico para separación de fases

Balance de energía completo del sistema

Visualización 3D de los recipientes del sistema

🧪 Componentes del Sistema
🔄 Unidades Principales
Generador: Donde se aplica calor para separar componentes

Evaporador: Donde ocurre la refrigeración

Absorbedor: Donde se reabsorbe el refrigerante

Intercambiadores de calor: Para recuperación energética

🧮 Especies Químicas
H₂O (agua) - Índice 0

NH₃ (amoníaco) - Índice 1

C₄H₁₀ (butano) - Índice 2

🛠️ Compilación y Ejecución
Requisitos del Sistema
bash
# Linux/Ubuntu
sudo apt-get install freeglut3 freeglut3-dev libglu1-mesa-dev

# macOS (usando Homebrew)
brew install freeglut

# Windows (MinGW)
# Instalar freeglut y enlazar con -lfreeglut -lopengl32 -lglu32
Compilación
bash
gcc -o sistema_absorcion main.c -lm -lGL -lGLU -lglut
Ejecución
bash
./sistema_absorcion
📊 Funcionalidades Principales
🔢 Cálculos Termodinámicos
Ecuación de estado Patel-Teja para mezclas ternarias

Coeficientes de fugacidad y equilibrio vapor-líquido

Propiedades residuales para entalpías

Flash isotérmico con método Rachford-Rice

⚖️ Balance de Energía
c
// Calor en el generador
double Q1 = calcular_Q_FL(T_gen, P_gen, F_gen, z_gen, x_gen, y_gen, &V_gen, &L_gen);

// Calor en el evaporador  
double Q2 = calcular_Q_FL(T_evap, P_evap, V_gen+F_but, z_evap, x_evap, y_evap, &V_evap, &L_evap);

// Calor en el absorbedor
double Q3 = calcular_Q_mezcla(T_absorb, P_absorb, F_absorb, x_gen, V_evap, y_evap, x_absorb, y_absorb, &V_absorb, &L_absorb);
📐 Cálculo de Geometría
Volúmenes de líquido en cada recipiente

Niveles basados en diferencias de presión

Masas de cada componente

🎨 Visualización 3D
Recipientes cilíndricos con niveles de líquido

Código de colores: azul para H₂O/NH₃, rojo para butano

Sistema de coordenadas de referencia

📈 Parámetros de Operación
🌡️ Condiciones Típicas
c
double T_gen = 393.15;      // 120°C - Generador
double T_absorb = 293.15;   // 20°C - Absorbedor  
double T_evap = 278.15;     // 5°C - Evaporador

double P_gen = 2.5e5;       // 2.5 bar - Presión inicial
🔄 Flujos Molares
c
double F_absorb = 0.017;    // mol/s - Flujo de absorbente
double F_gen = 0.017;       // mol/s - Flujo al generador
double F_but = 0.01;        // mol/s - Flujo de butano
📁 Estructura del Código
Funciones Principales
abc_pure() - Parámetros para componentes puros

P_patel_teja() - Ecuación de estado

flash_isotermico() - Cálculo de equilibrio VLE

entalpia_total() - Cálculo de entalpías

calcularAlturas() - Geometría del sistema

Funciones OpenGL - Visualización 3D

Variables Globales Importantes
c
float x1, x2, x3, x4;           // Niveles de líquido (cm)
double m_gen, m_absorb, m_evap;  // Masas (gramos)
double Vol_amo = 700.0;          // Volumen amoníaco (cm³)
double Vol_but = 400.0;          // Volumen butano (cm³)
📊 Salida del Programa
El programa genera:

Iteraciones del cálculo convergente

Calores Q1, Q2, Q3 en cada unidad

Presiones y temperaturas

Composiciones de todas las corrientes

Flujos vapor/líquido

Eficiencia del sistema (Q2/Q1)

Visualización 3D en tiempo real

🎯 Aplicaciones
Diseño de sistemas de refrigeración por absorción

Optimización de condiciones operativas

Estudio de mezclas refrigerante-absorbente

Análisis de eficiencia energética

⚠️ Notas Importantes
Unidades: El programa usa sistema CGS (cm, g, s) consistentemente

Convergencia: El método iterativo puede requerir ajustes para diferentes condiciones

Parámetros: Los coeficientes de interacción binaria son específicos para esta mezcla

Visualización: Requiere soporte OpenGL en el sistema

🔍 Personalización
Para modificar el sistema:

Ajustar z_gen[] para cambiar composiciones iniciales

Modificar temperaturas en T_gen, T_evap, T_absorb

Cambiar flujos en F_gen, F_absorb, F_but

Ajustar volúmenes en Vol_amo, Vol_but

Sistema desarrollado para simulación de ciclos de refrigeración por absorción 🚀