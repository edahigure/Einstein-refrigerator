#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <GLUT/glut.h>

#define LIQUID 1
#define VAPOR  0

// Variables globales
float x1, x2, x3, x4;
float h1 = 0.0, h2 = 5, h3 = 10.0;    //CM
float r = 5.0;      //CM
float H1,H2;
float rho1 = 0.89; // g/cm³
float rho2 = 0.58; // g/cm³
float L1 = 15;   //CM
float L2 = 15;   //CM
float L3 = 15;   //CM
double Vol_amo=700.0;
double Vol_but=400;
double m_gen ;
double m_absorb ;
double m_evap ;


void abc_pure(double T,double Tc,double Pc,double Zc,double F,double *sol);
double root_min_3p(double a, double b, double c);

double P_patel_teja(double T,double rho,double *xi)
{
    double sol[3];
    double am,bm,cm;
    double a_matrix[3][3],a[3],b[3],c[3],k[3][3];
    int n=3;
    double R=8.31446261815324;
    // Asumiendo el orden: 0=H2O, 1=NH3, 2=C4H10

    // Primero, inicializar todos a cero (auto-interacción)
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            k[i][j] = 0.0;
        }
    }
    // Ahora asignar los valores de la tabla (parámetros asimétricos)

    k[0][1] = -0.25054;  // k_ij H2O-NH3
    k[1][0] = -0.26238;    // k_ji NH3-H2O

    k[1][2] = 0.28337;  // k_ij NH3-C4H10
    k[2][1] = 0.13354;   // k_ji C4H10-NH3

    //H20
    abc_pure( T,647.1,220.64e5,0.26956,0.69365,sol);
    a[0]=sol[0];
    b[0]=sol[1];
    c[0]=sol[2];
    //NH3
    abc_pure( T,405.4,113.39e5,0.28289,0.64460,sol);
    a[1]=sol[0];
    b[1]=sol[1];
    c[1]=sol[2];
    //Butano
    abc_pure( T,425.12,37.907e5,0.31329,0.69851,sol);
    a[2]=sol[0];
    b[2]=sol[1];
    c[2]=sol[2];

    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        a_matrix[i][j]=(1-(k[i][j]+k[j][i])/2)*sqrt(a[i]*a[j]);
    }}

    am=0;
    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        am+= xi[i]*xi[j]*a_matrix[i][j];
    }}
    bm=0;
    for(int i=0;i<n;i++){
        bm+=xi[i]*b[i];
    }
    cm=0;
    for(int i=0;i<n;i++){
        cm+=xi[i]*c[i];
    }
    return  R*T*rho/(1.0 - bm*rho) - am*rho*rho/(1.0 + bm*rho + cm*rho*(1-bm*rho) );

}

void abc_pure(double T,double Tc,double Pc,double Zc,double F,double *sol)
{  
    double R=8.31446261815324;
   double Omega_a,Omega_b,Omega_c;
   Omega_c = 1-3*Zc;
   Omega_b = root_min_3p(2-3*Zc,3*pow(Zc,2),-pow(Zc,3));
   Omega_a = 3*pow(Zc,2) +  3* (1-2*Zc)*Omega_b+pow(Omega_b,2)+Omega_c;
   sol[0]=Omega_a*(pow(R*Tc,2)/Pc)*pow(1 + F*(1-sqrt(T/Tc)),2);
   sol[1]=Omega_b*(R*Tc/Pc);
   sol[2]=Omega_c*(R*Tc/Pc);
}

double root_min_3p(double a, double b, double c) {
    double EPS=1e-12;
    double p = b - a*a/3.0;
    double q = c - a*b/3.0 + 2.0*a*a*a/27.0;
    
    double discriminante = (q/2.0)*(q/2.0) + (p/3.0)*(p/3.0)*(p/3.0);
    
    // Caso 1: Una raíz real
    if (discriminante > EPS) {
        double term1 = -q/2.0 + sqrt(discriminante);
        double term2 = -q/2.0 - sqrt(discriminante);
        
        // Calcular raíz cúbica con signo correcto
        double raiz1 = (term1 >= 0) ? pow(term1, 1.0/3.0) : -pow(-term1, 1.0/3.0);
        double raiz2 = (term2 >= 0) ? pow(term2, 1.0/3.0) : -pow(-term2, 1.0/3.0);
        
        double x_real = raiz1 + raiz2 - a/3.0;
        
        // Verificar si es positiva
        if (x_real > EPS) {
            return x_real;
        } else {
            printf("no hay raiz positiva\n");
            return -1.0; // No hay raíz positiva
        }
    }
    // Caso 2: Tres raíces reales
    else {
        double r = sqrt(-p/3.0);
        double theta = acos(-q/(2.0*r*r*r));
        
        double min_positiva = -1.0; // -1 indica que no se encontró
        
        for (int k = 0; k < 3; k++) {
            double x_k = 2.0 * r * cos((theta - 2.0*M_PI*k)/3.0) - a/3.0;
            
            if (x_k > EPS) {
                if (min_positiva < 0 || x_k < min_positiva) {
                    min_positiva = x_k;
                }
            }
        }
        
        return (min_positiva > 0) ? min_positiva : -1.0;
    }
}

double root_max_3p(double a, double b, double c) {
    const double EPS = 1e-15;
    double p = b - a*a/3.0;
    double q = c - a*b/3.0 + 2.0*a*a*a/27.0;
    double disc = (q/2.0)*(q/2.0) + (p/3.0)*(p/3.0)*(p/3.0);

    double roots[3];
    int n_real = 0;

    if (disc > EPS) {
        double sd = sqrt(disc);
        double u = -q/2.0 + sd;
        double v = -q/2.0 - sd;
        double ru = (u >= 0) ? pow(u, 1.0/3.0) : -pow(-u, 1.0/3.0);
        double rv = (v >= 0) ? pow(v, 1.0/3.0) : -pow(-v, 1.0/3.0);
        roots[0] = ru + rv - a/3.0;
        n_real = 1;
    } else {
        double r = sqrt(-p/3.0);
        if (r < 1e-12) return -1.0;
        double theta = acos(-q / (2.0 * r*r*r));
        roots[0] = 2.0 * r * cos(theta / 3.0) - a/3.0;
        roots[1] = 2.0 * r * cos((theta + 2.0 * M_PI) / 3.0) - a/3.0;
        roots[2] = 2.0 * r * cos((theta + 4.0 * M_PI) / 3.0) - a/3.0;
        n_real = 3;
    }

    double max_pos = -1.0;
    for (int i = 0; i < n_real; i++) {
        if (roots[i] > 1e-10 && (max_pos < 0 || roots[i] > max_pos)) {
            max_pos = roots[i];
        }
    }
    return (max_pos > 0) ? max_pos : -1.0;
}

void solve_Z_cubic(double T, double P, double *xi, double *Z_L, double *Z_V) {
    double R = 8.31446261815324;
    double sol[3];
    double a_matrix[3][3],a[3],b[3],c[3],k[3][3];
    int n=3;

    // --- Calcular am, bm, cm (tu código actual) ---
    double am = 0, bm = 0, cm = 0;
    
     // Asumiendo el orden: 0=H2O, 1=NH3, 2=C4H10
  
    // Primero, inicializar todos a cero (auto-interacción)
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            k[i][j] = 0.0;
        }
    }
        // Ahora asignar los valores de la tabla (parámetros asimétricos)
        
    k[0][1] = -0.25054;  // k_ij H2O-NH3
    k[1][0] = -0.26238;    // k_ji NH3-H2O
        
    k[1][2] = 0.28337;  // k_ij NH3-C4H10
    k[2][1] = 0.13354;   // k_ji C4H10-NH3

    //H20
    abc_pure( T,647.1,220.64e5,0.26956,0.69365,sol);
    a[0]=sol[0];
    b[0]=sol[1];
    c[0]=sol[2];
    //NH3
    abc_pure( T,405.4,113.39e5,0.28289,0.64460,sol);
    a[1]=sol[0];
    b[1]=sol[1];
    c[1]=sol[2];
    //Butano
    abc_pure( T,425.12,37.907e5,0.31329,0.69851,sol);
    a[2]=sol[0];
    b[2]=sol[1];
    c[2]=sol[2];
    
    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        a_matrix[i][j]=(1-(k[i][j]+k[j][i])/2)*sqrt(a[i]*a[j]);
    }}

    am=0;
    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        am+= xi[i]*xi[j]*a_matrix[i][j];
    }}
    bm=0;
    for(int i=0;i<n;i++){
        bm+=xi[i]*b[i];
    }
    cm=0;
    for(int i=0;i<n;i++){
        cm+=xi[i]*c[i];
    }
    // --- Parámetros reducidos ---
    double RT = R * T;
    double A = am * P / (RT * RT);
    double B = bm * P / RT;
    double C = cm * P / RT;

    // --- Coeficientes Z-cúbica ---
    double coeff_A = C  - 1.0;
    double coeff_B = A - B - B*B - C - 2*B*C + B*C*C;  // + B*C*C
    double coeff_C = - A*B + B*C + B*B*C;    

    // --- Resolver ---
    double Z_min = root_min_3p(coeff_A, coeff_B, coeff_C);
    double Z_max = root_max_3p(coeff_A, coeff_B, coeff_C);

    double Z=Z_min;
   
    Z=Z_max; 

    *Z_L = (Z_min > 0) ? Z_min : -1.0;
    *Z_V = (Z_max > 0) ? Z_max : -1.0;

}

void fugacity_patel_teja(double T,double P, double Z,double *xi,double *fug)
{
    double sol[3];
    double am,bm,cm;
    double a_matrix[3][3],a[3],b[3],c[3],k[3][3];
    int n=3;
    double R=8.31446261815324;
    // Asumiendo el orden: 0=H2O, 1=NH3, 2=C4H10

    // Primero, inicializar todos a cero (auto-interacción)
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            k[i][j] = 0.0;
        }
    }
    // Ahora asignar los valores de la tabla (parámetros asimétricos)

    k[0][1] = -0.25054;  // k_ij H2O-NH3
    k[1][0] = -0.26238;    // k_ji NH3-H2O


    k[1][2] = 0.28337;  // k_ij NH3-C4H10
    k[2][1] = 0.13354;   // k_ji C4H10-NH3

    //H20
    abc_pure( T,647.1,220.64e5,0.26956,0.69365,sol);
    a[0]=sol[0];
    b[0]=sol[1];
    c[0]=sol[2];
    //NH3
    abc_pure( T,405.4,113.39e5,0.28289,0.64460,sol);
    a[1]=sol[0];
    b[1]=sol[1];
    c[1]=sol[2];
    //Butano
    abc_pure( T,425.12,37.907e5,0.31329,0.69851,sol);
    a[2]=sol[0];
    b[2]=sol[1];
    c[2]=sol[2];

    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
    a_matrix[i][j]=(1-(k[i][j]+k[j][i])/2)*sqrt(a[i]*a[j]);
    }}

    am=0;
    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
    am+= xi[i]*xi[j]*a_matrix[i][j];
    }}
    bm=0;
    for(int i=0;i<n;i++){
    bm+=xi[i]*b[i];
    }
    cm=0;
    for(int i=0;i<n;i++){
    cm+=xi[i]*c[i];
    }
  
    double V = Z*R*T/P;
    double sum_xa[3] ;

    double alpha = (bm+cm)/2;
    double beta = sqrt( bm*cm + (bm+cm)*(bm+cm)/4);
    double T1[3],T2[3],T3[3], T4[3] ;
    
    for(int i=0;i<3;i++){  

        sum_xa[i] =0;
        for(int j=0;j<n;j++){
            sum_xa[i]+=xi[j]*a_matrix[j][i];
        }  
        double gama_div_n2 =  - ( bm*c[i] + cm*b[i]   + (b[i]+c[i])*(bm+cm)/2  );

        double q = V*V +(bm+cm)*V -bm*cm;

        T1[i]= R*T*( log(V/(V-bm)) + b[i]/(V-bm)  );    
        T2[i] = (1/beta*(sum_xa[i]-am) + am/beta + am/(4*pow(beta,3))* gama_div_n2  )* log( (V+ alpha -beta)/(V+ alpha +beta)  );
        T3[i]= am/2  * (  b[i] + c[i] ) *1/q;
        

        T4[i] = (am / 2.0) * (gama_div_n2 / (beta*beta)) * (V + alpha) / ((V + alpha)*(V + alpha) - beta*beta);
        fug[i]=xi[i]*exp( 1/(R*T)*( T1[i]+ T2[i] +T3[i]+T4[i]- R*T*log(Z) ) );
       
    }
}
// Función de flash con flujo
void flash_isotermico(double T, double P, double F, double z[3], 
                   double *V, double *L, double x[3], double y[3]) {

    double Z_L, Z_V;
    double fug_L[3], fug_V[3];
    double K[3];
    double phi_L[3], phi_V[3];
    int iter = 0;
    double error = 1.0;

    double P_bar = P / 1e5;  // Pa → bar
    
    // 1. Guess inicial de K (basado en volatilidades típicas)
    K[0] = 0.3;   // H2O
    K[1] = 2.7;   // NH3
    K[2] = 1.7;   // C4H10
    
    double nu = 0.5; // Guess inicial: 50% vapor
    
    while(error > 1e-6 && iter < 100) {
        // 2. Resolver Rachford-Rice para nu (fracción vapor)
        double f_nu = 0.0;
        double df_nu = 0.0;
        
        for(int i=0; i<3; i++) {
            double denom = 1.0 + nu*(K[i] - 1.0);
            f_nu += z[i] * (K[i] - 1.0) / denom;
            df_nu += -z[i] * (K[i] - 1.0) * (K[i] - 1.0) / (denom * denom);
        }
        
        // Newton-Raphson
        double nu_new = nu - f_nu/df_nu;
        
        // Asegurar que nu está entre 0 y 1
        if(nu_new > 0.99) nu_new = 0.99;
        if(nu_new < 0.01) nu_new = 0.01;
        
        nu = nu_new;
        
        // 3. Calcular composiciones
        for(int i=0; i<3; i++) {
            x[i] = z[i] / (1.0 + nu*(K[i] - 1.0));
            y[i] = K[i] * x[i];
        }
        
        // 4. Normalizar composiciones
        double sum_x = 0.0, sum_y = 0.0;
        for(int i=0; i<3; i++) {
            sum_x += x[i];
            sum_y += y[i];
        }
        for(int i=0; i<3; i++) {
            x[i] /= sum_x;
            y[i] /= sum_y;
        }
        
        // 5. Calcular coeficientes de fugacidad para ambas fases
        solve_Z_cubic(T, P, x, &Z_L, &Z_V);
        fugacity_patel_teja(T, P, Z_L, x, fug_L);
        for(int i=0; i<3; i++) phi_L[i] = fug_L[i] / (x[i] * P_bar);
        
        solve_Z_cubic(T, P, y, &Z_L, &Z_V);
        fugacity_patel_teja(T, P, Z_V, y, fug_V);  
        for(int i=0; i<3; i++) phi_V[i] = fug_V[i] / (y[i] * P_bar);
        
        // 6. Actualizar K con atenuación
        error = 0.0;
        for(int i=0; i<3; i++) {
            double K_new = phi_L[i] / phi_V[i];
            error += fabs(K[i] - K_new);
            K[i] = 0.7 * K[i] + 0.3 * K_new; // Atenuación para estabilidad
        }
        
        iter++;
    }
    
    // 7. Convertir a flujos
    *V = F * nu;
    *L = F * (1.0 - nu);
}

double da_pure_dt(double T,double Tc,double Pc,double Zc,double F)
{  double R=8.31446261815324;
   double Omega_a,Omega_b,Omega_c;
   Omega_c = 1-3*Zc;
   Omega_b = root_min_3p(2-3*Zc,3*pow(Zc,2),-pow(Zc,3));
   Omega_a = 3*pow(Zc,2) +  3* (1-2*Zc)*Omega_b+pow(Omega_b,2)+Omega_c;
   return Omega_a*(pow(R*Tc,2)/Pc)*(1 + F*(1-sqrt(T/Tc)) ) *F*(-1/(sqrt(T/Tc))*1/Tc );
}

double entalpia_patel_teja_residual(double T,double Z,double *xi)
{
    double sol[3];
    double am,bm,cm;
    double a_matrix[3][3],a[3],b[3],c[3],k[3][3];
    int n=3;
    double R=8.31446261815324;
    double da_matrix_dT[3][3];
    double da_dT[3];
    double dam_dT;


    // Asumiendo el orden: 0=H2O, 1=NH3, 2=C4H10

    // Primero, inicializar todos a cero (auto-interacción)
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            k[i][j] = 0.0;
        }
    }
    // Ahora asignar los valores de la tabla (parámetros asimétricos)

    k[0][1] = -0.25054;  // k_ij H2O-NH3
    k[1][0] = -0.26238;    // k_ji NH3-H2O


    k[1][2] = 0.28337;  // k_ij NH3-C4H10
    k[2][1] = 0.13354;   // k_ji C4H10-NH3

    //H20
    abc_pure( T,647.1,220.64e5,0.26956,0.69365,sol);
    a[0]=sol[0];
    b[0]=sol[1];
    c[0]=sol[2];
    //NH3
    abc_pure( T,405.4,113.39e5,0.28289,0.64460,sol);
    a[1]=sol[0];
    b[1]=sol[1];
    c[1]=sol[2];
    //Butano
    abc_pure( T,425.12,37.907e5,0.31329,0.69851,sol);
    a[2]=sol[0];
    b[2]=sol[1];
    c[2]=sol[2];

    da_dT[0]=da_pure_dt( T,647.1,220.64e5,0.26956,0.69365);
    da_dT[1]=da_pure_dt( T,405.4,113.39e5,0.28289,0.64460);
    da_dT[2]=da_pure_dt( T,425.12,37.907e5,0.31329,0.69851);

    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
    a_matrix[i][j]=(1-(k[i][j]+k[j][i])/2)*sqrt(a[i]*a[j]);
    }}

    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
    da_matrix_dT[i][j]=(1-(k[i][j]+k[j][i])/2)*1/(2*sqrt(a[i]*a[j]))* (a[j] * da_dT[i]+a[i] * da_dT[j] ) ;
    }}

    am=0;
    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
    am+= xi[i]*xi[j]*a_matrix[i][j];
    }}

    dam_dT=0;
    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
    dam_dT+= xi[i]*xi[j]*da_matrix_dT[i][j];
    }}

    bm=0;
    for(int i=0;i<n;i++){
    bm+=xi[i]*b[i];
    }
    cm=0;
    for(int i=0;i<n;i++){
    cm+=xi[i]*c[i];
    }
    double alpha = (bm + cm)/2.0;
    double beta = sqrt(bm*cm + alpha*alpha);  // ¡ESTO ES EL TÉRMINO CORRECTO!
   
    return R*T*(Z -1)  + (T*dam_dT-am)/(2*beta) * log((Z +alpha +beta )/(Z +alpha -beta)  );

  
}

double entalpia_total(double T, double Z, double *xi, int phase) {
    double h_res = entalpia_patel_teja_residual(T, Z, xi);
    double Tref = 298.15;
        
    double Cp_liq[3] = {75.3, 80.8, 130.0};
    double Cp_ig[3] = {33.6, 37.1, 99.1};
    
    // ΔHvap a 25°C
    double Dh_vap[3] = {44020, 23350, 22400};
    
    if (phase == LIQUID) {
        // LÍQUIDO: referencia como líquido a 25°C
        double h_sens_liq = 0;
        for(int i=0; i<3; i++) h_sens_liq += xi[i] * Cp_liq[i] * (T - Tref);
        return h_sens_liq + h_res + 5000;
        
    } else {
        // VAPOR: referencia como GAS IDEAL a 25°C + ΔHvap
        double h_ig = 0;
        for(int i=0; i<3; i++) h_ig += xi[i] * Cp_ig[i] * (T - Tref);
        
        double total_Dh_vap = 0;
        for(int i=0; i<3; i++) total_Dh_vap += xi[i] * Dh_vap[i];
        
        return h_ig + total_Dh_vap + h_res+ 5000;
    }
}

double calcular_P(double P, double V_mol_s, double T, double *y) {
    double R = 8.31446261815324;
    double D = 0.0061;      // 6.1 mm
    double L = 0.5;         // 0.5 m
    
    // 1. Masa molar promedio (kg/mol)
    double PM = (y[0]*18.015 + y[1]*17.031 + y[2]*58.122) / 1000.0;  // kg/mol
    
    // 2. Flujo másico
    double m_dot = V_mol_s * PM;                                    // kg/s
    double area = M_PI * D * D / 4.0;
    double G = m_dot / area;                                        // kg/(m²·s)
    
    // 3. Reynolds y fricción
    double mu = 1.2e-5;                                             // Pa·s
    double Re = G * D / mu;
    double f = 0.316 / pow(Re, 0.25);                               // Blasius
      
    return sqrt(P * P -  2.0 * f * L * G * G * R * T / PM);
}

double calcular_Q_FL(double T, double P, 
                           double F, double *z,double *x,double *y,double *V,double *L) {
    
    // 1. Flash isotérmico en el generador
    //double V, L;  // moles de vapor y líquido
                   
    flash_isotermico(T, P, F, z, 
                    V, L, x, y);
    
    // 2. Factores de compresibilidad
    double Z_L, Z_V, Z_dum, ZF_L,ZF_V;
    solve_Z_cubic(T, P, x, &Z_L, &Z_dum);  // Z_L para líquido
    solve_Z_cubic(T, P, y, &Z_dum, &Z_V);  // Z_V para vapor
    solve_Z_cubic(T, P, z, &ZF_L, &ZF_V);  // Z_F
    
    // 3. Entalpías
    double h_F = entalpia_total(T, Z_L, z, LIQUID);
    double h_L = entalpia_total(T, Z_L, x, LIQUID);
    double h_V = entalpia_total(T, Z_V, y, VAPOR);
    
    // 4. Balance de energía
    return *V * h_V + *L * h_L - F * h_F;
       
}

double calcular_Q_FV(double T, double P, 
                           double F, double *z,double *x,double *y,double *V,double *L) {
    
    // 1. Flash isotérmico en el generador
    //double V, L;  // moles de vapor y líquido
                   
    flash_isotermico(T, P, F, z, 
                    V, L, x, y);
    
    // 2. Factores de compresibilidad
    double Z_L, Z_V, Z_dum, ZF_L,ZF_V;
    solve_Z_cubic(T, P, x, &Z_L, &Z_dum);  // Z_L para líquido
    solve_Z_cubic(T, P, y, &Z_dum, &Z_V);  // Z_V para vapor
    solve_Z_cubic(T, P, z, &ZF_L, &ZF_V);  // Z_F
    
    // 3. Entalpías
    double h_F = entalpia_total(T, Z_V, z, VAPOR);
    double h_L = entalpia_total(T, Z_L, x, LIQUID);
    double h_V = entalpia_total(T, Z_V, y, VAPOR);

    
    // 4. Balance de energía
    return *V * h_V + *L * h_L - F * h_F;
    
    
}

double calcular_Q_mezcla(double T, double P,
                        double F_liquido, double *x_F,  // Flujo y composición líquida
                        double F_vapor, double *y_F,    // Flujo y composición vapor  
                        double *x, double *y,           // Composiciones de salida del flash
                        double *V, double *L) {         // Flujos de salida del flash
    
    double F_total = F_liquido + F_vapor;
        
    double z[3] ;
    for (int i = 0; i < 3; i++) {
        z[i] = (F_liquido * x_F[i] + F_vapor * y_F[i]) / F_total;
    }
    
    // 2. Flash isotérmico con composición promedio
    flash_isotermico(T, P, F_total, z, V, L, x, y);
    
    // 3. Factores de compresibilidad
    double Z_L, Z_V, Z_F_liq, Z_F_vap, Z_dummy;
    solve_Z_cubic(T, P, x, &Z_L, &Z_dummy);           // Z_L para líquido salida
    solve_Z_cubic(T, P, y, &Z_dummy, &Z_V);           // Z_V para vapor salida
    solve_Z_cubic(T, P, x_F, &Z_F_liq, &Z_dummy);     // Z_F para líquido entrada
    solve_Z_cubic(T, P, y_F, &Z_dummy, &Z_F_vap);     // Z_F para vapor entrada
    
    // 4. Entalpías
    double h_L = entalpia_total(T, Z_L, x, LIQUID);
    double h_V = entalpia_total(T, Z_V, y, VAPOR);
    double h_F_liq = entalpia_total(T, Z_F_liq, x_F, LIQUID);
    double h_F_vap = entalpia_total(T, Z_F_vap, y_F, VAPOR);
    
    // 5. Balance de energía con ambas fases de entrada
    return (*V * h_V + *L * h_L) - (F_liquido * h_F_liq + F_vapor * h_F_vap);
}

void calcularAlturas(float P_gen,float P_evap,float P_abs) {
    //P en bares
    printf("\n=== CALCULANDO SISTEMA BUTANO ===\n");
    float pi = 4*atan(1);
    float g = 981.0;  // cm/s²
    
    float beta = rho2/rho1;
   
    H1 = Vol_amo/(pi*r*r)*2.0;
    H2 = Vol_but/(pi*r*r)*2.0;

    float m_but =  H2 *(pi*r*r) * rho2;
    float m_amonia =  H1 *(pi*r*r) * rho1;
    
    // ✅ TODAS LAS PRESIONES YA ESTÁN EN dyn/cm²
    
    float p1 = P_gen*1e6;
    float p2 = P_abs*1e6;
    float p3 = P_evap*1e6; 

    printf("p3: %f dyn/cm²\n", p3);
    printf("p1: %f dyn/cm²\n", p1);
    printf("p2: %f dyn/cm²\n", p2);
   
    // ✅ AHORA UNIDADES CONSISTENTES: (dyn/cm²) / ((g/cm³) × (cm/s²)) = cm
    float A = (p1-p2)/(rho1*g) + H1 - h2 ;
    float B = (p3-p2)/(rho2*g) + h3 - h2 + H2;
    
    x3 = (2*B-A)/(4-beta );
    x2 = B-2*x3;
    x1 = H1 - x2;
    x4 = H2 - x3;
    
    printf("x1=%.2f cm\nx2=%.2f cm\nx3=%.2f cm\nx4=%.2f cm\n", x1, x2, x3, x4);
    printf("m_but: %f gr\n", m_but);
    printf("m_amonia: %f gr\n", m_amonia);

    m_gen = rho1*x1*M_PI*r*r;
    m_absorb = (rho2*x3+rho1*x2)*M_PI*r*r;
    m_evap = rho2*x4*M_PI*r*r;

    printf("m_gen: %f gr\n", m_gen);
    printf("m_absorb: %f gr\n", m_absorb);
    printf("m_evap: %f gr\n", m_evap);
    
}

// Funciones de visualización OpenGL
void drawCylinder(float baseY, float height, float radius, float r, float g, float b) {
    glColor3f(r, g, b);
    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i <= 32; i++) {
        float angle = 2.0 * M_PI * i / 32;
        float x = cos(angle) * radius;
        float z = sin(angle) * radius;
        glVertex3f(x, baseY, z);
        glVertex3f(x, baseY + height, z);
    }
    glEnd();
}

void drawFloor() {
    glColor3f(0.3, 0.3, 0.3);
    glBegin(GL_QUADS);
    glVertex3f(-50, 0, -50);
    glVertex3f(50, 0, -50); 
    glVertex3f(50, 0, 50);
    glVertex3f(-50, 0, 50);
    glEnd();
    
    glColor3f(0.5, 0.5, 0.5);
    glBegin(GL_LINES);
    for(int i = -50; i <= 50; i += 5) {
        glVertex3f(-50, 0.1, i);
        glVertex3f(50, 0.1, i);
        glVertex3f(i, 0.1, -50);
        glVertex3f(i, 0.1, 50);
    }
    glEnd();
}

void drawEmptyContainer(float baseY, float height, float radius) {
    glColor3f(0.8, 0.8, 0.8);
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i <= 32; i++) {
        float angle = 2.0 * M_PI * i / 32;
        float x = cos(angle) * radius;
        float z = sin(angle) * radius;
        glVertex3f(x, baseY, z);
        glVertex3f(x, baseY + height, z);
    }
    glEnd();
    
    for (int ring = 0; ring <= 3; ring++) {
        float y = baseY + (height * ring / 3);
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i <= 32; i++) {
            float angle = 2.0 * M_PI * i / 32;
            float x = cos(angle) * radius;
            float z = sin(angle) * radius;
            glVertex3f(x, y, z);
        }
        glEnd();
    }
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void drawAxes() {
    glLineWidth(2.0);
    
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(20, 0, 0);
    
    glColor3f(0, 1, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 20, 0);
    
    glColor3f(0, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 20);
    glEnd();
    
    glLineWidth(1.0);
}

void display() {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt(0, 10, 100, 0, 0, 0, 0, 1, 0);
    drawFloor();  
    drawAxes(); 
    float pos1 = -15.0, pos2 = 0.0, pos3 = 15.0;
    
    // Recipiente 1 (solo líquido 1 - azul)
    glPushMatrix();
    glTranslatef(pos1, 0, 0);
    drawEmptyContainer(h1, L1, r+.1);    
    drawCylinder(h1, x1, r, 0.2, 0.4, 0.8);
    glPopMatrix();
    
    // Recipiente 2 (líquido 1 abajo + líquido 2 arriba)
    glPushMatrix();
    glTranslatef(pos2, 0, 0);
    drawEmptyContainer(h2, L2, r+.1);   
    drawCylinder(h2, x2, r, 0.2, 0.4, 0.8);      // Azul abajo
    drawCylinder(h2 + x2, x3, r, 0.8, 0.2, 0.2); // Rojo arriba
    glPopMatrix();
    
    // Recipiente 3 (solo líquido 2 - rojo)
    glPushMatrix();
    glTranslatef(pos3, 0, 0);
    drawEmptyContainer(h3, L3, r+.1);   
    drawCylinder(h3, x4, r, 0.8, 0.2, 0.2);
    glPopMatrix();
    
   
    glutSwapBuffers();
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, (float)w/h, 1, 200);
    glMatrixMode(GL_MODELVIEW);
}


int main(int argc, char** argv) {
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Sistema NH3 - Solución Físicamente Consistente");
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    
    glEnable(GL_DEPTH_TEST);

    double T_gen = 273.15+120;     // 90°C - temperatura típica
    double T_absorb = 273.15+20;  // 30°C
    double T_evap =  273.15+5; 


    printf("L1: %f cm\n",L1);  
    printf("L2: %f cm\n",L2);  
    printf("L3: %f cm\n",L3);  

    printf("h1: %f cm\n",h1);  
    printf("h2: %f cm\n",h2);  
    printf("h3: %f cm\n",h3);  

    printf("T_gen: %f\n",T_gen);
    printf("T_absorb: %f\n",T_absorb);
    printf("T_evap: %f\n",T_evap);
    
    
    double z_gen[3] = {0.8, 0.19, 0.01};  // % H2O, % NH3, % C4H10
    double z_evap[3];
    double z_absorb[3];

    double V_gen   , L_gen   , x_gen[3]   , y_gen[3];  // composiciones
    double V_absorb, L_absorb, x_absorb[3], y_absorb[3];
    double V_evap  , L_evap  , x_evap[3]  , y_evap[3];

    double Q1,Q2,Q3;
    double P_gen= 2.5e5,P_absorb,P_evap;

    double F_absorb = 0.017;
    double F_gen = 0.017;           
    double F_but = 0.01;           //  mol/s de alimentación
   
    
    calcularAlturas( P_gen/1e6, P_gen/1e6, P_gen/1e6);

    for(int iter=0;iter<1;iter++){

        printf("\niter:%d\n",iter);
        
        Q1 = calcular_Q_FL(T_gen, P_gen, F_gen, z_gen,x_gen,y_gen,&V_gen,&L_gen);

        printf("\nQ1: %f\n",Q1);
        printf("P_gen: %f\n",P_gen);
        printf("V_gen: %f\n",V_gen);
        printf("L_gen: %f\n",L_gen);    
        printf("z_gen: %f %f %f\n",z_gen[0],z_gen[1],z_gen[2]);
        printf("x_gen: %f %f %f\n",x_gen[0],x_gen[1],x_gen[2]);
        printf("y_gen: %f %f %f\n\n",y_gen[0],y_gen[1],y_gen[2]);

        P_evap=calcular_P(P_gen, V_gen, T_evap, y_gen);

        double masa_H2O_entra = V_gen * y_gen[0];
        double masa_NH3_entra = V_gen * y_gen[1]; 
        double masa_But_entra = V_gen * y_gen[2] + F_but;

        // 2. Masa total que entra
        double masa_total_entra = masa_H2O_entra + masa_NH3_entra + masa_But_entra;

        // 3. Composición de lo que ENTRA
        z_evap[0] = masa_H2O_entra / masa_total_entra;
        z_evap[1] = masa_NH3_entra / masa_total_entra;
        z_evap[2] = masa_But_entra / masa_total_entra;
    
        Q2 = calcular_Q_FL(T_evap, P_evap, V_gen+F_but , z_evap,x_evap,y_evap,&V_evap,&L_evap);

        printf("✅Q2: %f\n",Q2);
        printf("P_evap: %f\n",P_evap);
        printf("V_evap: %f\n",V_evap);
        printf("L_evap: %f\n",L_evap);
        printf("z_evap: %f %f %f\n",z_evap[0],z_evap[1],z_evap[2]);
        printf("x_evap: %f %f %f\n",x_evap[0],x_evap[1],x_evap[2]);
        printf("y_evap: %f %f %f\n\n",y_evap[0],y_evap[1],y_evap[2]);

        // Condiciones del absorbedor

        P_absorb=calcular_P(P_evap, V_evap, T_absorb, y_evap);

        // Composición global
        double F_total = V_evap + F_absorb ;

        for(int i=0; i<3; i++) {
            z_absorb[i] = (V_evap * y_evap[i] + F_absorb * x_gen[i]  )/ F_total ;
        }
                
        Q3 = calcular_Q_mezcla( T_absorb, P_absorb ,F_absorb,x_gen,V_evap,y_evap,x_absorb,y_absorb ,&V_absorb,&L_absorb) ;

        printf("Q3: %f\n",Q3);
        printf("P_absorb: %f\n",P_absorb); 
        printf("V_absorb: %f\n",V_absorb);
        printf("L_absorb: %f\n",L_absorb);    
        printf("z_absorb: %f %f %f\n",z_absorb[0],z_absorb[1],z_absorb[2]);
        printf("x_absorb: %f %f %f\n",x_absorb[0],x_absorb[1],x_absorb[2]);
        printf("y_absorb: %f %f %f\n\n",y_absorb[0],y_absorb[1],y_absorb[2]);
            
        F_gen = L_absorb * (x_absorb[0]+ x_absorb[1]) ;
        F_but = L_absorb * x_absorb[2] ;
        printf("F_gen:%f\n",F_gen);
        printf("F_but:%f\n",F_but);
        printf("F_absorb:%f\n",F_absorb);
        printf("Eff:%f\n",Q2/Q1*100);

        // Composición de lo que ENTRA al generador
        z_gen[0] = x_absorb[0] / (x_absorb[0] + x_absorb[1]);  // Fracción H₂O
        z_gen[1] = x_absorb[1] / (x_absorb[0] + x_absorb[1]);  // Fracción NH₃  
        z_gen[2] = 0.0;                                        // Sin butano
        
        calcularAlturas( P_gen/1e6, P_absorb/1e6, P_evap/1e6);
    }
    glutMainLoop();

    return 0;
}