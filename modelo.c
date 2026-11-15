#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// parâmetros do modelo 
double C_m = 1; //capacitância da membrana
double g_Na = 120; //condutância máxima do canal de sódio
double g_K = 36; //condutância máxima do canal de potássio
double g_L = 0.3; //condutância máxima do canal de vazamento
double E_Na = 50; //potencial de reversão do sódio
double E_K = -77; //potencial de reversão do potássio
double E_L = -54.4; //potencial de reversão do vazamento

// Parâmetros Sinápticos
double E_syn = 40.0; // Potencial de reversão sináptico
double tau_r = 0.4; // tempo de subida sináptica
double tau_d = 2.0; // tempo de decaimento sináptico

// Funções (5-10) v = V/mv
double alphaN(double V) {
    return (0.01 * (V + 55)) / (1 - exp(-(V + 55)/10));
}
double alphaM(double V) {
    return (0.1 * (V + 40)) / (1 - exp((V + 40)/10));
}
double alphaH(double V) {
    return 0.07 * exp(-(V + 65)/20);
}
double betaN(double V) {
    return 0.125 * exp(-(V + 65)/80);
}
double betaM(double V) {
    return 4 * exp(-(V + 65)/18);
}
double betaH(double V) {
    return 1 / (1 + exp(-(V + 35)/10));
}

// Funções (1-4) y = [V, n, m, h] I = corrente constante 

void derivadas(double y[4], double I, double dydt[4]){
    //variáveis do sistema 
    double V = y[0];
    double n = y[1];
    double m = y[2];
    double h = y[3];

    //correntes iônicas
    double I_Na = g_Na * m * m * m * h * (V - E_Na);
    double I_K = g_K * n * n * n * n * (V -E_K);
    double I_L = g_L * (V - E_L);

    //equações diferenciais
    dydt[0] = (-I_K - I_Na - I_L + I) / C_m;
    dydt[1] = alphaN(V) * (1 - n) - betaN(V) * n;
    dydt[2] = alphaM(V) * (1 - m) - betaM(V) * m;
    dydt[3] = alphaH(V) * (1 - h) - betaH(V) * h;
}

//RK4(Runge-Kutta 4 ordem)
double rk4(double y[4], double I, double dt, double resultado[4]) {
    double k1[4], k2[4], k3[4], k4[4], temp[4];
    //primeiro estágio  
    derivadas(y, I, k1);

    // segundo estágio
    for(int i = 0; i < 4; i++) {
        temp[i] = y[i] + 0.5 * dt * k1[i];
    }
    derivadas(temp, I, k2);

    //terceiro estágio
    for(int i = 0; i < 4; i++) {
        temp[i] = y[i] + 0.5 * dt * k2[i];
    }
    derivadas(temp, I, k3);

    //quarto estágio
    for(int i = 0; i < 4; i++) {
        temp[i] = y[i] + dt * k3[i];
    }
    derivadas(temp, I, k4);
    //atualização das variáveis 
    for(int i = 0; i < 4; i++) {
        resultado[i] = y[i] + (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
} 
