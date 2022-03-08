#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SEED 378520

#define N 200 // cantidad de particulas
#define m 1.0 // masa de las particulas

// parametros del potencial de Lennard-Jones: v(r)= 4*epsilon (- (sigma/r)^6 + (sigma/r)^12)
#define EPSILON 1
#define SIGMA 1
#define RC (2.5 * SIGMA) // radio de corte

// parametros de la dinamica de minimazion de energia
#define MAX_PASOS_minimizacion 10000
#define DT_min 0.0001 

// parametros de la dinamica con Verlet
#define MAX_PASOS_verlet 1000000
#define DT_verlet 0.002 

// parametro del termostato de Langevin
#define gamma 0.5 // coeficiente de friccion

// sampleo
#define dist_sample 500// pasos temporales entre dos sampleos para descorrelacionar

static double L; // lado de la caja
static double volumen; // volumen de la caja 
static double T; // temperatura

static double r[N][3]; // posiciones
static double v[N][3]; // velocidades
static double f[N][3]; // fuerzas

static double V;      // energia potencial
static double KE;     // energia cinetica
static double p;      // presion

// Subrutinas
void posicionesIniciales(int seed);
void velocidadesIniciales(double temperatura);
void mostrar();

void fuerzas_a_cero();
void calculo_fuerzas(double shift_v);
void lgv_force(double varianza);

void minimizacion_energia(double shift_v);

void verlet_positions();
void verlet_velocities();

double uni_rand();
double normal_rand();
double gausian_rand(double media, double sigma);

// input
FILE *input_file;

// output
FILE *file_positions;
FILE *file_verlet;


int main()
{
    // Abro archivo de input
    input_file = fopen("input.txt", "r");

    // Abro achivos de salida y les escribo el header
    file_verlet = fopen("MD2_verlet.txt", "w");
    fprintf(file_verlet, "# MD : SEED %i, N %d, MAX_PASOS_minimizacion %d, DT_min %f,MAX_PASOS_verlet %d, DT_verlet %f , masa %f\n", SEED, N, MAX_PASOS_minimizacion, DT_min, MAX_PASOS_verlet, DT_verlet, m);
    fprintf(file_verlet, "N,T,L,rho,paso,V,KE,p\n");

    // Calculo shift del potencial segun radio de corte RC
    double r_c6i = 1.0 / (pow(RC, 6));
    double shift_v = 4 * EPSILON * pow(SIGMA, 6) * (pow(SIGMA, 6) * r_c6i * r_c6i - r_c6i); 

    //Leo del input las densidades y temperaturas que quiero simular
    int cant; 
    fscanf(input_file, "%d", &cant); //cantidad de pares (rho,T)
    double rho[cant];
    double t[cant];
    int i;
    for (i=0; i<cant;i++){
        fscanf(input_file, "%lf %lf", &rho[i], &t[i]);
    }

    //Loop para distintas densidades y temperaturas
    for (i = 0; i < cant; i++)
    {
        //defino L y T
        L = pow((double)N/(rho[i]*pow(SIGMA,-3)),1.0/3.0);
        volumen = pow(L,3);
        T = t[i] * EPSILON;

        // Abro archivo de output para ver con VMD las configuraciones
        char nom[100];
        sprintf(nom,"MD2_movie_T%f_rho%f.vtf",T ,rho[i]);   
        file_positions=fopen(nom,"w");
        fprintf(file_positions, "atom 0:%d    radius 0.5 name Ar\n", N - 1);

        // Calculo varianza para langevin
        double varianza = pow((2 * m * gamma * T) / DT_verlet, 0.5); 
      
        // Inicializo
        posicionesIniciales(SEED);
        velocidadesIniciales(T);

        // Dinamica de minimizaxion de energia
        int paso;
        for (paso = 0; paso < MAX_PASOS_minimizacion; paso++)
        {
            minimizacion_energia(shift_v);
        }

        // Calculo fuerzas antes de arrancar el loop de velocity-Verlet
        fuerzas_a_cero();
        lgv_force(varianza);
        calculo_fuerzas(shift_v);
        
        // Dinamica de velocity-Verlet
        for (paso = 0; paso < MAX_PASOS_verlet; paso++)
        {
            verlet_positions();

            fuerzas_a_cero();
            lgv_force(varianza);
            calculo_fuerzas(shift_v); // calcula fuerzas, energia potencial V y presion p

            verlet_velocities(); // calcula velocidades y energia cinetica KE y temperatura instantanea

            // Escribo output
            if ((paso) % dist_sample == 0){

                printf("densidad: %f, temperatura: %f, paso verlet: %d \n",rho[i], T, paso);

                fprintf(file_verlet, "%d,%f,%f,%f,%d,%f,%f,%f\n", N, T, L,rho[i], paso, V, KE,p);

                fprintf(file_positions, "timestep\n\n");
                int i;
                for (i = 0; i < N; i++)
                {
                    fprintf(file_positions, "%f %f %f\n", r[i][0], r[i][1], r[i][2]);
                }
            }
        }
        fclose(file_positions);
    }

    fclose(file_verlet);
    fclose(input_file);
    
    
    return 0;
}


void verlet_velocities()
{
    KE = 0.0;
    for (int i = 0; i < N; i++)
    {
        v[i][0] += 0.5 * DT_verlet * f[i][0] / m;
        v[i][1] += 0.5 * DT_verlet * f[i][1] / m;
        v[i][2] += 0.5 * DT_verlet * f[i][2] / m;

        KE += v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
    }
    KE = 0.5 * m * KE;
}

void verlet_positions()
{
    for (int i = 0; i < N; i++)
    {
        // Actualizo posiciones de cada particula: r(t+DT_verlet)
        r[i][0] = r[i][0] + v[i][0] * DT_verlet + 0.5 * f[i][0] / m * DT_verlet * DT_verlet;
        r[i][1] = r[i][1] + v[i][1] * DT_verlet + 0.5 * f[i][1] / m * DT_verlet * DT_verlet;
        r[i][2] = r[i][2] + v[i][2] * DT_verlet + 0.5 * f[i][2] / m * DT_verlet * DT_verlet;

        // Ajusto según condiciones de borde periódicas
        r[i][0] = r[i][0] + (1 - ((int)((r[i][0] + L) / L))) * L;
        r[i][1] = r[i][1] + (1 - ((int)((r[i][1] + L) / L))) * L;
        r[i][2] = r[i][2] + (1 - ((int)((r[i][2] + L) / L))) * L;

        // Actualizo velocidades tiempo intermedio : v(t+ DT_verlet/2)
        v[i][0] += 0.5 * DT_verlet * f[i][0] / m;
        v[i][1] += 0.5 * DT_verlet * f[i][1] / m;
        v[i][2] += 0.5 * DT_verlet * f[i][2] / m;
    }
}

void fuerzas_a_cero()
{
    // Inicializo fuerzas en cero
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            f[i][j] = 0;
        }
    }
}

void lgv_force(double var)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            f[i][j] += -gamma * v[i][j] + gausian_rand(0.0, var);
        }
    }
}

void calculo_fuerzas(double shift_v)
{
    double dx, dy, dz, r2, r6i, fuerza;
    double sigma6 = pow(SIGMA, 6);
    
    // Inicializo energia V y presion p en cero
    V = 0.0;
    p = 0.0;
    
    // Analizo interacciones
    for (int i = 0; i < (N - 1); i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            // Calculo distancia entre partículas (i,j)
            dx = (r[i][0] - r[j][0]);
            dy = (r[i][1] - r[j][1]);
            dz = (r[i][2] - r[j][2]);

            // Ajusto según condiciones de borde periódicas
            dx = dx - ((int)((2 * dx) / L)) * L;
            dy = dy - ((int)((2 * dy) / L)) * L;
            dz = dz - ((int)((2 * dz) / L)) * L;

            r2 = dx * dx + dy * dy + dz * dz;

            // AGREGRO ACA CALCULO DE G(r)?

            // Calculo fuerzas según radio de corte RC
            if (r2 < (RC * RC))
            {
                r6i = 1.0 / (r2 * r2 * r2);

                V += 4 * EPSILON * sigma6 * (sigma6 * r6i * r6i - r6i) - shift_v;
                fuerza = 48 * EPSILON * sigma6 * (sigma6 * r6i * r6i - 0.5 * r6i);

                f[i][0] += (double)(dx * fuerza / r2);
                f[j][0] -= (double)(dx * fuerza / r2);

                f[i][1] += (double)(dy * fuerza / r2);
                f[j][1] -= (double)(dy * fuerza / r2);

                f[i][2] += (double)(dz * fuerza / r2);
                f[j][2] -= (double)(dz * fuerza / r2);
                
                // Calculo la presion 
                p += fuerza;
                                
            }
        }
    }
    //p = rho*T + p/(3*N*volumen); // ESTO LO HACEMOS EN EL POSTPROCESAMIENTO
}


void minimizacion_energia(double shift_v)
{

    // Calculo fuerzas
    fuerzas_a_cero();
    calculo_fuerzas(shift_v);

    // Actualizo posiciones de cada particula
    int i;
    for (i = 0; i < N; i++)
    {

        r[i][0] = r[i][0] + 0.5 * f[i][0] / m * DT_min * DT_min;
        r[i][1] = r[i][1] + 0.5 * f[i][1] / m * DT_min * DT_min;
        r[i][2] = r[i][2] + 0.5 * f[i][2] / m * DT_min * DT_min;

        // Ajusto según condiciones de borde periódicas (solo para cajas inmediatamente vecinas)
        r[i][0] = r[i][0] + (1 - ((int)((r[i][0] + L) / L))) * L;
        r[i][1] = r[i][1] + (1 - ((int)((r[i][1] + L) / L))) * L;
        r[i][2] = r[i][2] + (1 - ((int)((r[i][2] + L) / L))) * L;
    }
}

void posicionesIniciales(int seed)
{
    srand(seed);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            r[i][j] = uni_rand() * L;
        }
    }
}

void velocidadesIniciales(double temperatura)
{
    double sigma = sqrt(3 * temperatura / m); // ¿¿¿ VA UN 3 ???
    double media = 0.0;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            v[i][j] = gausian_rand(media, sigma);
        }
    }
}

void mostrar()
{
    printf("\n\n");

    for (int i = 0; i < N; i++)
    {
        // posiciones
        for (int j = 0; j < 3; j++)
        {
            printf("\t%f", r[i][j]);
        }
        printf("\t");

        // fuerzas
        for (int j = 0; j < 3; j++)
        {
            printf("\t%f", f[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

double gausian_rand(double media, double sigma)
{
    // distribucion  gaussiana (media, sigma)
    return media + sigma * normal_rand();
}

double normal_rand()
{
    // distribucion normal (media=0, sigma=1)
    return sqrt(-2 * log(uni_rand())) * cos(2 * M_PI * uni_rand());
}
double uni_rand()
{
    // distribucion uniforme (0,1]
    return (rand() + 1.0) / (RAND_MAX + 1.0);
}
