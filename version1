#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"

gsl_rng *tau;

//Definimos aparte las funciones que calculan w y theta porque del otro modo hay problemas de memoria
double w(int i, int j, int k, int l, int N, int p, int xi[][N][p], double a[]);
double th(int i, int j, int N, int p, int xi[][N][p], double a[]);
double solap(int mu,int N, int p, int k, int xi[][N][p], double a[], int s[][N]);

int main(void)
{
     //Comandos para el vector aleatorio
    extern gsl_rng *tau;
    int semilla=187465;
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);

    //Definimos el tamaño de la red.
    int N=50;
    //Definimos el número de patrones
    int p=1;
    //Definimos las matrices de todas las variables en juego
    int s[N][N], xi[N][N][p];
    double a[p];
    double T; //temperatura
    T=0.0001;

    //Los patrones están guardados en el fichero "p1.dat". Los metemos en xi
    FILE *fpatrones;
    //Haremos tres ficheros que guarden las configuraciones de la red: 
    //evols.txt guarda s en cada paso montecarlo
    //inic.txt guarda solo la configuracion de s inicial
    //fin.txt guarda solo la configuracion de s final
    FILE *fevols;
    FILE *finic;
    FILE *ffin;
    //En solap.txt guardaremos el solapamiento en cada paso
    FILE *fsolap;
    fevols=fopen("evols.txt","w");
    finic=fopen("inic.txt","w");
    ffin=fopen("fin.txt","w");
    fsolap=fopen("solap.txt","w");
    fpatrones=fopen("patrones.txt","r");

    //Índices que necesitaremos a lo largo del programa
    int i,j,k,l,n,m,mu;
    //Aprovechamos este bucle para calcular a, que es fijo y depende del patrón
    //Definimos el vector a, que tiene dimensión p (p=nºpatrones)
    for(mu=0;mu<p;mu++)
    {
        a[mu]=0;
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                fscanf(fpatrones,"%d,", &xi[i][j][mu]);
                a[mu]=a[mu]+xi[i][j][mu];
            }
            fscanf(fpatrones,"\n");
        }
        a[mu]=a[mu]/(N*N);
    }
    fclose(fpatrones);

    //CASO ALEATORIO: definimos el parámetro que dirá si la configuración inicial es o no aleatoria.
    int aleatoria;
    aleatoria=0;//La elección 1 implica aleatoria.
    
    //CASO DEFORMACION: definimos el número de iteraciones que hará para deformar el patrón
    int iter;
    int patron=0; //elige qué patrón deforma como inicial
    //Definimos el parámetro de deformación (se usará si aleatoria=0)
    double def; //Debe ser un número entre 0 (nada deformado) y 1 (completamente deformado).
    def=0.3; //Máximo 2 decimales.

    if(aleatoria==1)
    {
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                s[i][j]=gsl_rng_uniform_int(tau,2);
            }
        }
    }
    else
    {
        //Introducimos el patrón en s
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                s[i][j]=xi[i][j][patron];
            }
        }
        //Lo deformamos en función del parámetro def
        iter = def*N*N;
        for(int l=0;l<iter;l++)
        {
            //Tomamos una configuración al azar y la cambiamos
            i=gsl_rng_uniform_int(tau,N);
            j=gsl_rng_uniform_int(tau,N);
            s[i][j]=1-s[i][j];
        }
    }


    //Guardamos la configuracion inicial en evols e inic.
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            fprintf(fevols,"%i,",s[i][j]);
            fprintf(finic,"%i,",s[i][j]);
        }
        fprintf(fevols,"\n");
        fprintf(finic,"\n");
    }
    fprintf(fevols, "\n");
    fprintf(finic,"\n");


    int mont=30; //numero de pasos de montecarlo
    double delH,C,prob,x; //variables para el incremento de energía y los saltos de configuraciones
    double solapamiento;
    //Guardamos el solapamiento previo a hacer el algoritmo
    fprintf(fsolap,"%i\t",0);
    for(mu=0;mu<p;mu++)
    {
        solapamiento=solap(mu,N,p,k,xi,a,s);
        fprintf(fsolap,"%f\t",solapamiento);
    }
    fprintf(fsolap,"\n");

    for(k=0;k<mont;k++)
    {
        for(l=0;l<=N*N;l++)
        {
            n=gsl_rng_uniform_int(tau,N);
            m=gsl_rng_uniform_int(tau,N);

            
            delH=0;
            for(i=0;i<N;i++)
            {
                if(i!=n)
                {
                    for(j=0;j<N;j++)
                    {
                        if (j!=m)
                        {
                            delH=delH+w(i,j,n,m,N,p,xi,a)*s[i][j];
                        }               
                    }
                }
            }
            delH=(1-2*s[n][m])*(th(n,m,N,p,xi,a)-delH);
            C=exp(-delH/T);

            //determinamos cuanto vale p
            if(1<=C)
            {
                prob=1;
            }
            else prob=C;

            //Generamos un numero aleatorio uniformemente
            x=gsl_rng_uniform(tau); 
            if(x<prob)
            {
                s[n][m]=1-s[n][m];
            }
        }
        //Bucle para guardar en fichero en cada paso de montecarlo
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                fprintf(fevols,"%i,",s[i][j]);
            }
            fprintf(fevols,"\n");
        }
        fprintf(fevols, "\n");

        //Calculamos el solapamiento en cada paso de montecarlo;
        fprintf(fsolap,"%i\t",k+1);
        for(mu=0;mu<p;mu++)
        {
            solapamiento=solap(mu,N,p,k,xi,a,s);
            fprintf(fsolap,"%f\t",solapamiento);
        }
        fprintf(fsolap,"\n");

    }  
    fclose(fevols);
    fclose(ffin);  
    fclose(fsolap);
    
    return 0;
}



//FUNCIONES AUXILIARES

//Definimos la función que calcula w.
double w(int i, int j, int k, int l, int N, int p, int xi[][N][p], double a[])
{
    int mu;
    double suma=0;
    if((i==k)&&(j==l))
    {
        return 0; 
    }
    else
    {
        for(mu=0;mu<p;mu++)
        {
            suma=suma+(xi[i][j][mu]-a[mu])*(xi[k][l][mu]-a[mu]);
        }
        return suma/(N*N);
    }
}

//Definimos la función que calcula theta.
double th(int i, int j, int N, int p, int xi[][N][p], double a[])
{
    int k,l;
    double suma=0;
    for(k=0;k<N;k++)
    {
        for(l=0;l<N;l++)
        {
            suma=suma+w(i,j,k,l,N,p,xi,a);
        }
    }
    return 0.5*suma;
}

double solap(int mu, int N, int p, int k, int xi[][N][p], double a[], int s[][N])
{
    int i,j;
    double suma=0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            suma=suma+(xi[i][j][mu]-a[mu])*(s[i][j]-a[mu]);
        }
    }
    return suma/(N*N*a[mu]*(1-a[mu]));
}
