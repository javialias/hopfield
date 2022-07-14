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
    //Definimos el número de patrones que hay en nuestro archivo
    int p=5;
    //Definimos las matrices de todas las variables en juego
    int s[N][N], xi[N][N][p];
    double a[p];
    int t;
    double T;
    //Vector de temperaturas generado con escala logarítmica en base 10;
    double Temp[100]={1.00000000e-04,1.08984148e-04,1.18775445e-04,1.29446407e-04
,1.41076064e-04,1.53750546e-04,1.67563723e-04,1.82617896e-04
,1.99024558e-04,2.16905218e-04,2.36392304e-04,2.57630139e-04
,2.80776012e-04,3.06001344e-04,3.33492958e-04,3.63454458e-04
,3.96107745e-04,4.31694651e-04,4.70478737e-04,5.12747243e-04
,5.58813215e-04,6.09017821e-04,6.63732883e-04,7.23363628e-04
,7.88351687e-04,8.59178369e-04,9.36368225e-04,1.02049293e-03
,1.11217553e-03,1.21209502e-03,1.32099143e-03,1.43967126e-03
,1.56901346e-03,1.70997595e-03,1.86360272e-03,2.03103154e-03
,2.21350242e-03,2.41236676e-03,2.62909736e-03,2.86529935e-03
,3.12272209e-03,3.40327206e-03,3.70902706e-03,4.04225154e-03
,4.40541340e-03,4.80120226e-03,5.23254938e-03,5.70264936e-03
,6.21498382e-03,6.77334716e-03,7.38187469e-03,8.04507324e-03
,8.76785453e-03,9.55557156e-03,1.04140582e-02,1.13496727e-02
,1.23693440e-02,1.34806242e-02,1.46917434e-02,1.60116714e-02
,1.74501837e-02,1.90179340e-02,2.07265333e-02,2.25886358e-02
,2.46180322e-02,2.68297527e-02,2.92401774e-02,3.18671582e-02
,3.47301508e-02,3.78503590e-02,4.12508913e-02,4.49569324e-02
,4.89959297e-02,5.33977966e-02,5.81951336e-02,6.34234706e-02
,6.91215290e-02,7.53315095e-02,8.20994038e-02,8.94753358e-02
,9.75139324e-02,1.06274728e-01,1.15822607e-01,1.26228282e-01
,1.37568817e-01,1.49928203e-01,1.63397975e-01,1.78077891e-01
,1.94076672e-01,2.11512808e-01,2.30515432e-01,2.51225279e-01
,2.73795730e-01,2.98393944e-01,3.25202097e-01,3.54418735e-01
,3.86260238e-01,4.20962430e-01,4.58782318e-01,5.00000000e-01};
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
    for(t=0;t<100;t++)
    {
        T=Temp[t];
        fprintf(fsolap,"%f\t",T);
        fprintf(fevols,"%f\t",T);
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
        int patron=1; //elige qué patrón deforma como inicial (0,1,2,...)
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
    
        
        int mont=25; //numero de pasos de montecarlo
        double delH,C,prob,x; //variables para el incremento de energía y los saltos de configuraciones
        double solapamiento;
        double media[p];
        for(mu=0;mu<p;mu++)
        {
            media[mu]=0;
        }

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
            

            //Calculamos el solapamiento en cada paso de montecarlo
            //Aquí solamente nos interesa el último valor después de un cierto número de pasos de montecarlo
            
            
            int indic=0;
            if(k>19)
            {
                indic=1;
            }
            //Descomentar para que lo haga con la media
            //for(mu=0;mu<p;mu++)
            //{
            //    solapamiento=fabs(solap(mu,N,p,k,xi,a,s));//lo ponemos en valor absoluto por los estados espúreos
            //    media[mu]=media[mu]+indic*solapamiento/5;
            //}
            
        }
        //Descomentar para que lo haga con la media
        //for(mu=0;mu<p;mu++)
        //{
        //    fprintf(fsolap,"%f\t",media[mu]);
        //}
        //fprintf(fsolap,"\n");

        //Descomentar para que lo haga con el último paso MC
        for(mu=0;mu<p;mu++)
        {
            solapamiento=fabs(solap(mu,N,p,k,xi,a,s));
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
