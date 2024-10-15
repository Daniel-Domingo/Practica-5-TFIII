///     -------Practica 5 TF3, sim molecular------      ///

///LO DEJO POR AQUÍ PARA TENER LO QUE SABEMOS QUE FUNCIONA, PERO ESTA ES LA VERSION VIEJA DEL CÓDIGO

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define m 1     //masa
#define k 1     //cte elástica
#define K_b_T 1 //cte de Boltzmann por temperatura
#define PI 3.14159265358979323846

#define T 1000000 //Número de pasos de tiempo en cualquier algoritmo

#define N 10000 //Numero de datos gaussianos que se quieren, es solo para comprobaciones
#define Intervalo 100 //Numero de intervalos en el histograma

double Random_C ();
double num_aleatorio_gaussiano ();
void generador_vector_gaussiano (/*vector de salida*/ double *vector_numeros_gaussianos);
void histograma (double V[N]/*matriz de entrada*/, double *histograma_matriz /*matriz de salida*/);
double g_xn_pn (double posicion, double momento, double h, double nabla);
double f_pn (double momento);
double termino_estocastico_Z (double nabla, double h);
void Euler_Maruyama (double posicion, double momento, double h, double nabla);
void Runge_Kutta (double posicion, double momento, double h, double nabla);



///Practica 5 TF3, sim molecular

///Primera parte: Oscilador armónico
int main()
{
    //índices
    int i, j;

    //variables, aunque creo que no se toca su valor en ningún momento xd
    double /*coef de viscosidad*/nabla, /*paso de tiempo*/h=0.00001, /*Pto inicial en el espacio de fases*/ momento_inicial=0, posision_inicial=0;

    //cosas números aleatorios
    int seed=12131;
    srand(seed);

    for (i = 0; i < 4; i++){
        h = h*10;
        nabla=0.01;
        for (j = 0; j < 3; j++){
            nabla = nabla*10;
            Euler_Maruyama(posision_inicial, momento_inicial, h, nabla);
            Runge_Kutta (posision_inicial, momento_inicial, h, nabla);
        }
    }

    //aqui ejecuto las funciones que sacan los numeros gaussianos aleatorios
    //son solo de comprobación así que para ejecutar el programa de forma normal no se deben activar
    //double vector_numeros_gaussianos[N], histograma_matriz[Intervalo];
    //generador_vector_gaussiano(vector_numeros_gaussianos);
    //histograma(vector_numeros_gaussianos, histograma_matriz);

    return 0;
}



//genera un numero random en el intervalo [0,1)
//Mejor utilizar un Parisi-Rapuano
double Random_C (){
    double random;
    random=(rand()/((double)RAND_MAX+1));
    return random;
}



//Genera un número aleatorio con distribución gaussiana, a partir de un numero random en el intervalo [0,1)
double num_aleatorio_gaussiano (){
    double aleatorio_uniforme_1, aleatorio_uniforme_2, auxiliar_1, auxiliar_2, numero_gaussiano;

    ///Aquí hay un problema con el algoritmo: Cuando sale aleatorio_uniforme_1 = 0 ---> ln(0)=-inf; y da error
    //No se como se supone que deberíamos arreglarlo, de momento solo voy a poner una cláusula de que no sea igual a cero
    aleatorio_uniforme_1=Random_C ();
    while (aleatorio_uniforme_1 == 0.0)
        aleatorio_uniforme_1=Random_C ();
    aleatorio_uniforme_2=Random_C ();

    auxiliar_1=sqrt(-2*log(aleatorio_uniforme_1));
    auxiliar_2=2*PI*aleatorio_uniforme_2;

    //en la presentación dan como dos posibilidades, según las pruebas que he hecho es indistinto usar una u otra
    numero_gaussiano= auxiliar_1*cos(auxiliar_2);

    //Esta sería la segunda forma, solo cambia el cos por el sen
    //numero_gaussiano= auxiliar_1*sin(auxiliar_2);

    return numero_gaussiano;
}



//Esta función es para comprobar que la funcion de numeros gaussianos rulaba bien, probablemente se quitará en un futuro
void generador_vector_gaussiano (/*vector de salida*/ double *vector_numeros_gaussianos){
    int i;
    FILE* f;
    f = fopen("gaussiana.txt", "w");

    for (i = 0 ; i < N; i++ /*i=i+2*/){
        vector_numeros_gaussianos[i]= num_aleatorio_gaussiano();
        fprintf(f,"%f\n", vector_numeros_gaussianos[i]);
    }
    fclose(f);
}


//Lo mismo esta función es para comprobar que la funcion de numeros gaussianos rulaba bien, probablemente se quitará en un futuro
void histograma (double V[N]/*matriz de entrada*/, double *histograma_matriz /*matriz de salida*/){
    //Variables matriz histograma
    int asignacion, i;
    double anchura, min_desconocido, max_desconocido;

    max_desconocido=-100;
    min_desconocido=100;
    for (i=0;i<N;i++){
        if (V[i]>max_desconocido){
            max_desconocido=V[i];
        }

        if (V[i]<min_desconocido){
            min_desconocido=V[i];
        }
    }

    anchura=((max_desconocido-min_desconocido)/(double)Intervalo);

    for (i=0;i<Intervalo;i++){
        histograma_matriz[i]=0.0;
    }
    for (i=0;i<N;i++){
        asignacion=(int)((V[i]-min_desconocido)/anchura);
        if (asignacion==Intervalo){
                asignacion=asignacion-1;
        }
        histograma_matriz[asignacion]=histograma_matriz[asignacion]+1;
    }
    for (i=0;i<Intervalo;i++){
        histograma_matriz[i]=histograma_matriz[i]/N;
    }

//gráfico en gnuplot
    FILE* f;
    f = fopen("histo.txt", "w");
    int eje_x;
    eje_x=(int)(Intervalo/2);
    for (i=0;i<Intervalo;i++){
        fprintf(f,"%d %f\n", i-eje_x, histograma_matriz[i]);
    }
    fclose(f);

//comprobación de que las cosas vayan bien xd
    printf("%f\n", max_desconocido);
}



//ecuacion_oscilador_posicion_punto, en los algoritmos esta suele salir como f(p_n)
double f_pn (double momento){
    double x_punto;
    x_punto=momento/m;
    return x_punto;
}



//ecuacion_oscilador_momento_punto, calcula la parte NO ESTOCASTICA de p_punto, en los algoritmos esta suele salir como g(x_n, p_n)
double g_xn_pn (double posicion, double momento, double h, double nabla){
    ///IMPORTANTE, esto es (la derivada de) el potencial que tendremos que cambiar en la siguiente parte así que ojito
    double grad_Vx, p_punto;
    grad_Vx=-k*posicion;

    //aquí lo llamo p_punto, pero realmente no lo es porque le falta el término estocástico
    p_punto = -nabla*momento/m + grad_Vx;
    return p_punto;
}



//este es el término estocástico, que va incluido en la ecuación del momento
double termino_estocastico_Z (double nabla, double h){
    double eta, Z;
    eta=num_aleatorio_gaussiano();
    Z=sqrt(2*nabla*K_b_T*h)*eta;
    return Z;
}




///   ALGORITMOS   ///

//Primer algoritmo: Euler-Maruyama
void Euler_Maruyama (double posicion, double momento, double h, double nabla){
    int i;

    //Toda la parafernalia de los char es para poder sacar todos los archivos de distintas h y nabla en un solo bucle
    char nombre_archivo[1023]="Euler_Maruyama_h=", especificador_h[50], especificador_nabla[50], nombre_archivo_parte2[]="_nabla=", nombre_archivo_fin[]=".txt";
    sprintf(especificador_h, "%f", h);
    strncat(nombre_archivo, especificador_h, 1024);
    sprintf(especificador_nabla, "%f", nabla);
    strncat(nombre_archivo, nombre_archivo_parte2, 1024);
    strncat(nombre_archivo, especificador_nabla, 1024);
    strncat(nombre_archivo, nombre_archivo_fin, 1024);



    FILE *f;
    f = fopen(nombre_archivo, "w");

    for (i = 0 ; i < T; i++ /*i=i+2*/){
        posicion += f_pn(momento)*h;
        momento += g_xn_pn(posicion, momento, h, nabla)*h + termino_estocastico_Z(nabla, h);
        fprintf(f,"%f ", i*h); //Esto es el tiempo, pero lo mide en unidades de h, no se si ponerlo como i*h
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f\n", momento);
    }
    fclose(f);
}




void Runge_Kutta (double posicion, double momento, double h, double nabla){
    int i;
    //Uso la notación del power point de clase, las barras bajas se deben entender como "sub", p.ej f sub x1
    double f_x1, f_x2, g_p1, g_p2, Z, h_medios;

    //Toda la parafernalia de los char es para poder sacar todos los archivos de distintas h y nabla en un solo bucle
    char nombre_archivo[1023]="Runge_Kutta_h=", especificador_h[50], especificador_nabla[50], nombre_archivo_parte2[]="_nabla=", nombre_archivo_fin[]=".txt";
    sprintf(especificador_h, "%f", h);
    strncat(nombre_archivo, especificador_h, 1024);
    sprintf(especificador_nabla, "%f", nabla);
    strncat(nombre_archivo, nombre_archivo_parte2, 1024);
    strncat(nombre_archivo, especificador_nabla, 1024);
    strncat(nombre_archivo, nombre_archivo_fin, 1024);


    FILE *f;
    f = fopen(nombre_archivo, "w");

    h_medios=h*0.5;
    for (i = 0 ; i < T; i++ /*i=i+2*/){

        Z=termino_estocastico_Z(nabla, h);

        //Paso 1
        f_x1=f_pn(momento+Z);
        g_p1=g_xn_pn(posicion, momento+Z, h, nabla);

        //Paso 2
        f_x2=f_pn(momento+h*g_p1);
        g_p2=g_xn_pn(posicion+h*f_x1, momento+g_p1, h, nabla);

        //Final, calculo de nuevo punto en el espacio de fases
        posicion += h_medios*(f_x1+f_x2);
        momento += h_medios*(g_p1+g_p2)+Z;

        fprintf(f,"%f ", i*h); //Esto es el tiempo, pero lo mide en unidades de h, no se si ponerlo como i*h
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f\n", momento);
    }
    fclose(f);
}
