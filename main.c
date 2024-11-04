///     -------Practica 5 TF3, sim molecular------      ///
//actualizado el 04-11-24




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define m 1     //masa
#define k 1     //cte el�stica
#define A 1     //altura del doble pozo
#define K_b_T 0.2 //cte de Boltzmann por temperatura
#define PI 3.14159265358979323846

//para los numerros random, NO TOCAR
#define NormRANu (2.3283063671E-10F)

#define T 100000 //Tiempo en cualquier algoritmo
//El numero de pasos vendr� dado por T/h, as� pasa el mismo tiempo independientemente de h
#define h 0.01
const int pasos = (int)T/h;

#define Ng 1000 //Numero de datos gaussianos que se quieren, es solo para comprobaciones, SIEMPRE PAR
#define Intervalo 100 //Numero de intervalos en el histograma


///Estos �ltimos son para los ifdef:

//si se pone en "bucle" saca varios archivos de una para distintos h y nabla (en principio todos los que piden)
//si se pone en "unico" saca un solo archivo para un dado h y nabla
#define bucle //unico


//Calcular las energ�as durante el proceso, si (true) o no (false)
#define true

//Para elegir el tipo de potencial que se utiliza
#define doble_pozo //oscilador_armonico



//Numeros aleatorios uniformes
double Random_C();
void ini_ran(int SEMILLA);
double Random(void);
void num_aleatorio_gaussiano (double *dos_numeros_gaussianos);

void Histograma(int a, int b, int N, int divisiones,  double *array);

//Numeros aleatorios gaussianos
void generador_vector_gaussiano (/*vector de salida*/ double *t_estancia);
void histograma (double *V/*matriz de entrada*/, double *histograma_matriz /*matriz de salida*/);

//Ecuaciones oscilador
double fuerza (double posicion);
double damping (double momento, double nabla_dividido_m);
double g_xn_pn (double posicion, double momento, double nabla_dividido_m);
double f_pn (double momento);
void termino_estocastico_Z (double factor_estocastico, double *dos_terminos_estocasticos);

//Algoritmos
void verlet (double posicion, double momento, double nabla, double *array_posiciones);

//Datos para medir
double energia_cinetica (double momento);
double energia_potencial(double posicion);




//Variables globales para los numeros aleatorios
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;




/// -----Primera parte: OSCILADOR ARM�NICO----- ///
int main(){

    //variables
    double /*coef de viscosidad*/nabla, t_estancia_medio;
    double /*Pto inicial en el espacio de fases*/ momento_inicial=0, posicion_inicial=0;
    double t_estancia[500];
    int i, j, contador;
    FILE *f;
    // Inicializamos la rueda de n�meros random de Parisi-Rapuano
    ini_ran(123456789);

    for(i = 0; i < 500; i++){
        t_estancia[i] = 0;
    }
        

    ///EJECUTAMOS ALGORITMOS
    //Este ejecuta los algoritmos para un �nico valor de nabla
    #ifdef unico
        nabla = 1;
        verlet (posicion_inicial, momento_inicial, nabla, t_estancia);

        //comprobacion de t_estancia
        //for(i = 0; i < 500; i++){
        //printf("%f\n", t_estancia[i]);
        //}

        //contador es para saber cuantos valores validos (distintos de 0) tiene t_estancia
        contador = 0;

        for(i = 0; i < 500; i++){
            if(t_estancia[i] != 0){
                contador += 1;
            }
        }

        //comprobacion de contador
        //printf("contador = %i\n", contador);

        t_estancia_medio = 0;

        for(i = 0; i < contador; i++){
            t_estancia_medio += t_estancia[i];
        }

        t_estancia_medio = t_estancia_medio/contador;

        //comprobacion de t_estancia_medio
        //printf("t_estancia_medio = %f\n", t_estancia_medio);
        
        //nueva funcion histograma
        Histograma(0, 200000, contador, Intervalo, t_estancia);

    #endif // unico

    //Este bucle lo �nico que hace es que salgan varios archivos con distintos nabla de una, si se prefiere
    #ifdef bucle
        f = fopen("t_medio_frente_a_damping.txt", "w");
        nabla = 0.01;
        for (j = 0; j < 500; j++){
            nabla = nabla + 0.01*j;
            fprintf(f, "%f ", nabla);
            verlet (posicion_inicial, momento_inicial, nabla, t_estancia);

            contador = 0;

            for(i = 0; i < 500; i++){
                if(t_estancia[i] != 0){
                    contador += 1;
                }
            }

            //comprobacion de contador
            //printf("contador = %i\n", contador);

            t_estancia_medio = 0;

            for(i = 0; i < contador; i++){
                t_estancia_medio += t_estancia[i];
            }

            t_estancia_medio = t_estancia_medio/contador;

            //comprobacion de t_estancia_medio
            printf("t_estancia_medio = %f\n", t_estancia_medio);

            fprintf(f, "%f\n", t_estancia_medio);

            }
        fclose(f);
    #endif // bucle



    //aqui ejecuto las funciones que sacan los numeros gaussianos aleatorios
    //son solo de comprobaci�n as� que para ejecutar el programa de forma normal no se deben activar
    /*
    double vector_numeros_gaussianos[N], histograma_matriz[Intervalo];
    generador_vector_gaussiano(vector_numeros_gaussianos);
    histograma(vector_numeros_gaussianos, histograma_matriz);
    */

    return 0;
}


//Esta es la funcion histograma que use en compu. Creo que la normalizacion esta mal, pero por lo demas parece funcionar
void Histograma(int a, int b, int N, int divisiones,  double *array)
{

    int i,R,freq[divisiones];
    for (i=0; i<divisiones; i++)
    {
        freq[i]=0;
    }

    double w,anchura,Area;
    anchura=(b-a+0.0)/divisiones;
    Area=anchura*N;

    FILE*f;
    f=fopen ("Histograma.dat","w");
    if(f==NULL)
    {
        printf("ERROR");
        return 1;
    }
    for (i=0; i<N; i++)
    {
        w=array[i];
        R=(int)(w-a)/anchura;
        if (R >= divisiones){
             R = divisiones - 1;
        }
        freq[R]=freq[R]+1;
    }
    printf("La anchura es: %f El Area es: %f",anchura,Area);
    for (i=0; i<divisiones; i++)
    {
        fprintf(f,"%d %f\n",(int)(a+i*anchura),freq[i]/Area);
    }
    fclose(f);

    /*
    Entrada: L�mites intervalo (a y b, tipo int), n�mero de divisiones (divisiones tipo int), array con los valores de estudio (array tipo double) y su longitud (N, tipo int)
    Salida:Fichero con las frecuencias de cada intervalo.
    */
}








///N�meros aleatorios

//Generador de C de numeros random
//Mejor utilizar un Parisi-Rapuano
double Random_C (){
    double random;
    random=(rand()/((double)RAND_MAX+1));
    return random;

    //Inicializacion de numeros aleatorios de C, en principio no se usan
    //si se quieren usar COPIAR EN EL MAIN
    /*
    int seed=12131;
    srand(seed);
    */
}


//Parisi-Rapuano
//Esta funciones las guardo de cuando fisica computacional, ahora solo Dios sabe como chuchas funciona esta cosa
void ini_ran(int SEMILLA){
    int INI, FACTOR, SUM, i;
    srand(SEMILLA);
    INI = SEMILLA;
    FACTOR = 67397;
    SUM = 7364893;
    for (i = 0; i < 256; i++)
    {
        INI = (INI * FACTOR + SUM);
        irr[i] = INI;
    }
    ind_ran = ig1 = ig2 = ig3 = 0;
}

double Random(void){
    double r;
    ig1 = ind_ran - 24;
    ig2 = ind_ran - 55;
    ig3 = ind_ran - 61;
    irr[ind_ran] = irr[ig1] + irr[ig2];
    ir1 = (irr[ind_ran] ^ irr[ig3]);
    ind_ran++;
    r = ir1 * NormRANu;
    // printf("r=%f\n",r);
    return r;
}




//Genera un n�mero aleatorio con distribuci�n gaussiana, a partir de un numero random en el intervalo [0,1)
void num_aleatorio_gaussiano (double *dos_numeros_gaussianos){
    double aleatorio_uniforme_1, aleatorio_uniforme_2, auxiliar_1, auxiliar_2;

    ///Aqu� hay un problema con el algoritmo: Cuando sale aleatorio_uniforme_1 = 0 ---> ln(0)=-inf; y da error
    //No se como se supone que deber�amos arreglarlo, de momento solo voy a poner una cl�usula de que no sea igual a cero
    aleatorio_uniforme_1=Random ();
    while (aleatorio_uniforme_1 == 0.0)
        aleatorio_uniforme_1=Random ();
    aleatorio_uniforme_2=Random ();

    auxiliar_1=sqrt(-2*log(aleatorio_uniforme_1));
    auxiliar_2=2*PI*aleatorio_uniforme_2;

    //en la presentaci�n dan como dos posibilidades, seg�n las pruebas que he hecho es indistinto usar una u otra
    dos_numeros_gaussianos[0]= auxiliar_1*cos(auxiliar_2);

    //Esta ser�a la segunda forma, solo cambia el cos por el sen
    dos_numeros_gaussianos[1] = auxiliar_1*sin(auxiliar_2);
}



//Esta funci�n es para comprobar que la funcion de numeros gaussianos rulaba bien, probablemente se quitar� en un futuro
void generador_vector_gaussiano (/*vector de salida*/ double *vector_numeros_gaussianos){
    int i, N;

    double dos_numeros_gausianos[2];
    FILE* f;
    f = fopen("gaussiana.txt", "w");

    for (i = 0 ; i < N; i += 2){
        num_aleatorio_gaussiano(dos_numeros_gausianos);
        vector_numeros_gaussianos[i] = dos_numeros_gausianos[0];
        vector_numeros_gaussianos[i+1] = dos_numeros_gausianos[1];
        fprintf(f,"%f\n", vector_numeros_gaussianos[i]);
        fprintf(f,"%f\n", vector_numeros_gaussianos[i+1]);
    }
    fclose(f);
}


//Lo mismo, esta funci�n es para comprobar que la funcion de numeros gaussianos rulaba bien, probablemente se quitar� en un futuro
//Aunque, tambi�n nos sirve para lo de comprobar que la distribuci�n de posiciones y velocidades sea gaussiana
void histograma (double *V/*matriz de entrada*/, double *histograma_matriz /*matriz de salida*/){
    //Variables matriz histograma
    int asignacion, i, N;
    double anchura, min_desconocido, max_desconocido;

    N = sizeof(V);

    max_desconocido=-1000;
    min_desconocido=1000;
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

//gr�fico en gnuplot
    FILE* f;
    f = fopen("histo.txt", "w");
    int eje_x;
    eje_x=(int)(Intervalo/2);
    for (i=0;i<Intervalo;i++){
        fprintf(f,"%d %f\n", i-eje_x, histograma_matriz[i]);
    }
    fclose(f);
    //printf("%f\n", V[99]);
//comprobaci�n de que las cosas vayan bien xd
    //printf("%f\n", max_desconocido);
}






///     ECUACIONES OSCILADOR      ///
//ecuacion_oscilador_posicion_punto, en los algoritmos esta suele salir como f(p_n)
double f_pn (double momento){
    double x_punto;
    x_punto=momento/m;
    return x_punto;
}


#ifdef oscilador_armonico
///IMPORTANTE, esto es (la menos derivada de) el potencial que tendremos que cambiar en la siguiente parte as� que ojito
double fuerza (double posicion){
    double grad_Vx;
    grad_Vx=-k*posicion;
    return grad_Vx;
}
#endif // oscilador_armonico

#ifdef doble_pozo
//Fuerza del doble pozo, la he llamado solo fuerza para no tener que cambiar el resto de funciones
double fuerza (double posicion){
    double grad_Vx;
    grad_Vx=-A*4*posicion*(posicion*posicion-1);
    return grad_Vx;
}
#endif // doble_pozo

//T�rmino de rozamiento o damping
double damping (double momento, double nabla_dividido_m){
    double rozamiento;
    rozamiento = -nabla_dividido_m*momento;
    return rozamiento;
}



//ecuacion_oscilador_momento_punto, calcula la parte NO ESTOCASTICA de p_punto, en los algoritmos esta suele salir como g(x_n, p_n)
double g_xn_pn (double posicion, double momento, double nabla_dividido_m){
    double p_punto;  //aqu� lo llamo p_punto, pero realmente no lo es porque le falta el t�rmino estoc�stico
    p_punto = damping (momento, nabla_dividido_m) + fuerza (posicion);
    return p_punto;
}



//este es el t�rmino estoc�stico, que va incluido en la ecuaci�n del momento
//Saca 2 para aprovechar mejor la funcion generadora de numeros gaussianos
void termino_estocastico_Z (double factor_estocastico, double *dos_terminos_estocasticos){
    double dos_numeros_gaussianos[2];
    num_aleatorio_gaussiano(dos_numeros_gaussianos);
    //Para que no tenga que pararse cada vez a calcular, el factor = sqrt(2*nabla*K_b_T*h)
    dos_terminos_estocasticos[0] = factor_estocastico*dos_numeros_gaussianos[0];
    dos_terminos_estocasticos[1] = factor_estocastico*dos_numeros_gaussianos[1];
}






///     ALGORITMOS     ///


void verlet (double posicion, double momento, double nabla, double *t_estancia){
    int i, j, l, pasos, contador_pos, contador_neg;
    double mediapotencial,mediacinetica;
    //Uso la notaci�n del power point de clase, las barras bajas se deben entender como "sub", p.ej f sub x1
    double factor_estocastico, dos_terminos_estocasticos[2], a, b, factor_posicion, fuerza_ahora, h_medios, Z, posicion_anterior;

    //Toda la parafernalia de los char es para poder sacar todos los archivos de distintas h y nabla en un solo bucle
    char nombre_archivo[1023]="Verlet_h=", nombre_ocupacion[1023]="Ocupacion_Verlet_h=", especificador_h[50], especificador_nabla[50], nombre_archivo_parte2[]="_nabla=", nombre_archivo_fin[]=".txt";
    sprintf(especificador_h, "%f", h);
    strncat(nombre_archivo, especificador_h, 1024);
    sprintf(especificador_nabla, "%f", nabla);
    strncat(nombre_archivo, nombre_archivo_parte2, 1024);
    strncat(nombre_archivo, especificador_nabla, 1024);
    strncat(nombre_archivo, nombre_archivo_fin, 1024);

    strncat(nombre_ocupacion, especificador_h, 1024);
    strncat(nombre_ocupacion, nombre_archivo_parte2, 1024);
    strncat(nombre_ocupacion, especificador_nabla, 1024);
    strncat(nombre_ocupacion, nombre_archivo_fin, 1024);

    FILE *f, *g;
    f = fopen(nombre_archivo, "w");
    g = fopen(nombre_ocupacion, "w");

    pasos = (int)(T/h);
    factor_estocastico = sqrt(2*nabla*K_b_T*h);
    h_medios = h/2;
    b = 1/(1+nabla*h_medios/m);
    a = 2*b-1;
    factor_posicion = b*h_medios/m;
    mediacinetica=0;
    mediapotencial=0;
    contador_pos = 0;
    contador_neg = 0;
    l = 0;
    posicion_anterior = posicion;

    for (i = 0 ; i < pasos; i += 2){

        termino_estocastico_Z(factor_estocastico, dos_terminos_estocasticos);


        for(j=0;j<2;j++){
        //Paso 1 de tiempo
            fuerza_ahora = fuerza(posicion);
            ///AQUI NO TENGO CLARO SI HAY QUE MULTIPLICAR EL GAUSSIANO POR h
            Z = dos_terminos_estocasticos[j];

            posicion += factor_posicion*(2*momento+h*fuerza_ahora+Z);
            momento = a*momento+h_medios*(a*fuerza_ahora+fuerza(posicion))+b*Z;

            //tiempo de estancia (izquierda)
            if(posicion < 0){
                if(posicion_anterior > 0){
                    l = l+1;
                }
                t_estancia[l] += 1;
            }
            posicion_anterior = posicion;
            fprintf(f,"%f ", (i+j)*h); //Esto es el tiempo
            fprintf(f,"%f ", posicion);
            fprintf(f,"%f ", momento);

            #ifdef doble_pozo
                //ocupacion (normalizada) de la particula en el doble pozo, +1 para x>0 y -1 para x<0
                if(posicion > 0){
                    contador_pos += 1;
                    fprintf(g, "1 ");
                    fprintf(g, "%f\n", 1.0*contador_pos/pasos);
                }

                if(posicion < 0){
                    contador_neg += 1;
                    fprintf(g, "-1 ");
                    fprintf(g, "%f\n", 1.0*contador_neg/pasos);
                }
            #endif // doble_pozo


            #ifdef true
                mediacinetica+=energia_cinetica(momento);
                mediapotencial+=energia_potencial(posicion);
                fprintf(f,"%f ", mediacinetica/(i+j+1));
                fprintf(f,"%f ", mediapotencial/(i+j+1));
                fprintf(f,"%f\n", (mediapotencial+mediacinetica)/(i+j+1));

                /*fprintf(f,"%f ", energia_cinetica(momento));
                fprintf(f,"%f\n", energia_potencial(posicion));*/

            #else
                fprintf(f,"\n");

            #endif // true
        }


/*
        //Paso 2 de tiempo
        fuerza_ahora = fuerza(posicion);
        ///AQUI NO TENGO CLARO SI HAY QUE MULTIPLICAR EL GAUSSIANO POR h
        Z = dos_terminos_estocasticos[1];

        posicion += factor_posicion*(2*momento+h*fuerza_ahora+Z);
        momento = a*momento+h_medios*(a*fuerza_ahora+fuerza(posicion))+b*Z;

        fprintf(f,"%f ", (i+1)*h); //Esto es el tiempo
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);


        #ifdef true

            fprintf(f,"%f ", energia_cinetica(momento));
            fprintf(f,"%f\n", energia_potencial(posicion));

        #else
            fprintf(f,"\n");

        #endif // true*/

    }
    fclose(f);
    fclose(g);
}





///Datos a medir
double energia_cinetica (double momento){
    double Ec;
    Ec = momento*momento*0.5/m;
    return Ec;
}


double energia_potencial(double posicion){
    double Ep;
    Ep = 0.5*posicion*posicion*k;
    return Ep;
}
