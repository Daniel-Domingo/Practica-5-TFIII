///     -------Practica 5 TF3, sim molecular------      ///
       /// -----Segunda parte: DOBLE POZO----- ///


//actualizado el 11-11-24




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define m 1     //masa
#define k 1     //cte el�stica
#define K_b_T 0.2 //cte de Boltzmann por temperatura
#define PI 3.14159265358979323846


//para los numerros random, NO TOCAR
#define NormRANu (2.3283063671E-10F)

#define T 1000000 //Tiempo en cualquier algoritmo
#define h 0.001

#define Intervalo 100 //Numero de intervalos en el histograma


///Estos �ltimos son para los ifdef:

//si se pone en "bucle" saca varios archivos de una para distintos h y nabla (en principio todos los que piden)
//si se pone en "unico" saca un solo archivo para un dado h y nabla
#define bucle //unico


//Calcular las energ�as durante el proceso, si (energia) o no (cualquier otro)
#define false

//Pasar las posiciones y momentos a fichero, si(esp_de_fases) o no (cualquier otro)
#define no

//Para elegir el tipo de potencial que se utiliza (doble_pozo) o (pozo)
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
void termino_estocastico_Z (double factor_estocastico, double *dos_terminos_estocasticos);

//Algoritmos
void verlet (double posicion, double momento, double nabla, int *cuenta_saltos);

//Datos para medir
double energia_cinetica (double momento);
double energia_potencial(double posicion);









//Como la altura A va a cambiar es mejor tenerlo como variable global
double A = 1;


//Variables globales para los numeros aleatorios
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;




int main(){

    //variables
    double /*coef de viscosidad*/nabla;
    double /*Pto inicial en el espacio de fases*/ momento_inicial=0, posicion_inicial=0.00001, factor_bucle_nabla, tiempo_est_medio;
    int i, j, pasos_nabla = 10, cuenta_saltos;
    FILE *f;
    // Inicializamos la rueda de n�meros random de Parisi-Rapuano
    ini_ran(123456789);


    ///EJECUTAMOS ALGORITMOS
    //Este ejecuta los algoritmos para un �nico valor de nabla
    #ifdef unico ///LA OPCIÓN DE ÚNICO AHORA MISMO NO FUNCIONA, no por nada, simplemente no la he arreglado
        nabla = 1;

        for(i = 0; i < 500; i++){
            t_estancia[i] = 0;
        }


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
        fprintf(f, "#nabla    #Tiempo de estancia medio\n");

        for (i = 0; i < 2; i++){
            nabla = 0;
            if (i == 0){
                factor_bucle_nabla = 0.5;
                pasos_nabla = 10;
            }
            else{
                factor_bucle_nabla = 0.01;
                pasos_nabla = 99;
                }
            for (j = 0; j < pasos_nabla; j++){
                
                
                
                ///APAÑO MALO, QUITAR LUEGO
                if (j == 50)
                    j++;
                
                
                
                nabla += factor_bucle_nabla;

                //fprintf(f, "%f ", nabla);


                ///Ejecutamos Verlet
                verlet (posicion_inicial, momento_inicial, nabla, &cuenta_saltos);


                tiempo_est_medio = T/(cuenta_saltos+1);
                fprintf(f, "%f      %f\n", nabla, tiempo_est_medio);



                printf("nabla=%f: saltos=%d\n", nabla, cuenta_saltos);

            }
        }
        fclose(f);
    #endif // bucle




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
        return;
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

///IMPORTANTE, esto es (la menos derivada de) el potencial que tendremos que cambiar en la siguiente parte as� que ojito
double fuerza (double posicion){
    double grad_Vx;


    #ifdef oscilador_armonico
        grad_Vx=-k*posicion;
    #endif // oscilador_armonico

    #ifdef doble_pozo
        grad_Vx=-A*4*posicion*(posicion*posicion-1);
    #endif // doble_pozo


    return grad_Vx;
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






///     ALGORITMO VERLET     ///


void verlet (double posicion, double momento, double nabla, int *cuenta_saltos){
    int i, j, pasos, contador_pos, contador_neg, t_estancia;
    double factor_estocastico, dos_terminos_estocasticos[2], a, b, factor_posicion, fuerza_ahora, h_medios, Z, posicion_anterior;

    
    #ifdef energia
        double mediapotencial,mediacinetica;
        mediacinetica=0;
        mediapotencial=0;
    #endif // energia

/*
He cambiado el nombre de los archivos de salida
Como vamos a estar usando siempre Verlet y h=0.001 (igual la cambiamos pero en ese caso sería siempre esa) no tiene sentido ponerlo en el nombre del archivo. 
Lo que si tiene más sentido es poner A que vamos a tener que cambiarlo y T que hubiera tenido sentido ponerlo desde el pozo simple xd.
nabla/eta se queda como estaba. 
Ahora los archivos de posicion/momento/energías se llaman "Doble_pozo_(...)"; los de tiempo en cada estado "Ocupacion_(...)".
*/

    //Toda la parafernalia de los char es para poder sacar todos los archivos de distintas h y nabla en un solo bucle
    char nombre_archivo[1023]="Doble_pozo_A=", nombre_ocupacion[1023]="Ocupacion_A=", especificador_A[50], especificador_nabla[50], especificador_T[50], nombre_archivo_parte2[]="_nabla=", nombre_archivo_parte3[]="_T=", nombre_archivo_fin[]=".txt";
    sprintf(especificador_A, "%f", A);
    sprintf(especificador_nabla, "%f", nabla);
    sprintf(especificador_T, "%d", T);
    
    strncat(nombre_archivo, especificador_A, 1024);
    strncat(nombre_archivo, nombre_archivo_parte2, 1024);
    strncat(nombre_archivo, especificador_nabla, 1024);
    strncat(nombre_archivo, nombre_archivo_parte3, 1024);
    strncat(nombre_archivo, especificador_T, 1024);
    strncat(nombre_archivo, nombre_archivo_fin, 1024);

    strncat(nombre_ocupacion, especificador_A, 1024);
    strncat(nombre_ocupacion, nombre_archivo_parte2, 1024);
    strncat(nombre_ocupacion, especificador_nabla, 1024);
    strncat(nombre_archivo, nombre_archivo_parte3, 1024);
    strncat(nombre_archivo, especificador_T, 1024);
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


    contador_pos = 0;
    contador_neg = 0;
    posicion_anterior = posicion;
    t_estancia = 0;

    unsigned int control_t_estancia_columna;
    control_t_estancia_columna = 0;

    *cuenta_saltos = 0;


    int pasos_no_representados, cada_cuantos_pasos_representamos;
    cada_cuantos_pasos_representamos = 50; //Representamos uno de (cada_cuantos_pasos_representamos) pasos
    pasos_no_representados = 0; //Contador de los pasos a representar


    for (i = 0 ; i < pasos; i += 2){

        termino_estocastico_Z(factor_estocastico, dos_terminos_estocasticos);


        for(j=0;j<2;j++){
            fuerza_ahora = fuerza(posicion);
            Z = dos_terminos_estocasticos[j];

            posicion += factor_posicion*(2*momento+h*fuerza_ahora+Z);
            momento = a*momento+h_medios*(a*fuerza_ahora+fuerza(posicion))+b*Z ;




            //tiempo de estancia; como posicion_inicial > 0, la primera columna es el tiempo a la derecha (en x positivo)
            if(posicion*posicion_anterior > 0){
                t_estancia++;
            }
            if(posicion*posicion_anterior < 0){
                fprintf(g, "%f        ", t_estancia*h);
                t_estancia = 0;
                (*cuenta_saltos)++;
                control_t_estancia_columna++;
                if (control_t_estancia_columna == 2){
                    fprintf(g, "\n");
                    control_t_estancia_columna = 0;
                }
            }


            posicion_anterior = posicion;







            pasos_no_representados += 1;


            #ifdef esp_de_fases
                if (pasos_no_representados == cada_cuantos_pasos_representamos){
                    fprintf(f,"%f ", (i+j+1)*h); //Esto es el tiempo
                    fprintf(f,"%f ", posicion);
                    fprintf(f,"%f ", momento);
                }
            #endif // esp_de_fases

            //Esto es necesario para representar la energía y el espacio de fases a la vez


            #ifdef doble_pozo
                //ocupacion (normalizada) de la particula en el doble pozo, +1 para x>0 y -1 para x<0
                if(posicion > 0)
                    contador_pos += 1;
                if(posicion < 0)
                    contador_neg += 1;
            #endif // doble_pozo


            #ifdef energia
                mediacinetica += energia_cinetica(momento);
                mediapotencial += energia_potencial(posicion);

                if (pasos_no_representados == cada_cuantos_pasos_representamos){
                    #ifdef esp_de_fases

                    #else
                        fprintf(f,"%f ", (i+j)*h); //Esto es el tiempo
                    #endif // esp_de_fases
                    fprintf(f,"%f ", (i+j)*h); //Esto es el tiempo
                    fprintf(f,"%f ", mediacinetica/(i+j+1));
                    fprintf(f,"%f ", mediapotencial/(i+j+1));
                    fprintf(f,"%f\n", (mediapotencial+mediacinetica)/(i+j+1));
                }

            #else
                if (pasos_no_representados == cada_cuantos_pasos_representamos)
                    fprintf(f,"\n");

                if (pasos_no_representados == cada_cuantos_pasos_representamos)
                    pasos_no_representados = 0;
            #endif // energia


        }

    }


    #ifdef esp_de_fases
        fprintf(f,"%f ", (i+j)*h); //Esto es el tiempo
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);
    #endif // esp_de_fases

    #ifdef energia
        if (pasos_no_representados == cada_cuantos_pasos_representamos){
            fprintf(f,"%f ", (i+j+1)*h); //Esto es el tiempo
            fprintf(f,"%f ", mediacinetica/(i+j+1));
            fprintf(f,"%f ", mediapotencial/(i+j+1));
            fprintf(f,"%f\n", (mediapotencial+mediacinetica)/(i+j+1));
        }
        #endif // true



    fclose(g);
    fclose(f);
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
