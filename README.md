///     -------Practica 5 TF3, sim molecular------      ///





#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define m 1     //masa
#define k 1     //cte elástica
#define K_b_T 1 //cte de Boltzmann por temperatura
#define PI 3.14159265358979323846

//para los numerros random, NO TOCAR
#define NormRANu (2.3283063671E-10F)

#define T 500 //Tiempo en cualquier algoritmo
//El numero de pasos vendrá dado por T/h, así pasa el mismo tiempo independientemente de h

#define Ng 10000 //Numero de datos gaussianos que se quieren, es solo para comprobaciones, SIEMPRE PAR
#define Intervalo 100 //Numero de intervalos en el histograma


///Estos últimos son para los ifdef:

//si se pone en "bucle" saca varios archivos de una para distintos h y nabla (en principio todos los que piden)
//si se pone en "unico" saca un solo archivo para un dado h y nabla
#define unico //bucle


//Calcular las energías durante el proceso, si (true) o no (false)
#define true





//Numeros aleatorios uniformes
double Random_C();
void ini_ran(int SEMILLA);
double Random(void);
void num_aleatorio_gaussiano (double *dos_numeros_gaussianos);

//Numeros aleatorios gaussianos
void generador_vector_gaussiano (/*vector de salida*/ double *vector_numeros_gaussianos);
void histograma (double *V/*matriz de entrada*/, double *histograma_matriz /*matriz de salida*/);

//Ecuaciones oscilador
double fuerza (double posicion);
double damping (double momento, double nabla_dividido_m);
double g_xn_pn (double posicion, double momento, double nabla_dividido_m);
double f_pn (double momento);
void termino_estocastico_Z (double factor_estocastico, double *dos_terminos_estocasticos);

//Algoritmos
void Euler_Maruyama (double posicion, double momento, double h, double nabla);
void Runge_Kutta (double posicion, double momento, double h, double nabla);
void verlet (double posicion, double momento, double h, double nabla);

//Datos para medir
double energia_cinetica (double momento);
double energia_potencial(double posicion);




//Variables globales para los numeros aleatorios
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;











/// -----Primera parte: OSCILADOR ARMÓNICO----- ///
int main(){

    //variables
    double /*coef de viscosidad*/nabla, /*paso de tiempo*/h;
    double /*Pto inicial en el espacio de fases*/ momento_inicial=0, posision_inicial=0;



    // Inicializamos la rueda de números random de Parisi-Rapuano
    ini_ran(123456789);



    ///EJECUTAMOS ALGORITMOS
    //Este ejecuta los algoritmos para un único valor de h y nabla
    #ifdef unico
        h = 0.001;
        nabla = 0.1;
        Euler_Maruyama(posision_inicial, momento_inicial, h, nabla);
        Runge_Kutta (posision_inicial, momento_inicial, h, nabla);
        verlet (posision_inicial, momento_inicial, h, nabla);
    #endif // unico

    //Este bucle lo único que hace es que salgan varios archivos con distintos h y nabla de una, si se prefiere
    #ifdef bucle
        int i, j;
        h = 1;
        for (i = 0; i < 4; i++){
            h = h/10;
            nabla = 0.01;
            for (j = 0; j < 3; j++){
                nabla = nabla*10;
                Euler_Maruyama(posision_inicial, momento_inicial, h, nabla);
                Runge_Kutta (posision_inicial, momento_inicial, h, nabla);
                verlet (posision_inicial, momento_inicial, h, nabla);
            }
        }
    #endif // bucle


    ///Comprobacion de que la distribucion de posiciones y velocidades sea gaussiana
    //histograma()




    //aqui ejecuto las funciones que sacan los numeros gaussianos aleatorios
    //son solo de comprobación así que para ejecutar el programa de forma normal no se deben activar
    /*
    double vector_numeros_gaussianos[N], histograma_matriz[Intervalo];
    generador_vector_gaussiano(vector_numeros_gaussianos);
    histograma(vector_numeros_gaussianos, histograma_matriz);
    */

    return 0;
}











///Números aleatorios

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




//Genera un número aleatorio con distribución gaussiana, a partir de un numero random en el intervalo [0,1)
void num_aleatorio_gaussiano (double *dos_numeros_gaussianos){
    double aleatorio_uniforme_1, aleatorio_uniforme_2, auxiliar_1, auxiliar_2;

    ///Aquí hay un problema con el algoritmo: Cuando sale aleatorio_uniforme_1 = 0 ---> ln(0)=-inf; y da error
    //No se como se supone que deberíamos arreglarlo, de momento solo voy a poner una cláusula de que no sea igual a cero
    aleatorio_uniforme_1=Random ();
    while (aleatorio_uniforme_1 == 0.0)
        aleatorio_uniforme_1=Random ();
    aleatorio_uniforme_2=Random ();

    auxiliar_1=sqrt(-2*log(aleatorio_uniforme_1));
    auxiliar_2=2*PI*aleatorio_uniforme_2;

    //en la presentación dan como dos posibilidades, según las pruebas que he hecho es indistinto usar una u otra
    dos_numeros_gaussianos[0]= auxiliar_1*cos(auxiliar_2);

    //Esta sería la segunda forma, solo cambia el cos por el sen
    dos_numeros_gaussianos[1] = auxiliar_1*sin(auxiliar_2);
}



//Esta función es para comprobar que la funcion de numeros gaussianos rulaba bien, probablemente se quitará en un futuro
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


//Lo mismo, esta función es para comprobar que la funcion de numeros gaussianos rulaba bien, probablemente se quitará en un futuro
//Aunque, también nos sirve para lo de comprobar que la distribución de posiciones y velocidades sea gaussiana
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








///     ECUACIONES OSCILADOR      ///
//ecuacion_oscilador_posicion_punto, en los algoritmos esta suele salir como f(p_n)
double f_pn (double momento){
    double x_punto;
    x_punto=momento/m;
    return x_punto;
}



///IMPORTANTE, esto es (la menos derivada de) el potencial que tendremos que cambiar en la siguiente parte así que ojito
double fuerza (double posicion){
    double grad_Vx;
    grad_Vx=-k*posicion;
    return grad_Vx;
}


//Término de rozamiento o damping
double damping (double momento, double nabla_dividido_m){
    double rozamiento;
    rozamiento = -nabla_dividido_m*momento;
    return rozamiento;
}



//ecuacion_oscilador_momento_punto, calcula la parte NO ESTOCASTICA de p_punto, en los algoritmos esta suele salir como g(x_n, p_n)
double g_xn_pn (double posicion, double momento, double nabla_dividido_m){
    double p_punto;  //aquí lo llamo p_punto, pero realmente no lo es porque le falta el término estocástico
    p_punto = damping (momento, nabla_dividido_m) + fuerza (posicion);
    return p_punto;
}



//este es el término estocástico, que va incluido en la ecuación del momento
//Saca 2 para aprovechar mejor la funcion generadora de numeros gaussianos
void termino_estocastico_Z (double factor_estocastico, double *dos_terminos_estocasticos){
    double dos_numeros_gaussianos[2];
    num_aleatorio_gaussiano(dos_numeros_gaussianos);
    //Para que no tenga que pararse cada vez a calcular, el factor = sqrt(2*nabla*K_b_T*h)
    dos_terminos_estocasticos[0] = factor_estocastico*dos_numeros_gaussianos[0];
    dos_terminos_estocasticos[1] = factor_estocastico*dos_numeros_gaussianos[1];
}








///     ALGORITMOS     ///

//Primer algoritmo: Euler-Maruyama
void Euler_Maruyama (double posicion, double momento, double h, double nabla){
    int i, pasos;
    double factor_estocastico, dos_terminos_estocasticos[2], nabla_dividido_m;

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

    //El algoritmo en verdad son solo estas 5 líneas
    pasos = (int)(T/h);
    factor_estocastico = sqrt(2*nabla*K_b_T*h);
    nabla_dividido_m = nabla/m;

    for (i = 0 ; i < pasos; i += 2){
        termino_estocastico_Z(factor_estocastico, dos_terminos_estocasticos);
        //En cada paso salen 2 por lo de tener dos numeros gaussianos
        posicion += f_pn(momento)*h;
        momento += g_xn_pn(posicion, momento, nabla_dividido_m)*h + dos_terminos_estocasticos[0];

        fprintf(f,"%f ", i*h); //Esto es el tiempo
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);

        #ifdef true
            fprintf(f,"%f ", energia_cinetica(momento));
            fprintf(f,"%f\n", energia_potencial(posicion));

        #else
            fprintf(f,"\n");

        #endif // true


        //Segundo paso de tiempo
        posicion += f_pn(momento)*h;
        momento += g_xn_pn(posicion, momento, nabla_dividido_m)*h + dos_terminos_estocasticos[1];

        fprintf(f,"%f ", (i+1)*h); //Esto es el tiempo
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);

        #ifdef true
            fprintf(f,"%f ", energia_cinetica(momento));
            fprintf(f,"%f\n", energia_potencial(posicion));

        #else
            fprintf(f,"\n");

        #endif // true
    }
    fclose(f);
}




void Runge_Kutta (double posicion, double momento, double h, double nabla){
    int i, pasos;
    //Uso la notación del power point de clase, las barras bajas se deben entender como "sub", p.ej f sub x1
    double f_x1, f_x2, g_p1, g_p2, Z, h_medios, factor_estocastico, dos_terminos_estocasticos[2], nabla_dividido_m;

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

    pasos = (int)(T/h);
    h_medios = h*0.5;
    factor_estocastico = sqrt(2*nabla*K_b_T*h);
    nabla_dividido_m = nabla/m;

    for (i = 0 ; i < pasos; i += 2){

        termino_estocastico_Z(factor_estocastico, dos_terminos_estocasticos);



        //De nuevo, de cada paso de tiempo hay que hacer 2
        Z = dos_terminos_estocasticos[0];

        //Paso 1
        f_x1=f_pn(momento+Z);
        g_p1=g_xn_pn(posicion, momento+Z, nabla_dividido_m);

        //Paso 2
        f_x2=f_pn(momento+h*g_p1);
        g_p2=g_xn_pn(posicion+h*f_x1, momento+g_p1, nabla_dividido_m);

        //Final, calculo de nuevo punto en el espacio de fases
        posicion += h_medios*(f_x1+f_x2);
        momento += h_medios*(g_p1+g_p2)+Z;

        fprintf(f,"%f ", i*h); //Esto es el tiempo
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);

        #ifdef true
            fprintf(f,"%f ", energia_cinetica(momento));
            fprintf(f,"%f\n", energia_potencial(posicion));

        #else
            fprintf(f,"\n");

        #endif // true



        //Segundo paso de tiempo (no confundir con los pasos internos de RK)
        Z = dos_terminos_estocasticos[1];

        //Paso 1
        f_x1=f_pn(momento+Z);
        g_p1=g_xn_pn(posicion, momento+Z, nabla_dividido_m);

        //Paso 2
        f_x2=f_pn(momento+h*g_p1);
        g_p2=g_xn_pn(posicion+h*f_x1, momento+g_p1, nabla_dividido_m);

        //Final, calculo de nuevo punto en el espacio de fases
        posicion += h_medios*(f_x1+f_x2);
        momento += h_medios*(g_p1+g_p2)+Z;

        fprintf(f,"%f ", (i+1)*h); //Esto es el tiempo
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);

        #ifdef true
            fprintf(f,"%f ", energia_cinetica(momento));
            fprintf(f,"%f\n", energia_potencial(posicion));

        #else
            fprintf(f,"\n");

        #endif // true
    }
    fclose(f);
}




void verlet (double posicion, double momento, double h, double nabla){
    int i, pasos;
    //Uso la notación del power point de clase, las barras bajas se deben entender como "sub", p.ej f sub x1
    double factor_estocastico, dos_terminos_estocasticos[2], a, b, factor_posicion, fuerza_ahora, h_medios, Z;

    //Toda la parafernalia de los char es para poder sacar todos los archivos de distintas h y nabla en un solo bucle
    char nombre_archivo[1023]="Verlet_h=", especificador_h[50], especificador_nabla[50], nombre_archivo_parte2[]="_nabla=", nombre_archivo_fin[]=".txt";
    sprintf(especificador_h, "%f", h);
    strncat(nombre_archivo, especificador_h, 1024);
    sprintf(especificador_nabla, "%f", nabla);
    strncat(nombre_archivo, nombre_archivo_parte2, 1024);
    strncat(nombre_archivo, especificador_nabla, 1024);
    strncat(nombre_archivo, nombre_archivo_fin, 1024);


    FILE *f;
    f = fopen(nombre_archivo, "w");

    pasos = (int)(T/h);
    factor_estocastico = sqrt(2*nabla*K_b_T*h);
    h_medios = h/2;
    b = 1/(1+nabla*h_medios/m);
    a = 2*b-1;
    factor_posicion = b*h_medios/m;

    for (i = 0 ; i < pasos; i += 2){

        termino_estocastico_Z(factor_estocastico, dos_terminos_estocasticos);



        //Paso 1 de tiempo
        fuerza_ahora = fuerza(posicion);
        ///AQUI NO TENGO CLARO SI HAY QUE MULTIPLICAR EL GAUSSIANO POR h
        Z = dos_terminos_estocasticos[0];

        posicion += factor_posicion*(2*momento+h*fuerza_ahora+Z);
        momento = a*momento+h_medios*(a*fuerza_ahora+fuerza(posicion))+b*Z;

        fprintf(f,"%f ", i*h); //Esto es el tiempo
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);

        #ifdef true
            fprintf(f,"%f ", energia_cinetica(momento));
            fprintf(f,"%f\n", energia_potencial(posicion));

        #else
            fprintf(f,"\n");

        #endif // true



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

        #endif // true

    }
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


