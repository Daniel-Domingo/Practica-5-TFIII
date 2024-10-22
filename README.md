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

#define T 5000 //Tiempo en cualquier algoritmo
//El numero de pasos vendrá dado por T/h, así pasa el mismo tiempo independientemente de h

//Numero de datos gaussianos que se quieren, es solo para comprobaciones, SIEMPRE PAR
#define Ng 200000 //Para más de 250000 la memoria peta, habria que hacerlo con un malloc si se necesita

#define Intervalo 150 //Numero de intervalos en el histograma


///Estos últimos son para los ifdef:

//si se pone en "bucle" saca varios archivos de una para distintos h y nabla (en principio todos los que piden)
//si se pone en "unico" saca un solo archivo para un dado h y nabla
#define bucle //bucle

//Calcular las energías durante el proceso, si (energia) o no (cualquier otro)
#define energia

//Sacar histogramas de las distribuciones de posiciones y momentos, si (comprobar_dist) o no (cualquier otro)
#define otro

//Sacar una comprobacion del generador de numeros gaussianos si (gauss) o no (otro)
#define otro





//Numeros aleatorios uniformes
double Random_C();
void ini_ran(int SEMILLA);
double Random(void);

//Numeros aleatorios gaussianos
void num_aleatorio_gaussiano (double *dos_numeros_gaussianos);
void generador_vector_gaussiano (/*vector de salida*/ double *vector_numeros_gaussianos);
void calcular_media_desviacion(double *vector, int N, double *media, double *desviacion);
void histograma (double *V/*matriz de entrada*/, int num_elementos, const char *nombre_archivo);

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

#ifdef energia
    //Datos para medir
    double energia_cinetica (double momento);
    double energia_potencial(double posicion);
#endif // energia



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
    ini_ran(3456789);




    ///EJECUTAMOS ALGORITMOS
    //Este ejecuta los algoritmos para un único valor de h y nabla
    #ifdef unico
        h = 0.01;
        nabla = 1;
        Euler_Maruyama(posision_inicial, momento_inicial, h, nabla);
        Runge_Kutta (posision_inicial, momento_inicial, h, nabla);
        verlet (posision_inicial, momento_inicial, h, nabla);
    #endif // unico

    //Este bucle lo único que hace es que salgan varios archivos con distintos h y nabla de una, si se prefiere
    #ifdef bucle

        //Reiniciamos el fichero de las energías finales
        #ifdef energia
            FILE *f;
            f = fopen("tabla_energias", "w");
            fprintf(f, "");
        #endif // energia


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
            printf("finalizado h = %f\n", h);
        }
    #endif // bucle





    ///aqui ejecuto las funciones que sacan los numeros gaussianos aleatorios
    //son solo de comprobación así que para ejecutar el programa de forma normal no se deben activar
    #ifdef gauss
        double vector_numeros_gaussianos[Ng], media, desviacion;
        int N = Ng;
        generador_vector_gaussiano(vector_numeros_gaussianos);
        calcular_media_desviacion(vector_numeros_gaussianos, N, &media, &desviacion);
        histograma(vector_numeros_gaussianos, N, "comprobar gaussiana.txt");
    #endif // gauss



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




//Genera dos números aleatorios con distribución gaussiana, a partir de dos numeros random en el intervalo [0,1)
void num_aleatorio_gaussiano (double *dos_numeros_gaussianos){
    double aleatorio_uniforme_1, aleatorio_uniforme_2, auxiliar_1, auxiliar_2;

    ///Aquí hay un problema con el algoritmo: Cuando sale aleatorio_uniforme_1 = 0 ---> ln(0)=-inf; y da error
    //No se como se supone que deberíamos arreglarlo, de momento solo voy a poner una cláusula de que no sea igual a cero
    aleatorio_uniforme_1=Random ();
    while (aleatorio_uniforme_1 <= 1e-10)
        aleatorio_uniforme_1=Random ();
    aleatorio_uniforme_2=Random ();

    auxiliar_1=sqrt(-2*log(aleatorio_uniforme_1));
    auxiliar_2=2*PI*aleatorio_uniforme_2;

    dos_numeros_gaussianos[0]= auxiliar_1*cos(auxiliar_2);
    dos_numeros_gaussianos[1] = auxiliar_1*sin(auxiliar_2);
}



//Esta función es para comprobar que la funcion de numeros gaussianos rulaba bien, probablemente se quitará en un futuro
void generador_vector_gaussiano (/*vector de salida*/ double *vector_numeros_gaussianos){
    int i;

    double dos_numeros_gausianos[2];

    for (i = 0 ; i < Ng; i += 2){
        num_aleatorio_gaussiano(dos_numeros_gausianos);
        vector_numeros_gaussianos[i] = dos_numeros_gausianos[0];
        vector_numeros_gaussianos[i+1] = dos_numeros_gausianos[1];
    }
}


void calcular_media_desviacion(double *vector, int N, double *media, double *desviacion) {

    double suma = 0.0;
    double suma_cuadrada = 0.0;
    int i;

    for (i = 0; i < N; i++) {
        suma += vector[i];
        suma_cuadrada += vector[i] * vector[i];
    }

    *media = suma / N;
    *desviacion = sqrt((suma_cuadrada / N) - (*media) * (*media));

    printf("media = %f \n", *media);
    printf("desviacion = %f \n", *desviacion);

}



//Lo mismo, esta función es para comprobar que la funcion de numeros gaussianos rulaba bien, probablemente se quitará en un futuro
//Aunque, también nos sirve para lo de comprobar que la distribución de posiciones y velocidades sea gaussiana
void histograma(double *V, int num_elementos, const char *nombre_archivo) {
    // Variables
    int asignacion, i;
    double anchura, min_desconocido, max_desconocido;
    double histograma_matriz[Intervalo];  // Intervalo viene del #define

    // Inicializar mínimo y máximo
    max_desconocido = -1000;
    min_desconocido = 1000;

    // Encontrar el valor mínimo y máximo en V
    for (i = 0; i < num_elementos; i++) {
        if (V[i] > max_desconocido) {
            max_desconocido = V[i];
        }
        if (V[i] < min_desconocido) {
            min_desconocido = V[i];
        }
    }

    // Evitar división por cero si max == min
    if (max_desconocido == min_desconocido) {
        printf("Todos los valores son iguales. No se puede crear un histograma significativo.\n");
        return;
    }

    // Calcular la anchura de cada intervalo del histograma
    anchura = (max_desconocido - min_desconocido) / (double)Intervalo;

    // Inicializar el histograma
    for (i = 0; i < Intervalo; i++) {
        histograma_matriz[i] = 0.0;
    }

    // Llenar el histograma
    for (i = 0; i < num_elementos; i++) {
        asignacion = (int)floor((V[i] - min_desconocido) / anchura);

        // Asegurarse de que la asignación no exceda los límites
        if (asignacion >= Intervalo) {
            asignacion = Intervalo - 1;  // Ajuste para el valor máximo
        }
        histograma_matriz[asignacion] += 1;
    }

    // Normalizar el histograma (convertir a probabilidades)
    ///Ya no es una probabilidad, cambiar si quiere que sea el caso
    for (i = 0; i < Intervalo; i++) {
        histograma_matriz[i] = histograma_matriz[i] / (num_elementos * anchura);
    }

    // Escribir el histograma a un archivo para gnuplot
    FILE *f;
    f = fopen(nombre_archivo, "w");
    if (f == NULL) {
        printf("Error al abrir el archivo %s\n", nombre_archivo);
        return;
    }


    //int eje_x = (int)(Intervalo / 2); //La media queda en eje_x = 47
    for (i = 0; i < Intervalo; i++) {
        fprintf(f, "%f %f\n", (i+0.5) * anchura + min_desconocido, histograma_matriz[i]);
    }
    fclose(f);

    // Comprobación de que las cosas vayan bien
    printf("Max: %f, Min: %f\n", max_desconocido, min_desconocido);
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
    int i, pasos, j, w, pasos_representar, pasos_por_representado, pasos_ahora;
    double factor_estocastico, dos_terminos_estocasticos[2], nabla_dividido_m, tiempo;

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


    pasos = (int)(T/h);
    //Para evitar que los archivos que sacamos pesen gigas, hacemos que no se guarden todos; se guardan T*10 independientemente de h
    //Lo fácil sería añadir un if y rulando, pero para evitar cargar aún más el tiempo de computo, lo que hacemos es separar en 2 for
    pasos_representar = T*10;
    pasos_por_representado = (int)(1/(10*h) + 1);

    factor_estocastico = sqrt(2*nabla*K_b_T*h);
    nabla_dividido_m = nabla/m;





    //Vector para rellenar los histogramas para comprobar la distribucion de momentos y posiciones
    #ifdef comprobar_dist
        double *posicion_histo, *momento_histo;
        posicion_histo = malloc(pasos * sizeof(double));
        momento_histo = malloc(pasos * sizeof(double));

        //Comprobación de que la memoria se haya asignado correctamente
        if (posicion_histo == NULL || momento_histo == NULL) {
            fprintf(stderr, "Error: No se pudo asignar memoria.\n");
            exit(1);
        }

    #endif // comprobar_dist_gaussiana


    //Guardamos la posicion inicial
    fprintf(f,"%d ", 0); //Esto es el tiempo
    fprintf(f,"%f ", posicion);
    fprintf(f,"%f ", momento);

    #ifdef energia
        fprintf(f,"%f ", 0);
        fprintf(f,"%f ", 0);
        fprintf(f,"%f\n", 0);
    #else
        fprintf(f,"\n");
    #endif



    #ifdef energia
        double mediapotencial,mediacinetica;
        mediacinetica=0;
        mediapotencial=0;
    #endif // energia

    pasos_ahora = 0;

    for (i = 0 ; i < pasos_representar; i++){
        for (w = 0 ; w < pasos_por_representado; w += 2){

            termino_estocastico_Z(factor_estocastico, dos_terminos_estocasticos);

            for(j=0 ; j<2;j++){
                //En cada paso salen 2 por lo de tener dos numeros gaussianos
                posicion += f_pn(momento)*h;
                momento += g_xn_pn(posicion, momento, nabla_dividido_m)*h + dos_terminos_estocasticos[j];


                pasos_ahora++;

                //Rellenamos los vectores histogramas para comprobar la distribucion de momentos y posicion
                #ifdef comprobar_dist
                    posicion_histo[i+j] = posicion;
                    momento_histo[i+j] = momento;
                #endif // comprobar_dist_gaussiana


                #ifdef energia
                    mediacinetica+=energia_cinetica(momento);
                    mediapotencial+=energia_potencial(posicion);
                    //fprintf(f,"%f ", energia_cinetica(momento));
                    //fprintf(f,"%f\n", energia_potencial(posicion));
                #endif // energia

            }
        }

        tiempo = (i+1)/10;

        fprintf(f,"%f ", tiempo);
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);

        #ifdef energia
            fprintf(f,"%f ", mediacinetica/pasos_ahora);
            fprintf(f,"%f ", mediapotencial/pasos_ahora);
            fprintf(f,"%f\n", (mediapotencial+mediacinetica)/pasos_ahora);
        #else
            fprintf(f, "\n");
        #endif // energia

    }
    fclose(f);



    //Pasamos a un fichero las energías finales
    #ifdef energia
        f = fopen("tabla_energias", "a");
        char nombre_archivo_cin[1024], nombre_archivo_pot[1024], nombre_archivo_tot[1024];


        strcpy(nombre_archivo_cin, nombre_archivo);
        strcpy(nombre_archivo_pot, nombre_archivo);
        strcpy(nombre_archivo_tot, nombre_archivo);
        strncat(nombre_archivo_cin, " cinetica = %f\n", 1024);
        strncat(nombre_archivo_pot, " potencial = %f\n", 1024);
        strncat(nombre_archivo_tot, " total = %f\n", 1024);

        fprintf(f, nombre_archivo_cin, mediacinetica/pasos_ahora);
        fprintf(f, nombre_archivo_pot, mediapotencial/pasos_ahora);
        fprintf(f, nombre_archivo_tot, (mediacinetica+mediapotencial)/pasos_ahora);

        fclose(f);
    #endif // energia




    //Llamamos a la funcion histograma para comprobar la distribucion de momentos y posicion
    #ifdef comprobar_dist
        int media, desviacion;
        calcular_media_desviacion(posicion_histo, pasos, &media, &desviacion);
        printf(nombre_archivo" posicion: media = %f, desviacion = %f", media, desviacion);

        calcular_media_desviacion(posicion_histo, pasos, &media, &desviacion);
        printf(nombre_archivo" momento: media = %f, desviacion = %f", media, desviacion);


        strcpy(nombre_archivo, "Histo_posicion_");
        strncat(nombre_archivo, "Euler_Maruyama_h=", 1024);
        strncat(nombre_archivo, especificador_h, 1024);
        strncat(nombre_archivo, nombre_archivo_parte2, 1024);
        strncat(nombre_archivo, especificador_nabla, 1024);
        strncat(nombre_archivo, nombre_archivo_fin, 1024);

        histograma(posicion_histo, pasos, nombre_archivo);


        strcpy(nombre_archivo, "Histo_momento_");
        strncat(nombre_archivo, "Euler_Maruyama_h=", 1024);
        strncat(nombre_archivo, especificador_h, 1024);
        strncat(nombre_archivo, nombre_archivo_parte2, 1024);
        strncat(nombre_archivo, especificador_nabla, 1024);
        strncat(nombre_archivo, nombre_archivo_fin, 1024);

        histograma(momento_histo, pasos, nombre_archivo);



        //Importante liberar memoria
        free(momento_histo);
        momento_histo = NULL;

        free(posicion_histo);
        posicion_histo = NULL;


    #endif // comprobar_dist

}




void Runge_Kutta (double posicion, double momento, double h, double nabla){
    int i, pasos, j, w, pasos_representar, pasos_por_representado, pasos_ahora;
    //Uso la notación del power point de clase, las barras bajas se deben entender como "sub", p.ej f sub x1
    double f_x1, f_x2, g_p1, g_p2, Z, h_medios, factor_estocastico, dos_terminos_estocasticos[2], nabla_dividido_m, tiempo;

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
    //Para evitar que los archivos que sacamos pesen gigas, hacemos que no se guarden todos; se guardan T*10 independientemente de h
    //Lo fácil sería añadir un if y rulando, pero para evitar cargar aún más el tiempo de computo, lo que hacemos es separar en 2 for
    pasos_representar = T*10;
    pasos_por_representado = (int)(pasos/pasos_representar + 1);

    //printf("pasos_por_representado = %d\n", pasos_por_representado);

    h_medios = h*0.5;
    factor_estocastico = sqrt(2*nabla*K_b_T*h);
    nabla_dividido_m = nabla/m;


    //Vector para rellenar los histogramas para comprobar la distribucion de momentos y posiciones
    #ifdef comprobar_dist
        double *posicion_histo, *momento_histo;
        posicion_histo = malloc(pasos * sizeof(double));
        momento_histo = malloc(pasos * sizeof(double));

        //Comprobación de que la memoria se haya asignado correctamente
        if (posicion_histo == NULL || momento_histo == NULL) {
            fprintf(stderr, "Error: No se pudo asignar memoria.\n");
            exit(1);
        }

    #endif // comprobar_dist_gaussiana


    #ifdef energia
        double mediapotencial,mediacinetica;
        mediacinetica=0;
        mediapotencial=0;
    #endif // energia

    //Guardamos posicion y momento inicial
    fprintf(f,"%d ", 0); //Esto es el tiempo
    fprintf(f,"%f ", posicion);
    fprintf(f,"%f ", momento);

    #ifdef energia
        fprintf(f,"%f ", 0);
        fprintf(f,"%f ", 0);
        fprintf(f,"%f\n", 0);
    #else
        fprintf(f,"\n");
    #endif


    pasos_ahora = 0;

    for (i = 0 ; i < pasos_representar; i++){
        for (w = 0 ; w < pasos_por_representado; w += 2){

            termino_estocastico_Z(factor_estocastico, dos_terminos_estocasticos);

            for(j = 0; j < 2; j++){
                //De nuevo, de cada paso de tiempo hay que hacer 2
                Z = dos_terminos_estocasticos[j];

                //Paso 1
                f_x1=f_pn(momento+Z);
                g_p1=g_xn_pn(posicion, momento+Z, nabla_dividido_m);

                //Paso 2
                f_x2=f_pn(momento+h*g_p1);
                g_p2=g_xn_pn(posicion+h*f_x1, momento+h*g_p1, nabla_dividido_m);

                //Final, calculo de nuevo punto en el espacio de fases
                posicion += h_medios*(f_x1+f_x2);
                momento += h_medios*(g_p1+g_p2)+Z;

                pasos_ahora++;

                //Rellenamos los vectores histogramas para comprobar la distribucion de momentos y posicion
                #ifdef comprobar_dist
                    posicion_histo[i+j] = posicion;
                    momento_histo[i+j] = momento;
                #endif // comprobar_dist_gaussiana


                #ifdef energia
                    mediacinetica+=energia_cinetica(momento);
                    mediapotencial+=energia_potencial(posicion);
                #endif // energia
            }
        }

        tiempo = (i+1)/10;


        fprintf(f,"%f ", tiempo);
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);

        #ifdef energia
            fprintf(f,"%f ", mediacinetica/pasos_ahora);
            fprintf(f,"%f ", mediapotencial/pasos_ahora);
            fprintf(f,"%f\n", (mediapotencial+mediacinetica)/pasos_ahora);
        #else
            fprintf(f, "\n");
        #endif // energia

    }

    fclose(f);



    //Pasamos a un fichero las energías finales
    #ifdef energia
        f = fopen("tabla_energias", "a");
        char nombre_archivo_cin[1024], nombre_archivo_pot[1024], nombre_archivo_tot[1024];


        strcpy(nombre_archivo_cin, nombre_archivo);
        strcpy(nombre_archivo_pot, nombre_archivo);
        strcpy(nombre_archivo_tot, nombre_archivo);
        strncat(nombre_archivo_cin, " cinetica = %f\n", 1024);
        strncat(nombre_archivo_pot, " potencial = %f\n", 1024);
        strncat(nombre_archivo_tot, " total = %f\n", 1024);

        fprintf(f, nombre_archivo_cin, mediacinetica/pasos_ahora);
        fprintf(f, nombre_archivo_pot, mediapotencial/pasos_ahora);
        fprintf(f, nombre_archivo_tot, (mediacinetica+mediapotencial)/pasos_ahora);

        fclose(f);
    #endif // energia



    //Llamamos a la funcion histograma para comprobar la distribucion de momentos y posicion
    #ifdef comprobar_dist
        int media, desviacion;
        calcular_media_desviacion(posicion_histo, pasos, &media, &desviacion);
        printf(nombre_archivo" posicion: media = %f, desviacion = %f", media, desviacion);

        calcular_media_desviacion(posicion_histo, pasos, &media, &desviacion);
        printf(nombre_archivo" momento: media = %f, desviacion = %f", media, desviacion);



        strcpy(nombre_archivo, "Histo_posicion_");
        strncat(nombre_archivo, "Runge_Kutta_h=", 1024);
        strncat(nombre_archivo, especificador_h, 1024);
        strncat(nombre_archivo, nombre_archivo_parte2, 1024);
        strncat(nombre_archivo, especificador_nabla, 1024);
        strncat(nombre_archivo, nombre_archivo_fin, 1024);

        histograma(posicion_histo, pasos, nombre_archivo);



        strcpy(nombre_archivo, "Histo_momento_");
        strncat(nombre_archivo, "Runge_Kutta_h=", 1024);
        strncat(nombre_archivo, especificador_h, 1024);
        strncat(nombre_archivo, nombre_archivo_parte2, 1024);
        strncat(nombre_archivo, especificador_nabla, 1024);
        strncat(nombre_archivo, nombre_archivo_fin, 1024);

        histograma(momento_histo, pasos, nombre_archivo);



        //Importante liberar memoria
        free(momento_histo);
        momento_histo = NULL;

        free(posicion_histo);
        posicion_histo = NULL;


    #endif // comprobar_dist
}




void verlet (double posicion, double momento, double h, double nabla){
    int i, j, pasos, w, pasos_representar, pasos_por_representado, pasos_ahora;
    //Uso la notación del power point de clase, las barras bajas se deben entender como "sub", p.ej f sub x1
    double factor_estocastico, dos_terminos_estocasticos[2], a, b, factor_posicion, fuerza_ahora, h_medios, Z, tiempo;

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



    pasos = (int)(T/h + 0.9999);
    //Para evitar que los archivos que sacamos pesen gigas, hacemos que no se guarden todos; se guardan T*10 independientemente de h
    //Lo fácil sería añadir un if y rulando, pero para evitar cargar aún más el tiempo de computo, lo que hacemos es separar en 2 for
    pasos_representar = T*10;
    pasos_por_representado = (int)(pasos/pasos_representar + 0.999);



    factor_estocastico = sqrt(2*nabla*K_b_T*h);
    h_medios = h/2;
    b = 1/(1+nabla*h_medios/m);
    a = 2*b-1;
    factor_posicion = b*h_medios/m;

    #ifdef energia
        double mediapotencial,mediacinetica;
        mediacinetica=0;
        mediapotencial=0;
    #endif // energia


    //Vector para rellenar los histogramas para comprobar la distribucion de momentos y posiciones
    #ifdef comprobar_dist
        double *posicion_histo, *momento_histo;
        posicion_histo = malloc(pasos * sizeof(double));
        momento_histo = malloc(pasos * sizeof(double));

        //Comprobación de que la memoria se haya asignado correctamente
        if (posicion_histo == NULL || momento_histo == NULL) {
            fprintf(stderr, "Error: No se pudo asignar memoria.\n");
            exit(1);
        }

    #endif // comprobar_dist_gaussiana



    //Guardamos posicion y momento inicial
    fprintf(f,"%d ", 0); //Esto es el tiempo
    fprintf(f,"%f ", posicion);
    fprintf(f,"%f ", momento);

    #ifdef energia
        fprintf(f,"%f ", 0);
        fprintf(f,"%f ", 0);
        fprintf(f,"%f\n", 0);
    #else
        fprintf(f,"\n");
    #endif


    pasos_ahora = 0;


    for (i = 0 ; i < pasos_representar; i++){
        for (w = 0 ; w < pasos_por_representado; w += 2){

            termino_estocastico_Z(factor_estocastico, dos_terminos_estocasticos);

            for(j=0;j<2;j++){
            //Paso 1 de tiempo
                fuerza_ahora = fuerza(posicion);
                Z = dos_terminos_estocasticos[j];

                posicion += factor_posicion*(2*momento+h*fuerza_ahora+Z);
                momento = a*momento+h_medios*(a*fuerza_ahora+fuerza(posicion))+b*Z;


                pasos_ahora++;

                //Rellenamos los vectores histogramas para comprobar la distribucion de momentos y posicion
                #ifdef comprobar_dist
                    posicion_histo[i+j] = posicion;
                    momento_histo[i+j] = momento;
                #endif // comprobar_dist_gaussiana


                #ifdef energia
                    mediacinetica+=energia_cinetica(momento);
                    mediapotencial+=energia_potencial(posicion);

                    /*fprintf(f,"%f ", energia_cinetica(momento));
                    fprintf(f,"%f\n", energia_potencial(posicion));*/

                #else
                    fprintf(f,"\n");

                #endif // energia
            }
        }

        tiempo = (i+1)/10;

        //printf("pasos_ahora = %d\n", pasos_ahora);

        fprintf(f,"%f ", tiempo);
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);

        #ifdef energia
            fprintf(f,"%f ", mediacinetica/pasos_ahora);
            fprintf(f,"%f ", mediapotencial/pasos_ahora);
            fprintf(f,"%f\n", (mediapotencial+mediacinetica)/pasos_ahora);
        #else
            fprintf(f, "\n");
        #endif // energia


    }
    fclose(f);



    //Pasamos a un fichero las energías finales
    #ifdef energia
        f = fopen("tabla_energias", "a");
        char nombre_archivo_cin[1024], nombre_archivo_pot[1024], nombre_archivo_tot[1024];


        strcpy(nombre_archivo_cin, nombre_archivo);
        strcpy(nombre_archivo_pot, nombre_archivo);
        strcpy(nombre_archivo_tot, nombre_archivo);
        strncat(nombre_archivo_cin, " cinetica = %f\n", 1024);
        strncat(nombre_archivo_pot, " potencial = %f\n", 1024);
        strncat(nombre_archivo_tot, " total = %f\n", 1024);

        fprintf(f, nombre_archivo_cin, mediacinetica/pasos_ahora);
        fprintf(f, nombre_archivo_pot, mediapotencial/pasos_ahora);
        fprintf(f, nombre_archivo_tot, (mediacinetica+mediapotencial)/pasos_ahora);

        fclose(f);
    #endif // energia



    #ifdef comprobar_dist
        int media, desviacion;
        calcular_media_desviacion(posicion_histo, pasos, &media, &desviacion);
        printf(nombre_archivo" posicion: media = %f, desviacion = %f", media, desviacion);

        calcular_media_desviacion(posicion_histo, pasos, &media, &desviacion);
        printf(nombre_archivo" momento: media = %f, desviacion = %f", media, desviacion);


        strcpy(nombre_archivo, "Histo_posicion_");
        strncat(nombre_archivo, "Verlet_h=", 1024);
        strncat(nombre_archivo, especificador_h, 1024);
        strncat(nombre_archivo, nombre_archivo_parte2, 1024);
        strncat(nombre_archivo, especificador_nabla, 1024);
        strncat(nombre_archivo, nombre_archivo_fin, 1024);

        histograma(posicion_histo, pasos, nombre_archivo);


        strcpy(nombre_archivo, "Histo_momento_");
        strncat(nombre_archivo, "Verlet_h=", 1024);
        strncat(nombre_archivo, especificador_h, 1024);
        strncat(nombre_archivo, nombre_archivo_parte2, 1024);
        strncat(nombre_archivo, especificador_nabla, 1024);
        strncat(nombre_archivo, nombre_archivo_fin, 1024);

        histograma(momento_histo, pasos, nombre_archivo);



        //Importante liberar memoria
        free(momento_histo);
        momento_histo = NULL;

        free(posicion_histo);
        posicion_histo = NULL;


    #endif // comprobar_dist
}




#ifdef energia
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

#endif //energia
