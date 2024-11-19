///     -------Practica 5 TF3, sim molecular------      ///
       /// -----Segunda parte: DOBLE POZO----- ///


//actualizado el 13-11-24




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define m 1     //masa
#define k 1     //cte el stica
#define K_b_T 0.2 //cte de Boltzmann por temperatura
#define PI 3.14159265358979323846


//para los numerros random, NO TOCAR
#define NormRANu (2.3283063671E-10F)

#define T 500000 //Tiempo en cualquier algoritmo
#define h 0.001  //Paso de tiempo

#define Intervalo 100 //Numero de intervalos en el histograma



///Estos ultimos son para los ifdef:


//He quitado la opción de bucle, de todas formas se puede poner un bucle que vaya de 0 a 1 y hace lo mismo que la opcion de unico, ensuciaba el codigo

//Calcular las energ as durante el proceso (energia)
    #define energia

//Pasar las posiciones y momentos a fichero (esp_de_fases)
    //#define esp_de_fases

//Activar (ocupacion) o desactivar todo lo de la ocupación y tal
    //#define ocupacion

//Comprobar la distribucion de posiciones y momentos usando un array (comprobar_dist)
    //#define comprobar_dist


//Para elegir el tipo de potencial que se utiliza (doble_pozo) o (fuerza_constante)
    #define doble_pozo //oscilador_armonico




//Numeros aleatorios uniformes
void ini_ran(int SEMILLA);
double Random(void);
void num_aleatorio_gaussiano (double *dos_numeros_gaussianos);

//Histogramas, en principio su funcion es la misma, se que son distintas pero no se exactamente en que
void Histograma(int a, int b, int N, int divisiones,  double *array);
void histograma (double *V/*matriz de entrada*/,int num_elementos, const char *nombre_archivo);

//Ecuaciones oscilador
double fuerza (double posicion);
void termino_estocastico_Z (double factor_estocastico, double *dos_terminos_estocasticos);

//Algoritmo
void verlet (double posicion, double momento, double nabla, int *cuenta_saltos);

//Datos para medir
double energia_cinetica (double momento);
double energia_potencial(double posicion);








//Como la altura A va a cambiar es mejor tenerlo como variable global
double A = 0;


//Variables globales para los numeros aleatorios
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;




int main(){

    //Variables necesarias siempre
    double /*coef de viscosidad*/nabla;
    double /*Pto inicial en el espacio de fases*/ momento_inicial=0, posicion_inicial=0.00001;


    //Variables bucle
    int i;


    //variables ocupación
    #ifdef ocupacion
        double tiempo_est_medio, T_aux = T;
    #endif // ocupacion
    //Como verlet necesita que le pasemos la variable pues hace falta que esté aunque no se use, se podría hacer pero que palo
    int cuenta_saltos = 0;


    //Archivo tabla energias finales, debe resetearse en cada ejecución porque el resto son apends
    #ifdef energia
        char nombre_tabla_energias[]="Tabla_energias_finales_T=", especificador_T[20], nombre_archivo_fin[]=".txt";

        sprintf(especificador_T, "%d", T);

        strncat(nombre_tabla_energias, especificador_T, 1024);
        strncat(nombre_tabla_energias, nombre_archivo_fin, 1024);

        FILE *e;
        e = fopen(nombre_tabla_energias, "w");
        fprintf(e, "#Dampin\t Altura maximo\t Energia cinetica\t Energia potencial\t Energia total\n");
        fclose(e);
    #endif // energia


    // Inicializamos la rueda de n meros random de Parisi-Rapuano
    ini_ran(123456789);





    ///EJECUTAMOS ALGORITMOS

///LA OPCIÓN DE ÚNICO AHORA MISMO NO FUNCIONA, no por nada, simplemente no la he arreglado
///    ¡¡¡¡LA OPCIÓN DE BUCLE SE PUEDE USAR COMO UNICO!!!!  solo hay que poner un bucle de una iteracion

    //Este bucle solo sirve para que salgan varios archivos con distintos valores de nabla, altura... De una, si se prefiere

    //Abrimos el archivo de ocupación frente a damping
    #ifdef ocupacion
        FILE *f;
        f = fopen("t_medio_frente_a_altura.txt", "w");
        fprintf(f, "#A    #Tiempo de estancia medio\n");
    #endif // ocupacion



    //Esto es un apaño un poco cutre, para sacar más puntos por etas entre 2 y 3 para la ocupacion, si se hace otra cosa leer mas abajo
    int j;

    nabla = 0;
    for (i = 0; i < 3; i++){
        //Esto es para representar distintos nablas o A de una, las cosas de la ocupacion principalmente
        /*
        if (i == 0)
            A = 0.5;
        if (i == 1)
            A = 1;
        if (i == 2)
            A = 2;
        */
        nabla += 0.5;
        A = 0;
        for (j = 0; j < 10; j++){   //Acordarse de descomentar el final del bucle cuando se use
            //nabla += factor_bucle_nabla;

            A += 0.5;


            ///Ejecutamos Verlet
            verlet (posicion_inicial, momento_inicial, nabla, &cuenta_saltos);



            //Sacamos aquí la ocupación
            #ifdef ocupacion
                tiempo_est_medio = T_aux/(cuenta_saltos+1.0);
                fprintf(f, "%f      %f\n", A, tiempo_est_medio);
                //Comprobación y bueno, para ver como va mientras se ejecuta porque si no me pienso que no funciona xd
            #endif // ocupacion

            printf("nabla=%f; A=%f\n", nabla, A);
        }
    }

    //Cerramos el archivo de ocupación frente a damping
    #ifdef ocupacion
        fclose(f);
    #endif // ocupacion





    return 0;
}







///Números aleatorios
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




///Numero aleatorio gaussiano
//Genera un n mero aleatorio con distribuci n gaussiana, a partir de un numero random en el intervalo [0,1)
void num_aleatorio_gaussiano (double *dos_numeros_gaussianos){
    double aleatorio_uniforme_1, aleatorio_uniforme_2, auxiliar_1, auxiliar_2;

    ///Aqui hay un problema con el algoritmo: Cuando sale aleatorio_uniforme_1 = 0 ---> ln(0)=-inf; y da error
    //No se como se supone que deber amos arreglarlo, de momento solo voy a poner una clausula de que no sea igual a cero
    aleatorio_uniforme_1=Random ();
    while (aleatorio_uniforme_1 == 0.0)
        aleatorio_uniforme_1=Random ();
    aleatorio_uniforme_2=Random ();

    auxiliar_1=sqrt(-2*log(aleatorio_uniforme_1));
    auxiliar_2=2*PI*aleatorio_uniforme_2;

    //en la presentaci n dan como dos posibilidades, seg n las pruebas que he hecho es indistinto usar una u otra
    dos_numeros_gaussianos[0]= auxiliar_1*cos(auxiliar_2);

    //Esta ser a la segunda forma, solo cambia el cos por el sen
    dos_numeros_gaussianos[1] = auxiliar_1*sin(auxiliar_2);
}




///Histogramas
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
    f=fopen ("Histograma.txt","w");
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
    Entrada: Limites intervalo (a y b, tipo int), numero de divisiones (divisiones tipo int), array con los valores de estudio (array tipo double) y su longitud (N, tipo int)
    Salida:Fichero con las frecuencias de cada intervalo.
    */
}




//Nos sirve para lo de comprobar que la distribuci n de posiciones y velocidades sea gaussiana; y quizá otras cosas
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

///IMPORTANTE, esto es (la menos derivada de) el potencial que tendremos que cambiar en la siguiente parte as  que ojito
double fuerza (double posicion){
    double grad_Vx;


    #ifdef fuerza_constante
        grad_Vx=-A*4*posicion*(posicion*posicion-1) + constante_fuerza;
    #endif // fuerza_constante

    #ifdef doble_pozo
        grad_Vx=-A*4*posicion*(posicion*posicion-1);
    #endif // doble_pozo


    return grad_Vx;
}




//este es el t rmino estoc stico, que va incluido en la ecuaci n del momento
//Saca 2 para aprovechar mejor la funcion generadora de numeros gaussianos
void termino_estocastico_Z (double factor_estocastico, double *dos_terminos_estocasticos){
    double dos_numeros_gaussianos[2];
    num_aleatorio_gaussiano(dos_numeros_gaussianos);
    //Para que no tenga que pararse cada vez a calcular, el factor = sqrt(2*nabla*K_b_T*h)
    dos_terminos_estocasticos[0] = factor_estocastico*dos_numeros_gaussianos[0];
    dos_terminos_estocasticos[1] = factor_estocastico*dos_numeros_gaussianos[1];
}






///     -------------------EJECUCION ALGORITMO VERLET---------------------     ///
//Además del propio algoritmo, contiene un monton de cosas mas
void verlet (double posicion, double momento, double nabla, int *cuenta_saltos){


    //Variables necesarias para ejecutar Verlet en cualquier caso
    int i, j, pasos;
    double factor_estocastico, dos_terminos_estocasticos[2], a, b, factor_posicion, fuerza_ahora, h_medios, Z;

    //Calculamos las variables necesarias para Verlet en cualquier caso, que se mantendrán constantes durante el bucle
    pasos = (int)(T/h);
    factor_estocastico = sqrt(2*nabla*K_b_T*h);
    h_medios = h/2;
    b = 1/(1+nabla*h_medios/m);
    a = 2*b-1;
    factor_posicion = b*h_medios/m;



    //Variables ocupación de cada estado
    #ifdef ocupacion
        unsigned int t_estancia, lado_actual;

        //Inicializamos las variables necesarias para la ocupación
        t_estancia = 0;
        *cuenta_saltos = 0;
        ///Lado actual depende de como se ponga el x inicial, 1 si x inicial positivo y 0 si es negativo
        if (posicion > 0)
            lado_actual = 1;
        if (posicion < 0)
            lado_actual = 0;
        if (posicion == 0){
            printf("La posicion inicial no puede ser 0");
            return;
        }
    #endif // ocupacion



    //Variables calcular energías promedio
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

    ///NOMBRES ARCHIVOS, ojo que el orden en el que está puestoimporta y mucho da fallo de memoria que se yo porque
    //Toda la parafernalia de los char es para poder sacar todos los archivos de distintas A y nabla en un solo bucle

    //TITULOS, dependiendo de que ifdef pongamos y lo que queramos sacar basicamente
    #ifdef esp_de_fases
        char nombre_archivo[]="Espacio_fases_dp_A=";
    #else
        #ifdef energia
            char nombre_archivo[]="Energias_dp_A=", nombre_tabla_energias[]="Tabla_energias_finales_T=";
        #endif
    #endif // esp_de_fases


    #ifdef ocupacion
        char nombre_ocupacion[]="Ocupacion_A=";
    #endif // ocupacion


    #ifdef comprobar_dist
        char archivo_posiciones[]="Distribucion_posiciones_A=";
        char archivo_momentos[]="Distribucion_momentos_A=";
    #endif // comprobar_dist




    //ESPECIFICADORES, se usan siempre y vienen a poner en el nombre del archivo que valores de T, nabla y A tenemos
    char especificador_A[1024], especificador_nabla[1024], nombre_archivo_parte2[]="_nabla=";
    char especificador_T[20], nombre_archivo_parte3[]="_T=", nombre_archivo_fin[]=".txt";
    sprintf(especificador_A, "%f", A);
    sprintf(especificador_nabla, "%f", nabla);
    sprintf(especificador_T, "%d", T);




    //UNIMOS ESPECIFICADORES, pues añadimos los especificadores necesarios a cada archivo
    //Archivo ocupacion
    #ifdef ocupacion
        strncat(nombre_ocupacion, especificador_A, 1024);
        strncat(nombre_ocupacion, nombre_archivo_parte2, 1024);
        strncat(nombre_ocupacion, especificador_nabla, 1024);
        strncat(nombre_ocupacion, nombre_archivo_parte3, 1024);
        strncat(nombre_ocupacion, especificador_T, 1024);
        strncat(nombre_ocupacion, nombre_archivo_fin, 1024);

        //Archivo ocupación de cada estado
        FILE *g;
        g = fopen(nombre_ocupacion, "w");
    #endif // ocupacion



    //Archivo energias finales
    #ifdef energia
        strncat(nombre_tabla_energias, especificador_T, 1024);
        strncat(nombre_tabla_energias, nombre_archivo_fin, 1024);

        printf(nombre_tabla_energias);
        FILE *e;
        e = fopen(nombre_tabla_energias, "a");
    #endif // energia



    //Archivo trayectorias y/o energías promedio
    #ifdef esp_de_fases
        strncat(nombre_archivo, especificador_A, 1024);
        strncat(nombre_archivo, nombre_archivo_parte2, 1024);
        strncat(nombre_archivo, especificador_nabla, 1024);
        strncat(nombre_archivo, nombre_archivo_parte3, 1024);
        strncat(nombre_archivo, especificador_T, 1024);
        strncat(nombre_archivo, nombre_archivo_fin, 1024);

        FILE *f;
        f = fopen(nombre_archivo, "w");

        int pasos_no_representados, cada_cuantos_pasos_representamos;
        cada_cuantos_pasos_representamos = 1; //Representamos uno de (cada_cuantos_pasos_representamos) pasos
        pasos_no_representados = 0; //Contador de los pasos a representar
    #else
        #ifdef energia
            strncat(nombre_archivo, especificador_A, 1024);
            strncat(nombre_archivo, nombre_archivo_parte2, 1024);
            strncat(nombre_archivo, especificador_nabla, 1024);
            strncat(nombre_archivo, nombre_archivo_parte3, 1024);
            strncat(nombre_archivo, especificador_T, 1024);
            strncat(nombre_archivo, nombre_archivo_fin, 1024);

            FILE *f;
            f = fopen(nombre_archivo, "w");


            int pasos_no_representados, cada_cuantos_pasos_representamos;
            //Representamos uno de (cada_cuantos_pasos_representamos) pasos
            cada_cuantos_pasos_representamos = (int)(pasos/1000 + 0.999);
            pasos_no_representados = 0; //Contador de los pasos a representar
        #endif
    #endif // esp_de_fases


    //Aqui para comprobar la distribucion de las posiciones y velocidades
    #ifdef comprobar_dist
        double (*array_posiciones) = malloc(pasos*sizeof(double));
        double (*array_momentos) = malloc(pasos*sizeof(double));
        for (i = 0 ; i < pasos; i++){
            array_posiciones[i] = 0;
            array_momentos[i] = 0;
        }
    #endif // comprobar_dist












    ///Empieza el bucle con los pasos del algoritmo
;    for (i = 0 ; i < pasos; i += 2){

        termino_estocastico_Z(factor_estocastico, dos_terminos_estocasticos);


        for(j=0;j<2;j++){
            fuerza_ahora = fuerza(posicion);
            Z = dos_terminos_estocasticos[j];

            posicion += factor_posicion*(2*momento+h*fuerza_ahora+Z);
            momento = a*momento+h_medios*(a*fuerza_ahora+fuerza(posicion))+b*Z ;

            ///El algoritmo en si es hasta aquí, ahora cálculos


            //Aqui va todo lo necesario para la ocupacion
            #ifdef ocupacion
                //tiempo de estancia; como posicion_inicial > 0, la primera columna es el tiempo a la derecha (en x positivo)
                //Para que salga bien hay que eliminar las posiciones cercanas a 0
                if (lado_actual == 1){
                    if(posicion < -0.2){
                        lado_actual = 0;
                        (*cuenta_saltos)++;

                        fprintf(g, "%f ", t_estancia*h);
                        t_estancia = 0;
                    }
                    else
                        t_estancia++;
                }

                else{
                    if(posicion > 0.2){
                        lado_actual = 1;

                        fprintf(g, "%f ", t_estancia*h);
                        t_estancia = 0;

                        (*cuenta_saltos)++;
                        //Si se pone que empiece en x < 0, habría que cambiar esto
                        fprintf(g, "\n");
                    }
                    else
                        t_estancia++;
                }
            #endif // ocupacion




            #ifdef esp_de_fases
                if (pasos_no_representados == cada_cuantos_pasos_representamos){
                    fprintf(f,"%f ", (i+j+1)*h); //Esto es el tiempo
                    fprintf(f,"%f ", posicion);
                    fprintf(f,"%f ", momento);
                }
            #endif // esp_de_fases




            //Este putisimo lio es necesario para representar la energía y el espacio de fases a la vez:
            #ifdef energia
                mediacinetica += energia_cinetica(momento);
                mediapotencial += energia_potencial(posicion);

                if (pasos_no_representados == cada_cuantos_pasos_representamos){
                    //Ojala existiera esto mas facil, porque esto es un "si no esta definido esp_fases:"
                    #ifdef esp_de_fases
                    #else
                        fprintf(f,"%f ", (i+j+1)*h); //Esto es el tiempo
                    #endif // esp_de_fases
                    fprintf(f,"%f ", mediacinetica/(i+j+1));
                    fprintf(f,"%f ", mediapotencial/(i+j+1));
                    fprintf(f,"%f\n", (mediapotencial+mediacinetica)/(i+j+1));
                    pasos_no_representados = 0;
                }
                pasos_no_representados ++;

            #else
                #ifdef esp_de_fases
                    if (pasos_no_representados == cada_cuantos_pasos_representamos){
                        fprintf(f,"\n");
                        pasos_no_representados = 0;
                    }
                    pasos_no_representados ++;
                #endif
            #endif // esp_de_fases



            //Rellenamos los arrays para comprobar la distribucion de posiciones y momentos
            #ifdef comprobar_dist
                array_posiciones[i+j] = posicion;
                array_momentos[i+j] = momento;
            #endif // comprobar_dist
        }

    }


    #ifdef esp_de_fases
        fprintf(f,"%f ", (i+j)*h); //Esto es el tiempo
        fprintf(f,"%f ", posicion);
        fprintf(f,"%f ", momento);
    #endif // esp_de_fases

    #ifdef energia
        //energias promedio al final del archivo
        fprintf(f,"%f ", (i+j+1)*h); //Esto es el tiempo
        fprintf(f,"%f ", mediacinetica/(i+j+1));
        fprintf(f,"%f ", mediapotencial/(i+j+1));
        fprintf(f,"%f\n", (mediapotencial+mediacinetica)/(i+j+1));

        //Pasamos las energias promedio finales en a tabla
        fprintf(e, "%f\t %f\t ", nabla, A);
        fprintf(e,"%f\t ", mediacinetica/(i+j+1));
        fprintf(e,"%f\t ", mediapotencial/(i+j+1));
        fprintf(e,"%f\n", (mediapotencial+mediacinetica)/(i+j+1));
    #endif // true




    //Cerramos el archivo de ocupación
    #ifdef ocupacion
        fclose(g);
    #endif // ocupacion


    //Cerramos el archivo de energias finales
    #ifdef energia
        fclose(e);
    #endif // energia


    //Cerramos el archivo de trayectorias y/o energias promedio
    #ifdef esp_de_fases
        fclose(f);
    #else
        #ifdef energia
            fclose(f);
        #endif
    #endif // esp_de_fases


    //Pasamos a histograma y limpiamos memoria de la comprobacion de distribucion de posiciones y velocidades
    #ifdef comprobar_dist
        strncat(archivo_posiciones, especificador_A, 1024);
        strncat(archivo_posiciones, nombre_archivo_parte2, 1024);
        strncat(archivo_posiciones, especificador_nabla, 1024);
        strncat(archivo_posiciones, nombre_archivo_parte3, 1024);
        strncat(archivo_posiciones, especificador_T, 1024);
        strncat(archivo_posiciones, nombre_archivo_fin, 1024);

        histograma(array_posiciones, pasos, archivo_posiciones);


        strncat(archivo_momentos, especificador_A, 1024);
        strncat(archivo_momentos, nombre_archivo_parte2, 1024);
        strncat(archivo_momentos, especificador_nabla, 1024);
        strncat(archivo_momentos, nombre_archivo_parte3, 1024);
        strncat(archivo_momentos, especificador_T, 1024);
        strncat(archivo_momentos, nombre_archivo_fin, 1024);

        histograma(array_momentos, pasos, archivo_momentos);


        //Liberamos memoria
        free(array_posiciones);
        free(array_momentos);
    #endif // comprobar_dist
}
///     -------------------HASTA AQUI VERLET---------------------       ///








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
