///     -------Practica 5 TFIII, Sim molecular------      ///
   /// -----Tercera parte: POLIMERO (con fuerza)----- ///


//actualizado el 23-11-24


///DUDAS:
//valor de nabla?



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>





#define m 1     //Ahora es la masa de un monómero
#define K_b_T 1 //cte de Boltzmann por temperatura
#define PI 3.14159265358979323846


//para los numerros random, NO TOCAR
#define NormRANu (2.3283063671E-10F)

#define T 500000 //Tiempo en cualquier algoritmo
#define h 0.001  //Paso de tiempo
#define nabla 1

#define Intervalo 100 //Numero de intervalos (cajitas) en el histograma

//No se si nos hará falta, pero igual para comprobar es bueno empezar en 2D, de momento lo dejo en 3D:
#define dim 3

//Los archivos no pueden contener todos los pasos, si no mi ordenador y gnuplot explotan, por eso escribimos aquí un quiero uno de cada cuantos pasos:
#define pasos_representar T //Así pues, si ponemos T = 10000, h = 0.001, hay 10000000 pasos, pero solo estamos guardando uno de cada 1000, es decir, 10000


///Datos nuevos para esta parte:
#define k_e 100 //Cte elástica, una cte multiplicativa al potencial
#define b 1 //Dist de equilibrio entre los muelles, en el doble pozo teníamos (x^2-1)^2, pues ahora será (x^2-b)^2



//Numero de monomeros en la cadena, luego se puede poner con un #define, a partir del punto vii
//int N = 2;
#define N 16






//Variables globales para los numeros aleatorios
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;




//Estos ultimos son para los ifdef:
//De momento ninguno xd




///Predeclaración de funciones
//Numeros aleatorios uniformes
void ini_ran(int SEMILLA);
double Random(void);
void num_aleatorio_gaussiano (double *dos_numeros_gaussianos);


//Histogramas, en principio su funcion es la misma, se que son distintas pero no se exactamente en que
void histograma (double *V/*matriz de entrada*/,int num_elementos, const char *nombre_archivo);

//Calcular distancias entre monomeros consecutivos
void distancias(double **posiciones, double **diferencia_de_posiciones, double *modulo_distancias);

//Ecuaciones oscilador
void termino_estocastico_Z (double factor_estocastico, double *dos_terminos_estocasticos);
void sacar_array_numeros_Z(double **array_Z_aux, double factor_estocastico);

void funcion_fuerzas(double **fuerzas, double **diferencia_de_posiciones, double *modulo_distancias, double fuerza_extremo);


//Algoritmo
void paso_verlet_posiciones(double **posiciones, double **momentos, double **fuerza_old, double **array_Z_aux, double factor_posicion);
void paso_verlet_momentos(double **momentos, double **fuerza, double **fuerza_old, double **array_Z_aux, double a, double b_verlet, double h_medios);
void verlet (double **posiciones, double **momentos, double fuerza_extremo, double *vector_extremo_extremo_medio);


//Datos para medir
double energia_cinetica_cadena (double **momentos);
double energia_potencial_cadena (double *modulo_distancias);
double energia_cinetica_monomero_1D (double momento);
double energia_potencial_monomero_1D (double distancia);

double distancia_extremo_extremo(double **posiciones, double *vector_estremo_estremo);

//Todas menos esta se podrían eliminar
double radio_de_giro(double **posiciones, double *CDM);


//Cosas auxiliares
void config_inicial_polimero(double **posiciones, double **momentos);
void calcular_CDM_y_momentum(double **posiciones, double **momentos, double *CDM, double *momentum_CDM);
void posiciones_respecto_al_CDM(double **posiciones, double **momentos, double *CDM, double *momentum_CDM);









///     --------MAIN--------       ///
int main(){
    int i;

    // Inicializamos la rueda de n meros random de Parisi-Rapuano
    ini_ran(time(NULL));


    ///Arrays de las coordenadas generalizadas
    /*
    En esta parte, vamos a tener N partículas con 3 posiciones, momentos cada una; 6N variables.
    Se podría hacer todo en un solo array, pero para no hacer un lio con los indices vamos a hacer 2 arrays de 2 indices, p. ej: posiciones[N][3]
    El primer indice nos indica a que monomero de la cadena nos referimos y el segundo la posicion en los ejes x,y,z que identificamos con 0,1,2.
    */
    double **posiciones = malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
        posiciones[i] = malloc(dim * sizeof(double));

    double **momentos = malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
        momentos[i] = malloc(dim * sizeof(double));

    //Y les damos una configuracion inicial
    config_inicial_polimero(posiciones, momentos);




    ///Ejecutamos Verlet y sacamos absolutamente todo
    //Vamos a hacer un barrido en fuerzas (adimensionalizada):
    double fuerza_extremo, vector_extremo_extremo_medio[3], modulo_vector_extremo_extremo_medio = 0, suma_vector_extremo_extremo_medio = 0;
    char nombre_archivo_fuerza_vs_extension[1024];
    sprintf(nombre_archivo_fuerza_vs_extension, "fuerza_vs_extension_T=%d_N=%d", T, N);
    FILE *archivo_fuerza_vs_extension;
    archivo_fuerza_vs_extension = fopen(nombre_archivo_fuerza_vs_extension, "w")
    fprintf(archivo_fuerza_vs_extension, "#Fuerza\t Extension media en:\t x\t y\t z\t x+y+z\t modulo\t");

    for (fuerza_extremo = 0.2; fuerza_extremo < 1.1; fuerza_extremo += 0.2){
        for (j = 0; j < 3; j++)
            vector_extremo_extremo_medio[j] = 0;
        verlet (posiciones, momentos, fuerza_extremo, vector_extremo_extremo_medio);
        for (j = 0; j < 3; j++){
            modulo_vector_extremo_extremo_medio += vector_extremo_extremo_medio[j]*vector_extremo_extremo_medio[j];
            suma_vector_extremo_extremo_medio += vector_extremo_extremo_medio[j];
        }

        fprintf(archivo_fuerza_vs_extension, "%f\t %f\t %f\t %f\t %f\t %f\t ", vector_extremo_extremo_medio[0], vector_extremo_extremo_medio[1], vector_extremo_extremo_medio[2], suma_vector_extremo_extremo_medio, modulo_vector_extremo_extremo_medio);
    }

    fclose(archivo_fuerza_vs_extension);





    //Liberamos memoria
    for (i = 0; i < N; i++)
        free(posiciones[i]);
    free(posiciones);

    for (i = 0; i < N; i++)
        free(momentos[i]);
    free(momentos);
}
///     ------FIN DEL MAIN------     ///










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
    aleatorio_uniforme_1 = Random ();
    while (aleatorio_uniforme_1 == 0.0)
        aleatorio_uniforme_1 = Random ();
    aleatorio_uniforme_2 = Random ();

    auxiliar_1=sqrt(-2*log(aleatorio_uniforme_1));
    auxiliar_2=2*PI*aleatorio_uniforme_2;

    //en la presentaci n dan como dos posibilidades, seg n las pruebas que he hecho es indistinto usar una u otra
    dos_numeros_gaussianos[0]= auxiliar_1*cos(auxiliar_2);

    //Esta ser a la segunda forma, solo cambia el cos por el sen
    dos_numeros_gaussianos[1] = auxiliar_1*sin(auxiliar_2);
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





///Calcular distancias
//Calculamos la distancia entre monómeros consecutivos, tanto en el mismo eje como en total (distacia euclidea), se usa tanto en la función fuerzas como en la energía potencial
void distancias(double **posiciones, double **diferencia_de_posiciones, double *modulo_distancias){

    int i, j;

    for (i = 0; i < (N - 1); i++){
        for (j = 0; j < dim; j++)
            diferencia_de_posiciones[i][j] = posiciones[i][j] - posiciones[i+1][j];
    }

    for (i = 0; i < (N - 1); i++){
        modulo_distancias[i] = 0;
        for (j = 0; j < dim; j++)
            modulo_distancias[i] += diferencia_de_posiciones[i][j]*diferencia_de_posiciones[i][j];
        modulo_distancias[i] = sqrt(modulo_distancias[i]);
    }
}





///     ECUACIONES OSCILADOR      ///
//este es el termino estocastico, que va incluido en la ecuacion del momento
//Saca 2 para aprovechar mejor la funcion generadora de numeros gaussianos
void termino_estocastico_Z (double factor_estocastico, double *dos_terminos_estocasticos){
    double dos_numeros_gaussianos[2];
    num_aleatorio_gaussiano(dos_numeros_gaussianos);
    //Para que no tenga que pararse cada vez a calcular, el factor = sqrt(2*nabla*K_b_T*h)
    dos_terminos_estocasticos[0] = factor_estocastico*dos_numeros_gaussianos[0];
    dos_terminos_estocasticos[1] = factor_estocastico*dos_numeros_gaussianos[1];
}


//Esta funcion rellena el array auiliar de términos estocásticos
void sacar_array_numeros_Z(double **array_Z_aux, double factor_estocastico){
    double dos_terminos_estocasticos[2];

    for (int i = 0; i < N; i+= 2){
        for (int j = 0; j < dim; j++){
            termino_estocastico_Z (factor_estocastico, dos_terminos_estocasticos);
            array_Z_aux[i][j] = dos_terminos_estocasticos[0];
            array_Z_aux[i+1][j] = dos_terminos_estocasticos[1];
        }
    }
}


//Aqui rellenamos el array auxiliar de fuerzas, sea el viejo o el nuevo
void funcion_fuerzas(double **fuerzas, double **diferencia_de_posiciones, double *modulo_distancias, double fuerza_extremo){
    int i, j;
    double constante = -k_e*1.0;


    //En la cadena, a cada polimero (excepto los dos extrmos) le actuan dos fuerzas. Pero, bien sabemos que varias se repiten pues la que hace el
    //polimero 1 sobre el 2 es igual a menos la del 2 sobre el 1. Asi que para variar, otro array auxiliar, considerando solo las fuerza en un sentido:
    double **fuerzas_un_sentido = malloc(N * sizeof(double *));
    for (i = 0; i < (N - 1); i++){
        fuerzas_un_sentido[i] = malloc(dim * sizeof(double));
        for (j = 0; j < dim; j++)
            fuerzas_un_sentido[i][j] = -constante * diferencia_de_posiciones[i][j] * (1 - (1.0*b)/modulo_distancias[i]);
    }

    //Ahora si por fin, sacamos la fuerza final sobre un monomero:
    for (j = 0; j < dim; j++){
        // Fuerza sobre el primer monómero (solo tiene conexión con el segundo), PERO ya que ahora está fijo pues vale 0, aunque daría igual porque su posición no se actualiza
        fuerzas[0][j] = 0;     //-fuerzas_un_sentido[0][j];
        // Fuerza sobre el último monómero (solo tiene conexión con el penúltimo); y ahora también le añadimos una fuerza en la diracción x, que viene fuera del bucle
        fuerzas[N-1][j] = fuerzas_un_sentido[N-2][j];

        // Fuerza sobre los monómeros intermedios
        for (i = 1; i < (N - 1); i++)
            fuerzas[i][j] = fuerzas_un_sentido[i-1][j] - fuerzas_un_sentido[i][j];
    }
    ///Añadimos una fuerza cte en la dirección x:
    fuerzas[N-1][0] += fuerza_extremo;

    //Liberamos memoria
    for (i = 0; i < (N - 1); i++)
        free(fuerzas_un_sentido[i]);
    free(fuerzas_un_sentido);
}





///     -------------------EJECUCION ALGORITMO VERLET---------------------     ///
///Quizas sería bueno hacer un cambio de sistema de referencia, bien sea en torno al CDM o a un monomero concreto pàra evitar overflow en las posiciones, pero sería un jaleo tambien
void paso_verlet_posiciones(double **posiciones, double **momentos, double **fuerza_old, double **array_Z_aux, double factor_posicion){
    for (int i = 1; i < N; i++){
        for (int j = 0; j < dim; j++)
            posiciones[i][j] += factor_posicion*(2*momentos[i][j]+h*fuerza_old[i][j]+array_Z_aux[i][j]);
    }
}


void paso_verlet_momentos(double **momentos, double **fuerza, double **fuerza_old, double **array_Z_aux, double a, double b_verlet, double h_medios){
    for (int i = 1; i < N; i++){
        for (int j = 0; j < dim; j++)
            momentos[i][j] = a * momentos[i][j] + h_medios * (a * fuerza_old[i][j] + fuerza[i][j]) + b_verlet * array_Z_aux[i][j];
    }
}




void verlet (double **posiciones, double **momentos, double fuerza_extremo, double *vector_extremo_extremo_medio){

    int t, i, j, pasos;
    double a, b_verlet, factor_posicion, h_medios, factor_estocastico;

    //Calculamos las variables necesarias para Verlet en cualquier caso, que se mantendrán constantes durante el bucle
    pasos = (int)(T/h);

    h_medios = h/2;
    b_verlet = 1/(1+nabla*h_medios/m);
    a = 2*b_verlet-1;
    factor_posicion = b_verlet*h_medios/m;
    factor_estocastico = sqrt(2*nabla*K_b_T*h);




    ///Arrays auxiliares
    //Es conveniente tambien crear dos arrayses de fuerzas para simplificar la implementacion, uno con las viejas y otro las nuevas
    double **fuerzas = malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
        fuerzas[i] = malloc(dim * sizeof(double));

    double **fuerzas_old = malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
        fuerzas_old[i] = malloc(dim * sizeof(double));


    //Tambien es conveniente sacar las distancias entre monómeros consecutivos, se usa tanto en la función fuerzas como en la energía potencial
    double **diferencia_de_posiciones = malloc(N * sizeof(double *));
    for (i = 0; i < (N - 1); i++)
        diferencia_de_posiciones[i] = malloc(dim * sizeof(double));

    double modulo_distancias[N-1];


    //Por ultimo (creo), un array de numeros aleatorios, ya que los que se aplican en cada particula y eje deben ser distintos
    double **array_Z_aux = malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
        array_Z_aux[i] = malloc(dim * sizeof(double));







    ///Datos a medir
     //DISTANCIA EXTREMO A EXTREMO (promedio)
    double distancia_extremo_extremo_media = 0, vector_extremo_extremo[3], CDM[3], momentum_CDM[3];

    char nombre_archivo_dist_ext_ext[1024];
    sprintf(nombre_archivo_dist_ext_ext, "Distancia_extremo_extremo_N=%d_T=%d.dat", N, T);
    FILE *archivo_dist_ext_ext;
    archivo_dist_ext_ext = fopen(nombre_archivo_dist_ext_ext, "w");

    fprintf(archivo_dist_ext_ext, "#Tiempo T\t Dist extremo-extremo: \t x\t y\t z\t modulo\n");






    ///Empieza el bucle con los pasos del algoritmo

    //Inicializamos fuerzas y distancias
    distancias(posiciones, diferencia_de_posiciones, modulo_distancias);
    funcion_fuerzas(fuerzas, diferencia_de_posiciones, modulo_distancias, fuerza_extremo);


    for (t = 0; t < pasos; t++) {
        // Guardar las fuerzas viejas
        for (i = 0; i < N; i++) {
            for (j = 0; j < dim; j++) {
                fuerzas_old[i][j] = fuerzas[i][j];
            }
        }


        // Generar términos estocásticos
        sacar_array_numeros_Z(array_Z_aux, factor_estocastico);

        // Actualizar posiciones
        paso_verlet_posiciones(posiciones, momentos, fuerzas_old, array_Z_aux, factor_posicion);

        //Actualizamos las distancias entre monomeros consecutivos
        distancias(posiciones, diferencia_de_posiciones, modulo_distancias);

        // Calcular nuevas fuerzas
        funcion_fuerzas(fuerzas, diferencia_de_posiciones, modulo_distancias, fuerza_extremo);

        // Actualizar momentos
        paso_verlet_momentos(momentos, fuerzas, fuerzas_old, array_Z_aux, a, b_verlet, h_medios);

        //Por último, actualizamos el CDM
        calcular_CDM_y_momentum(posiciones, momentos, CDM, momentum_CDM);





        ///Datos a medir
        //DISTANCIA EXTREMO A EXTREMO, de nuevo un promedio, pero en este caso calculamos tanto la media vectorial como de distancia
        distancia_extremo_extremo_media += distancia_extremo_extremo(posiciones, vector_extremo_extremo);
        for (j = 0; j < dim; j++)
            vector_extremo_extremo_medio[j] += vector_extremo_extremo[j];


        // Escribir en el archivo cada `pasos_representar` pasos
        if (t % pasos_representar == 0) {
            fprintf(archivo_Rg, "%f\t%lf\n", t * h, Rg_media / (t + 1));
            fprintf(archivo_dist_ext_ext, "%f\t%lf\t%lf\t%lf\t%lf\n", t * h, vector_extremo_extremo_medio[0] / (t + 1), vector_extremo_extremo_medio[1] / (t + 1), vector_extremo_extremo_medio[2] / (t + 1), distancia_extremo_extremo_media / (t+1));

            //Y ya que estamos, renormalizamos las posiciones respecto del CDM
            //void posiciones_respecto_al_CDM(posiciones, momentos, CDM, momentum_CDM);
        }

        // Comprobar errores en la escritura
        if (ferror(archivo_dist_ext_ext)) {
            printf("Error escribiendo en el archivo\n");
            fclose(archivo_dist_ext_ext);
            return;
        }



    }
    ///Fin del bucle principal


    //Cerramos archivos
    fclose(archivo_dist_ext_ext);
    printf("Archivo de dist ext-ext generado correctamente: %s\n", nombre_archivo_dist_ext_ext);


    //Liberamos memoria
    for (i = 0; i < N; i++){
        free(fuerzas[i]);
        free(fuerzas_old[i]);
        free(array_Z_aux[i]);
    }
    free(fuerzas);
    free(fuerzas_old);
    free(array_Z_aux);

    for (i = 0; i < (N - 1); i++)
        free(diferencia_de_posiciones[i]);
    free(diferencia_de_posiciones);

}
///     -------------------HASTA AQUI VERLET---------------------       ///






///Datos a medir
double energia_cinetica_cadena (double **momentos){
    int i, j;
    double Ecin_cadena = 0;

    for (i = 0; i < N; i++){
        for (j = 0; j < dim; j++){
            Ecin_cadena += energia_cinetica_monomero_1D(momentos[i][j]);
            //if ((j == 0) && (i == 0))
                //printf("%f\n",momentos[i][j]);
                //printf("%f\n",Ecin_cadena);
        }
    }

    return Ecin_cadena;
}

double energia_potencial_cadena(double *modulo_distancias) {
    int i;
    double Epot_cadena = 0;

    //A diferencia de la cietica, la potencial va asociada a un "muelle" (o enlace) y no a la partícula; por eso el bucle es hasta N-1
    for (i = 0; i < N - 1; i++)
        Epot_cadena += energia_potencial_monomero_1D(modulo_distancias[i]);

    return Epot_cadena;
}

double energia_cinetica_monomero_1D (double momento){
    double Ec;
    Ec = momento*momento*0.5/m;
    return Ec;
}

double energia_potencial_monomero_1D (double distancia){
    double Ep, parentesis;
    parentesis = distancia - b;
    Ep = 0.5*k_e*parentesis*parentesis;
    return Ep;
}


double radio_de_giro(double **posiciones, double *CDM){
    double dist_al_CDM, Rg = 0, parentesis;

    for (int i = 0; i < N; i++){
        dist_al_CDM = 0;
        for (int j = 0; j < dim; j++){
            parentesis = posiciones[i][j] - CDM[j];
            dist_al_CDM += parentesis*parentesis;
        }
        Rg += dist_al_CDM;
    }
    Rg = sqrt(Rg/N);
    return Rg;
}


//Si se quiere, se pueden eliminar todas las funciones de esta categoría menos esta, que es la única que vamos a utilizar para esta parte
double distancia_extremo_extremo(double **posiciones, double *vector_estremo_estremo){
    double dist_ext_ext = 0;
    for (int j = 0; j < dim; j++){
        vector_estremo_estremo[j] = posiciones[0][j] - posiciones[N-1][j];
        dist_ext_ext += vector_estremo_estremo[j] * vector_estremo_estremo[j];
    }
    dist_ext_ext = sqrt(dist_ext_ext);
    return dist_ext_ext;
}







///Otras cosas que no van en ninguna categoria
//Inicializamos las posiciones y momentos
void config_inicial_polimero(double **posiciones, double **momentos){

    int i, j, direccion;


    //Ponemos simplemente que los momentos iniciales sean 0
    for (i = 0; i < N; i++){
        for (j = 0; j < dim; j++)
            momentos[i][j] = 0;
    }


    //En cuanto a las posiciones, ponemos el primero en 0, y las sucesivas las ponemos en uno de los ejes aleatoriamente a distancia b del anterior:
    for (j = 0; j < dim; j++)
        posiciones[0][j] = 0;

    for (i = 1; i < N; i++){
        for (j = 0; j < dim; j++){
            posiciones[i][j] = posiciones[i-1][j];
        }
        direccion = (int)(Random() * 6);
        if (direccion < 3)
            posiciones[i][direccion] += b;
        else{
            direccion = direccion - 3;
            posiciones[i][direccion] -= b;
        }
    }
}



//Lo necesitamos para calcular el radio de giro, y en el futuro para representar el polímero
void calcular_CDM_y_momentum(double **posiciones, double **momentos, double *CDM, double *momentum_CDM) {
    for (int j = 0; j < dim; j++) {
        CDM[j] = 0.0;
        momentum_CDM[j] = 0.0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) {
            CDM[j] += posiciones[i][j] / N;
            momentum_CDM[j] += momentos[i][j] / N;
        }
    }
}


//Renormaliza las distancias y momentos respecto del Centro De Masas
void posiciones_respecto_al_CDM(double **posiciones, double **momentos, double *CDM, double *momentum_CDM) {
    // Ajustamos posiciones respecto al CDM
    for (int i = 0; i < N; i++) { // Ahora incluye todas las partículas
        for (int j = 0; j < dim; j++) {
            posiciones[i][j] -= CDM[j];
            momentos[i][j] -= momentum_CDM[j]; // Ajustamos momentos también
        }
    }
}

