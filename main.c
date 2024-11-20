#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define Intervalo 100  // Número de intervalos (cajas) en el histograma


void histograma(const char *nombre_archivo_entrada, int columna, int num_elementos, const char *nombre_archivo_salida);


int main() {
    char nombre_archivo_entrada[1024];
    char nombre_archivo_salida[1024];
    int num_elementos = -1;  // Número de elementos, -1 para contar automáticamente
    int columna = 1;

    histograma(nombre_archivo_entrada, columna, num_elementos, nombre_archivo_salida);
}

void histograma(const char *nombre_archivo_entrada, int columna, int num_elementos, const char *nombre_archivo_salida) {
    int asignacion, i;
    double anchura, min_desconocido, max_desconocido;
    double variable_histo;
    double histograma_matriz[Intervalo];  // Inicializar el histograma
    for (i = 0; i < Intervalo; i++) histograma_matriz[i] = 0;

    FILE *g = fopen(nombre_archivo_entrada, "r");
    if (g == NULL) {
        printf("Error al abrir el archivo de entrada %s\n", nombre_archivo_entrada);
        return;
    }

    // Determinar el número de elementos si no se proporciona
    if (num_elementos < 0) {
        num_elementos = 0;
        char linea[1024];
        while (fgets(linea, sizeof(linea), g)) {
            num_elementos++;
        }
        rewind(g);  // Regresar al inicio del archivo
    }

    // Inicializar valores máximos y mínimos
    max_desconocido = -INFINITY;
    min_desconocido = INFINITY;

    // Determinar los valores mínimo y máximo en la columna especificada
    for (i = 0; i < num_elementos; i++) {
        double valores[100];  // Ajustar según el número esperado de columnas
        int col_actual = 0;
        while (fscanf(g, "%lf", &valores[col_actual]) == 1) {
            col_actual++;
        }

        if (columna > col_actual) {
            printf("Columna especificada (%d) fuera de rango.\n", columna);
            fclose(g);
            return;
        }

        variable_histo = valores[columna - 1];
        if (variable_histo > max_desconocido)
            max_desconocido = variable_histo;
        if (variable_histo < min_desconocido)
            min_desconocido = variable_histo;
    }

    fclose(g);

    // Evitar división por cero si todos los valores son iguales
    if (max_desconocido == min_desconocido) {
        printf("Todos los valores son iguales. No se puede crear un histograma significativo.\n");
        return;
    }

    anchura = (max_desconocido - min_desconocido) / Intervalo;

    // Rellenar el histograma
    g = fopen(nombre_archivo_entrada, "r");
    for (i = 0; i < num_elementos; i++) {
        double valores[100];
        int col_actual = 0;
        while (fscanf(g, "%lf", &valores[col_actual]) == 1) {
            col_actual++;
        }

        variable_histo = valores[columna - 1];
        asignacion = (int)floor((variable_histo - min_desconocido) / anchura);
        if (asignacion >= Intervalo)
            asignacion = Intervalo - 1;

        histograma_matriz[asignacion]++;
    }
    fclose(g);

    // Normalizar el histograma
    for (i = 0; i < Intervalo; i++) {
        histograma_matriz[i] /= (num_elementos * anchura);
    }

    // Guardar el histograma en el archivo de salida
    FILE *f = fopen(nombre_archivo_salida, "w");
    if (f == NULL) {
        printf("Error al abrir el archivo de salida %s\n", nombre_archivo_salida);
        return;
    }

    for (i = 0; i < Intervalo; i++) {
        fprintf(f, "%f %f\n", (i + 0.5) * anchura + min_desconocido, histograma_matriz[i]);
    }
    fclose(f);

    printf("Histograma generado: %s\n", nombre_archivo_salida);
}
