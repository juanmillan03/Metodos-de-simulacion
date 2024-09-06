#!/bin/bash

# Lista de números de threads que quieres probar
threads_list=(2 4 6 8 10 12 16 18)

# Nombre del ejecutable
executable="./LBM"

# Bucle para ejecutar el programa con diferentes números de threads
for threads in "${threads_list[@]}"; do
    export OMP_NUM_THREADS=$threads
    echo "Ejecutando con $threads threads..."
    # Medir el tiempo de ejecución
    time $executable
done
