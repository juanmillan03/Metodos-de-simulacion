# Configuración de PM3D y tamaño de la imagen
set pm3d map
set size ratio 1

# Comando para graficar los datos
splot "Espejo.dat"

# Pausa para mostrar el gráfico antes de continuar
pause -1 "Presione Enter para continuar..."

# Configuración del terminal y salida
set terminal jpeg enhanced
set output "Espejo.jpg"

# Comando para graficar los datos nuevamente y guardar la imagen
splot "Espejo.dat"

# Limpiar configuraciones para futuros usos
reset
