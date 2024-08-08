# Configuración de PM3D y tamaño de la imagen
set pm3d map
set size ratio 1

# Comando para graficar los datos
splot "Ondas.dat"

# Pausa para mostrar el gráfico antes de continuar
pause -1 "Presione Enter para continuar..."

# Configuración del terminal y salida
set terminal jpeg enhanced
set output "Ondas.jpg"

# Comando para graficar los datos nuevamente y guardar la imagen
splot "Ondas.dat"

# Limpiar configuraciones para futuros usos
reset
