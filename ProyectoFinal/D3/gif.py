import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import glob
import os
import re

# Lista de directorios donde están los archivos
directories = ["D3/8/", "D3/16/", "D3/12/"]
# Directorio donde se guardarán los GIFs
gif_directory = "D3/"

# Función para extraer el número de tiempo del nombre del archivo
def extract_time_from_filename(filename):
    match = re.search(r"(\d+)\.txt$", filename)
    if match:
        return int(match.group(1))
    return float('inf')  # Para asegurar que los archivos sin número se ordenen al final

# Loop sobre cada directorio
for directory in directories:
    print(f"Procesando directorio: {directory}")

    # Lista de archivos .txt en el directorio actual ordenados por número de tiempo
    file_list = sorted(glob.glob(directory + "Ondas_*.txt"), key=extract_time_from_filename)

    # Verificar que se encontraron archivos
    if not file_list:
        print(f"No se encontraron archivos .txt en el directorio {directory}.")
        continue

    # Crear una lista para almacenar las imágenes
    images = []

    # Definir límites de los ejes (ajusta según tus datos)
    x_min, x_max = 0, 10  # Rango para el eje x
    y_min, y_max = 0, 10  # Rango para el eje y
    presion_min, presion_max = -1,1  # Rango para la presión (eje de color

    # Loop sobre cada archivo y generar un gráfico
    for filename in file_list:
        try:
            # Leer los datos
            data = np.loadtxt(filename)
            
            # Verificar que los datos tienen la forma correcta
            if data.shape[1] < 3:
                print(f"Advertencia: El archivo {filename} no tiene suficientes columnas de datos.")
                continue

            # Asignar las columnas a variables
            x = data[:, 0]
            y = data[:, 1]
            presion = data[:, 2]
            
            # Crear la figura
            plt.figure(figsize=(10, 6))
            scatter = plt.scatter(x, y, c=presion, cmap='magma')
            #scatter = plt.scatter(x, y, c=presion, cmap='magma')
            plt.colorbar(scatter, label='Presión')
            plt.xlabel('x(m)')
            plt.ylabel('y(m)')
            plt.title(f'Evolución del sistema - {os.path.basename(filename)}')
            plt.grid(True)
            
            # Guardar la imagen temporalmente
            temp_filename = f"{directory}temp_{os.path.splitext(os.path.basename(filename))[0]}.png"
            plt.savefig(temp_filename)
            plt.close()
            
            # Verificar si el archivo de imagen se ha guardado correctamente
            if not os.path.isfile(temp_filename):
                print(f"Error: No se pudo guardar la imagen temporal {temp_filename}.")
                continue
            
            # Leer la imagen y agregarla a la lista
            images.append(imageio.imread(temp_filename))
        
        except Exception as e:
            print(f"Error procesando el archivo {filename}: {e}")

    # Verificar que se generaron imágenes
    if not images:
        print(f"No se generaron imágenes para el directorio {directory}.")
        continue

    # Crear el GIF
    try:
        # Definir el nombre del GIF basado en el nombre del directorio
        gif_filename = f"{gif_directory}evolucion_sistema_{os.path.basename(directory.rstrip('/'))}.gif"
        imageio.mimsave(gif_filename, images, duration=1, loop=0)
        print(f"GIF creado exitosamente como '{gif_filename}'")
    except Exception as e:
        print(f"Error al crear el GIF para el directorio {directory}: {e}")

    # Eliminar archivos temporales
    for temp_file in glob.glob(f"{directory}temp_*.png"):
        os.remove(temp_file)
