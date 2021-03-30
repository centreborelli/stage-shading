import iio                   # bibliothèque d'entrée/sortie de fichiers image
x = iio.read("fuji.npy")     # lecture du tableau d'hauteurs sur l'image x
y = x[1:,] - x[:-1,]         # dérivée partielle par différences finies
Y = 127 + 2 * y              # cadrage du rang dans [0,255] (à la louche)
iio.write("fuji_dy.png", Y)  # écriture de la dérivée partielle comme image
