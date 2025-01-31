import numpy as np
import matplotlib.pyplot as plt

# Dati forniti
#x_values = np.array([ 50, 100, 200, 400, 800])
#x_values = np.array([ 10, 20, 40, 80])
x_values = np.array([ 20, 40, 80])
y_values = np.array([0.00405746, 0.00103989, 0.000261453])
#y_values = np.array([0.0050279, 0.00040189, 1.6584953e-05, 1.510174e-06, 1.364273137e-07])
#y_values = np.array([0.0050279, 0.00040189, 1.6584953e-05, 1.510174e-06, 1.364273137e-07])
#y_values = np.array([0.0099723, 0.00127229, 0.00014061, 1.65480080e-05, 2.0395803e-06 ])

# Funzione 1/x
#x_line = np.linspace(49, 800, 500)  # Range per la linea 1/x
x_line = np.linspace(min(x_values), max(x_values), 500)
y_line = 5e1 / (x_line**(3))
z_line = 1e5 / (x_line**(4))
w_line = 1e5 / (x_line**(5))
k_line = 1 / (x_line**(2))

# Creazione del grafico log-log
plt.figure(figsize=(8, 6))
plt.loglog(x_values, y_values, 'o-', label='Erreur', markersize=8)
#plt.loglog(x_line, y_line, '--', label='$1/h^{3}$', color='red')
#plt.loglog(x_line, z_line, '--', label='$1/h^4$', color='green')
#plt.loglog(x_line, w_line, '--', label='$1/h^5$', color='red')
plt.loglog(x_line, k_line, '--', label='$1/h^2$', color='red')


# Personalizzazione del grafico
plt.title('Erreur par rapport au passé spatial SSPRK3-WENO5', fontsize=14)
plt.xlabel('$dx$', fontsize=12)
plt.ylabel('$E$', fontsize=12)
plt.legend(fontsize=12)
plt.grid(which="both", linestyle="--", linewidth=0.5)

# Mostra il grafico
plt.show()
