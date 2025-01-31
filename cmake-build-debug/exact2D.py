import numpy as np
import matplotlib.pyplot as plt

def extract_diagonal_from_file(file_path):
    try:
        # Leggi la matrice dal file
        with open(file_path, 'r') as file:
            matrix = []
            for line in file:
                row = list(map(float, line.split()))  # Converte la riga in numeri float
                matrix.append(row)

        # Verifica che la matrice sia quadrata
        num_rows = len(matrix)
        for row in matrix:
            if len(row) != num_rows:
                raise ValueError("La matrice nel file non è quadrata.")

        # Estrai la diagonale principale
        diagonal = [matrix[i][i] for i in range(num_rows)]

        return diagonal

    except FileNotFoundError:
        print("Errore: File non trovato.")
    except ValueError as e:
        print(f"Errore: {e}")
    except Exception as e:
        print(f"Si è verificato un errore: {e}")


def plot_diagonal_with_function(diagonal, known_function):
    n = len(diagonal)
    x = np.linspace(0,1,n)  # Genera gli indici per la diagonale
    y_known = known_function(x)  # Calcola la funzione nota

    # Plot della diagonale e della funzione nota
    plt.figure(figsize=(10, 6))
    plt.plot(x, diagonal, 'o-', label='Diagonale estratta', markersize=8)
    plt.plot(x, y_known, 'r--', label='Funzione conosciuta', linewidth=2)
    plt.xlabel('Indice')
    plt.ylabel('Valore')
    plt.title('Confronto tra diagonale e funzione conosciuta')
    plt.legend()
    plt.grid()
    plt.show()


# Esempio di utilizzo
file_path = "rho.txt"  # Sostituisci con il percorso del tuo file
diagonal = extract_diagonal_from_file(file_path)

#rho_ref = np.loadtxt('rho.txt')

def known_function(x):
        return 1+ 0.2*np.sin(3.14159*(2*(x) - .1*(.5)))

    # Plot della diagonale e confronto con la funzione
plot_diagonal_with_function(diagonal, known_function)

# Calcola la norma L1 della differenza (escludendo estremi)
def calculate_l1_norm(data1, data2):
    diff = np.abs(data1 - data2)
    return np.sum(diff)

n = len(diagonal)
x = np.linspace(0,1,n)
# Valori della funzione conosciuta
y_known = known_function(x)

# Calcola la norma L1 della differenza escludendo i primi e gli ultimi 3 elementi
l1_norm = calculate_l1_norm(diagonal, y_known)

# Stampa il risultato della norma L1
print(f"Norma L1 della differenza (escludendo 3 elementi iniziali e finali): {l1_norm:.5f}")