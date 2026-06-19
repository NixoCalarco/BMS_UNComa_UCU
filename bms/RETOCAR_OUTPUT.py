import csv

# Inicializamos contenedores
Vh, Se0, XXe0, TT0, TE0 = [], [], [], [], []
Vhat, Se1, XXe1, TT1, TE1 = [], [], [], [], []

with open("data.csv", "r") as f:
    reader = csv.reader(f)

    for row in reader:
        # Ejemplo de row:
        # ['Kalman', ' V', ' 3.200991 Se', ' 0.100000 XXe', ' 0.228819 TT', ' 0.000000 TE', ' 0.000000']

        # Normalizar quitando espacios
        tokens = [t.strip() for t in row]

        # Determinar si es KF o EKF
        mode = tokens[0]
        isKF  = mode == "Kalman"
        isEKF = mode == "Kalman Extended"

        # Extraer valores según posiciones
        V   = float(tokens[2].split()[0])
        Se  = float(tokens[3].split()[0])
        XXe = float(tokens[4].split()[0])
        TT  = float(tokens[5].split()[0])
        TE  = float(tokens[6].split()[0])

        if isKF:
            Vh.append(V)
            Se0.append(Se)
            XXe0.append(XXe)
            TT0.append(TT)
            TE0.append(TE)

        elif isEKF:
            Vhat.append(V)
            Se1.append(Se)
            XXe1.append(XXe)
            TT1.append(TT)
            TE1.append(TE)

