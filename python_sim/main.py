import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator, interp1d

# data
lambert = np.loadtxt("./archivos/Lamb.txt", dtype=float)
emf = np.loadtxt("./archivos/EMF_flipdim.txt", dtype=float)
_, _, corriente, _, _, _, tension = np.loadtxt(
    "./archivos/Ciclado_pack-2_vc_2.txt", dtype=float, unpack=True
)

corriente *= -1
tension /= 1_000

muestreo_original = 1
muestreo_deseado = 1
corriente = corriente[:: int(muestreo_original / muestreo_deseado)]
tension = tension[:: int(muestreo_original / muestreo_deseado)]

muestreo_deseado /= 3600
cant_muestras = int(len(corriente) / 4)

# batería
bat_ci = 0.1
capacidad = 7.7  # en A/h

# parámetros que no entiendo
a = 9.25e-02  # inversa del cero del filtro RCE
b = 3.8125e-02  # inversa del polo del filtro RCE
r = 5.125e-02  # resistencia interna

A = np.array(
    [[1, 0], [1 - np.exp(-muestreo_deseado / b), np.exp(-muestreo_deseado / b)]]
)

B = np.array([-1, (b - a) * (1 - np.exp(-muestreo_deseado / b)) - 1]) * (
    muestreo_deseado / capacidad
)

tension_min = 3.05  # tensión mínima de descarga

# observador
# más parámetros que no entiendo
R2 = 1
Pb = np.eye(2)  # matriz identidad
Xb = np.array([bat_ci, bat_ci])

# derivada EMF
dV = emf[1:, 0] - emf[:-1, 0]
dX = emf[1:, 1] - emf[:-1, 1]
Df = dV / dX
Df = np.column_stack((emf[1:, 1], Df))

# gráfica de la derivada
# plt.figure()
# plt.plot(Df[:, 0], Df[:, 1])
# plt.grid(True)
# plt.show()

# en python primero se crean interpoladores,
# dsp se les pasan parámetros
interp_emf_pchip = PchipInterpolator(emf[:, 1], emf[:, 0])
interp_emf_linear = interp1d(emf[:, 0], emf[:, 1], kind="linear")
interp_emf_cubic = interp1d(
    emf[:, 1], emf[:, 0], kind="cubic", fill_value="extrapolate"
)
interm_emf_cubic_2 = interp1d(
    emf[:, 0], emf[:, 1], kind="cubic", fill_value="extrapolate"
)
interp_df_cubic = interp1d(Df[:, 0], Df[:, 1], kind="cubic", fill_value="extrapolate")

jk = 1
Vhat = np.zeros(cant_muestras)
Vh = np.zeros(cant_muestras)
Vhat[0] = interp_emf_pchip(bat_ci)

for EKF in range(2):
    xxmin = np.zeros(cant_muestras)
    TT = np.zeros(cant_muestras)
    TE = np.zeros(cant_muestras)
    ne = []
    # ne = np.zeros(cant_muestras)
    Se = np.zeros(cant_muestras)
    XXe = np.zeros(cant_muestras)
    Ps = np.zeros(cant_muestras)
    Px = np.zeros(cant_muestras)
    J = np.zeros(2)
    e = []

    for i in range(1, cant_muestras):
        Xm = interp_emf_linear(tension[i] + corriente[i] * r)
        df = interp_df_cubic(Xb[1])

        if not EKF:
            R1 = np.array([[1, 0], [0, 10]])
            C = np.array([0, 1])
            K = Pb @ C.T * (1.0 / (R2 + C @ Pb @ C.T)) *0
            # K = Pb @ C.T @ np.linalg.inv(1.0 / (R2 + C @ Pb @ C.T))
            Xe = Xb + K * (Xm - Xb[1])
            P = np.linalg.inv(np.linalg.inv(Pb) + C.T @ C / R2)
            Xb = A @ Xe + B * corriente[i]
            Pb = A @ P @ A.T + R1

            Vh[i] = interp_emf_cubic(Xb[1]) - corriente[i] * r
        else:
            R1 = np.array([[0.1, 0], [0, 10]])
            C = np.array([0, df])
            K = Pb @ C.T * 1.0 / (R2 + C @ Pb @ C.T) *0
            # K = Pb @ C.T @ np.linalg.inv(R2 + C @ Pb @ C.T)
            Vhat[i] = interp_emf_pchip(Xb[1]) - corriente[i] * r
            Xe = Xb + K * (tension[i] - Vhat[i])
            P = np.linalg.inv(np.linalg.inv(Pb) + C.T @ C / R2)
            Ps[i] = np.sqrt(P[0, 0])
            Px[i] = np.sqrt(P[1, 1])
            Xb = A @ Xe + B * corriente[i]
            Pb = A @ P @ A.T + R1

        Xb[0] = np.clip(Xb[0], 0, 1)
        Xb[1] = np.clip(Xb[1], 0, 1)

        # Guardar estados
        Se[i] = Xb[0]
        XXe[i] = Xb[1]

        # --- Tiempo remanente real ---
        Te, j = 0, 1
        while (
            corriente[i] > 0
            and (i + j) < cant_muestras
            and tension[i + j] > tension_min
        ):
            Te += 1
            j += 1

        if corriente[i] > 0 and tension[i] > tension_min:
            TE[i] = muestreo_deseado * Te
            xmin = interm_emf_cubic_2(tension_min + corriente[i] * r)
            xxmin[i] = xmin

            # Lambert
            t1 = (xmin - Se[i]) * capacidad / corriente[i] - (b - a)
            t2 = (Se[i] - XXe[i]) * capacidad / corriente[i] + (b - a)
            arg = -(t2 / b) * np.exp(t1 / b)
            m2 = np.argmin(np.abs(lambert[:, 0] - arg))
            w = lambert[m2, 1]
            TT[i] = w * b - t1
            # if TT[i] < 0:
            #     TT[i] = 0

            e.append(TE[i] - TT[i])
            jk += 1

        if i > 0 and corriente[i] < 0 and corriente[i - 1] > 0:
            ne.append(100 * np.linalg.norm(e) / np.sqrt(len(e)))
            e = []
            jk = 1

    # Métrica J
    ne = np.array(ne)
    J[EKF] = np.median(ne[ne < 60])

    # Plot
    plt.figure()
    plt.plot(TE, "--", label="Real")
    plt.plot(TT, "--", label="Prediccion")
    plt.title(f"EKF = {EKF}")
    plt.ylabel("Tiempo remanente")
    plt.legend()
    plt.grid(True)

# --- Resultados finales ---
minJ = np.min(J)
idx = np.argmin(J)

plt.figure(9)
plt.plot(corriente)
plt.title("Registro de corriente del ensayo")

plt.figure(31)
plt.plot(tension[:-1], label="Real")
plt.plot(Vhat, label="EKF")
plt.plot(Vh, label="KF")
plt.legend()
plt.title("Tensiones estimadas")

plt.figure(41)
plt.plot(TE, "--", label="Real")
plt.plot(TT, "--", label="Prediccion")
np.set_printoptions(threshold=sys.maxsize)
print(TT)
plt.title("Tiempo remanente")
plt.legend()
plt.grid(True)

plt.figure(51)
plt.plot(Se, "--", label="Se")
plt.plot(XXe, "--", label="XXe")
plt.title("Tiempo remanente")
plt.legend()
plt.grid(True)

plt.show()

def main():
    pass


if __name__ == "__main__":
    main()
