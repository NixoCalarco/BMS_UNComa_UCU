import matplotlib.pyplot as plt
from output import (
    Corriente,
    TE0,
    TT0,
    TE1,
    TT1,
    Se0,
    XXe0,
    Se1,
    XXe1,
    Vh,
    Vhat,
    tension,
)

fig, axs = plt.subplots(2, 3)

# Plot
axs[0, 0].plot(TE0, "--", label="Real")
axs[0, 0].plot(TT0, "--", label="Prediccion")
axs[0, 0].set_title("KF")
axs[0, 0].set_ylabel("Tiempo remanente")
axs[0, 0].legend()
axs[0, 0].grid(True)

axs[0, 1].plot(TE1, "--", label="Real")
axs[0, 1].plot(TT1, "--", label="Prediccion")
axs[0, 1].set_title("EKF")
axs[0, 1].set_ylabel("Tiempo remanente")
axs[0, 1].legend()
axs[0, 1].grid(True)

axs[0, 2].plot(Corriente)
axs[0, 2].set_title("Registro de corriente del ensayo")

axs[1, 0].plot(tension[:-1], label="Real")
axs[1, 0].plot(Vhat, label="EKF")
axs[1, 0].plot(Vh, label="KF")
axs[1, 0].legend()
axs[1, 0].set_title("Tensiones estimadas")

axs[1, 1].plot(Se0, "--", label="Se")
axs[1, 1].plot(XXe0, "--", label="XXe")
axs[1, 1].set_title("KF")
axs[1, 1].legend()
axs[1, 1].grid(True)


axs[1, 2].plot(Se1, "--", label="Se")
axs[1, 2].plot(XXe1, "--", label="XXe")
axs[1, 2].set_title("EKF")
axs[1, 2].legend()
axs[1, 2].grid(True)

plt.show()
