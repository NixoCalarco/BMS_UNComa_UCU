#ifndef KALMAN_H
#define KALMAN_H

// estructuras para pasar resultados entre tasks
enum kalman_alg { NORMAL, EXTENDED };

typedef struct {
    enum kalman_alg algorithm;
    double V, Se, XXe, TT, TE;
} KalmanVals;

void kalman(void (*send_values)(KalmanVals vals));

#endif
