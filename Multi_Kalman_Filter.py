import numpy as np
import serial
from SOC_OCV_curve_fit import *
from matplotlib import pyplot as plt

arduino = serial.Serial(port='COM9', baudrate=115200, timeout=0.1)

SOC_initial = 1

Qn = 3 * 3600
delta_t = 5.0

SOC_OCV = best_fit()
dSOC_OCV = np.polyder(SOC_OCV)

# R0 = 1.5581e-5
# R1 = 0.0031
# C1 = 1.8565e4

R0 = 1.611e-4
R1 = 0.116
C1 = 10.386

tau_1 = R1 * C1
a1 = np.exp(-delta_t / tau_1)
b1 = R1 * (1 - np.exp(-delta_t / tau_1))

A = np.array([[1, 0], [0, a1]])
B = np.array([[-1 * 0.9 / (Qn * 3600)], [b1]])
D = np.array([-R0])

fig = plt.figure(figsize=(12, 10))
plot = fig.add_subplot(111)

plt.ylim([0, 100])


def read_data(ser):
    line = ser.readline()
    line = line.decode("utf-8")
    data = line.strip()
    return data


class Extended_Kalman_Filter:
    def __init__(self):
        self.X = np.array([[SOC_initial], [0]])
        self.Q = np.array([[1e-5, 0], [0, 1e-4]])
        self.R = np.array([0.004]).reshape([1, 1])
        self.P = np.array([[0.025, 0], [0, 0.01]])

    def update(self, i, v):
        SOC = self.X[0][0]
        V1 = self.X[1][0]
        OCV = np.polyval(SOC_OCV, SOC)
        Vt = OCV - V1 - (i * R0)
        dOCV = np.polyval(dSOC_OCV, SOC)
        C = np.array([dOCV, -1]).reshape([1, 2])
        error = v - Vt
        SOC_hat = self.X[0][0]
        self.X = np.dot(A, self.X) + np.dot(B, i)
        self.P = (A @ self.P @ A.T) + self.Q
        KalmanGain = (self.P @ C.T) / ((C @ self.P @ C.T) + self.R)
        self.X = self.X + (KalmanGain * error)
        self.P = (np.ones(shape=[2, 2]) - (KalmanGain @ C)) @ self.P

        current_data.append(SOC_hat)
        plot.scatter([t], [SOC_hat * 100], color="blue", s=20)


class Adaptive_Extended_Kalman_Filter:
    def __init__(self):
        self.X = np.array([[SOC_initial], [0]])
        self.Q = np.array([[1e-5, 0], [0, 1e-4]])
        self.R = np.array([0.004]).reshape([1, 1])
        self.P = np.array([[0.025, 0], [0, 0.01]])

    def update(self, i, v):
        SOC = self.X[0][0]
        V1 = self.X[1][0]
        OCV = np.polyval(SOC_OCV, SOC)
        Vt = OCV - V1 - (i * R0)
        dOCV = np.polyval(dSOC_OCV, SOC)
        C = np.array([dOCV, -1]).reshape([1, 2])
        error = v - Vt
        SOC_hat = self.X[0][0]
        self.X = np.dot(A, self.X) + np.dot(B, i)
        self.P = (A @ self.P @ A.T) - self.Q
        KalmanGain = (self.P @ C.T) / ((C @ self.P @ C.T) + self.R)
        self.X = self.X + (KalmanGain * error)
        self.P = (np.ones(shape=[2, 2]) - (KalmanGain @ C)) @ self.P
        error_matrix = np.array([error]).reshape([1, 1])
        self.Q = KalmanGain @ error_matrix @ KalmanGain.T

        current_data.append(SOC_hat)
        plot.scatter([t], [SOC_hat * 100], color="green", s=20)


class Iterative_Extended_Kalman_Filter:
    def __init__(self):
        self.X = np.array([[SOC_initial], [0]])
        self.Q = np.array([[1e-5, 0], [0, 1e-4]])
        self.R = np.array([0.004]).reshape([1, 1])
        self.P = np.array([[0.025, 0], [0, 0.01]])
        self.X_j = []
        self.N = 10
        self.X_initial = None

    def update(self, i, v):
        for j in range(1, self.N+1):
            SOC = self.X[0][0]
            V1 = self.X[1][0]
            OCV = np.polyval(SOC_OCV, SOC)
            Vt = OCV - V1 - (i * R0)
            dOCV = np.polyval(dSOC_OCV, SOC)
            C = np.array([dOCV, -1]).reshape([1, 2])

            self.X = np.dot(A, self.X) + np.dot(B, i)
            self.P = (A @ self.P @ A.T) - self.Q

            if j == 1:
                self.X_initial = self.X
                self.X_j.append(self.X_initial)

            error = v - Vt - np.dot(np.array([np.polyval(dSOC_OCV, SOC), -1]), self.X_initial - self.X)
            KalmanGain = (self.P @ C.T) / ((C @ self.P @ C.T) + self.R)
            self.X = self.X_initial + (KalmanGain * error)
            self.P = (np.ones(shape=[2, 2]) - (KalmanGain @ C)) @ self.P

            self.X_j.append(self.X)

            if abs(self.X_j[-1][0][0] - self.X_j[-2][0][0]) < 0.01:
                self.X = self.X_j[-1]
            else:
                self.X = self.X_j[-2]
                break

        SOC_hat = self.X[0][0]
        current_data.append(SOC_hat)
        # plot.scatter([t], [SOC_hat * 100], color="red", s=20)


class Q_Additive_Kalman_Filter:
    def __init__(self):
        self.X = np.array([[SOC_initial], [0]])
        self.Q = np.array([[1e-5, 0], [0, 1e-4]])
        self.R = np.array([0.004]).reshape([1, 1])
        self.P = np.array([[0.025, 0], [0, 0.01]])
        self.N = 5
        self.eta = []
        self.Cn = None

    def update(self, i, v):
        self.Cn = 0
        SOC = self.X[0][0]
        V1 = self.X[1][0]
        OCV = np.polyval(SOC_OCV, SOC)
        Vt = OCV - V1 - (i * R0)
        dOCV = np.polyval(dSOC_OCV, SOC)
        C = np.array([dOCV, -1]).reshape([1, 2])
        error = v - Vt
        SOC_hat = self.X[0][0]
        self.X = np.dot(A, self.X) + np.dot(B, i)
        self.P = (A @ self.P @ A.T) - self.Q
        KalmanGain = (self.P @ C.T) / ((C @ self.P @ C.T) + self.R)
        self.X = self.X + (KalmanGain * error)
        self.P = (np.ones(shape=[2, 2]) - (KalmanGain @ C)) @ self.P
        error_matrix = np.array([error]).reshape([1, 1])
        self.eta.append(error)

        if int(len(self.eta)) < 5:
            self.Q = KalmanGain @ error_matrix @ KalmanGain.T

        else:
            self.Cn = self.Cn + sum(self.eta[-5:])
            self.Cn = self.Cn/self.N
            self.Q = KalmanGain @ np.array([self.Cn]).reshape([1, 1]) @ KalmanGain.T

        current_data.append(SOC_hat)
        plot.scatter([t], [SOC_hat * 100], color="orange", s=20)


t = 0
header = ["time", "SOC_EKF", "SOC_AEKF", "SOC_IEKF", "SOC_Q_Additive", "Current"]
# plotting_data = np.array([header, [t, SOC_initial, SOC_initial, SOC_initial, SOC_initial, 0]])


EKF = Extended_Kalman_Filter()
AEKF = Adaptive_Extended_Kalman_Filter()
IEKF = Iterative_Extended_Kalman_Filter()
Q_adapt = Q_Additive_Kalman_Filter()


while True:
    raw_data = read_data(arduino)
    if raw_data:
        current_data = [t]
        t += 5
        U = float(raw_data.split("$")[0])
        V = float(raw_data.split("$")[1])
        EKF.update(U, V)
        AEKF.update(U, V)
        IEKF.update(U, V)
        Q_adapt.update(U, V)

        current_data.append(U)
        # np.append(plotting_data, np.array([[current_data]]))
        '''
        if t % 1800 == 0:
            with open("Multi_Kalman_Filter_data.csv", "w") as f:
                writer = csv.writer(f)
                writer.writerows(plotting_data)
        '''

        print(f"Time: {t} , SOC_EKF: {current_data[1]} , SOC_AEKF: {current_data[2]} , SOC_IEKF: {current_data[3]} , SOC_QAdpt: {current_data[4]} , Current: {U} , Voltage: {V}")
        # print(f"Time: {t} , SOC_EKF: {current_data[1]} , Current: {U} , Voltage: {V}")
        plt.pause(0.001)

plt.show()
