import numpy as np
import serial
from SOC_OCV_curve_fit import *
from matplotlib import pyplot as plt

arduino = serial.Serial(port='COM9', baudrate=115200, timeout=0.1)

fig = plt.figure(figsize=(12, 10))
plot = fig.add_subplot(111)

plt.ylim([0, 110])


def read_data(ser):
    line = ser.readline()
    line = line.decode("utf-8")
    data = line.strip()
    return data


t = 0
header = ["time", "Current", "Voltage", "SOC"]
plotting_data = [[0, 0, 0, 100]]
SOC_Idt = 1

with open("Kalman_Filter_Arduino1.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(header)

while True:
    raw_data = read_data(arduino)
    if raw_data:
        current_data = [t]
        t += 5
        temp = raw_data.split("$")
        U = float(temp[0])
        V = float(temp[1])
        SOC = float(temp[2])

        plotting_data.append([t, U, V, SOC])

        if t % 900 == 0:
            with open("Kalman_Filter_Arduino1.csv", "a") as f:
                writer = csv.writer(f)
                writer.writerows(plotting_data)
            plotting_data = []

        if t == 10:
            SOC_Idt = SOC
        elif t > 10:
            SOC_Idt -= 0.9*U*5.00/10800

        plot.scatter([t], [SOC * 100], color="orange", s=50)
        plot.scatter([t], [SOC_Idt * 100], color="blue", s=50)
        print(f"t : {t}    SOC : {SOC}   Current : {U}   Voltage : {V}")
        plt.pause(0.001)

plt.show()
