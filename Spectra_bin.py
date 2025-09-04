import matplotlib.pyplot as plt
import dolfyn
import numpy as np
from scipy.signal import find_peaks
from scipy.fft import fft, fftfreq


def nearest(lst, target):
    return min(lst, key=lambda x: abs(x - target))

def spectra_bin(file_name_slave, file_name_master):
    df_slave = dolfyn.read(file_name_slave)
    df_master = dolfyn.read(file_name_master)
    time_slave = df_slave.time.values
    time_master = df_master.time.values

    time_diff_s = [time_slave[i] - time_slave[i - 1] for i in range(1, len(time_slave))]
    time_diff_m = [time_master[i] - time_master[i - 1] for i in range(1, len(time_master))]
    for x in [10, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8, 10**9, 10**10]:
            if len(find_peaks(time_diff_s, x)[0]) <= 20:
                peaks_slave = find_peaks(time_diff_s, x)
                break
    for x in [10, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8, 10**9, 10**10]:
            if len(find_peaks(time_diff_m, x)[0]) <= 20:
                peaks_master = find_peaks(time_diff_m, x)
                break

    N_s = peaks_slave[0][0]
    N_m = peaks_master[0][0]
    T = 1.0 / 16.0

    for x in range(len(peaks_slave[0])):
        if x != len(peaks_slave[0]) - 1:
            fig, ax = plt.subplots(2, 3, figsize=(15, 10))
            plt.title('Slave ' + str(x) + 'bin')
            yf = np.log(fft(df_slave.vel[0][peaks_slave[0][x]:peaks_slave[0][x + 1]].values))
            xf = np.log(fftfreq(N_s, T)[:N_s // 2])
            dy_dx = np.gradient(2.0 / N_s * np.abs(yf[0:N_s // 2]), xf)

            tangent_value_1 = nearest(dy_dx[2000:], -5 / 3)

            yf = np.log(fft(df_slave.vel[1][peaks_slave[0][x]:peaks_slave[0][x + 1]].values))
            xf = np.log(fftfreq(N_s, T)[:N_s // 2])
            dy_dx = np.gradient(2.0 / N_s * np.abs(yf[0:N_s // 2]), xf)

            tangent_value_2 = nearest(dy_dx[2000:], -5 / 3)

            yf = np.log(fft(df_slave.vel[2][peaks_slave[0][x]:peaks_slave[0][x + 1]].values))
            xf = np.log(fftfreq(N_s, T)[:N_s // 2])
            dy_dx = np.gradient(2.0 / N_s * np.abs(yf[0:N_s // 2]), xf)

            tangent_value_3 = nearest(dy_dx[2000:], -5 / 3)
            #x_tangent = xf[slope.index(tangent_value)]
            #y_tangent = 2.0 / N_s * np.abs(yf[0:N_s // 2])[slope.index(tangent_value)]
            #tangent_line = lambda x_1: tangent_value * (x_1 - x_tangent) + y_tangent
            #x_values = np.linspace(xf[slope.index(tangent_value)] - 0.1, xf[slope.index(tangent_value)] + 0.1, 13919)
            #y_values = tangent_line(x_values)
            u_s = np.mean(np.array(
                [g - np.mean(df_slave.vel.values[0][peaks_slave[0][x]:peaks_slave[0][x + 1]]) for g in
                 df_slave.vel.values[0][peaks_slave[0][x]:peaks_slave[0][x + 1]]]) ** 2)
            v_s = np.mean(np.array(
                [g - np.mean(df_slave.vel.values[1][peaks_slave[0][x]:peaks_slave[0][x + 1]]) for g in
                 df_slave.vel.values[1][peaks_slave[0][x]:peaks_slave[0][x + 1]]]) ** 2)
            w_s = np.mean(np.array(
                [g - np.mean(df_slave.vel.values[2][peaks_slave[0][x]:peaks_slave[0][x + 1]]) for g in
                 df_slave.vel.values[2][peaks_slave[0][x]:peaks_slave[0][x + 1]]]) ** 2)
            TKE_s = 1 / 2 * (u_s + v_s + w_s)

            ax[0,0].plot(xf, 2.0 / N_s * np.abs(yf[0:N_s // 2]))
            yf = np.log(fft(df_slave.vel[1][peaks_slave[0][x]:peaks_slave[0][x + 1]].values))
            xf = np.log(fftfreq(N_s, T)[:N_s // 2])
            ax[0,1].plot(xf, 2.0 / N_s * np.abs(yf[0:N_s // 2]))
            yf = np.log(fft(df_slave.vel[2][peaks_slave[0][x]:peaks_slave[0][x + 1]].values))
            xf = np.log(fftfreq(N_s, T)[:N_s // 2])
            ax[0,2].plot(xf, 2.0 / N_s * np.abs(yf[0:N_s // 2]))
            #ax[0,0].axvline(xf[dy_dx.tolist().index(tangent_value_3)], ymin=0, ymax=0.1, color='red')
            ax[1, 0].plot(df_slave.vel[0].values[peaks_slave[0][x]:peaks_slave[0][x + 1]])
            ax[1, 0].plot(df_slave.vel[1].values[peaks_slave[0][x]:peaks_slave[0][x + 1]])
            ax[1, 0].plot(df_slave.vel[2].values[peaks_slave[0][x]:peaks_slave[0][x + 1]])
            ax[1,1].plot(df_slave.vel[0].values)

            ax[1,1].plot(df_slave.vel[1].values)
            ax[1, 1].plot(df_slave.vel[1].values)
            list_ticks = [0, int(len(df_slave.time.to_dataframe())/4), int(len(df_slave.time.to_dataframe())/2), len(df_slave.time.to_dataframe())]
            ax[1, 1].set_xticks(list_ticks)
            ax[1, 1].set_xticklabels([str(df_slave.time.to_dataframe().iloc[0][0].hour) + ':' + str(df_slave.time.to_dataframe().iloc[0][0].minute), str(df_slave.time.to_dataframe().iloc[list_ticks[1]][0].hour) + ':' + str(df_slave.time.to_dataframe().iloc[list_ticks[1]][0].minute),str(df_slave.time.to_dataframe().iloc[list_ticks[2]][0].hour) + ':' + str(df_slave.time.to_dataframe().iloc[list_ticks[2]][0].minute), str(df_slave.time.to_dataframe().iloc[list_ticks[3] - 1][0].hour) + ':' + str(df_slave.time.to_dataframe().iloc[list_ticks[3] - 1][0].minute)], rotation=90)

            ax[1, 1].plot([peaks_slave[0][x], peaks_slave[0][x+1]], [3, 3], 'k-', lw=2)
            ax[1, 2].text(0.05, 0.5, "TKE: " + str(TKE_s), fontsize=12)
            ax[1, 2].text(0.05, 0.6, "Tangent X: " + str(tangent_value_1), fontsize=12)
            ax[1, 2].text(0.05, 0.7, "Tangent Y: " + str(tangent_value_2), fontsize=12)
            ax[1, 2].text(0.05, 0.8, "Tangent Z: " + str(tangent_value_3), fontsize=12)

    for x in range(len(peaks_master[0])):
        if x != len(peaks_master[0]) - 1:
            fig, ax = plt.subplots(2, 3, figsize=(15, 10))
            plt.title('Slave ' + str(x) + 'bin')
            yf = np.log(fft(df_master.vel[0][peaks_master[0][x]:peaks_master[0][x + 1]].values))
            xf = np.log(fftfreq(N_m, T)[:N_m // 2])
            dy_dx = np.gradient(2.0 / N_m * np.abs(yf[0:N_m // 2]), xf)

            tangent_value_1 = nearest(dy_dx[2000:], -5 / 3)

            yf = np.log(fft(df_master.vel[1][peaks_master[0][x]:peaks_master[0][x + 1]].values))
            xf = np.log(fftfreq(N_m, T)[:N_m // 2])
            dy_dx = np.gradient(2.0 / N_m * np.abs(yf[0:N_m // 2]), xf)

            tangent_value_2 = nearest(dy_dx[2000:], -5 / 3)

            yf = np.log(fft(df_master.vel[2][peaks_master[0][x]:peaks_master[0][x + 1]].values))
            xf = np.log(fftfreq(N_m, T)[:N_m // 2])
            dy_dx = np.gradient(2.0 / N_m * np.abs(yf[0:N_m // 2]), xf)

            tangent_value_3 = nearest(dy_dx[2000:], -5 / 3)
            # x_tangent = xf[slope.index(tangent_value)]
            # y_tangent = 2.0 / N_s * np.abs(yf[0:N_s // 2])[slope.index(tangent_value)]
            # tangent_line = lambda x_1: tangent_value * (x_1 - x_tangent) + y_tangent
            # x_values = np.linspace(xf[slope.index(tangent_value)] - 0.1, xf[slope.index(tangent_value)] + 0.1, 13919)
            # y_values = tangent_line(x_values)
            u_s = np.mean(np.array(
                [g - np.mean(df_master.vel.values[0][peaks_master[0][x]:peaks_master[0][x + 1]]) for g in
                 df_master.vel.values[0][peaks_master[0][x]:peaks_master[0][x + 1]]]) ** 2)
            v_s = np.mean(np.array(
                [g - np.mean(df_master.vel.values[1][peaks_master[0][x]:peaks_master[0][x + 1]]) for g in
                 df_master.vel.values[1][peaks_master[0][x]:peaks_master[0][x + 1]]]) ** 2)
            w_s = np.mean(np.array(
                [g - np.mean(df_master.vel.values[2][peaks_master[0][x]:peaks_master[0][x + 1]]) for g in
                 df_master.vel.values[2][peaks_master[0][x]:peaks_master[0][x + 1]]]) ** 2)
            TKE_s = 1 / 2 * (u_s + v_s + w_s)

            ax[0, 0].plot(xf, 2.0 / N_m * np.abs(yf[0:N_m // 2]))
            yf = np.log(fft(df_master.vel[1][peaks_master[0][x]:peaks_master[0][x + 1]].values))
            xf = np.log(fftfreq(N_m, T)[:N_m // 2])
            ax[0, 1].plot(xf, 2.0 / N_m * np.abs(yf[0:N_m // 2]))
            yf = np.log(fft(df_master.vel[2][peaks_master[0][x]:peaks_master[0][x + 1]].values))
            xf = np.log(fftfreq(N_m, T)[:N_m // 2])
            ax[0, 2].plot(xf, 2.0 / N_m * np.abs(yf[0:N_m // 2]))
            #ax[0,0].axvline(xf[dy_dx.tolist().index(tangent_value_3)], ymin=0, ymax=0.1, color='red')
            ax[1, 0].plot(df_master.vel[0].values[peaks_master[0][x]:peaks_master[0][x + 1]])
            ax[1, 0].plot(df_master.vel[1].values[peaks_master[0][x]:peaks_master[0][x + 1]])
            ax[1, 0].plot(df_master.vel[2].values[peaks_master[0][x]:peaks_master[0][x + 1]])
            ax[1, 1].plot(df_master.vel[0].values)
            ax[1, 1].plot(df_master.vel[1].values)
            ax[1, 1].plot(df_master.vel[1].values)
            ax[1, 1].plot([peaks_master[0][x], peaks_master[0][x + 1]], [3, 3], 'k-', lw=2)
            list_ticks = [0, int(len(df_master.time.to_dataframe()) / 4), int(len(df_master.time.to_dataframe()) / 2), len(df_master.time.to_dataframe())]
            ax[1, 1].set_xticks(list_ticks)
            ax[1, 1].set_xticklabels([str(df_master.time.to_dataframe().iloc[0][0].hour) + ':' + str(df_master.time.to_dataframe().iloc[0][0].minute), str(df_master.time.to_dataframe().iloc[list_ticks[1]][0].hour) + ':' + str(df_master.time.to_dataframe().iloc[list_ticks[1]][0].minute),str(df_master.time.to_dataframe().iloc[list_ticks[2]][0].hour) + ':' + str(df_master.time.to_dataframe().iloc[list_ticks[2]][0].minute), str(df_master.time.to_dataframe().iloc[list_ticks[3] - 1][0].hour) + ':' + str(df_master.time.to_dataframe().iloc[list_ticks[3] - 1][0].minute)], rotation=90)
            ax[1, 2].text(0.05, 0.5, "TKE: " + str(TKE_s), fontsize=12)
            ax[1, 2].text(0.05, 0.6, "Tangent X: " + str(tangent_value_1), fontsize=12)
            ax[1, 2].text(0.05, 0.7, "Tangent Y: " + str(tangent_value_2), fontsize=12)
            ax[1, 2].text(0.05, 0.8, "Tangent Z: " + str(tangent_value_3), fontsize=12)




