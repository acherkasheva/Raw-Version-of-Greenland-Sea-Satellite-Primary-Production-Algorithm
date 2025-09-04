import os
import pandas as pd
import numpy as np
import xarray as xr
import statistics
import math
import scipy.interpolate
import scipy.io
import matplotlib.pyplot as plt

def PP_sum_func(chl=None, par=None, gaussian=None, spectra=None, depth=None):
    PP_sum = []
    PP_mean_list = []
    for y in range(25): #25 all years, 0 - 1998 №10
        for x in range(6): #6 all months, 0 - april
            os.chdir('D:\\')
            chl_file = xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc')['CHL'][3 + 12*y + x].to_dataframe().pivot_table(index='lat', columns='lon', values='CHL', dropna=False) #111+x
            chl_file[chl_file>=5] = 5
            date = np.datetime_as_string(xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc')['CHL'][3 + 12*y + x]['time'].values)[:7]
            os.chdir('D:\\DATABASE\\CHL_PAR')
            chl3_list = [i for i in os.listdir() if 'AV' not in i]
            chl3 = ''.join([i for i in chl3_list if i[4:10] == date[:4] + date[5:]])
            chl3_data = xr.open_dataset(chl3).to_dataframe().pivot_table(index='lat', columns='lon', values='CHL1_mean', dropna=False)
            chl3_data.drop(chl3_data.columns[[0]], axis = 1, inplace = True)
            chl3_data[chl3_data>=5] = 5
            par3_list = [i for i in os.listdir() if '_PAR' in i]
            par3 = ''.join([i for i in par3_list if i[4:10] == date[:4] + date[5:]])
            par3_data = xr.open_dataset(par3).to_dataframe().pivot_table(index='lat', columns='lon', values='PAR_mean', dropna=False)
            homepage()
            reanalysis = xr.open_dataset('dswrf.sfc.mon.mean.nc')
            reanalysis_time = [np.datetime_as_string(reanalysis['dswrf']['time'].values[i])[:7] for i in range(len(reanalysis['dswrf']['time']))]
            time_ind = reanalysis_time.index(date[:4] + '-' + date[5:])
            reanalysis_data = reanalysis['dswrf'][time_ind, :15, 72:105]
            flux_lat = list(reanalysis['lat'].values)[:15]
            flux_lon = list(reanalysis['lon'].values)[72:105]
            flux_lon_normalized = list(np.array(flux_lon) - 180)
            flux_pd = reanalysis_data.to_dataframe().pivot_table(index='lat', columns='lon', values='dswrf', dropna=False)
            x_, y_ = np.meshgrid(flux_lon_normalized, flux_lat)  # Initiating the grid of original coordinates
            os.chdir('D:\\')
            xq, yq = np.meshgrid(xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc')['lon'].values.tolist(), xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc')['lat'].values.tolist())  # Initiating the grid of desired coordinates
            f = scipy.interpolate.griddata((y_.ravel(), x_.ravel()), flux_pd.values.ravel(), (yq.ravel(), xq.ravel()))  # Interpolation
            f.resize(xq.shape)  # Restoring the interpolated values into the grid
            f = pd.DataFrame(f)/4.6
            f.index = chl_file.index.to_list()
            f.columns = chl_file.columns.to_list()
            int_file = pd.DataFrame(columns=range(chl_file.shape[1]), index=range(chl_file.shape[0]))
            int_file.index = chl_file.index.to_list()
            int_file.columns = chl_file.columns.to_list()
            os.chdir('D:\\matlabcode')
            PAR = scipy.io.loadmat('PAR_2.mat')
            if x == 0:
                PAR = pd.DataFrame(PAR['PAR_itog_res_april_2022'])
                Zeu = Zeu_final_april
                for i in range(chl_file.shape[0]):
                    for j in range(chl_file.shape[1]):
                        if chl == 4:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_april[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_april[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_april_2[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_april_2[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)), 1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j], 1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)), 1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j], 1) * 10) - 1] else 0
                        if chl == 3:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_april[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_april[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_april_2[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_april_2[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.iloc[i, j] < 0.06:
                                    chl3_data.iloc[i, j] = 0.06  #Or 0?
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)), 1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j], 1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)), 1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j], 1) * 10) - 1] else 0
            if x == 1:
                PAR = pd.DataFrame(PAR['PAR_itog_res_may_2022'])
                Zeu = Zeu_final_may
                for i in range(chl_file.shape[0]):
                    for j in range(chl_file.shape[1]):
                        if chl == 4:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_may[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_may[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_may_2[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_may_2[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)),1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j],1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)),1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j],1) * 10) - 1] else 0
                        if chl == 3:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_may[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_may[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_may_2[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_may_2[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.iloc[i, j] < 0.06:
                                    chl3_data.iloc[i, j] = 0.06  # Or 0?
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)),1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j],1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)),1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j], 1) * 10) - 1] else 0
            if x == 2:
                PAR = pd.DataFrame(PAR['PAR_itog_res_june_2022'])
                Zeu = Zeu_final_june
                for i in range(chl_file.shape[0]):
                    for j in range(chl_file.shape[1]):
                        if chl == 4:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_june[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_june[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_june_2[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_june_2[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)),1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j],1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)),1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j],1) * 10) - 1] else 0
                        if chl == 3:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_june[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_june[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_june_2[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_june_2[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.iloc[i, j] < 0.06:
                                    chl3_data.iloc[i, j] = 0.06  # Or 0?
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)),1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j],1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)),1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j],1) * 10) - 1] else 0
            if x == 3:
                PAR = pd.DataFrame(PAR['PAR_itog_res_july_2022'])
                Zeu = Zeu_final_july
                for i in range(chl_file.shape[0]):
                    for j in range(chl_file.shape[1]):
                        if chl == 4:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_july[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_july[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_july_2[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_july_2[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)),1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j],1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)),1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j],1) * 10) - 1] else 0
                        if chl == 3:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_july[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_july[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_july_2[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_july_2[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.iloc[i, j] < 0.06:
                                    chl3_data.iloc[i, j] = 0.06  # Or 0?
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)),1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j],1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)),1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j],1) * 10) - 1] else 0
            if x == 4:
                PAR = pd.DataFrame(PAR['PAR_itog_res_august_2022'])
                Zeu = Zeu_final_august
                for i in range(chl_file.shape[0]):
                    for j in range(chl_file.shape[1]):
                        if chl == 4:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_august[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_august[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_august_2[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_august_2[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)),1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j],1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)),1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j],1) * 10) - 1] else 0
                        if chl == 3:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_august[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_august[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_august_2[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_august_2[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.iloc[i, j] < 0.06:
                                    chl3_data.iloc[i, j] = 0.06  # Or 0?
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)),1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j],1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)),1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j],1) * 10) - 1] else 0
            if x == 5:
                PAR = pd.DataFrame(PAR['PAR_itog_res_september_2022'])
                Zeu = Zeu_final_september
                for i in range(chl_file.shape[0]):
                    for j in range(chl_file.shape[1]):
                        if chl == 4:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_september[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_september[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_september_2[int(chl_file.iloc[i, j].round(1) * 10)]) if math.isnan(chl_file.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_september_2[int(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl_file.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)), 1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j], 1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl_file.iloc[i, j], 1)),1)].values) if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] >= 70 else 38 * chl_file.iloc[i, j] ** 0.423 if not math.isnan(chl_file.iloc[i, j]) and chl_file.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl_file.iloc[i, j], 1) * 10) - 1] else 0
                        if chl == 3:
                            if gaussian == True:
                                if depth == 0:
                                    int_file.iloc[i, j] = float(CHL_int_september[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_september[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                                if depth == 1:
                                    int_file.iloc[i, j] = float(CHL_int_september_2[int(chl3_data.iloc[i, j].round(1) * 10)]) if math.isnan(chl3_data.iloc[i, j]) is False else 0
                                    if int_file.iloc[i, j] == 0:
                                        int_file.iloc[i, j] = float(CHL_int_september_2[int(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean() * 10) if math.isnan(chl3_data.iloc[i - 1:i + 1, j - 1:j + 1].mean().mean()) is False else 0])
                            if gaussian == False:
                                if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.iloc[i, j] < 0.06:
                                    chl3_data.iloc[i, j] = 0.06  # Or 0?
                                if depth == 0:
                                    int_file.iloc[i, j] = float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)), 1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j],1) * 10) - 1] else 0
                                if depth == 1:
                                    int_file.iloc[i, j] = float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(chl3_data.iloc[i, j], 1)),1)].values) if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] >= 70 else 38 * chl3_data.iloc[i, j] ** 0.423 if not math.isnan(chl3_data.iloc[i, j]) and chl3_data.index.tolist()[i] < 70 and mixed['Mix'][mixed_grid.iloc[i, j]] < Zeu[int(round(chl3_data.iloc[i, j],1) * 10) - 1] else 0
            PAR = PAR.T
            PAR = PAR * 0.0864
            os.chdir('D:\\')
            PAR.index = xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc')['lat'].values
            PAR.columns = xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc')['lon'].values
            par3_data = par3_data.iloc[:, 1:]
            par3_data.index = xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc')['lat'].values
            par3_data.columns = xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc')['lon'].values
            if par == 2 and spectra == 'aph':
                PP = (12/4.6) * int_file * PAR * 0.4068584935483871 * 16 * 0.01
            if par == 2 and spectra == 'od':
                PP = (12/4.6) * int_file * PAR * 0.37380660122121206 * 16 * 0.01
            if par == 2 and spectra == 0:
                PP = (12 / 4.6) * int_file * PAR * 16 * 0.01
            if par == 3 and spectra == 'aph':
                PP = (12 / 4.6) * int_file * par3_data * 0.4068584935483871 * 16 * 0.01
            if par == 3 and spectra == 'od':
                PP = (12 / 4.6) * int_file * par3_data * OD_grid / 100 * 16 * 0.01
            if par == 3 and spectra == 0:
                PP = (12 / 4.6) * int_file * par3_data * 16 * 0.01
            if par == 1 and spectra == 'aph':
                PP = (12 / 4.6) * int_file * f * 0.4068584935483871 * 16 * 0.01
            if par == 1 and spectra == 'od':
                PP = (12 / 4.6) * int_file * f * OD_grid / 100 * 16 * 0.01
            if par == 1 and spectra == 0:
                PP = (12 / 4.6) * int_file * f * 16 * 0.01

            for i in range(PP.shape[0]):
                for j in range(PP.shape[1]):
                    if area_sum_2.iloc[i,j] == -999:
                        PP.iloc[i,j] = 0.0

            if x == 0:
                PP_sum.append((PP * 16 * 10**3 * 30).sum().sum())
                PP_mean_list.append(PP.mean(skipna=True).mean(skipna=True))
            if x == 1:
                PP_sum.append((PP * 16 * 10**3 * 31).sum().sum())
                PP_mean_list.append(PP.mean(skipna=True).mean(skipna=True))
            if x == 2:
                PP_sum.append((PP * 16 * 10**3 * 30).sum().sum())
                PP_mean_list.append(PP.mean(skipna=True).mean(skipna=True))
            if x == 3:
                PP_sum.append((PP * 16 * 10**3 * 31).sum().sum())
                PP_mean_list.append(PP.mean(skipna=True).mean(skipna=True))
            if x == 4:
                PP_sum.append((PP * 16 * 10**3 * 31).sum().sum())
                PP_mean_list.append(PP.mean(skipna=True).mean(skipna=True))
            if x == 5:
                PP_sum.append((PP * 16 * 10**3 * 30).sum().sum())
                PP_mean_list.append(PP.mean(skipna=True).mean(skipna=True))
            c = PP
            name = str([x, y]) + '_new_chl4_par3_without_masking' + '.csv'
            c.to_csv(name)
    a = PP_sum
    b = PP_mean_list
    return a, b