import os
import pandas as pd
import numpy as np
import xarray as xr
import statistics
import math
import scipy.interpolate
import scipy.io
import matplotlib.pyplot as plt

def PP_model(year_input, month_input, chl_level, par_type, gaussian=True, spectra=None, depth = None):

    def nearest(lst, target):
        return min(lst, key=lambda x: abs(x - target))

    year_n = [int(x) for x in year_input]
    month_n = [int(x) for x in month_input]
    os.chdir('C:\\Users\\iopan\\PycharmProjects\\untitled')
    df = xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_1671189154637.nc')
    chl = df['CHL']
    chl_time = pd.DataFrame(df['time'].values)

    df_wider = xr.open_dataset('cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_1678284524388.nc')
    chl_wider = df_wider['CHL']
    chl_time_wider = pd.DataFrame(df_wider['time'].values)
    lon_wider = chl_wider['lon'].values.tolist()
    lat_wider = chl_wider['lat'].values.tolist()

    df_pp = xr.open_dataset('cmems_obs-oc_glo_bgc-pp_my_l4-multi-4km_P1M_1671595284205.nc')
    pp = df_pp['PP']
    pp_time = pd.DataFrame(df_pp['time'].values)

    df_pp_wider = xr.open_dataset('cmems_obs-oc_glo_bgc-pp_my_l4-multi-4km_P1M_1678467079124.nc')
    pp_wider = df_pp_wider['PP']
    pp_time_wider = pd.DataFrame(df_pp_wider['time'].values)

    reanalysis = xr.open_dataset('dswrf.sfc.mon.mean.nc') # Importing the Flux data
    reanalysis_time = [np.datetime_as_string(reanalysis['dswrf']['time'].values[i])[:7] for i in range(len(reanalysis['dswrf']['time']))] # Creating a list of dates of the Flux data
    y_data = pd.DataFrame(year_to_check) # Converting a list of years of the field points
    m_data = pd.DataFrame(month_to_check) # Converting a list of months of the field points
    df = pd.DataFrame({'PP Field': []}) # Pre-locating a database for the output
    for n in range(len(year_n)): # For every year of input
        for m in range(len(month_n)): # For every month of input
            if month_n[m] in date_dict[year_n[n]]: # Condition for a month to be connected with a year in the field points  P.S. update date_dict
                y_l = list(y_data[y_data[0] == year_n[n]].index) # Looking for indexes of years in the field points
                m_l = list(m_data[m_data[0] == month_n[m]].index) # Looking for indexes of months in the field points
                common_index = list(set(y_l) & set(m_l)) # Searching for mutual indexes of months and years
                lat_point = [lat_to_check_normalized[x] for x in common_index] # Storing the latitudes of common indexes
                lon_point = [lon_to_check_normalized[x] for x in common_index] # Storing the longitudes of common indexes
                Zeu = [Zeu_final_april if month_n[m] == 4 else Zeu_final_may if month_n[m] == 5 else Zeu_final_june if month_n[m] == 6 else Zeu_final_july if month_n[m] == 7 else Zeu_final_august if month_n[m] == 8 else Zeu_final_september if month_n[m] == 9 else 0][0] # Storing the Euphotic depth profile for the current year and month
                closest_mixed_index = [[abs(lat_point[y] - mixed['Lat'][x]) + abs(lon_point[y] - mixed['Lon'][x]) for x in range(mixed['Lon'].shape[0])].index(min([abs(lat_point[y] - mixed['Lat'][x]) + abs(lon_point[y] - mixed['Lon'][x]) for x in range(mixed['Lon'].shape[0])])) for y in range(len(common_index))] # Searching for a closest coordinate in mixed-layer depth data
                closest_OD_index = [[abs(lat_point[y] - OD_mean['Lat'][x]) + abs(lon_point[y] - OD_mean['Lon'][x]) for x in range(OD_mean['Lon'].shape[0])].index(min([abs(lat_point[y] - OD_mean['Lat'][x]) + abs(lon_point[y] - OD_mean['Lon'][x]) for x in range(OD_mean['Lon'].shape[0])])) for y in range(len(common_index))] # Searching for a closest coordinate in mixed-layer depth data

                PAR_hermes = xr.open_dataset('L3m_' + str(year_n[n]) + '0' + str(month_n[m]) + '01_PAR.nc')['PAR_mean'].to_dataframe().pivot_table(index='lat', columns='lon', values='PAR_mean', dropna=False)  # Storing a PAR Hermes data
                PAR_hermes.index = lat_sat  # Latitude calibration
                PAR_hermes.columns = lon_sat

                date = chl_time[chl_time[0] == str(year_n[n]) + '-0' + str(month_n[m]) + '-01'].index[0] # Saving the date of a CHL Copernic for further processes
                chl_df = chl[date].to_dataframe().pivot_table(index='lat', columns='lon', values='CHL', dropna=False) # Storing a CHL Copernic database
                chl_df[chl_df > 5] = 5 # Filtering the CHL database as there is no values higher than 5 in gaussian constants
                chl_df.index = lat_sat # Latitude calibration
                chl_df.columns = lon_sat # Longitude calibration

                date_wider = chl_time_wider[chl_time_wider[0] == str(year_n[n]) + '-0' + str(month_n[m]) + '-01'].index[0]  # Saving the date of a CHL Copernic for further processes
                chl_df_wider = chl_wider[date_wider].to_dataframe().pivot_table(index='lat', columns='lon', values='CHL', dropna=False)  # Storing a CHL Copernic database
                chl_df_wider[chl_df_wider > 5] = 5  # Filtering the CHL database as there is no values higher than 5 in gaussian constants

                chl_to_int = [chl_df.loc[lat_sat[lat_sat.index(lat_to_check_normalized[x]) - 1] : lat_sat[lat_sat.index(lat_to_check_normalized[x]) + 1],lon_sat[lon_sat.index(lon_to_check_normalized[x]) - 1] : lon_sat[lon_sat.index(lon_to_check_normalized[x]) + 1]].values if np.isnan(chl_df.loc[lat_to_check_normalized[x]][lon_to_check_normalized[x]]) else chl_df.loc[lat_to_check_normalized[x]][lon_to_check_normalized[x]] if lon_to_check_normalized[x] >= -18 and lon_to_check_normalized[x] <= 15 else chl_df_wider.loc[nearest(lat_wider, lat_to_check_normalized[x])][nearest(lon_wider, lon_to_check_normalized[x])] for x in common_index] # Either storing surrounding cells of Chl values if current one is Nan, or storing an actual Chl value
                chl_to_int = [np.nanmean(chl_to_int[x]) for x in range(len(common_index))] # Calculating Chl mean values of surrounding cells

                if gaussian == True: # If gaussian is chosen
                    if depth == 0:
                        CHL_int = [float(CHL_int_june[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 6 else float(CHL_int_july[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 7 else float(CHL_int_august[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 8 else float(CHL_int_september[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 9 else float(CHL_int_april[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 4 else float(CHL_int_may[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 5 else 0 for x in chl_to_int] # Searching for a value in gaussian profile
                    if depth == 1:
                        CHL_int = [float(CHL_int_june_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 6 else float(CHL_int_july_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 7 else float(CHL_int_august_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 8 else float(CHL_int_september_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 9 else float(CHL_int_april_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 4 else float(CHL_int_may_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 5 else 0 for x in chl_to_int]
                if gaussian == False: # If linear is chosen
                    if depth == 0:
                        CHL_int = [float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(x, 1)), 1)].values) if not math.isnan(x) and lat_point[ind] >= 70 else 38 * x ** 0.423 if not math.isnan(x) and lat_point[ind] < 70 and mixed['Mixed'][closest_mixed_index[ind]] < Zeu[int(round(x, 1) * 10) - 1] else 0 for ind, x in enumerate(chl_to_int)] # Either using a linear profile or an equation for low latitudes in case if Euphotic depth is higher than mixed-layer depth
                    if depth == 1:
                        CHL_int = [float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(x, 1)), 1)].values) if not math.isnan(x) and lat_point[ind] >= 70 else 38 * x ** 0.423 if not math.isnan(x) and lat_point[ind] < 70 and mixed['Mixed'][closest_mixed_index[ind]] < Zeu[int(round(x, 1) * 10) - 1] else 0 for ind, x in enumerate(chl_to_int)] # Either using a linear profile or an equation for low latitudes in case if Euphotic depth is higher than mixed-layer depth

                os.chdir('D:\\DATABASE\\HermesCHL') # Changing the path to the folder where all CHL Hermes files are located
                chl3_list = [i for i in os.listdir() if 'AV' not in i] # Filtering CHL Hermes files
                chl3 = ''.join([i for i in chl3_list if i[4:10] == str(year_n[n]) + '0' + str(month_n[m])]) # Storing the name of the CHL Hermes file depending on a date
                chl3_data = xr.open_dataset(chl3) # Storing a CHL Hermes data
                chl3_data = chl3_data.to_dataframe().pivot_table(index='lat', columns='lon', values='CHL1_mean', dropna=False) # Converting the CHL Hermes data into the grid
                chl3_data.index = lat_sat  # Latitude calibration
                chl3_data.columns = lon_sat  # Longitude calibration
                chl3_data[chl3_data < 0.06] = 0.06

                os.chdir('D:\\DATABASE\\HermesCHL\\Wider Net')
                chl3_wider_list = [i for i in os.listdir() if 'AV' not in i]  # Filtering CHL Hermes files
                chl3_wider = ''.join([i for i in chl3_wider_list if i[4:10] == str(year_n[n]) + '0' + str(month_n[m])])  # Storing the name of the CHL Hermes file depending on a date
                chl3_wider_data = xr.open_dataset(chl3_wider)  # Storing a CHL Hermes data
                chl3_wider_data = chl3_wider_data.to_dataframe().pivot_table(index='lat', columns='lon', values='CHL1_mean', dropna=False)  # Converting the CHL Hermes data into the grid
                chl3_wider_data.index = lat_wider
                chl3_wider_data.columns = lon_wider

                os.chdir('C:\\Users\\iopan\\PycharmProjects\\untitled') # Changing back the folder path
                chl3_data[chl3_data>5] = 5
                chl3_wider_data[chl3_wider_data > 5] = 5
                chl3_to_int = [chl3_data.loc[lat_to_check_normalized[x]][lon_to_check_normalized[x]] if lon_to_check_normalized[x] >= -18.0 and lon_to_check_normalized[x] <= 15.0 else chl3_wider_data.loc[nearest(lat_wider, lat_to_check_normalized[x])][nearest(lon_wider, lon_to_check_normalized[x])] for x in common_index] # Searching for CHL Hermes values
                if gaussian == True:  # If gaussian is chosen
                    if depth == 0:
                        CHL3_int = [float(CHL_int_june[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 6 else float(CHL_int_july[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 7 else float(CHL_int_august[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 8 else float(CHL_int_september[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 9 else float(CHL_int_april[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 4 else float(CHL_int_may[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 5 else 0 for x in chl3_to_int] # Searching for a value in the gaussian profile
                    if depth == 1:
                        CHL3_int = [float(CHL_int_june_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 6 else float(CHL_int_july_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 7 else float(CHL_int_august_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 8 else float(CHL_int_september_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 9 else float(CHL_int_april_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 4 else float(CHL_int_may_2[int(x.round(1) * 10)]) if math.isnan(x) == False and month_n[m] == 5 else 0 for x in chl3_to_int]
                else:  # If linear is chosen
                    if depth == 0:
                        CHL3_int = [float(universal_pd['Int'][round(universal_pd['ind'], 1) == round(float(round(x, 1)), 1)].values) if not math.isnan(x) else 0 for x in chl3_to_int] # Using the linear profile to store the CHL value
                    if depth == 1:
                        CHL3_int = [float(universal_pd_2['Int'][round(universal_pd_2['ind'], 1) == round(float(round(x, 1)), 1)].values) if not math.isnan(x) else 0 for x in chl3_to_int]
                time_ind = reanalysis_time.index(str(year_n[n]) + '-0' + str(month_n[m])) # Storing the date for a Flux data
                reanalysis_data = reanalysis['dswrf'][time_ind,2:13,82:106] # Storing the Flux data 2:13,86:105
                flux_lat = list(reanalysis['lat'].values)  # Storing the original latitude values
                flux_lon = list(reanalysis['lon'].values)  # Storing the original longitude values
                flux_lat = flux_lat[2:13]  # Filtering stored latitude values by area of interest
                flux_lon = flux_lon[82:106]  # Filtering stored longitude values by area of interest
                flux_pd = reanalysis_data.to_dataframe()  # Preparation to conversion obtained data frame into the grid format
                flux_pd = flux_pd.pivot_table(index='lat', columns='lon', values='dswrf', dropna=False)  # Conversion obtained data frame into the grid format
                flux_lon_normalized = list(np.array(flux_lon)-180)
                x_, y_ = np.meshgrid(flux_lon_normalized, flux_lat)  # Initiating the grid of original coordinates
                xq, yq = np.meshgrid(lon_sat, lat_sat)  # Initiating the grid of desired coordinates
                f = scipy.interpolate.griddata((y_.ravel(), x_.ravel()), flux_pd.values.ravel(), (yq.ravel(), xq.ravel()))  # Interpolation
                f.resize(xq.shape)  # Restoring the interpolated values into the grid
                f = pd.DataFrame(f) # Converting obtained values into the grid
                f.index = lat_sat # Latitude calibration
                f.columns = lon_sat # Longitude calibration

                pp_chl4_par3 = [(12/4.6)*CHL_int[ind]*PAR_hermes.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 for ind, x in enumerate(common_index)] # Calculating PP using CHL Copernicus and PAR Hermes
                pp_chl3_par3 = [(12/4.6)*CHL3_int[ind]*PAR_hermes.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 for ind, x in enumerate(common_index)] # Calculating PP using CHL Hermes and PAR Hermes
                pp_chl4_par2 = [(12/4.6)*CHL_int[ind]*PAR_april.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 4 else (12/4.6)*CHL_int[ind]*PAR_may.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 5 else (12/4.6)*CHL_int[ind]*PAR_june.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 6 else (12/4.6)*CHL_int[ind]*PAR_july.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 7 else (12/4.6)*CHL_int[ind]*PAR_august.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 8 else (12/4.6)*CHL_int[ind]*PAR_september.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 for ind, x in enumerate(common_index)] # Calculating PP using CHL Copernicus and PAR Eumetsat
                pp_chl3_par2 = [(12/4.6)*CHL3_int[ind]*PAR_april.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 4 else (12/4.6)*CHL3_int[ind]*PAR_may.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 5 else (12/4.6)*CHL3_int[ind]*PAR_june.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 6 else (12/4.6)*CHL3_int[ind]*PAR_july.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 7 else (12/4.6)*CHL3_int[ind]*PAR_august.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 8 else (12/4.6)*CHL3_int[ind]*PAR_september.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*16*0.01 if month_n[m] == 9 else 0 for ind, x in enumerate(common_index)] # Calculating PP Using CHL Hermes and PAR Eumetsat
                pp_chl4_f = [(12/4.6)*CHL_int[ind]*f.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]/4.6*16*0.01*0.43 for ind, x in enumerate(common_index)] # Calculating PP using CHL Copernicus and Flux
                pp_chl3_f = [(12/4.6)*CHL3_int[ind]*f.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]/4.6*16*0.01*0.43 for ind, x in enumerate(common_index)] # Calculating PP using CHL Hermes and Flux

                pp_chl4_par2_multi_by_aph = [(12/4.6)*CHL_int[ind]*PAR_april.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 4 else (12/4.6)*CHL_int[ind]*PAR_may.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 5 else (12/4.6)*CHL_int[ind]*PAR_june.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 6 else (12/4.6)*CHL_int[ind]*PAR_july.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 7 else (12/4.6)*CHL_int[ind]*PAR_august.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 8 else (12/4.6)*CHL_int[ind]*PAR_september.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 for ind, x in enumerate(common_index)] # Calculating PP using CHL Copernicus and PAR Eumetsat
                pp_chl4_par3_multi_by_aph = [(12/4.6)*CHL_int[ind]*PAR_hermes.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 for ind, x in enumerate(common_index)] # Calculating PP using CHL Copernicus and PAR Hermes
                pp_chl4_f_multi_by_aph = [(12/4.6)*CHL_int[ind]*f.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]/4.6*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01*0.43 for ind, x in enumerate(common_index)] # Calculating PP using CHL Copernicus and Flux
                pp_chl3_par2_multi_by_aph = [(12/4.6)*CHL3_int[ind]*PAR_april.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 4 else (12/4.6)*CHL3_int[ind]*PAR_may.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 5 else (12/4.6)*CHL3_int[ind]*PAR_june.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 6 else (12/4.6)*CHL3_int[ind]*PAR_july.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 7 else (12/4.6)*CHL3_int[ind]*PAR_august.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 8 else (12/4.6)*CHL3_int[ind]*PAR_september.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 if month_n[m] == 9 else 0 for ind, x in enumerate(common_index)] # Calculating PP Using CHL Hermes and PAR Eumetsat
                pp_chl3_par3_multi_by_aph = [(12/4.6)*CHL3_int[ind]*PAR_hermes.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01 for ind, x in enumerate(common_index)] # Calculating PP using CHL Hermes and PAR Hermes
                pp_chl3_f_multi_by_aph = [(12/4.6)*CHL3_int[ind]*f.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]/4.6*OD_mean['Mean'][closest_OD_index[ind]]/100*16*0.01*0.43 for ind, x in enumerate(common_index)] # Calculating PP using CHL Hermes and Flux

                pp_chl4_par2_spectra_mean = [(12 / 4.6)*CHL_int[ind]*PAR_april.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 4 else (12/4.6)*CHL_int[ind]*PAR_may.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 5 else (12 / 4.6)*CHL_int[ind]*PAR_june.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 6 else (12/4.6)*CHL_int[ind]*PAR_july.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 7 else (12/4.6)*CHL_int[ind]*PAR_august.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 8 else (12/4.6)*CHL_int[ind]*PAR_september.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 for ind, x in enumerate(common_index)]
                pp_chl4_par3_spectra_mean = [(12/4.6)*CHL_int[ind]*PAR_hermes.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 for ind, x in enumerate(common_index)] # Calculating PP using CHL Copernicus and PAR Hermes
                pp_chl4_f_spectra_mean = [(12/4.6)*CHL_int[ind]*f.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]/4.6*0.4068584935483871*16*0.01*0.43 for ind, x in enumerate(common_index)] # Calculating PP using CHL Copernicus and Flux
                pp_chl3_par2_spectra_mean = [(12/4.6)*CHL3_int[ind]*PAR_april.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 4 else (12/4.6)*CHL3_int[ind]*PAR_may.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 5 else (12/4.6)*CHL3_int[ind]*PAR_june.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 6 else (12/4.6)*CHL3_int[ind]*PAR_july.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 7 else (12/4.6)*CHL3_int[ind]*PAR_august.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 8 else (12/4.6)*CHL3_int[ind]*PAR_september.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 if month_n[m] == 9 else 0 for ind, x in enumerate(common_index)] # Calculating PP Using CHL Hermes and PAR Eumetsat
                pp_chl3_par3_spectra_mean = [(12/4.6)*CHL3_int[ind]*PAR_hermes.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]*0.4068584935483871*16*0.01 for ind, x in enumerate(common_index)] # Calculating PP using CHL Hermes and PAR Hermes
                pp_chl3_f_spectra_mean = [(12/4.6)*CHL3_int[ind]*f.loc[nearest(lat_sat, lat_to_check_normalized[x])][nearest(lon_sat, lon_to_check_normalized[x])]/4.6*0.4068584935483871*16*0.01*0.43 for ind, x in enumerate(common_index)] # Calculating PP using CHL Hermes and Flux

                pp_behrenfeld_chl4 = [10 ** (0.559 * math.log10(chl_to_int[ind]) + 2.793) if not np.isnan(CHL_int[ind]) and CHL_int[ind] != 0 else 0 for ind, x in enumerate(common_index)]
                pp_behrenfeld_chl3 = [10 ** (0.559 * math.log10(chl3_to_int[ind]) + 2.793) if CHL3_int[ind] != 0 else 0 for ind, x in enumerate(common_index)]
                date_pp = pp_time[pp_time[0] == str(year_n[n]) + '-0' + str(month_n[m]) + '-01'].index[0] # Storing the date for Copernicus PP Model
                pp_df = pp[date_pp].to_dataframe().pivot_table(index='lat', columns='lon', values='PP', dropna=False) # Storing and converting Copernicus PP Model data
                pp_df.index = lat_sat # Latitude calibration
                pp_df.columns = lon_sat # Longitude calibration

                date_pp_wider = pp_time_wider[pp_time_wider[0] == str(year_n[n]) + '-0' + str(month_n[m]) + '-01'].index[0]  # Storing the date for Copernicus PP Model
                pp_df_wider = pp_wider[date_pp_wider].to_dataframe().pivot_table(index='lat', columns='lon', values='PP', dropna=False)  # Storing and converting Copernicus PP Model data
                pp_df_wider.index = lat_wider  # Latitude calibration
                pp_df_wider.columns = lon_wider  # Longitude calibration

                pp_mod = [pp_df.loc[lat_to_check_normalized[x]][lon_to_check_normalized[x]] if lon_to_check_normalized[x] <=15 and lon_to_check_normalized[x] >=-18 else pp_df_wider.loc[nearest(lat_wider, lat_to_check_normalized[x])][nearest(lon_wider, lon_to_check_normalized[x])] for x in common_index] # Searching for Copernicus PP Model values

                chl3_out = [x for x in chl3_to_int]

                if par_type == 2 and chl_level == 4: # Output for CHL Copernicus and PAR Eumetsat
                    if not spectra:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'PAR Eumetsat CHL4' : pp_chl4_par2}), ignore_index=True)
                    if spectra == 1:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'PUR Eumetsat CHL4':pp_chl4_par2_multi_by_aph}), ignore_index=True)
                    if spectra == 2 :
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Spectra Eumetsat CHL4':pp_chl4_par2_spectra_mean}), ignore_index=True)


                if par_type == 3 and chl_level == 4: # Output for CHL Copernicus and PAR Hermes
                    if not spectra:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'PAR Hermes CHL4' : pp_chl4_par3}), ignore_index=True)
                    if spectra == 1:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'PUR Hermes CHL4' : pp_chl4_par3_multi_by_aph}), ignore_index=True)
                    if spectra == 2:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Spectra Hermes CHL4' : pp_chl4_par3_spectra_mean}), ignore_index=True)


                if par_type == 1 and chl_level == 4: # Output for CHL Copernicus and Flux
                    if not spectra:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Solar R. Flux CHL4' : pp_chl4_f}), ignore_index=True)
                    if spectra == 1 :
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Solar R. PUR CHL4' : pp_chl4_f_multi_by_aph}), ignore_index=True)
                    if spectra == 2:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Solar R. Spectra CHL4' : pp_chl4_f_spectra_mean}), ignore_index=True)


                if par_type == 2 and chl_level == 3: # Output for CHL Hermes and PAR Eumetsat
                    if not spectra:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'PAR Eumetsat CHL3' : pp_chl3_par2}), ignore_index=True)
                    if spectra == 1:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'PUR Eumetsat CHL3' : pp_chl3_par2_multi_by_aph}), ignore_index=True)
                    if spectra == 2:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Spectra Eumetsat CHL3' : pp_chl3_par2_spectra_mean}), ignore_index=True)


                if par_type == 3 and chl_level == 3: # Output for CHL Hermes and PAR Hermes
                    if not spectra:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'PAR Hermes CHL3' : pp_chl3_par3, 'CHL3': chl3_out}), ignore_index=True)  # Chl3
                    if spectra == 1:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'PUR Hermes CHL3' : pp_chl3_par3_multi_by_aph}), ignore_index=True)
                    if spectra == 2:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Spectra Hermes CHL3' : pp_chl3_par3_spectra_mean}), ignore_index=True)


                if par_type == 1 and chl_level == 3: # Output for CHL Hermes and FLux
                    if not spectra:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Solar R. Flux CHL3' : pp_chl3_f}), ignore_index=True)
                    if spectra == 1:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Solar R. PUR CHL3' : pp_chl3_f_multi_by_aph}), ignore_index=True)
                    if spectra == 2:
                        df = df.append(pd.DataFrame( {'PP Field' : [PP_int_from_paper[x] for x in common_index], 'Solar R. Spectra CHL3' : pp_chl3_f_spectra_mean}), ignore_index=True)


                if par_type == 0 and chl_level == 0: # Output for Copernicus Model
                    df = df.append(pd.DataFrame({'PP Field': [PP_int_from_paper[x] for x in common_index], 'Copernicus Model': pp_mod}), ignore_index=True)

                if par_type == -1 and chl_level == 4: # Output for Copernicus Model
                    df = df.append(pd.DataFrame({'PP Field': [PP_int_from_paper[x] for x in common_index], 'PP Behrenfeld 4': pp_behrenfeld_chl4}), ignore_index=True)

                if par_type == -1 and chl_level == 3: # Output for Copernicus Model
                    df = df.append(pd.DataFrame({'PP Field': [PP_int_from_paper[x] for x in common_index], 'PP Behrenfeld 3': pp_behrenfeld_chl3}), ignore_index=True)

    df = df[df.iloc[:,1]!=0]
    df.dropna(inplace=True)
    df.reset_index(inplace=True)
    del df['index']

    rmsd = [(math.log10(df.iloc[x,1]) - math.log10(df['PP Field'][x]))**2 for x in range(df.shape[0]) if not np.isnan(df.iloc[x,1])]
    rmsd = np.sqrt(sum(rmsd)/df.shape[0])
    df['RMSD'] = rmsd

    pp_field_log = [math.log10(df['PP Field'][x]) for x in range(df.shape[0]) if not np.isnan(df.iloc[x,1])]
    pp_mod_log = [math.log10(df.iloc[x,1]) for x in range(df.shape[0]) if not np.isnan(df.iloc[x,1])]
    bias = np.mean(pp_mod_log) - np.mean(pp_field_log)
    bias_norm =bias/statistics.stdev(pp_field_log)
    df['BIAS'] = bias
    df['BIAS norm'] = bias_norm

    df['uRMSD'] = np.sqrt(sum([((math.log10(df.iloc[x,1]) - np.mean(pp_mod_log)) - (math.log10(df['PP Field'][x]) - np.mean(pp_field_log)))**2 for x in range(df.shape[0]) if not np.isnan(df.iloc[x,1])])/df.shape[0])
    df['uRMSD'] = df['uRMSD'] * -1 if (statistics.stdev(pp_mod_log) - statistics.stdev(pp_field_log)) < 0 else df['uRMSD']
    df['uRMSD norm'] = df['uRMSD']/statistics.stdev(pp_field_log)

    df['Corr coef'] = np.corrcoef(df.iloc[:, 0], df.iloc[:, 1])[0][1]

    return df