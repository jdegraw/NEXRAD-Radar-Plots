from radar_plots_with_gps import *

basefolder = r"C:\Users\jonat\MetEnv\code\NEXRAD-Radar-Plots"
custom_animation(zoom_step = 0.2, station = 'KLOT', gps_data=basefolder+r"\GPS_data\GPS_data_combined_Sort_Merge.shp", 
                 tornado_tracks = {basefolder+r'\Tornadoes\south_kirkland.shp':{'start':datetime(2021,8,9,21,44),'end':datetime(2021,8,9,21,47)},
                 basefolder+r"\Tornadoes\burlington.shp":{'start':datetime(2021,8,9,22,24),'end':datetime(2021,8,9,22,39)},
                 basefolder+r"\Tornadoes\sycamore_1.shp":{'start':datetime(2021,8,9,23,13),'end':datetime(2021,8,9,23,18)},
                 basefolder+r"\Tornadoes\sycamore_2.shp":{'start':datetime(2021,8,9,23,26),'end':datetime(2021,8,9,23,27)},
                 basefolder+r"\Tornadoes\sycamore_3.shp":{'start':datetime(2021,8,9,23,28),'end':datetime(2021,8,9,23,37)}})