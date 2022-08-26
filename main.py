from radar_plots_with_gps import *


custom_animation(zoom_step = 0.2, station = 'KLOT', gps_data=r"C:\Users\Jonathan DeGraw\NEXRAD Radar Code\GPS_data_combined_sorted\GPS_data_combined_Sort_Merge.shp", 
                 tornado_tracks = {r'C:\Users\Jonathan DeGraw\NEXRAD Radar Code\Tornadoes\south_kirkland.shp':{'start':datetime(2021,8,9,21,44),'end':datetime(2021,8,9,21,47)},
                 r"C:\Users\Jonathan DeGraw\NEXRAD Radar Code\Tornadoes\burlington.shp":{'start':datetime(2021,8,9,22,24),'end':datetime(2021,8,9,22,39)},
                 r"C:\Users\Jonathan DeGraw\NEXRAD Radar Code\Tornadoes\sycamore_1.shp":{'start':datetime(2021,8,9,23,13),'end':datetime(2021,8,9,23,18)},
                 r"C:\Users\Jonathan DeGraw\NEXRAD Radar Code\Tornadoes\sycamore_2.shp":{'start':datetime(2021,8,9,23,26),'end':datetime(2021,8,9,23,27)},
                 r"C:\Users\Jonathan DeGraw\NEXRAD Radar Code\Tornadoes\sycamore_3.shp":{'start':datetime(2021,8,9,23,28),'end':datetime(2021,8,9,23,37)}})