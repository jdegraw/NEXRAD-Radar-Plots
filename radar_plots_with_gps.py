'''
Import the nessesary Packages/ Libraries
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta # read/ display the time
import cartopy.crs as ccrs # Map projections
import cartopy.feature as cfeature # Map features
import cartopy.io.shapereader as shpreader # Helps read shapefiles; can be used to add additional features to the map
from metpy.plots import ctables, USCOUNTIES, add_timestamp # Color tables, US County lines, and function to add time to plots
from metpy.io import Level2File # Function that reads level 2 NEXRAD Data
from pyproj import Geod # Used to change cartesian/ polar coordinates to latitude and longitude
import boto3 # Accessing radar data
import botocore
from botocore.client import Config
import geopandas
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import shapely.geometry as sgeom
from shapely.ops import substring


'''
The following function accesses level 2 nexrad data files. For more information
search docs.opendata.aws/noaa-nexrad/readme.html

The function takes two arguments, the date and the station. The date MUST
be a datetime object. The station must be a four character string, 
the station ID.
'''

def access_datafiles(d, station):
    s3 = boto3.resource('s3', config=Config(signature_version=botocore.UNSIGNED, 
                                            user_agent_extra='Resource'))
    bucket = s3.Bucket('noaa-nexrad-level2')
    prefix = f'{d:%Y}/{d:%m}/{d:%d}/{station}/{station}{d:%Y%m%d_%H}'
    #print(prefix)
    objects = []
    for obj in bucket.objects.filter(Prefix = prefix):
        print(obj.key)
        objects.append(obj)
    return objects
    
'''
Function that will be called to plot the radar data. Accepts the output of Level2File() 
and a matplotlib.axes.Axes object as arguements.
'''

def plot_radar(f):
    sweep = 0

    az = np.array([ray[0].az_angle for ray in f.sweeps[sweep]])

    ref_hdr = f.sweeps[sweep][0][4][b'REF'][0]
    ref_range = np.arange(ref_hdr.num_gates) * ref_hdr.gate_width + ref_hdr.first_gate

    ref = np.array([ray[4][b'REF'][1] for ray in f.sweeps[sweep]])
    #radarfilter = f.clutter_filter_bypass_map['data'][0]

    truedata = []
    for data_ray in ref:
        rows = []
        for data_point in data_ray:
            if data_point < 10:
                data_point = np.nan
                rows.append(data_point)
            else:
                rows.append(data_point)
        truedata.append(rows)

    truedata = np.array(truedata)
    data = np.ma.array(truedata)
    data[np.isnan(data)] = np.ma.masked
    
    lat = f.sweeps[sweep][0][1].lat
    lon = f.sweeps[sweep][0][1].lon
    
    g = Geod(ellps='clrk66')
    
    center_lat = np.ones([len(az),len(ref_range)])*lat    
    center_lon = np.ones([len(az),len(ref_range)])*lon
    az2D = np.ones_like(center_lat)*az[:,None]
    rng2D = np.ones_like(center_lat)*np.transpose(ref_range[:,None])*1000
    lon,lat,back=g.fwd(center_lon,center_lat,az2D,rng2D)
    del center_lon
    del center_lat
    del az2D
    del rng2D

    norm, cmap = ctables.registry.get_with_steps('NWSStormClearReflectivity',-20,0.5)
    
    return lon, lat, data, cmap, norm
    
'''
Function that plots city names. Accepts an extent (a list of latitudes and longitudes
that define the edges of the plot) and a matplotlib.axes.Axes object as arguments.
'''

def plot_cities(extent, ax):
    pop = cfeature.NaturalEarthFeature(
            category = 'cultural',
            name = 'populated_places',
            scale = '10m',
            facecolor = 'none')
    pop = shpreader.natural_earth(resolution='10m', category='cultural', name='populated_places')
    shp = shpreader.Reader(pop)

    lats = []
    lons = []
    USCities = []
    populations = []
    for record in shp.records():
        country = record.attributes['SOV0NAME']
        lat = record.attributes['LATITUDE']
        lon = record.attributes['LONGITUDE']
        population = record.attributes['POP_MAX']
        if extent[0]<lon<extent[1] and extent[2]<lat<extent[3]:
            USCities.append(record.attributes['NAME'])
            lats.append(lat)
            lons.append(lon)
            populations.append(population)
    if len(populations) > 100:
        Sorted = sorted(populations, reverse = True)[0:29]
        maximum = max(Sorted)
        minimum = min(Sorted)
    elif len(populations) > 0:
        maximum = max(populations)
        minimum = min(populations)
    else:
        return

    transform = ccrs.PlateCarree()._as_mpl_transform(ax)
    for i in range(len(USCities)):
        if minimum <= populations[i] <= maximum:
            ax.annotate(USCities[i], (lons[i],lats[i]), 
                        c = 'lightyellow', xycoords = transform, size = 24, 
                        fontfamily = 'Calibri')
                        
def create_bucketlist(datelist, station, frequency):
    bucketlist = []
    for d in datelist:
#         try:
        subbucketlist = access_datafiles(d, station)
        if frequency == 'hourly':
            if subbucketlist[0] not in bucketlist:
                bucketlist.append(subbucketlist[0])
        elif frequency == 'all_scans':
            for subbucket in subbucketlist:
                if subbucket not in bucketlist:
                    bucketlist.append(subbucket)
                
#         except:
#             print('Failed!')
#             continue
    del subbucketlist
    return bucketlist
    
    
def create_figure():
    
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize = (20,20))
    ax = plt.subplot(projection = proj)
    ax.set_facecolor('black')

    ax.add_feature(USCOUNTIES.with_scale('20m'), edgecolor = 'white',
                   linewidth = 2, zorder = 0)
    roads = cfeature.NaturalEarthFeature(
            category = 'cultural',
            name = 'roads',
            scale = '10m',
            facecolor = 'none'
            )
    ax.add_feature(roads, edgecolor = 'lightblue', linewidth = 2)
    
    return fig, ax
    
def custom_animation(zoom_step = 0.5, station = 'KIND', **kwargs):
    
    from matplotlib.animation import ArtistAnimation #import function to do the animation
    
#     proj = ccrs.PlateCarree()
#     fig = plt.figure(figsize = (20,20))
#     ax = plt.subplot(projection = proj)
#     ax.set_facecolor('black')

#     ax.add_feature(USCOUNTIES.with_scale('20m'), edgecolor = 'white',
#                    linewidth = 0.5, zorder = 0)
#     roads = cfeature.NaturalEarthFeature(
#             category = 'cultural',
#             name = 'roads',
#             scale = '10m',
#             facecolor = 'none'
#             )
#     ax.add_feature(roads, edgecolor = 'lightblue', linewidth = 0.5)
    
    if 'lat' in kwargs and 'lon' in kwargs:
        if 'time' not in kwargs:
            raise ValueError("Invalid Kwargs. Missing 'time' kwarg.")
        lat = kwargs['lat']
        lon = kwargs['lon']
        extent = [lat-zoom_step,lat+zoom_step,lon-zoom_step,lon+zoom_step]
        ax.set_extent(extent)
        plot_cities(extent, ax)
        time = kwargs['time']
        
        bucketlist = create_bucketlist(time, station)
        
        artists = []
        for file in bucketlist[0:len(bucketlist)-1]:
            f = Level2File(file.get()['Body'])
            try:
                mesh, time_stamp = plot_radar(f, ax)
                artists.append([mesh,time_stamp])
            except:
                continue
        
    elif 'gps_data' in kwargs:
        gps_data = geopandas.read_file(kwargs['gps_data'])
        
        gps_data['real_time'] = [datetime.strptime(point,'%Y-%m-%dT%H:%M:%SZ').replace(month = 8)+timedelta(hours=5) for point in gps_data['DateTimeS']]
        gps_data['real_geometry'] = [sgeom.Point(point.coords.xy[0][0],point.coords.xy[1][0]) for point in gps_data['geometry']]
        hours = set([date.replace(minute=0, second=0, microsecond=0) for date in gps_data['real_time']])
        day = gps_data['DateTime'][0]
        
        tornado = False
        if 'tornado_tracks' in kwargs:
            tor_list = []
            for tornado_shpfile in kwargs['tornado_tracks']:
                
                tor_data = geopandas.read_file(tornado_shpfile)
                dt = (kwargs['tornado_tracks'][tornado_shpfile]['end']-kwargs['tornado_tracks'][tornado_shpfile]['start'])
               
                mp = sgeom.MultiPoint()
                line = tor_data['geometry'][0]

                diff = line.length/dt.total_seconds()
                for i in np.arange(0, line.length, diff):
                    s = substring(line, i, i+diff)
                    mp = mp.union(s.boundary)
                pointlist = []
                timelist = []
                time_counter = kwargs['tornado_tracks'][tornado_shpfile]['start']
                for point in mp:
                    pointlist.append(point)
                    timelist.append(time_counter)
                    time_counter = time_counter + timedelta(seconds = 1)
                
                tor_df = pd.DataFrame({'geometry':pointlist,'real_time':timelist})
                tor_gdf = geopandas.GeoDataFrame(tor_df)
                tor_list.append(tor_df)
                print("Read Tornado file: " + tornado_shpfile)
            del tor_data
            del dt
            del line
            del diff
            del s

            del tor_df
            del pointlist
            del timelist
            del time_counter
            
                    
                    
            tornado = True
            print("----All Tornado Files Read Successfully----")
        bucketlist = create_bucketlist(list(hours), station, 'all_scans')
        first_frame = True
        nframe = 0
        for date in gps_data['real_time']:
            
            print(f'Time of Frame: {date}')
            
            i = gps_data['geometry'].index[gps_data['real_time']==date].tolist()[0]
                    
            lat = gps_data['geometry'][i].coords.xy[0][0]
            lon = gps_data['geometry'][i].coords.xy[1][0]


            for bucket in bucketlist:

                if bucket.key[20:35] == date.strftime("%Y%m%d_%H%M%S"):
                    f = Level2File(bucket.get()['Body'])
                    first_frame = False
                    radar_lon, radar_lat, data, cmap, norm = plot_radar(f)
                    del f
                    bucketlist.remove(bucket)
                    break

            if first_frame:
                continue
            fig, ax = create_figure()
            ax.set_extent([lat-zoom_step,lat+zoom_step,lon-zoom_step,lon+zoom_step])
            
            mesh = ax.pcolormesh(radar_lon, radar_lat, data, cmap=cmap, norm=norm, shading = 'auto')
            if i <= 750:
                j = 0
            else:
                j = i-750
            geopandas.GeoSeries(gps_data['real_geometry'][j:i]).plot(ax=ax, marker='o', color='blue', markersize=40, transform = ccrs.PlateCarree())
            image = plt.imread("location_indicator.png")
            imagebox = OffsetImage(image, zoom=.12)
            imagebox.image.axes = ax

            ab = AnnotationBbox(imagebox, (lat, lon), xycoords='data', frameon=False)
            ax.add_artist(ab)

            if tornado:
                for tor_df in tor_list:
                    for t in tor_df['real_time']:
                        if t == date:
                            pt = tor_df.loc[tor_df['real_time'] == date,'geometry'].item()
                            x = pt.coords.xy[0][0]
                            y = pt.coords.xy[1][0]
                            image = plt.imread("tornado_icon.png")
                            imagebox = OffsetImage(image, zoom=.07)
                            imagebox.image.axes = ax
                            ab = AnnotationBbox(imagebox, (x, y), xycoords='data', frameon=False, box_alignment = (0.6,0))
                            ax.add_artist(ab)
                            
                            geopandas.GeoSeries(tor_df['geometry']).plot(ax=ax, marker='o', color='gray', markersize=40, transform = ccrs.PlateCarree())
                            geopandas.GeoSeries(tor_df['geometry'][0:tor_df['geometry'][tor_df['real_time'] == date].index.tolist()[0]]).plot(ax=ax, marker='o', color='purple', markersize=40, transform = ccrs.PlateCarree())      
            
            
            time_stamp = add_timestamp(ax, time=date)
            time_stamp.set_font('Calibri')
            time_stamp.set_fontsize(36)
            fname = r"C:\Users\Jonathan DeGraw\NEXRAD Radar Code\Frames" + '\\' + str(nframe) + "_Animation_Frame_.png"
            plt.savefig(fname)
            print(f"Saved File: {fname}")
            plt.close(fig)
            nframe += 1
            
            del fig
            del ax
            del mesh
            del lat
            del lon
        
    
            
    else:
        raise ValueError("Invalid Kwargs. Please use either the 'lat' and 'lon' kwargs or enter the GPS data as a shpfile under the kwarg 'gps_data'.")
        