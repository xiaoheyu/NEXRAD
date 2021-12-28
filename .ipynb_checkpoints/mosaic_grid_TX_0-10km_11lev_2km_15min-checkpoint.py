import sys
import numpy as np
import datetime
import pyart
import boto
import os
import tempfile
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import time
#suppress deprecation warnings
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)
import moviepy.editor as mpy
#import dill
import pickle
from pathlib import Path
import glob


###############################################################################################################################

def plot_grid(grid,plotVariable = 'reflectivity',level = 13):
               
    fig = plt.figure(figsize=[15, 8])
    ax = fig.add_subplot(111)
    displaygrid = pyart.graph.GridMapDisplayBasemap(grid)
    displaygrid.plot_basemap(lat_lines=(25,30,35,40,45,50), lon_lines=(-120,-110,-100,-90,-80,-70))
    displaygrid.plot_grid(plotVariable, level, vmin=-20, vmax=40,
                 cmap = pyart.graph.cm.NWSRef)
    fig.savefig(os.path.join("./plot/", grid.time["units"] +str(level) + ".png"), format='png', dpi=200)
    print("Map saved")
    plt.close()
#     grid.write('Mosaic_grid05_US', format='NETCDF4', arm_time_variables=False, arm_alt_lat_lon_variables=False)
#     print("Grid saved")
    
#     pyart.io.write_grid_geotiff(grid, "Mosaic_grid_US", 'reflectivity', rgb=False, level=None, cmap='viridis', vmin=0, vmax=75, color_levels=None, warp=False, sld=False, use_doublequotes=False)
#     print("Geotiff saved")
###############################################################################################################################
def Diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1]
    return li_dif

#get all sites' name
def get_allsites():
    all_sites = []
    locs = pyart.io.nexrad_common.NEXRAD_LOCATIONS
    for key in locs:
        all_sites.append(key)
    return all_sites
###############################################################################################################################

def customized_download(date='2020/05/16',time=190000,all_sites=['KFWS'],download = "NO"):

    #create a datetime object for the current time in UTC and use the
    # year, month, and day to drill down into the NEXRAD directory structure.
    #now = datetime.datetime.utcnow()
    #pirnt(now)
    #date = ("{:4d}".format(now.year) + '/' + "{:02d}".format(now.month) + '/' +
    #        "{:2d}".format(now.day) + '/')

    #----------------------------------------------------------
#     date = '2020/01/28/'
#     time = 220000
#     #all_sites = get_allsites()

#     # assign a target site for testing purpose
#     all_sites = ['KFWS']
    #---------------------------------------------------------

    print("Search Time %s %s" %(date, time))
    print(all_sites)

    radars = []
    i = 0
    count = 1
    awscount = 0
    Fail_lst =[]
    awslst = []
    #get the bucket list for the selected date

    #use boto to connect to the AWS nexrad holdings directory
    s3conn = boto.connect_s3()
    bucket = s3conn.get_bucket('noaa-nexrad-level2')

    #Note: this returns a list of all of the radar sites with data for
    # the selected date
    ls = bucket.list(prefix=date + '/',delimiter='/')
    for item in ls:
        awslst.append(item.name.split('/')[-2])

    #Find the Missing sites from AWS lst
    li3 = Diff(awslst, all_sites)
    print("Missing sites : %s " %(li3))

    for key in ls:
        #print(key.name)
        awscount+=1
    print("%d sites are selected, total %d sites returned from AWS at %s %s" %(len(all_sites)-len(li3),awscount-1,date,time))
    print("===================================================================================")
    print('\n')

    for site in all_sites:
        for key in ls:
            #only pull the data and save the arrays for the site we want
            if site in key.name.split('/')[-2]:
                print("%s has been found in AWS return list (%d/%d)" %(site,count,awscount-1))
                #set up the path to the NEXRAD files
                path = date +'/' + site + '/' + site
                #grab the last file in the file list
                keys = bucket.get_all_keys(prefix=path)

                #find the file closest to myTime
                temp = 0
                minimum = 20000
                for key in keys:
                    radarTime = int(key.name.split('/')[-1].split('_')[1])
                    temp = abs(int(time) - radarTime)
                    #print(temp)
                    if temp < minimum:
                        minimum = temp
                        fname = key
                try:
                    print(fname)

                    #get the file
                    s3key = bucket.get_key(fname)
                    #save a temporary file to the local host
                    localfile = tempfile.NamedTemporaryFile()
                    #write the contents of the NEXRAD file to the temporary file
                    s3key.get_contents_to_filename(localfile.name)
                    print(localfile.name)
                    #use the read_nexrad_archive function from PyART to read in NEXRAD file

                    radars.append(pyart.io.read_nexrad_archive(localfile.name))
                    #get the date and time from the radar file for plot enhancement
                    fileTime = radars[i].time['units'].split(' ')[-1].split('T')
                    print(site + ' ' + 'has been read and added to radars list' + ': ' + fileTime[0] + ' at ' + fileTime[1] )
                    i+=1
                    count+=1
                    print('\n')
#                     localfile.close()
#                     os.remove(localfile.name)
                except:
                    print("%s not read sucsseful <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" %(site))
                    Fail_lst.append(site)
                    print('\n')
    print("Failed reading sites %s" %(Fail_lst))
    
     # check all radar objects are valid
    print("%d radars objects are saved in radar list " %len(radars))
    field_names = ['reflectivity','velocity','spectrum_width',
            'differential_phase','differential_reflectivity','cross_correlation_ratio']
    remove_list = []
    
    for i in range(len(radars)):
        try:
            for field_name in field_names:
                radars[i].check_field_exists(field_name)
        except:
            remove_list.append(i)
            

    for index in sorted(remove_list, reverse=True):
        del radars[index]
        
    print("%d radar obejcts are valid with all six variables " % len(radars))
    
    
    Path(os.path.join(date.replace("/","-")[:7],date.replace("/","-"),time)).mkdir(parents = True, exist_ok=True)
    if download == "YES":
        print("Saving radars object ......")
        
        pickle.dump(radars, open(os.path.join(date.replace("/","-"),time,date.replace("/","-") + "_" + time +"_" + 'radars.pickle'), 'wb'))
    return radars
    
    
###############################################################################################################################
    
    
def mosaic_radar(radars,
#                grid_shape =(11, 140, 135),
#                grid_limits=((0, 10001), (-800000, 600000),(-850000, 500000)),
               
               grid_shape=(4, 200, 250), # 2km resolution
               grid_limits=((0, 3001), (-200000, 200000), (-250000, 250000)),
               grid_origin = (32.57861, -97.303611),
               date= 0,
               time = 0
                
#                grid_shape =(11, 110, 100),
#                grid_limits=((0, 10001),(-750000, 350000),(-400000, 600000)),
#                grid_origin = (39.00, -120.00),
#                date= 0,
#                time = 0
                ):
    
    
    
    
    print("Converting radar to grid......")
    grid = pyart.map.grid_from_radars(
    radars,
    grid_shape=grid_shape,
    grid_limits=grid_limits,
    grid_origin = grid_origin,
    gridding_algo='map_gates_to_grid',
    fields=['reflectivity','velocity','spectrum_width',
            'differential_phase','differential_reflectivity','cross_correlation_ratio'])
    
    print("Mosaic grid created")
    
#     print("Plottting radar map......")
#     fig = plt.figure(figsize=[15, 8])
#     ax = fig.add_subplot(111)
    
#     displaygrid = pyart.graph.GridMapDisplayBasemap(grid)
#     displaygrid.plot_basemap(lat_lines=(25,30,35,40,45,50), lon_lines=(-120,-110,-100,-90,-80,-70))
#     displaygrid.plot_grid(plotVariable, level=0, vmin=-20, vmax=40,
#                  cmap = pyart.graph.cm.NWSRef)
    
    
    
    
    
#     fig.savefig(os.path.join("./plot/"+str(date), grid.time["units"] + ".png"), format='png', dpi=200)
#     print("Map saved")
    grid.write(os.path.join(date.replace("/","-")[:7],date.replace("/","-"), time,date.replace("/","-") + "_" + time + ".nc"), format='NETCDF4', arm_time_variables=False, arm_alt_lat_lon_variables=False)
    print("Grid saved")
    
#     pyart.io.write_grid_geotiff(grid, os.path.join(date.replace("/","-")[:7],date.replace("/","-"), time, date.replace("/","-") +"_" + time), 'reflectivity', rgb=False, level=None, cmap='viridis', vmin=0, vmax=75, color_levels=None, warp=False, sld=False, use_doublequotes=False)
#     print("Geotiff saved")
    
    return grid
    
    
    
###############################################################################################################################
    
def plot_grid_elevation(grid,date,time,plotVariable = 'reflectivity'):
        # Plot radars at each elevation and create the gif movie

    #images = []
    for i in range(4):
        fig = plt.figure(figsize=[15, 8])
        #ax = fig.add_subplot(111)

        displaygrid = pyart.graph.GridMapDisplayBasemap(grid)
        displaygrid.plot_basemap(lat_lines=(26,28,30,32,34,36,38), lon_lines=(-106,-104,-102,-100,-98,-96,-94), auto_range=True, min_lon=-92, max_lon=-86, min_lat=25, max_lat=44) # extent for dallas
        
#         displaygrid.plot_basemap(lat_lines=(25,30,35,40,45,50), lon_lines=(-120,-110,-100,-90,-80,-70), auto_range=True, min_lon=-92, max_lon=-86, min_lat=25, max_lat=44) # extent for US
            
        displaygrid.plot_grid(plotVariable, level=i, vmin=-20, vmax=40,
                         cmap = pyart.graph.cm.NWSRef
        ,title = "TEXAS Mosaic at %2.1f km %s %s %s"%(i, date, time, plotVariable))
        
        plt.savefig(os.path.join(date.replace("/","-")[:7],date.replace("/","-"),time, date.replace("/","-")+ "_" + time +  "_" + "level" + "_" + str(f"{i:02d}")+ ".png") )
        #images.append(plt)
        print("level %d image saved" %(i))
        plt.close()
        
        
        
#     for i in range(20,101,10):
#         fig = plt.figure(figsize=[15, 8])
#         #ax = fig.add_subplot(111)

#         displaygrid = pyart.graph.GridMapDisplayBasemap(grid)
#         displaygrid.plot_basemap(lat_lines=(25,30,35,40,45,50), lon_lines=(-120,-110,-100,-90,-80,-70), auto_range=True, min_lon=-92, max_lon=-86, min_lat=25, max_lat=44) # extent for dallas
        
# #         displaygrid.plot_basemap(lat_lines=(25,30,35,40,45,50), lon_lines=(-120,-110,-100,-90,-80,-70), auto_range=True, min_lon=-92, max_lon=-86, min_lat=25, max_lat=44) # extent for US
            
#         displaygrid.plot_grid(plotVariable, level=i, vmin=-20, vmax=40,
#                          cmap = pyart.graph.cm.NWSRef)
#         #,title = "U.S Mosaic at %2.1f km %s %s %s"%(i*0.1, date, time, plotVariable) 
        
#         plt.savefig(os.path.join(date.replace("/","-"),time, date.replace("/","-")+ "_" + time +  "_" + "level" + "_" + str(f"{i:02d}")+ ".png") )
#         #images.append(plt)
#         print("level %d image saved" %(i))


###############################################################################################################################

def make_gif(date,time,fps=5, gif_name = 'KFWS_matchPlot'):
    folder=os.path.join(date.replace("/","-")[:7],date.replace("/","-"),time)
    file_list = glob.glob(folder+'/*.png')
    print(file_list)
    list.sort(file_list, key=lambda x: int(x.split('_')[3].split('.png')[0])) # Sort the images by #, this may need to be tweaked for your use case
    clip = mpy.ImageSequenceClip(file_list, fps=fps)
    clip.write_gif(os.path.join(date.replace("/","-")[:7],date.replace("/","-"),time,'{}.gif'.format(gif_name)), fps=fps)

   

########################################################################################################
def daterange(date1, date2):
    for n in range(int ((date2 - date1).days)+1):
        yield date1 + datetime.timedelta(n)
        
#########################################################################################################################
def implement_mosaic(by,bm,bd,ey,em,ed):
    start_dt = datetime.date(by,bm,bd)  # custommized your starting date and end date here
    end_dt = datetime.date(ey,em,ed)
    TimeRange=['000000','001500','003000','004500','010000','011500','013000','014500','020000','090000','091500','093000','094500','100000','101500','103000','104500','110000','111500','113000','114500','120000','121500','123000','124500','130000','131500','133000','134500','140000','141500','143000','144500','150000','151500','153000','154500','160000','161500','163000','164500','170000','171500','173000','174500','180000','181500','183000','184500','190000','191500','193000','194500','200000','201500','203000','204500','210000','211500','213000','214500','220000','221500','223000','224500','230000','231500','233000','234500']


    for dt in daterange(start_dt, end_dt):
        for myTime in TimeRange:
            start_time = time.time()
            #check if target mosaic exists, if yes then plot directly
            if Path(os.path.join(dt.strftime("%Y-%m"),dt.strftime("%Y-%m-%d"), myTime,dt.strftime("%Y-%m-%d") + "_" + myTime + ".gif")).exists():
                print("%s exists" %(os.path.join(dt.strftime("%Y-%m"),dt.strftime("%Y-%m-%d"), myTime,dt.strftime("%Y-%m-%d") + "_" + myTime + ".gif")))


            elif Path(os.path.join(dt.strftime("%Y-%m"),dt.strftime("%Y-%m-%d"), myTime,dt.strftime("%Y-%m-%d") + "_" + myTime + ".nc")).exists():

                print("%s exists, start plotting now......" %(os.path.join(dt.strftime("%Y-%m"),dt.strftime("%Y-%m-%d"), myTime,dt.strftime("%Y-%m-%d") + "_" + myTime + ".nc")))

#                 try:
#                     grid = pyart.io.read_grid(os.path.join(dt.strftime("%Y-%m"),dt.strftime("%Y-%m-%d"), myTime,dt.strftime("%Y-%m-%d") + "_" + myTime + ".nc"))
#                     plot_grid_elevation(grid,dt.strftime("%Y/%m/%d"),myTime)

#                     make_gif(dt.strftime("%Y/%m/%d"), myTime, gif_name = dt.strftime("%Y-%m-%d") +"_"+ myTime)
#                     print("--- %s seconds ---" % (time.time() - start_time))
#                 except:
#                     print("Reading grid or plotting failed")



            else:
                radars =customized_download(date=dt.strftime("%Y/%m/%d"),time=myTime,all_sites=all_sites,download = 'NO')


                try:
                    grid = mosaic_radar(radars,date=dt.strftime("%Y/%m/%d"),time = myTime)
#                     plot_grid_elevation(grid,dt.strftime("%Y/%m/%d"),myTime)

#                     make_gif(dt.strftime("%Y/%m/%d"), myTime, gif_name = dt.strftime("%Y-%m-%d") +"_"+ myTime)
                    print("--- %s seconds ---" % (time.time() - start_time))
                    del grid, radars
                    gc.collect()
                except:
                    print("convert failed")

#####################################################################################################
if __name__ == "__main__":
    
#     now = datetime.datetime.utcnow()
#     myDate = ("{:4d}".format(now.year) + '/' + "{:02d}".format(now.month) + '/' +"{:02d}".format(now.day))
#     myTime = ("{:02d}".format(now.hour)  + "{:02d}".format(now.minute)  +"{:02d}".format(now.second) )
#     outPut_file = os.path.join(myDate,myTime,"nohup.txt")

    # TEXAS
    #all_sites=['KAMA','KLBB','KMAF','KSJT','KDFX','KFDR','KDYX','KFWS','KGRK','KEWX','KCRP','KBRO','KSHV','KHGX','KROE','KLCH']
#     all_sites=['KMAX','KBHX','KBBX','KDAX','KMVX','KHNX','KVBX','KVTX','KEYX','KSOX','KNKX','KYUX']
# specify the nexrad site name, staring and ending date range
    all_sites=['KFWS']
    by = int(sys.argv[1])
    bm = int(sys.argv[2])
    bd = int(sys.argv[3])
    ey = int(sys.argv[4])
    em = int(sys.argv[5])
    ed = int(sys.argv[6])

    implement_mosaic(by,bm,bd,ey,em,ed)

    