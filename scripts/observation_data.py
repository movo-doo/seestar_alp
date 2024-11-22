"""
Working Test script for implementation extracts fit header info and moon position and other attributes from date/time of observation.
Files data in csv file.  Cannot be run effectively without setting up environment (files and directories) for it to work with. See code below
This is a few snippets from a larger python program provides pipeline workflow from seestar files to anaylsing and stacking and backing up files.
Not neccessarily part of seestar alp however code to gather data from fit header and moon ephemeris data may be useful to project.
I am a novice python programmer at best so use best programming practices have not been implemented. As well coding performed in such a way
to allow myself to follow my own workflow.  See cls_target.py accompanying code in seperate file
"""

import sys
import os
import shutil
import glob
import ephem
from pathlib import Path
from datetime import datetime, timezone
import csv
import math
import astropy as ap
from astropy.io import fits
from cls_target import Target                     # LOCAL PY FILE STORED IN SAME DIRECTORY AS THIS PY PROGRAM


def setup_observer(location_lat, location_lon, location_elevation, date_time):
    observer = ephem.Observer()
    observer.lat = str(location_lat)
    observer.long = str(location_lon)
    observer.elevation = location_elevation
    observer.date = date_time

    return observer


def get_moon_illumination(observer):

    # Get the moon object
    moon = ephem.Moon(observer)

    # Calculate moon phase percentage
    moon_phase = int(moon.phase)

    return moon_phase


def calculate_angular_separation(observer, target):

    # Get the Moon's position at Observer.lat,long,height,datetime from Image Header
    moon = ephem.Moon(observer)
    # Convert Moon's RA and Dec to degrees
    moon_ra = math.degrees(moon.ra)
    moon_dec = math.degrees(moon.dec)

    # Use the spherical law of cosines to calculate angular separation
    ra_diff = math.radians(moon_ra - target.ra)
    moon_dec_rad = math.radians(moon_dec)
    target_dec_rad = math.radians(target.dec)

    # Formula for angular separation in degrees
    angle = math.acos(
        math.sin(moon_dec_rad) * math.sin(target_dec_rad) +
        math.cos(moon_dec_rad) * math.cos(target_dec_rad) * math.cos(ra_diff)
    )

    # Convert angle from radians to degrees
    angular_distance = "{:.4f}".format(math.degrees(angle))

    return angular_distance


def is_moon_above_horizon(observer):

    if observer.date is None:
        print('appears no date was supplied')
        observer.date = datetime.datetime.now(timezone.utc)

    # Get the Moon's current position
    moon = ephem.Moon(observer)

    # Get the Moon's altitude (in degrees)
    moon_altitude = moon.alt

    # Check if the Moon is above the horizon
    is_above = moon_altitude > 0  # True if altitude is positive (above horizon)
    degrees_total = math.degrees(moon_altitude)

    return is_above, degrees_total


def format_time_difference(seconds):
    # Convert time difference into days, hours, minutes, and seconds
    days, remainder = divmod(seconds, 86400)  # 86400 seconds in a day
    hours, remainder = divmod(remainder, 3600)  # 3600 seconds in an hour
    minutes, seconds = divmod(remainder, 60)  # 60 seconds in a minute

    # Build a list of time components
    time_parts = []

    if days > 0:
        time_parts.append(f"{int(days)} days")

    if hours > 1:
        time_parts.append(f"{int(hours)} hours")
    elif hours == 1:
        time_parts.append(f"{int(hours)} hour")

    if minutes > 1:
        time_parts.append(f"{int(minutes)} minutes")
    elif minutes == 1:
        time_parts.append(f"{int(minutes)} minute")

    if seconds > 1:
        time_parts.append(f"{int(seconds)} seconds")
    elif seconds == 1:
        time_parts.append(f"{int(seconds)} second")

    time_string = ', '.join(time_parts)
    # Join the time components with commas and return the result
    return time_string


# ASSIGNMENT USUALLY FOR DRAG AND DROP SEESTAR DEVICE MAIN TARGET FOLDER INTO PY PROG DESKTOP ICON
# PY PROGRAM CURRENTLY IN FOLDER AT SAME LEVEL AS TARGET MAIN FOLDERS
# SEESTAR MAIN TARGET FOLDERS CREATED LOCALLY IF NOT EXIST AND FILES COPIES MINUS THNS
# SEESTAR TARGET SUBS ANALYSED FOR THEIR DT TM AND LOCAL (SESSION) FOLDER CREATED WITH SAME DATATIME AND TARGET ID UNDER TARGET MAIN FOLDER
# SEESTAR FIT SUBS ONLY MOVED OVER TO LOCAL SESSION FOLDER
# CODE BELOW DOES NOT INCLUDE THIS CODE. LEAVING USER TO DECIDE WHERE SUBS SUBSIDE (NO PUN INTENTED)
# SOURCE BELOW CONSIDERED SEESTAR FOLDERS...DESTINATION CONSIDERED LOCAL STORED/ARCHIVED LOCATION


if len(sys.argv) > 1:
    # input("There is sys.arg avail... assigning to source")
    source = sys.argv[1]
else:
    # input("There is no sys.arg avail going to assign source")
    source_directory = 'IC 4587'
    # source_subdirectory = 'IC 4587 sub 20241118-183742'


base_directory_path = os.getcwd()

destination_directory = source_directory
destination_subdirectory = destination_directory + ' sub 20241118-183742'
destination_directory_path = base_directory_path + '\\' + destination_directory
destination_subdirectory_path = destination_directory_path + '\\' + destination_subdirectory

print(f'destination_subdirectory_path : {destination_subdirectory_path}')
input("press any key to continue..")

fit_files = '*.fit'
file_location = destination_subdirectory_path + '\\' + fit_files
num_sess_subs = 0
for f in glob.glob(file_location):
    num_sess_subs += 1
files = [f for f in glob.glob(file_location)]
files.sort(key=lambda x: os.path.getmtime(os.path.join(destination_subdirectory_path, x)), reverse=True)

if files:
    no_files = len(files)
    earliest_file = files[no_files-1]
    earliest_file_mod_time = os.path.getmtime(os.path.join(destination_subdirectory_path, earliest_file))
    earliest_file_mod_time_formatted = datetime.fromtimestamp(earliest_file_mod_time).strftime('%Y-%m-%d %H:%M:%S')
    middle_file = files[int(no_files/2)]
    middle_file_mod_time = os.path.getmtime(os.path.join(destination_subdirectory_path, middle_file))
    middle_file_mod_time_formatted = datetime.fromtimestamp(middle_file_mod_time).strftime('%Y-%m-%d %H:%M:%S')
    latest_file = files[0]
    latest_file_mod_time = os.path.getmtime(os.path.join(destination_subdirectory_path, latest_file))
    latest_file_mod_time_formatted = datetime.fromtimestamp(latest_file_mod_time).strftime('%Y-%m-%d %H:%M:%S')
else:
    print("No files found in the directory.")


time_difference_in_seconds = round(latest_file_mod_time - earliest_file_mod_time)
formatted_time = format_time_difference(time_difference_in_seconds)
print(f"Time difference : {time_difference_in_seconds} seconds or {formatted_time}")

# DEFINE WHICH FILE TO HAVE HEADER INFO GATHERED FROM

measured_file = middle_file
average_fwhm = 2.73948            # THIS VALUE IS TAKEN FROM ANALYSIS OF MIDDLE FIT SUB IMAGE. CODE NOT IN THIS SNIPPET ADDITIONAL METRICS NEAR FUTURE

# GET FITS HEADER INFORMATION FOR TARGET LOG RECORD

# input("about to gather fits header info...")
# get information from the fits header of the latest file
fits_header = fits.getheader(measured_file)

ave_fwhm = "{:.2f}".format(average_fwhm)
efficiency = ((num_sess_subs*fits_header['EXPOSURE'])/time_difference_in_seconds) * 100
efficiency = "{:.1f}".format(efficiency)
print("efficiency is: " + efficiency)
fits_header['EXPOSURE'] = "{:.0f}".format(fits_header['EXPOSURE'])

print("EXPOSURE is: ", fits_header['EXPOSURE'])
print("object was : ", fits_header['OBJECT'])

fits_header['CCD-TEMP'] = "{:.1f}".format(fits_header['CCD-TEMP'])
print("CCD-TEMP is: ", fits_header['CCD-TEMP'])
print("filter was : ", fits_header['FILTER'])
print("gain was : ", fits_header['GAIN'])
print("focuspos was : ", fits_header['FOCUSPOS'])
print("SiteLong was : ", fits_header['SITELONG'])
print("SiteLat was : ", fits_header['SITELAT'])
print("Observation Date was : ", fits_header['DATE-OBS'])


# GATHER OBSERVER INFORMATION AND MOON POSITION, PHASE AND ANGULAR DISTANCE FROM TARGET 

print("datetime.now(timezone.utc) :", {datetime.now(timezone.utc)})
datetime_object = datetime.strptime(fits_header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')

# following gathers image centre properties 

observer_properties = setup_observer(fits_header['SITELAT'], fits_header['SITELONG'], 100, datetime_object)   # all values from fit file header
print(f'observer_properties .lat: {observer_properties.lat} .lon: {observer_properties.lon} .elevation: {observer_properties.elevation} .date: {observer_properties.date}')

image_centre = Target(fits_header['OBJECT'], fits_header['RA'], fits_header['DEC'])  # all values from fit file header
print(f'image_centre properties .name: {image_centre.name} .ra: {image_centre.ra} .dec: {image_centre.dec}')


# GET MOON STATUS AND ALTITUDE ON OBSERVATION DATE

moon_status, moon_altitude_degrees = is_moon_above_horizon(observer_properties)


# GET LUNAR PERCENT ILLUMINATION AT DATE AND TIME OF IMAGE

illumination = get_moon_illumination(observer_properties)


# GET ANGULAR SEPARATION BETWEEN MOON AND IMAGE_CENTRE AT DATE AND TIME OF IMAGE

angular_separation = calculate_angular_separation(observer_properties, image_centre)


# FORMAT DATA FOR CSV FILE

fits_header['RA'] = "{:.4f}".format(fits_header['RA'])
fits_header['DEC'] = "{:.4f}".format(fits_header['DEC'])

passes = 1      # passes CURRENTLY HARDCODED : Which Siril Script is later run via Siril. Passes to be calculated in future versions based on observing parameters

# STRUCTURE DATA FOR TARGET LOG FILE

header = [
    ("Object Directory", "Passes", "Earliest File", "Latest File", "Num Sess Subs", "Exp (sec)", "Sens Temp °C", "Filter", "Gain", "Focus", "RA", "DEC", "EFF %", "MOON POS", "MOON ILL %", "OBJ SEP °", "AVE FWHM")
]


data = [
    (destination_subdirectory+" Labelled", passes, earliest_file_mod_time_formatted, latest_file_mod_time_formatted, num_sess_subs, fits_header['EXPOSURE'], fits_header['CCD-TEMP'], fits_header['FILTER'], fits_header['GAIN'], fits_header['FOCUSPOS'], fits_header['RA'], fits_header['DEC'], efficiency, moon_status, illumination, angular_separation, ave_fwhm)
]


base_path = r"S:"  # WHEN USING COMPUTER IN MAIN ROOM
# base_path = r"F:"    # WHEN USING MAIN PROCESSING COMPUTER IE SHARED LOG FILE DUE TO LOCAL MAPPING
file_name = "/folderinfo_results_TBD.csv"
file_name_w_path = base_path + file_name
file_exists = os.path.isfile(file_name_w_path)
print(file_exists)
print(file_name_w_path, " to be appended!")
# input("opening output filename...")
with open(file_name_w_path, mode='a', newline='') as file:
    writer = csv.writer(file, delimiter='|')
    if not file_exists:
        # input("file doesnt exist so will create header...")
        writer.writerows(header) 
   
    writer.writerows(data)

print(file_name_w_path, " appended successfully!")
